use crate::component_meta_graph::ComponentMetaGraph;
use crate::global_demand_coordinator::GlobalDemandCoordinator;
use crate::graph::Graph;
use crate::macro_splicer::{splice_macro_tour, verify_tour_on_raw_graph};
use crate::pinpointed_strip_solver::PinpointedStripSolver;
use crate::two_tier_decomposer::decompose_graph;
use std::collections::{HashMap, HashSet};
use std::path::Path;
use std::time::Instant;

pub struct TwoTierSolverOptions {
    pub timeout_secs: f64,
    pub max_iterations: usize,
    pub enable_patching: bool,
    pub output_path: Option<String>,
}

impl Default for TwoTierSolverOptions {
    fn default() -> Self {
        Self {
            timeout_secs: 1800.0,
            max_iterations: 50_000,
            enable_patching: true,
            output_path: None,
        }
    }
}

#[derive(Debug, Clone)]
pub struct TwoTierOptions {
    pub timeout_secs: f64,
    pub output_tour: Option<String>,
    pub cycle_repair: bool,
}

impl Default for TwoTierOptions {
    fn default() -> Self {
        Self {
            timeout_secs: 1800.0,
            output_tour: None,
            cycle_repair: true,
        }
    }
}

pub struct TwoTierOrchestrator;

impl TwoTierOrchestrator {
    pub fn solve(g: &Graph, options: &TwoTierOptions) -> Option<Vec<i32>> {
        let opt = TwoTierSolverOptions {
            timeout_secs: options.timeout_secs,
            max_iterations: 50_000,
            enable_patching: options.cycle_repair,
            output_path: options.output_tour.clone(),
        };
        solve_graph_two_tier(g, &opt)
    }
}

/// Writes a certified tour in standard TSPLIB/HCP format.
pub fn write_hcp_tour(tour: &[i32], output_path: &str) -> std::io::Result<()> {
    let name = Path::new(output_path)
        .file_stem()
        .and_then(|s| s.to_str())
        .unwrap_or("tour");
    crate::tour_verifier::TourVerifier::write_tsplib_hcp(tour, name, output_path)
}

/// Top-level solver function integrating:
/// 1. decompose_graph(&g)
/// 2. GlobalDemandCoordinator
/// 3. PinpointedStripSolver
/// 4. splice_macro_tour + patch_cycles_2opt
/// 5. Indicator Cut-Crossing SECs on subtours
/// 6. Independent tour certification & output to HCP format
pub fn solve_graph_two_tier(g: &Graph, options: &TwoTierSolverOptions) -> Option<Vec<i32>> {
    let start_time = Instant::now();
    println!("=== Starting Two-Tier Demand-Coordinated Solver in Rust ===");

    let decomp = decompose_graph(g);
    let large_strips_count = decomp.strips.iter().filter(|s| s.len() >= 10).count();
    println!(
        "Decomposition: {} hubs ({} S, {} B, {} M), {} strips ({} large)",
        decomp.all_hubs.len(),
        decomp.s_hubs.len(),
        decomp.b_hubs.len(),
        decomp.m_hubs.len(),
        decomp.strips.len(),
        large_strips_count
    );

    let enable_mtz = decomp.all_hubs.len() >= 2;
    if enable_mtz {
        println!(
            "Active Chinese Remainder Encoding (CRT [2, 3, 7]) enabled on {} hubs",
            decomp.all_hubs.len()
        );
    }

    let mut coordinator = GlobalDemandCoordinator::new_with_mtz(g, &decomp, enable_mtz);
    let mut strip_solver = PinpointedStripSolver::new(g, &decomp);
    let mut strip_cache: HashMap<(usize, Vec<(i32, usize)>, usize), Result<Vec<Vec<i32>>, Vec<i32>>> = HashMap::new();

    let mut outer_it = 0;
    while outer_it < options.max_iterations {
        outer_it += 1;
        let elapsed = start_time.elapsed().as_secs_f64();
        if elapsed > options.timeout_secs {
            println!(
                "[TIMEOUT] Reached global {:.1}s limit at iteration {}",
                options.timeout_secs, outer_it
            );
            return None;
        }

        if outer_it % 10 == 1 || outer_it <= 5 {
            println!("\n--- Outer Iteration {} ({:.1}s) ---", outer_it, elapsed);
        }

        let assignment = coordinator.solve_assignment();
        if assignment.is_none() {
            println!("Coordinator returned UNSAT: search space exhausted.");
            return None;
        }

        let (hh_edges, strip_demands) = assignment.unwrap();
        if outer_it % 10 == 1 || outer_it <= 5 {
            println!(
                "Coordinator assigned {} Hub-Hub edges across {} strips",
                hh_edges.len(),
                decomp.strips.len()
            );
        }

        let mut all_strips_sat = true;
        let mut strip_paths = HashMap::new();

        for (si, s) in decomp.strips.iter().enumerate() {
            if start_time.elapsed().as_secs_f64() > options.timeout_secs {
                println!("[TIMEOUT] Global timeout reached during strip solving.");
                return None;
            }

            let dem = strip_demands.get(&si).cloned().unwrap_or_default();
            let tot_d: usize = dem.values().sum();
            let k = if tot_d >= 2 {
                tot_d / 2
            } else if s.len() < 10 {
                1
            } else {
                4
            };

            let s_hub = decomp.strip_adj_hubs.get(&si).and_then(|adj| {
                decomp.s_hubs.iter().find(|h| adj.contains(h)).copied()
            });
            let b_hub = decomp.strip_adj_hubs.get(&si).and_then(|adj| {
                decomp.b_hubs.iter().find(|h| adj.contains(h)).copied()
            });

            let mut dem_vec: Vec<(i32, usize)> = dem.iter().map(|(&h, &d)| (h, d)).filter(|&(_, d)| d > 0).collect();
            dem_vec.sort_unstable();
            let cache_key = (si, dem_vec, k);

            let res = if let Some(cached) = strip_cache.get(&cache_key) {
                cached.clone()
            } else {
                let r = strip_solver.solve_strip(si, &dem, s_hub, b_hub, k);
                strip_cache.insert(cache_key, r.clone());
                r
            };

            match res {
                Ok(paths) => {
                    strip_paths.insert(si, paths);
                }
                Err(failed_core) => {
                    all_strips_sat = false;
                    coordinator.add_conflict_clause(si, &dem, &failed_core);
                    if outer_it % 10 == 1 || outer_it <= 5 {
                        println!(
                            "  Strip {:2} ({}v) UNSAT with core {:?} -> conflict learned",
                            si,
                            s.len(),
                            failed_core
                        );
                    }
                    break;
                }
            }
        }

        if !all_strips_sat {
            continue;
        }

        println!(
            "All {} strips SATISFIED at iter {} ({:.1}s)! Splicing full tour...",
            decomp.strips.len(),
            outer_it,
            start_time.elapsed().as_secs_f64()
        );

        let (is_single, cycles) = splice_macro_tour(
            g,
            &decomp,
            &hh_edges,
            &strip_paths,
            &strip_demands,
            options.enable_patching,
        );

        if is_single && cycles.len() == 1 {
            let tour = &cycles[0];
            println!(
                "SUCCESS! Single Hamiltonian tour formed with {} vertices!",
                tour.len()
            );
            if verify_tour_on_raw_graph(tour, g) {
                println!("CERTIFICATION PASSED: Verified tour independently on raw graph G!");
                let out_path = options
                    .output_path
                    .as_deref()
                    .unwrap_or("scratch/graph950/found_tour_rust.hcp");
                if let Err(e) = write_hcp_tour(tour, out_path) {
                    eprintln!("Warning: failed to write tour to {}: {}", out_path, e);
                } else {
                    println!("Wrote certified tour to {}", out_path);
                }
                println!(
                    "Total time elapsed: {:.2}s",
                    start_time.elapsed().as_secs_f64()
                );
                return Some(tour.clone());
            } else {
                eprintln!("Verification on raw graph failed.");
                return None;
            }
        } else {
            if !cycles.is_empty() {
                println!(
                    "Splicer detected {} disconnected subtours -> adding macro cut clauses",
                    cycles.len()
                );
                if cycles.len() > 1 {
                    let meta_graph = ComponentMetaGraph::build(&cycles, g);
                    if meta_graph.meta_components.len() > 1 {
                        println!(
                            "Meta-graph partitioned into {} disconnected components -> generating multi-component SEC cuts",
                            meta_graph.meta_components.len()
                        );
                        coordinator.add_meta_component_cuts(meta_graph.get_meta_components(), &cycles);
                    }
                }
                coordinator.add_exact_subtour_block(&hh_edges, &strip_demands, &cycles);
                for cyc in &cycles {
                    let cyc_verts: HashSet<i32> = cyc.iter().copied().collect();
                    coordinator.add_macro_cut(&cyc_verts);
                    coordinator.inject_component_mtz(cyc);
                }
            } else {
                println!("Splicer failed boundary matching -> learning conflict clauses on strips");
                for (si, _) in decomp.strips.iter().enumerate() {
                    if let Some(dem) = strip_demands.get(&si) {
                        let failed_hubs: Vec<i32> = dem.keys().copied().collect();
                        coordinator.add_conflict_clause(si, dem, &failed_hubs);
                    }
                }
            }
        }
    }

    None
}
