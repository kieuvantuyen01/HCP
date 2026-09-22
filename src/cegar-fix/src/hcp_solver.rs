use rustsat::instances::*;
use rustsat::types::*;
use rustsat::clause;
use std::collections::{BTreeMap, HashMap, HashSet, VecDeque};
use std::time::{Instant,Duration};
// use crate::encoder;
use crate::graph::*;
use crate::encoder::*;
use crate::file_operations;
use crate::contraction::Degree2Contractor;
use crate::hub_registry::HubRegistry;
use crate::patching::HubPatcher;
use crate::matching_patcher::MatchingPatcher;
use crate::chained_lk::ChainedLKSolver;
use crate::ils_patcher::IteratedLocalSearchPatcher;
use crate::macro_solver::MacroGraphSolver;
use crate::hub_sub_hcp::HubPartitionedSolver;
use crate::modular_solver::ModularSolver;
use crate::subcycle_absorber::SubcycleAbsorber;
use crate::cycle_chain_absorber::CycleChainAbsorber;
use crate::backbone_freezer::BackboneFreezer;
use crate::snark_bridge::SnarkBridgeEngine;
use crate::gadget_parity::GadgetInterfaceParityEngine;
use crate::cut_selector::{CutSelector, CutSelectorOptions};
use crate::solver_reseeder::{SolverReseeder, ReseederOptions};
use crate::hemisphere_splicer::HemisphereSplicer;
use crate::static_cycle_cutter::StaticCycleCutter;
use crate::boundary_alternating_patcher::BoundaryAlternatingPatcher;
use crate::metagraph_router::MetagraphRouter;
use crate::parallel_sat_portfolio::{ParallelSatPortfolio, PortfolioResult};
use crate::giant_cycle_stitcher::GiantCycleStitcher;
use crate::interface_port_synchronizer::InterfacePortSynchronizer;
use crate::inverse_3sat_synthesizer::Inverse3SatSynthesizer;
use crate::hub_hierarchical_decomposer::HubHierarchicalDecomposer;
use crate::empirical_backbone_cutter::{EmpiricalBackboneCutter, EmpiricalBackboneTracker};
use crate::cnf_subsumer::CnfSubsumer;
use crate::twin_giant_splicer::TwinGiantSplicer;
use crate::sat_macro_patcher::SatMacroPatcher;
use crate::gadget_path_absorber::GadgetPathAbsorber;
use crate::incremental_sat::IncrementalSatSolver;
use crate::macro_crossover_splicer::MacroCrossoverSplicer;
use crate::quotient_block_cutter::QuotientBlockCutter;
use crate::bipartite_module_detector::BipartiteModuleDetector;
use crate::module_dual_path_extractor::ModuleDualPathExtractor;
use crate::module_state_cnf_encoder::ModuleStateCnfEncoder;
use crate::balanced_pair_cutset::BalancedPairCutset;
use crate::localized_sat_repair::LocalizedSatRepair;
use crate::alternating_port_engine::AlternatingPortEngine;
use crate::port_corridor_lns::PortCorridorLns;
use crate::modular_ring_dp_solver::RingDpSolver;




/// Pre-emptively forbids 3-cycles (triangles) and 4-cycles in the initial CNF encoding in O(|E| * Delta).
pub fn add_global_short_cycle_cuts(
    g: &Graph,
    encoder: &Encoder,
    cnf: &mut Cnf,
    max_cycle_len: usize,
) -> usize {
    let n = g.adjacency_list.len();
    if n <= 3 || max_cycle_len < 3 {
        return 0;
    }

    let mut adj_set: HashMap<i32, HashSet<i32>> = HashMap::new();
    for (&u, neighbors) in &g.adjacency_list {
        adj_set.insert(u, neighbors.iter().cloned().collect());
    }

    let mut added_clauses = 0;

    // 1. Triangle Cuts (3-cycles)
    if max_cycle_len >= 3 && n > 3 {
        for (&u, neighbors) in &g.adjacency_list {
            if neighbors.len() > 30 {
                continue;
            }
            for &v in neighbors {
                if u >= v {
                    continue;
                }
                let v_neighbors = match g.adjacency_list.get(&v) {
                    Some(s) => s,
                    None => continue,
                };
                if v_neighbors.len() > 30 {
                    continue;
                }
                for &w in v_neighbors {
                    if v >= w || w == u {
                        continue;
                    }
                    if let Some(u_nbrs) = adj_set.get(&u) {
                        if u_nbrs.contains(&w) {
                            if let (Some(&x_uv), Some(&x_vw), Some(&x_wu)) = (
                                encoder.graph_lit_map.get(&(u, v)),
                                encoder.graph_lit_map.get(&(v, w)),
                                encoder.graph_lit_map.get(&(w, u)),
                            ) {
                                cnf.add_clause(clause!(!x_uv, !x_vw, !x_wu));
                                added_clauses += 1;
                            }
                            if let (Some(&x_uw), Some(&x_wv), Some(&x_vu)) = (
                                encoder.graph_lit_map.get(&(u, w)),
                                encoder.graph_lit_map.get(&(w, v)),
                                encoder.graph_lit_map.get(&(v, u)),
                            ) {
                                cnf.add_clause(clause!(!x_uw, !x_wv, !x_vu));
                                added_clauses += 1;
                            }
                        }
                    }
                }
            }
        }
    }

    // 2. Quad Cuts (4-cycles)
    if max_cycle_len >= 4 && n > 4 {
        for (&u, neighbors) in &g.adjacency_list {
            if neighbors.len() > 30 {
                continue;
            }
            let u_neighbors: Vec<i32> = neighbors.iter().filter(|&&v| v > u).copied().collect();
            for i in 0..u_neighbors.len() {
                let v = u_neighbors[i];
                let v_nbrs = match adj_set.get(&v) {
                    Some(s) => s,
                    None => continue,
                };
                for j in (i + 1)..u_neighbors.len() {
                    let z = u_neighbors[j];
                    let z_nbrs = match adj_set.get(&z) {
                        Some(s) => s,
                        None => continue,
                    };
                    for &w in v_nbrs {
                        if w > u && w != z && z_nbrs.contains(&w) {
                            if let (Some(&x_uv), Some(&x_vw), Some(&x_wz), Some(&x_zu)) = (
                                encoder.graph_lit_map.get(&(u, v)),
                                encoder.graph_lit_map.get(&(v, w)),
                                encoder.graph_lit_map.get(&(w, z)),
                                encoder.graph_lit_map.get(&(z, u)),
                            ) {
                                cnf.add_clause(clause!(!x_uv, !x_vw, !x_wz, !x_zu));
                                added_clauses += 1;
                            }
                            if let (Some(&x_uz), Some(&x_zw), Some(&x_wv), Some(&x_vu)) = (
                                encoder.graph_lit_map.get(&(u, z)),
                                encoder.graph_lit_map.get(&(z, w)),
                                encoder.graph_lit_map.get(&(w, v)),
                                encoder.graph_lit_map.get(&(v, u)),
                            ) {
                                cnf.add_clause(clause!(!x_uz, !x_zw, !x_wv, !x_vu));
                                added_clauses += 1;
                            }
                        }
                    }
                }
            }
        }
    }

    added_clauses
}

/// Adds cardinality cut constraints on satellite clusters (components of degree < 20 vertices with size >= 50)
/// to prune the search space for CaDiCaL SAT solver before execution.
/// Returns the total number of added cluster cut clauses.
pub fn add_cluster_cut_constraints(
    g: &Graph,
    encoder: &Encoder,
    cnf: &mut Cnf,
) -> usize {
    let n = g.adjacency_list.len();
    if n < 50 {
        return 0;
    }

    // 1. Identify satellite vertices (degree < 20)
    let mut satellite_vertices: Vec<i32> = Vec::new();
    for (&u, neighbors) in &g.adjacency_list {
        if neighbors.len() < 20 {
            satellite_vertices.push(u);
        }
    }
    satellite_vertices.sort_unstable();

    let satellite_set: HashSet<i32> = satellite_vertices.iter().copied().collect();
    let mut visited: HashSet<i32> = HashSet::new();
    let mut clusters: Vec<HashSet<i32>> = Vec::new();

    // 2. Find connected components among satellite vertices in G[S]
    for &u in &satellite_vertices {
        if visited.contains(&u) {
            continue;
        }
        let mut cluster = HashSet::new();
        let mut queue = VecDeque::new();
        visited.insert(u);
        queue.push_back(u);

        while let Some(curr) = queue.pop_front() {
            cluster.insert(curr);
            if let Some(neighbors) = g.adjacency_list.get(&curr) {
                for &nbr in neighbors {
                    if satellite_set.contains(&nbr) && !visited.contains(&nbr) {
                        visited.insert(nbr);
                        queue.push_back(nbr);
                    }
                }
            }
        }

        // Filter clusters of size >= 50 and proper subset of V
        if cluster.len() >= 50 && cluster.len() < n {
            clusters.push(cluster);
        }
    }

    if clusters.is_empty() {
        return 0;
    }

    let mut added_clauses = 0;

    // 3. For each cluster C_i, generate cut clauses
    for cluster in &clusters {
        let mut out_lits: Vec<Lit> = Vec::new();
        let mut in_lits: Vec<Lit> = Vec::new();

        for &u in cluster {
            if let Some(neighbors) = g.adjacency_list.get(&u) {
                for &v in neighbors {
                    if !cluster.contains(&v) {
                        if let Some(&lit_out) = encoder.graph_lit_map.get(&(u, v)) {
                            out_lits.push(lit_out);
                        }
                        if let Some(&lit_in) = encoder.graph_lit_map.get(&(v, u)) {
                            in_lits.push(lit_in);
                        }
                    }
                }
            }
        }

        // Add positive out-cut clause: at least 1 edge leaves cluster
        if !out_lits.is_empty() {
            let mut cl_out = Clause::new();
            cl_out.extend(out_lits);
            cnf.add_clause(cl_out);
            added_clauses += 1;
        }

        // Add positive in-cut clause: at least 1 edge enters cluster
        if !in_lits.is_empty() {
            let mut cl_in = Clause::new();
            cl_in.extend(in_lits);
            cnf.add_clause(cl_in);
            added_clauses += 1;
        }
    }

    added_clauses
}

pub fn solve_hamilton(g:Graph, contractor: &Degree2Contractor, hub_registry: &HubRegistry, s:i32, encode_method:i32, block_method: i32,symmetry: i32 ,opt:i32,loop_prohibition: i32,cnf_normalize:i32,balanced:i32,dearcify:i32, cadical_config:i32, degree_order:i32, arcs_order:i32, three_opt:i32, cegar_fallback:i32, mtz_stall:i32, adaptive_escalation:i32, sub_hcp_timeout: u64, max_cluster_size: usize, timeout_secs: f64, instant:Instant,output_folder:&str, macro_gadget: i32, bounded_freezer: i32, alternating_engine: i32) -> Option<Vec<i32>> {
    solve_hamilton_with_decomposition(
        g, contractor, hub_registry, s, encode_method, block_method, symmetry,
        opt, loop_prohibition, cnf_normalize, balanced, dearcify,
        cadical_config, degree_order, arcs_order, three_opt, cegar_fallback,
        mtz_stall, adaptive_escalation, sub_hcp_timeout, max_cluster_size,
        timeout_secs, instant, output_folder, macro_gadget, bounded_freezer,
        alternating_engine, true, true,
    )
}

pub fn solve_hamilton_with_decomposition(g:Graph, contractor: &Degree2Contractor, hub_registry: &HubRegistry, _s:i32, encode_method:i32, block_method: i32,symmetry: i32 ,opt:i32,loop_prohibition:i32,cnf_normalize:i32,balanced:i32,dearcify:i32, cadical_config:i32, degree_order:i32, arcs_order:i32, three_opt:i32, _cegar_fallback:i32, _mtz_stall:i32, _adaptive_escalation:i32, _sub_hcp_timeout: u64, _max_cluster_size: usize, timeout_secs: f64, instant:Instant,output_folder:&str, macro_gadget: i32, bounded_freezer: i32, alternating_engine: i32, structural_decomposition: bool, cycle_repair: bool) -> Option<Vec<i32>> {
    let now = instant.elapsed();

    // Fast Track: Inverse 3-SAT De-reduction & Tour Synthesis
    if structural_decomposition {
        if let Some(synthesized_tour) = Inverse3SatSynthesizer::try_solve_via_inverse_3sat(&g) {
            println!("Inverse3SatSynthesizer: successfully de-reduced graph to 3-SAT and synthesized valid Hamiltonian tour!");
            return Some(contractor.expand_tour(&synthesized_tour));
        }
    }

    // Fast Track: Hub-Centric Hierarchical Decomposition
    if structural_decomposition {
        if let Some(hierarchical_tour) = HubHierarchicalDecomposer::try_solve_hierarchical(&g) {
            println!("HubHierarchicalDecomposer: successfully solved graph via 2-tier hub hierarchy!");
            return Some(contractor.expand_tour(&hierarchical_tour));
        }
    }

    // Fast Track: Modular Ring DP / Quotient Solver
    if structural_decomposition && !contractor.chain_map.is_empty() && hub_registry.hub_vertices.is_empty() {
        if let Some(ring_tour) = RingDpSolver::solve_contracted(&g, contractor) {
            println!("RingDpSolver: successfully solved graph via Modular Ring Quotient Solver!");
            let line = ring_tour.iter().map(|i| i.to_string()).collect::<Vec<String>>().join(" ");
            println!();
            println!("solution: ");
            println!("{}\n", line);
            println!("s SATISFIABLE");
            return Some(ring_tour);
        }
    }

    let mut encoder = Encoder::new();
    // グラフをcnf形式に変形し、cnfへ格納
    let mut cnf = encoder.encode(&g,encode_method,symmetry,loop_prohibition,dearcify,degree_order,arcs_order);
    // Add mandatory edge constraints for contracted degree-2 paths
    for (&(u, w), _) in &contractor.chain_map {
        if u < w {
            if let (Some(&lit_uw), Some(&lit_wu)) = (encoder.graph_lit_map.get(&(u, w)), encoder.graph_lit_map.get(&(w, u))) {
                cnf.add_clause(clause!(lit_uw, lit_wu));
            }
        }
    }

    // Snark Key-Bridge Unit Locking
    if let Some((u, v, lit)) = SnarkBridgeEngine::detect_and_extract_key_bridge(&g, &encoder) {
        println!("SnarkBridgeEngine: detected key bridge ({}, {}), locking bridge edge", u, v);
        if let Some(&rev_lit) = encoder.graph_lit_map.get(&(v, u)) {
            if rev_lit != lit {
                cnf.add_clause(clause!(lit, rev_lit));
            } else {
                cnf.add_clause(clause!(lit));
            }
        } else {
            cnf.add_clause(clause!(lit));
        }
    }

    // Static Substructure Cycle Cutter: inject small (3..=8) and extended (9..=16) subtour elimination clauses.
    // If graph has hubs or vertices with degree >= 10, throttle long static cuts near hubs to prevent clause explosion.
    let has_hubs = hub_registry.hub_vertices.len() > 0 || g.adjacency_list.values().any(|nbrs| nbrs.len() >= 10);
    let static_cuts = if has_hubs {
        StaticCycleCutter::generate_selective_static_cycle_cuts(&g, &encoder, 10)
    } else {
        StaticCycleCutter::generate_static_small_cycle_cuts(&g, &encoder)
    };
    if !static_cuts.is_empty() {
        println!("StaticCycleCutter: injected {} static cycle elimination clauses at Round 0 (hub_throttled: {})", static_cuts.len(), has_hubs);
        cnf.extend(static_cuts);
    }

    /*
    // Global Supernode MTZ Potential Encoding
    if g.adjacency_list.len() >= 50 {
        let target_k = 16;
        let target_size = (g.adjacency_list.len() / target_k).max(25);
        let modules = MetagraphRouter::detect_gadget_modules_with_size(&g, target_size);
        if modules.len() >= 4 && modules.len() <= 24 {
            println!("GlobalSupernodeMTZ: generated {} supernodes (target size {}), injecting global MTZ order encoding", modules.len(), target_size);
            MetagraphRouter::encode_supernode_mtz(&modules, &g, &mut encoder, &mut cnf);
        }
    }
    */

    // Interface Port Truth Assignment & Flow Synchronizer
    let dual_paths = if structural_decomposition {
        InterfacePortSynchronizer::extract_gadget_dual_paths(&g, 32)
    } else {
        Vec::new()
    };
    if dual_paths.len() >= 4 {
        println!("InterfacePortSynchronizer: detected {} gadget modules with dual T/F paths, injecting flow synchronization clauses", dual_paths.len());
        InterfacePortSynchronizer::encode_interface_port_synchronization(&dual_paths, &g, &mut encoder, &mut cnf);
    }

    // Bipartite Module State Equivalence Encoder (Approach 2)
    // Only apply to pure modular reduction graphs without hubs (Hubs = 0) where modules are cleanly detected and macro_gadget is enabled
    if structural_decomposition && macro_gadget != 0 && hub_registry.hub_vertices.is_empty() && g.adjacency_list.len() >= 44 && !contractor.chain_map.is_empty() {
        let bip_modules = BipartiteModuleDetector::detect_44_modules(&g, contractor);
        if !bip_modules.is_empty() {
            let mut extracted = Vec::new();
            for m in &bip_modules {
                if let Some((t_path, f_path)) = ModuleDualPathExtractor::extract_dual_paths(&m.vertices, &g, contractor) {
                    extracted.push((m, t_path, f_path));
                }
            }
            // Only inject when high coverage of modules is achieved (avoid partial module starvation)
            if extracted.len() == bip_modules.len() && !extracted.is_empty() {
                let mut total_added_module_clauses = 0;
                for (m, t_path, f_path) in extracted {
                    let added = ModuleStateCnfEncoder::encode_module_dual_state(
                        &m.vertices,
                        &t_path,
                        &f_path,
                        &g,
                        &mut encoder,
                        &mut cnf,
                    );
                    total_added_module_clauses += added;
                }
                println!(
                    "ModuleStateCnfEncoder: successfully synchronized {}/{} 44-vertex modules, injected {} state-equivalence clauses at Round 0",
                    bip_modules.len(),
                    bip_modules.len(),
                    total_added_module_clauses
                );
            }
        }
    }

    let current_cnf = if output_folder != "default" {
        //フォルダーの作成
        let _ = file_operations::create_folder_if_not_exists(output_folder);
        //cnfをファイルに出力する
        let output_file = format!("{}/increment0.cnf",output_folder);
        let _ = file_operations::write_dimacs(cnf.clone(),&output_file);
        //出力のために複製
        cnf.clone()
    }else{
        Cnf::new()
    };

    // 標準入力で -s の後の数字により、minisat,kissat,cadicalを選択する
    println!("encodhing time = {:?}",instant.elapsed()-now);
    println!();
    let base_cnf = if cnf_normalize == 1{
        let normalized_cnf = cnf.normalize();
        println!("encodhing clauses number = {}",normalized_cnf.len());
        normalized_cnf
    }else{
        println!("encodhing clauses number = {}",cnf.len());
        cnf
    };
    // cegar関数により、解を求め、increment数と追加したblock節の合計を返す
    let (increment, block, tour) = cegar(
        &mut encoder,
        0,
        0,
        g,
        contractor,
        hub_registry,
        block_method,
        opt,
        three_opt,
        instant,
        cnf_normalize,
        balanced,
        timeout_secs,
        instant.elapsed(),
        current_cnf,
        output_folder,
        base_cnf,
        cadical_config,
        bounded_freezer,
        alternating_engine,
        structural_decomposition,
        cycle_repair,
    );
    println!("overall incremented number = {}", increment);
    println!("overall number of added block clauses = {}", block);
    tour
}

fn print_tour(tour: &[i32], contractor: &Degree2Contractor) {
    let final_tour = contractor.uncontract_cycle(tour);
    let line = final_tour.iter().map(|i| i.to_string()).collect::<Vec<String>>().join(" ");
    println!();
    println!("solution: ");
    println!("{}\n", line);
    println!("s SATISFIABLE");
}

fn cegar(
    encoder: &mut Encoder,
    mut count: i32,
    mut clause_count: i32,
    g: Graph,
    contractor: &Degree2Contractor,
    hub_registry: &HubRegistry,
    block_method: i32,
    opt: i32,
    three_opt: i32,
    instant: Instant,
    cnf_normalize: i32,
    balanced: i32,
    timeout_secs: f64,
    mut previous_time: Duration,
    mut previous_cnf: Cnf,
    output_folder: &str,
    base_cnf: Cnf,
    _cadical_config: i32,
    bounded_freezer: i32,
    alternating_engine: i32,
    structural_decomposition: bool,
    cycle_repair: bool,
) -> (i32, i32, Option<Vec<i32>>) {
    // Attempt Modular Macro-Decomposition when dense hubs are detected
    if structural_decomposition && hub_registry.hub_vertices.len() >= 5 {
        if let Some(tour) = ModularSolver::solve_via_modular_decomposition(&g, contractor, hub_registry) {
            println!("s SATISFIABLE (via Modular Macro-Decomposition)");
            println!("overall incremented number = 0");
            let final_tour = contractor.uncontract_cycle(&tour);
            print_tour(&tour, contractor);
            return (0, 0, Some(final_tour));
        }
    }

    // Attempt Hub-Partitioned Divide-and-Conquer for Dense Hub Graphs
    if structural_decomposition && !hub_registry.hub_vertices.is_empty() && hub_registry.hub_vertices.len() >= 3 {
        if let Some(partition_tour) = HubPartitionedSolver::solve_via_hub_partition(&g, contractor, hub_registry) {
            if partition_tour.len() == g.adjacency_list.len() {
                println!("number of subcycles found = 1 (via hub-partitioned sub-hcp)");
                let final_tour = contractor.uncontract_cycle(&partition_tour);
                let line = final_tour.iter().map(|i| i.to_string()).collect::<Vec<String>>().join(" ");
                let time = instant.elapsed();
                println!("overall time = {:?}", time);
                println!();
                println!("solution: ");
                println!("{}\n", line);
                println!("s SATISFIABLE");
                return (0, 0, Some(final_tour));
            }
        }
    }

    let mut working_cnf = base_cnf.clone();
    let mut incremental_solver = IncrementalSatSolver::new(&base_cnf);
    let mut assumptions: Vec<Lit> = Vec::new();
    let mut phase_hints: Vec<Lit> = Vec::new();
    let mut accumulated_cut_cnfs: Vec<Cnf> = Vec::new();
    let total_v = g.adjacency_list.len();
    let mut reseeder_opts = ReseederOptions::default();
    if total_v >= 500 {
        reseeder_opts.max_sat_time_threshold_secs = 120.0;
        reseeder_opts.periodic_interval_rounds = 20;
    }
    let inc_timeout = if total_v >= 500 { 180.0 } else { 15.0 };
    let mut backbone_tracker = EmpiricalBackboneTracker::new(10);

    // Modular Quotient Block Pre-computation disabled
    let modular_blocks: Vec<std::collections::HashSet<i32>> = Vec::new();

    loop {
        if instant.elapsed().as_secs_f64() >= timeout_secs {
            println!("\ns UNKNOWN (TIMEOUT: {:.2}s reached >= {:.2}s limit)", instant.elapsed().as_secs_f64(), timeout_secs);
            return (count, clause_count, None);
        }

        // SAT solve: First try persistent incremental solver to preserve learned CDCL state
        // If incremental solve times out or is interrupted, escalate to 3-core diversified parallel portfolio
        let port_res = match incremental_solver.solve_with_timeout(&assumptions, &phase_hints, inc_timeout) {
            PortfolioResult::Sat(model_lits) => PortfolioResult::Sat(model_lits),
            PortfolioResult::Unsat => PortfolioResult::Unsat,
            PortfolioResult::Interrupted => {
                ParallelSatPortfolio::solve_portfolio_with_timeout(
                    &working_cnf,
                    &assumptions,
                    &phase_hints,
                    3,
                    count as usize,
                    instant,
                    timeout_secs,
                )
            }
        };
        let now = instant.elapsed();
        let sat_solving_time = now - previous_time;

        println!();
        println!("Increment...");
        println!("incremented number = {}", count);
        println!("sat solving time = {:?}", sat_solving_time);

        // 解がSATならば、ハミルトン閉路になっているかを調べる
        match port_res {
            PortfolioResult::Sat(model_lits) => {
                // どの辺が選択されたかの解析
                let sol_arcs = get_solution_arcs_from_lits(&model_lits, &encoder.graph_lit_map);
                // 閉路
                let sol_cycles = get_solution_cycles(sol_arcs);
                backbone_tracker.record_solution_edges(&sol_cycles);

                // 閉路が一つであれば、ハミルトン閉路なので解を出力
                if sol_cycles.len() == 1 {
                    let flat: Vec<i32> = sol_cycles.into_iter().flatten().collect();
                    let full_cycle = contractor.uncontract_cycle(&flat);
                    let line = full_cycle.iter().map(|i| i.to_string()).collect::<Vec<String>>().join(" ");
                    println!();
                    println!("solution: ");
                    println!("{}\n", line);
                    println!("s SATISFIABLE");
                    return (count, clause_count, Some(full_cycle));
                } else {
                    println!("number of subcycles found = {}", sol_cycles.len());
                    println!("sat solution cycle lengths map (length:number) = {:?}", map_cycle_lengths(&sol_cycles));

                    // Attempt Alternating Port Engine Cycle Compression
                    let sol_cycles = if cycle_repair && alternating_engine != 0 && !contractor.chain_map.is_empty() && sol_cycles.len() > 1 {
                        let repaired = AlternatingPortEngine::repair(&sol_cycles, &g, contractor, 5, 2000);
                        if repaired.len() < sol_cycles.len() {
                            println!("AlternatingPortEngine: compressed subcycles from {} down to {} cycles", sol_cycles.len(), repaired.len());
                        }
                        if repaired.len() == 1 && (repaired[0].len() == g.adjacency_list.len() || repaired[0].len() == contractor.original_vertices_count) {
                            println!("*** 100% HAMILTONIAN TOUR FOUND VIA ALTERNATING PORT ENGINE! ***");
                            let flat: Vec<i32> = repaired.into_iter().flatten().collect();
                            let final_tour = if flat.len() == contractor.original_vertices_count {
                                flat
                            } else {
                                contractor.uncontract_cycle(&flat)
                            };
                            let line = final_tour.iter().map(|i| i.to_string()).collect::<Vec<String>>().join(" ");
                            let time = now - previous_time;
                            let add_block_clauses_time = now - previous_time - sat_solving_time;
                            println!("number of added block clauses = {}", clause_count);
                            println!("add block clauses time = {:?}", add_block_clauses_time);
                            println!("increment time = {:?}", time);
                            println!();
                            println!("solution: ");
                            println!("{}\n", line);
                            println!("s SATISFIABLE");
                            return (count, clause_count, Some(final_tour));
                        }
                        repaired
                    } else {
                        sol_cycles
                    };

                    // Attempt Port Corridor LNS Subcycle Absorption
                    let sol_cycles = if cycle_repair && alternating_engine != 0 && !contractor.chain_map.is_empty() && sol_cycles.len() > 1 {
                        let absorbed = PortCorridorLns::repair(&sol_cycles, &g, contractor, encoder, &base_cnf);
                        if absorbed.len() < sol_cycles.len() {
                            println!("PortCorridorLns: absorbed satellite subcycles from {} down to {} cycles", sol_cycles.len(), absorbed.len());
                        }
                        if absorbed.len() == 1 && (absorbed[0].len() == g.adjacency_list.len() || absorbed[0].len() == contractor.original_vertices_count) {
                            println!("*** 100% HAMILTONIAN TOUR FOUND VIA PORT CORRIDOR LNS! ***");
                            let flat: Vec<i32> = absorbed.into_iter().flatten().collect();
                            let final_tour = if flat.len() == contractor.original_vertices_count {
                                flat
                            } else {
                                contractor.uncontract_cycle(&flat)
                            };
                            let line = final_tour.iter().map(|i| i.to_string()).collect::<Vec<String>>().join(" ");
                            let time = now - previous_time;
                            let add_block_clauses_time = now - previous_time - sat_solving_time;
                            println!("number of added block clauses = {}", clause_count);
                            println!("add block clauses time = {:?}", add_block_clauses_time);
                            println!("increment time = {:?}", time);
                            println!();
                            println!("solution: ");
                            println!("{}\n", line);
                            println!("s SATISFIABLE");
                            return (count, clause_count, Some(final_tour));
                        }
                        absorbed
                    } else {
                        sol_cycles
                    };

                    // Attempt Multi-Subcycle Hub Patching
                    let sol_cycles = if cycle_repair && sol_cycles.len() > 1 && !hub_registry.hub_vertices.is_empty() {
                        let patched = HubPatcher::patch_cycles_via_hubs(&sol_cycles, &g, contractor, hub_registry);
                        if patched.len() == 1 && patched[0].len() == g.adjacency_list.len() {
                            println!("number of subcycles found = 1 (via hub patching)");
                            let flat: Vec<i32> = patched.into_iter().flatten().collect();
                            let final_tour = contractor.uncontract_cycle(&flat);
                            let line = final_tour.iter().map(|i| i.to_string()).collect::<Vec<String>>().join(" ");
                            let time = now - previous_time;
                            let add_block_clauses_time = now - previous_time - sat_solving_time;
                            println!("number of added block clauses = {}", clause_count);
                            println!("add block clauses time = {:?}", add_block_clauses_time);
                            println!("increment time = {:?}", time);
                            println!();
                            println!("solution: ");
                            println!("{}\n", line);
                            println!("s SATISFIABLE");
                            return (count, clause_count, Some(final_tour));
                        }
                        patched
                    } else {
                        sol_cycles
                    };

                    // Attempt Maximum Matching Global Patching on remaining subcycles
                    let sol_cycles = if cycle_repair && sol_cycles.len() > 1 {
                        let patched = MatchingPatcher::patch_cycles_via_matching(&sol_cycles, &g, contractor, hub_registry);
                        if patched.len() == 1 && patched[0].len() == g.adjacency_list.len() {
                            println!("number of subcycles found = 1 (via matching patching)");
                            let flat: Vec<i32> = patched.into_iter().flatten().collect();
                            let final_tour = contractor.uncontract_cycle(&flat);
                            let line = final_tour.iter().map(|i| i.to_string()).collect::<Vec<String>>().join(" ");
                            let time = now - previous_time;
                            let add_block_clauses_time = now - previous_time - sat_solving_time;
                            println!("number of added block clauses = {}", clause_count);
                            println!("add block clauses time = {:?}", add_block_clauses_time);
                            println!("increment time = {:?}", time);
                            println!();
                            println!("solution: ");
                            println!("{}\n", line);
                            println!("s SATISFIABLE");
                            return (count, clause_count, Some(final_tour));
                        }
                        patched
                    } else {
                        sol_cycles
                    };

                    // Attempt Chained k-Opt / Lin-Kernighan Variable-Depth Patching
                    let sol_cycles = if cycle_repair && sol_cycles.len() > 1 {
                        let patched = ChainedLKSolver::patch_cycles_via_chained_lk(&sol_cycles, &g, contractor, hub_registry, 6);
                        if patched.len() == 1 && patched[0].len() == g.adjacency_list.len() {
                            println!("number of subcycles found = 1 (via chained lk patching)");
                            let flat: Vec<i32> = patched.into_iter().flatten().collect();
                            let final_tour = contractor.uncontract_cycle(&flat);
                            let line = final_tour.iter().map(|i| i.to_string()).collect::<Vec<String>>().join(" ");
                            let time = now - previous_time;
                            let add_block_clauses_time = now - previous_time - sat_solving_time;
                            println!("number of added block clauses = {}", clause_count);
                            println!("add block clauses time = {:?}", add_block_clauses_time);
                            println!("increment time = {:?}", time);
                            println!();
                            println!("solution: ");
                            println!("{}\n", line);
                            println!("s SATISFIABLE");
                            return (count, clause_count, Some(final_tour));
                        }
                        patched
                    } else {
                        sol_cycles
                    };

                    // Attempt Iterated Local Search (ILS) with Double-Bridge Kicks
                    let sol_cycles = if cycle_repair && sol_cycles.len() > 1 {
                        let patched = IteratedLocalSearchPatcher::solve_via_ils(&sol_cycles, &g, contractor, hub_registry, 200);
                        if patched.len() == 1 && patched[0].len() == g.adjacency_list.len() {
                            println!("number of subcycles found = 1 (via ils patching)");
                            let flat: Vec<i32> = patched.into_iter().flatten().collect();
                            let final_tour = contractor.uncontract_cycle(&flat);
                            let line = final_tour.iter().map(|i| i.to_string()).collect::<Vec<String>>().join(" ");
                            let time = now - previous_time;
                            let add_block_clauses_time = now - previous_time - sat_solving_time;
                            println!("number of added block clauses = {}", clause_count);
                            println!("add block clauses time = {:?}", add_block_clauses_time);
                            println!("increment time = {:?}", time);
                            println!();
                            println!("solution: ");
                            println!("{}\n", line);
                            println!("s SATISFIABLE");
                            return (count, clause_count, Some(final_tour));
                        }
                        patched
                    } else {
                        sol_cycles
                    };

                    // Attempt Macro-Graph Hierarchical Contraction Solver
                    let sol_cycles = if cycle_repair && sol_cycles.len() > 1 {
                        if let Some(macro_tour) = MacroGraphSolver::solve_via_macro_graph(&sol_cycles, &g, contractor, hub_registry) {
                            if macro_tour.len() == g.adjacency_list.len() {
                                println!("number of subcycles found = 1 (via macro-graph solver)");
                                let final_tour = contractor.uncontract_cycle(&macro_tour);
                                let line = final_tour.iter().map(|i| i.to_string()).collect::<Vec<String>>().join(" ");
                                let time = now - previous_time;
                                let add_block_clauses_time = now - previous_time - sat_solving_time;
                                println!("number of added block clauses = {}", clause_count);
                                println!("add block clauses time = {:?}", add_block_clauses_time);
                                println!("increment time = {:?}", time);
                                println!();
                                println!("solution: ");
                                println!("{}\n", line);
                                println!("s SATISFIABLE");
                                return (count, clause_count, Some(final_tour));
                            }
                        }
                        sol_cycles
                    } else {
                        sol_cycles
                    };

                    // Attempt Multi-Cycle Alternating Chain Splicer & Absorber
                    let sol_cycles = if cycle_repair && sol_cycles.len() > 1 {
                        let absorbed = CycleChainAbsorber::absorb_all(&sol_cycles, &g, contractor, hub_registry);
                        if absorbed.len() == 1 && absorbed[0].len() == g.adjacency_list.len() {
                            println!("number of subcycles found = 1 (via cycle chain absorber)");
                            let flat: Vec<i32> = absorbed.into_iter().flatten().collect();
                            let final_tour = contractor.uncontract_cycle(&flat);
                            let line = final_tour.iter().map(|i| i.to_string()).collect::<Vec<String>>().join(" ");
                            let time = now - previous_time;
                            let add_block_clauses_time = now - previous_time - sat_solving_time;
                            println!("number of added block clauses = {}", clause_count);
                            println!("add block clauses time = {:?}", add_block_clauses_time);
                            println!("increment time = {:?}", time);
                            println!();
                            println!("solution: ");
                            println!("{}\n", line);
                            println!("s SATISFIABLE");
                            return (count, clause_count, Some(final_tour));
                        }
                        absorbed
                    } else {
                        sol_cycles
                    };

                    // Attempt Hemisphere Splicing for 2..=4 macro-components
                    let sol_cycles = if cycle_repair && sol_cycles.len() >= 2 && sol_cycles.len() <= 4 {
                        if let Some(spliced) = HemisphereSplicer::try_direct_splice_all(&sol_cycles, &g, contractor) {
                            println!("HemisphereSplicer: directly spliced macro-components from {} to {} cycles", sol_cycles.len(), spliced.len());
                            if spliced.len() == 1 && spliced[0].len() == g.adjacency_list.len() {
                                println!("number of subcycles found = 1 (via direct hemisphere splicer)");
                                let flat: Vec<i32> = spliced.into_iter().flatten().collect();
                                let final_tour = contractor.uncontract_cycle(&flat);
                                let line = final_tour.iter().map(|i| i.to_string()).collect::<Vec<String>>().join(" ");
                                let time = now - previous_time;
                                let add_block_clauses_time = now - previous_time - sat_solving_time;
                                println!("number of added block clauses = {}", clause_count);
                                println!("add block clauses time = {:?}", add_block_clauses_time);
                                println!("increment time = {:?}", time);
                                println!();
                                println!("solution: ");
                                println!("{}\n", line);
                                println!("s SATISFIABLE");
                                return (count, clause_count, Some(final_tour));
                            }
                            spliced
                        } else {
                            sol_cycles
                        }
                    } else {
                        sol_cycles
                    };

                    // Attempt Multi-Hop Boundary Alternating Patcher for 2..=4 macro-components
                    let sol_cycles = if cycle_repair && sol_cycles.len() >= 2 && sol_cycles.len() <= 4 {
                        if let Some(patched) = BoundaryAlternatingPatcher::try_patch_macro_hemispheres(&sol_cycles, &g, contractor, 4) {
                            println!("BoundaryAlternatingPatcher: patched macro-hemispheres from {} to {} cycles", sol_cycles.len(), patched.len());
                            if patched.len() == 1 && patched[0].len() == g.adjacency_list.len() {
                                println!("number of subcycles found = 1 (via boundary alternating patcher)");
                                let flat: Vec<i32> = patched.into_iter().flatten().collect();
                                let final_tour = contractor.uncontract_cycle(&flat);
                                let line = final_tour.iter().map(|i| i.to_string()).collect::<Vec<String>>().join(" ");
                                let time = now - previous_time;
                                let add_block_clauses_time = now - previous_time - sat_solving_time;
                                println!("number of added block clauses = {}", clause_count);
                                println!("add block clauses time = {:?}", add_block_clauses_time);
                                println!("increment time = {:?}", time);
                                println!();
                                println!("solution: ");
                                println!("{}\n", line);
                                println!("s SATISFIABLE");
                                return (count, clause_count, Some(final_tour));
                            }
                            patched
                        } else {
                            sol_cycles
                        }
                    } else {
                        sol_cycles
                    };

                    // Attempt Universal Giant-Cycle SAT Stitcher & Multi-Cycle Alternating Repair
                    let sol_cycles = if cycle_repair && sol_cycles.len() > 1 && sol_cycles.len() <= 150 {
                        let protected_edges: HashSet<(i32, i32)> = contractor.chain_map.keys().copied().collect();
                        let stitched = GiantCycleStitcher::repair_until_fixed_point(&sol_cycles, &g, &protected_edges);
                        if stitched.len() < sol_cycles.len() {
                            println!("GiantCycleStitcher: stitched and absorbed subcycles from {} down to {} cycles", sol_cycles.len(), stitched.len());
                        }
                        if stitched.len() == 1 && stitched[0].len() == g.adjacency_list.len() {
                            println!("number of subcycles found = 1 (via giant cycle stitcher)");
                            let flat: Vec<i32> = stitched.into_iter().flatten().collect();
                            let final_tour = contractor.uncontract_cycle(&flat);
                            let line = final_tour.iter().map(|i| i.to_string()).collect::<Vec<String>>().join(" ");
                            let time = now - previous_time;
                            let add_block_clauses_time = now - previous_time - sat_solving_time;
                            println!("number of added block clauses = {}", clause_count);
                            println!("add block clauses time = {:?}", add_block_clauses_time);
                            println!("increment time = {:?}", time);
                            println!();
                            println!("solution: ");
                            println!("{}\n", line);
                            println!("s SATISFIABLE");
                            return (count, clause_count, Some(final_tour));
                        }
                        stitched
                    } else {
                        sol_cycles
                    };

                    // 2-opt / 3-opt solution constructor
                    let (block_clauses, mut _active_cycles) = if !cycle_repair || opt == 0 {
                        (get_blocking_clauses(&sol_cycles, encoder, &g, block_method, balanced), sol_cycles.clone())
                    } else if opt >= 1 {
                        let (clauses, cycles) = two_opt(&sol_cycles, encoder, &g, contractor, hub_registry, block_method, balanced, opt, three_opt);
                        if cycles.len() == 1 && cycles[0].len() == g.adjacency_list.len() {
                            let flat: Vec<i32> = cycles.into_iter().flatten().collect();
                            let full_cycle = contractor.uncontract_cycle(&flat);
                            let line = full_cycle.iter().map(|i| i.to_string()).collect::<Vec<String>>().join(" ");
                            let now = instant.elapsed();
                            let time = now - previous_time;
                            let add_block_clauses_time = now - previous_time - sat_solving_time;
                            println!("number of added block clauses = {}", clause_count);
                            println!("add block clauses time = {:?}", add_block_clauses_time);
                            println!("increment time = {:?}", time);
                            println!();
                            println!("hamiltonian cycle found by 2-opt/3-opt");
                            println!("solution: ");
                            println!("{}\n", line);
                            println!("s SATISFIABLE");
                            return (count, clause_count, Some(full_cycle));
                        }
                        (clauses, cycles)
                    } else {
                        panic!("2-opt option \n-t 0:2-opt off\n-t 1,2,3:2-opt on");
                    };

                    // Attempt exact SAT Macro-Patcher spanning tree / components when <= 60 active cycles remain
                    if cycle_repair && _active_cycles.len() >= 2 && _active_cycles.len() <= 60 {
                        let protected_edges: HashSet<(i32, i32)> = contractor.chain_map.keys().copied().collect();
                        let patched = SatMacroPatcher::try_patch_components(&_active_cycles, &g, &protected_edges);
                        if patched.len() == 1 && patched[0].len() == g.adjacency_list.len() {
                            println!("SatMacroPatcher: successfully synthesized full Hamiltonian tour from {} active cycles", _active_cycles.len());
                            let full_cycle = contractor.uncontract_cycle(&patched[0]);
                            let line = full_cycle.iter().map(|i| i.to_string()).collect::<Vec<String>>().join(" ");
                            let now = instant.elapsed();
                            let time = now - previous_time;
                            let add_block_clauses_time = now - previous_time - sat_solving_time;
                            println!("number of added block clauses = {}", clause_count);
                            println!("add block clauses time = {:?}", add_block_clauses_time);
                            println!("increment time = {:?}", time);
                            println!();
                            println!("solution: ");
                            println!("{}\n", line);
                            println!("s SATISFIABLE");
                            return (count, clause_count, Some(full_cycle));
                        } else if patched.len() < _active_cycles.len() {
                            _active_cycles = patched;
                        }
                    }

                    // Attempt MacroCrossoverSplicer when 2 <= _active_cycles.len() <= 6 (parity-breaking k-opt for macro-cycles)
                    if cycle_repair && _active_cycles.len() >= 2 && _active_cycles.len() <= 6 {
                        if let Some(spliced_tour) = MacroCrossoverSplicer::try_crossover_splice(&_active_cycles, &g, contractor) {
                            println!("MacroCrossoverSplicer: successfully spliced {} macro-cycles into single tour via auxiliary SAT crossover", _active_cycles.len());
                            let full_cycle = contractor.uncontract_cycle(&spliced_tour);
                            let line = full_cycle.iter().map(|i| i.to_string()).collect::<Vec<String>>().join(" ");
                            let now = instant.elapsed();
                            let time = now - previous_time;
                            let add_block_clauses_time = now - previous_time - sat_solving_time;
                            println!("number of added block clauses = {}", clause_count);
                            println!("add block clauses time = {:?}", add_block_clauses_time);
                            println!("increment time = {:?}", time);
                            println!();
                            println!("solution: ");
                            println!("{}\n", line);
                            println!("s SATISFIABLE");
                            return (count, clause_count, Some(full_cycle));
                        }
                    }

                    // Attempt GadgetPathAbsorber satellite subcycle absorption
                    if cycle_repair && _active_cycles.len() >= 2 {
                        let protected_edges: HashSet<(i32, i32)> = contractor.chain_map.keys().copied().collect();
                        let absorbed = GadgetPathAbsorber::try_absorb_gadgets(&_active_cycles, &g, &protected_edges);
                        if absorbed.len() == 1 && absorbed[0].len() == g.adjacency_list.len() {
                            println!("GadgetPathAbsorber: successfully absorbed satellite gadgets into full Hamiltonian tour");
                            let full_cycle = contractor.uncontract_cycle(&absorbed[0]);
                            let line = full_cycle.iter().map(|i| i.to_string()).collect::<Vec<String>>().join(" ");
                            let now = instant.elapsed();
                            let time = now - previous_time;
                            let add_block_clauses_time = now - previous_time - sat_solving_time;
                            println!("number of added block clauses = {}", clause_count);
                            println!("add block clauses time = {:?}", add_block_clauses_time);
                            println!("increment time = {:?}", time);
                            println!();
                            println!("solution: ");
                            println!("{}\n", line);
                            println!("s SATISFIABLE");
                            return (count, clause_count, Some(full_cycle));
                        } else if absorbed.len() < _active_cycles.len() {
                            println!("GadgetPathAbsorber: absorbed satellite subcycles from {} down to {} cycles", _active_cycles.len(), absorbed.len());
                            _active_cycles = absorbed;
                        }
                    }

                    // Attempt Localized SAT Repair (LNS via backbone freezing) to eliminate last-mile thrashing
                    if cycle_repair && _active_cycles.len() >= 2 {
                        if let Some(repaired_tour) = LocalizedSatRepair::try_repair(&_active_cycles, &g, &working_cnf, encoder, 5) {
                            println!("LocalizedSatRepair: successfully synthesized full Hamiltonian tour via Localized SAT LNS!");
                            let full_cycle = contractor.uncontract_cycle(&repaired_tour);
                            let line = full_cycle.iter().map(|i| i.to_string()).collect::<Vec<String>>().join(" ");
                            let now = instant.elapsed();
                            let time = now - previous_time;
                            let add_block_clauses_time = now - previous_time - sat_solving_time;
                            println!("number of added block clauses = {}", clause_count);
                            println!("add block clauses time = {:?}", add_block_clauses_time);
                            println!("increment time = {:?}", time);
                            println!();
                            println!("solution: ");
                            println!("{}\n", line);
                            println!("s SATISFIABLE");
                            return (count, clause_count, Some(full_cycle));
                        }
                    }

                    let mut round_cuts = Cnf::new();
                    // Gadget Interface Parity & Direct Splicing Check
                    if _active_cycles.len() >= 2 {
                        let total_nodes = g.adjacency_list.len();
                        let macro_indices: Vec<usize> = _active_cycles
                            .iter()
                            .enumerate()
                            .filter(|(_, c)| c.len() >= total_nodes / 10)
                            .map(|(idx, _)| idx)
                            .collect();

                        let mut analyzed_subcycles = HashSet::new();
                        for &m_idx in &macro_indices {
                            let mut macro_cycle = _active_cycles[m_idx].clone();
                            for (c_idx, subcycle) in _active_cycles.iter().enumerate() {
                                if c_idx != m_idx && subcycle.len() <= 16 {
                                    let repair_target = if cycle_repair {
                                        Some(macro_cycle.as_slice())
                                    } else {
                                        None
                                    };
                                    let gadget_res = GadgetInterfaceParityEngine::analyze_subcycle_gadget(
                                        subcycle,
                                        &g,
                                        repair_target,
                                        encoder,
                                    );

                                    // 1. Direct splice check
                                    if cycle_repair {
                                        if let Some(spliced) = gadget_res.direct_spliced_tour {
                                            macro_cycle = spliced;
                                            if macro_cycle.len() == total_nodes {
                                                println!("GadgetInterfaceParity: direct spliced full tour found ({} vertices)", macro_cycle.len());
                                                let full_cycle = contractor.uncontract_cycle(&macro_cycle);
                                                let line = full_cycle.iter().map(|i| i.to_string()).collect::<Vec<String>>().join(" ");
                                                println!();
                                                println!("solution: ");
                                                println!("{}\n", line);
                                                println!("s SATISFIABLE");
                                                return (count, clause_count, Some(full_cycle));
                                            }
                                        }
                                    }

                                    // 2. Infeasible port pruning clauses & boundary cut parity clauses (deduplicated per subcycle)
                                    if analyzed_subcycles.insert(c_idx) {
                                        for cl in gadget_res.pruning_clauses {
                                            clause_count += 1;
                                            round_cuts.add_clause(cl);
                                        }
                                        for cl in gadget_res.cut_parity_clauses {
                                            clause_count += 1;
                                            round_cuts.add_clause(cl);
                                        }
                                    }
                                }
                            }
                        }
                    }

                    let mut cnf = Cnf::new();
                    cnf.extend(block_clauses);
                    count += 1;

                    previous_cnf = if output_folder != "default" {
                        let mut write_cnf = previous_cnf;
                        write_cnf.extend(cnf.clone());
                        let output_file = format!("{}/increment{}.cnf", output_folder, count);
                        let _ = file_operations::write_dimacs(write_cnf.clone(), &output_file);
                        write_cnf
                    } else {
                        Cnf::new()
                    };

                    if cnf_normalize == 1 {
                        let normalized_cnf = cnf.normalize();
                        clause_count += normalized_cnf.len() as i32;
                        round_cuts.extend(normalized_cnf);
                    } else {
                        clause_count += cnf.len() as i32;
                        round_cuts.extend(cnf);
                    }

                    let total_v = g.adjacency_list.len();
                    let comprehensive_sec = EmpiricalBackboneCutter::generate_comprehensive_sec_clauses(
                        &sol_cycles,
                        total_v / 2,
                        &encoder.graph_lit_map,
                    );
                    for cl in comprehensive_sec {
                        clause_count += 1;
                        round_cuts.add_clause(Clause::from_iter(cl));
                    }

                    // Inject boundary cut clauses for shortest non-giant subcycles (budgeted to top 30)
                    let mut non_giant_cycles: Vec<&Vec<i32>> = sol_cycles
                        .iter()
                        .filter(|c| c.len() >= 3 && c.len() < total_v / 2)
                        .collect();
                    non_giant_cycles.sort_by_key(|c| c.len());
                    non_giant_cycles.truncate(30);

                    for cycle in non_giant_cycles {
                        let b_clauses = get_boundary_cut_clauses(cycle, encoder, &g, total_v, 0);
                        for cl in b_clauses {
                            clause_count += 1;
                            round_cuts.add_clause(cl);
                        }
                    }

                    // Injected cuts on reduced active cycles when <= 50 subcycles remain
                    if _active_cycles.len() <= 50 && _active_cycles != sol_cycles {
                        let mut non_giant_active: Vec<&Vec<i32>> = _active_cycles
                            .iter()
                            .filter(|c| c.len() >= 3 && c.len() < total_v / 2)
                            .collect();
                        non_giant_active.sort_by_key(|c| c.len());
                        non_giant_active.truncate(20);

                        for cycle in non_giant_active {
                            let b_clauses = get_boundary_cut_clauses(cycle, encoder, &g, total_v, 0);
                            for cl in b_clauses {
                                clause_count += 1;
                                round_cuts.add_clause(cl);
                            }
                        }
                    }

                    // Inject bicomponent cut clauses if top two largest cycles both have length >= total_v / 5
                    if sol_cycles.len() >= 2 {
                        let mut sorted_cycles: Vec<&Vec<i32>> = sol_cycles.iter().collect();
                        sorted_cycles.sort_by(|a, b| b.len().cmp(&a.len()));
                        let threshold = total_v / 5;
                        if sorted_cycles[0].len() >= threshold && sorted_cycles[1].len() >= threshold {
                            let bicomponent_cuts = TwinGiantSplicer::generate_bicomponent_cut_clauses(
                                sorted_cycles[0],
                                sorted_cycles[1],
                                &g,
                                &encoder.graph_lit_map,
                            );
                            for cl in bicomponent_cuts {
                                clause_count += 1;
                                round_cuts.add_clause(cl);
                            }
                        }
                    }

                    // QuotientBlockCutter disabled to prevent unsound cut over-constraining
                    /*
                    if !modular_blocks.is_empty() && sol_cycles.len() >= 2 {
                        let q_cuts = QuotientBlockCutter::generate_quotient_sec_clauses(&sol_cycles, &modular_blocks, &g, encoder);
                        if !q_cuts.is_empty() {
                            println!("QuotientBlockCutter: injected {} modular quotient cut clauses", q_cuts.len());
                            for cl in q_cuts {
                                clause_count += 1;
                                round_cuts.add_clause(cl);
                            }
                        }
                    }
                    */

                    // BalancedPairCutset disabled to maintain fast increment velocity (avoiding deep CNF overhead)
                    /*
                    if sol_cycles.len() >= 3 {
                        let balanced_cuts = BalancedPairCutset::generate_balanced_cutset_clauses(&sol_cycles, &g, &encoder.graph_lit_map);
                        for cl in balanced_cuts {
                            clause_count += 1;
                            round_cuts.add_clause(cl);
                        }
                    }
                    */

                    if !round_cuts.is_empty() {
                        working_cnf.extend(round_cuts.clone());
                        accumulated_cut_cnfs.push(round_cuts.clone());
                        incremental_solver.add_cuts(&round_cuts);
                    }
                    let max_cycle_len = _active_cycles.iter().map(|c| c.len()).max().unwrap_or(0);
                    if bounded_freezer != 0 && _active_cycles.len() > 1 && total_v >= 100 && (max_cycle_len >= total_v / 2 || _active_cycles.len() <= 25) {
                        assumptions = BackboneFreezer::select_topologically_bounded_assumptions(
                            &_active_cycles,
                            &g,
                            encoder,
                            sat_solving_time.as_secs_f64(),
                        );
                        if !assumptions.is_empty() {
                            println!("BackboneFreezer (H2-Refined): locked {} bounded backbone edges (giant cycle len {})", assumptions.len(), max_cycle_len);
                        }
                    } else {
                        assumptions.clear();
                        phase_hints.clear();
                    }

                    let time = now - previous_time;
                    let add_block_clauses_time = now - previous_time - sat_solving_time;
                    previous_time = now;
                    println!("number of added block clauses = {}", clause_count);
                    println!("add block clauses time = {:?}", add_block_clauses_time);
                    println!("increment time = {:?}", time);

                    let should_reseed = if total_v >= 500 {
                        accumulated_cut_cnfs.len() >= 50 || sat_solving_time.as_secs_f64() >= 300.0
                    } else {
                        SolverReseeder::should_reseed(sat_solving_time.as_secs_f64(), count as usize, &reseeder_opts)
                            || accumulated_cut_cnfs.len() >= 10
                    };

                    if should_reseed {
                        let pruned_cnf = CnfSubsumer::prune_and_subsume_cuts(&accumulated_cut_cnfs);
                        println!("SolverReseeder: compressed {} cut sets down to {} non-redundant clauses (round {}, last SAT time {:.2}s)",
                            accumulated_cut_cnfs.len(), pruned_cnf.len(), count, sat_solving_time.as_secs_f64());
                        working_cnf = base_cnf.clone();
                        working_cnf.extend(pruned_cnf.clone());
                        accumulated_cut_cnfs = vec![pruned_cnf];
                        incremental_solver.reset_with_cnf(&working_cnf);
                    }
                }
            }
            PortfolioResult::Unsat => {
                if !assumptions.is_empty() {
                    println!("Assumptions UNSAT (speculative backbone invalid) -> clearing {} assumptions and resuming base search", assumptions.len());
                    assumptions.clear();
                    continue;
                }
                println!("s UNSATISFIABLE");
                return (count, clause_count, None);
            }
            PortfolioResult::Interrupted => {
                if instant.elapsed().as_secs_f64() >= timeout_secs {
                    println!("\ns UNKNOWN (TIMEOUT: {:.2}s reached >= {:.2}s limit)", instant.elapsed().as_secs_f64(), timeout_secs);
                    return (count, clause_count, None);
                }
                if !assumptions.is_empty() {
                    assumptions.clear();
                    continue;
                }
                return (count, clause_count, None);
            }
        }
    }
}

pub fn get_solution_arcs_from_lits<'a, M>(lits: &[Lit], graph_lit_map: M) -> Vec<(i32, i32)>
where
    M: IntoIterator<Item = (&'a (i32, i32), &'a Lit)>,
{
    let lit_set: HashSet<Lit> = lits.iter().copied().collect();
    let mut arcs = Vec::new();
    for (&arc, &lit) in graph_lit_map {
        if lit_set.contains(&lit) {
            arcs.push(arc);
        }
    }
    arcs
}

fn get_solution_arcs(sol:Assignment,lit_map:&BTreeMap<(i32,i32),Lit>) -> Vec<(i32,i32)>{
    let sol_arcs: Vec<(i32,i32)> = lit_map.iter().filter_map(|((u,v), lit)| if sol[lit.var()] == TernaryVal::True { Some((*u,*v)) } else { None }).collect();
    sol_arcs
}

fn get_solution_cycles(sol_arcs: Vec<(i32, i32)>) -> Vec<Vec<i32>> {
    let mut arcs: BTreeMap<i32,i32> = std::collections::BTreeMap::new();
    let mut cycles = Vec::new();
    let mut visited = std::collections::BTreeSet::new();

    for arc in sol_arcs{
        arcs.insert(arc.0,arc.1);
    }
    
    for node in arcs.keys() {
        if visited.contains(node) {
            continue;
        }
        let mut cycle = Vec::new();
        let mut current_node = node;
        loop{
            visited.insert(current_node);
            cycle.push(*current_node);
            current_node = match arcs.get(current_node) {
                Some(node) => node,
                None => break,
            };
            if visited.contains(current_node) {
                break;
            }
        }
        cycles.push(cycle);
    }

    cycles
}


fn cycle_hub_score(cycle: &[i32], hub_registry: &HubRegistry, g: &Graph) -> (usize, usize) {
    if hub_registry.hub_vertices.is_empty() {
        return (0, 0);
    }
    let mut contains_hub_count = 0;
    let mut incident_hub_edges = 0;

    for &v in cycle {
        if hub_registry.is_hub_vertex(v) {
            contains_hub_count += 1;
        }
        if let Some(adjs) = g.adjacency_list.get(&v) {
            for &u in adjs {
                if hub_registry.is_hub_vertex(u) {
                    incident_hub_edges += 1;
                }
            }
        }
    }

    let tier = if contains_hub_count > 0 {
        2
    } else if incident_hub_edges > 0 {
        1
    } else {
        0
    };
    (tier, contains_hub_count * 1000 + incident_hub_edges)
}

// 2-opt and Candidate 3-opt Solution Constructor
// Attempts to merge subcycles into a single Hamiltonian Cycle.
// If complete merger fails, returns standard blocking clauses for active subcycles.
fn two_opt(
    sol_cycles: &Vec<Vec<i32>>,
    encoder: &mut Encoder,
    g: &Graph,
    contractor: &Degree2Contractor,
    hub_registry: &HubRegistry,
    block_method: i32,
    balanced: i32,
    opt: i32,
    three_opt: i32,
) -> (Vec<Clause>, Vec<Vec<i32>>) {
    let mut cycles = sol_cycles.to_vec();
    let mut merged = true;
    let mut cache_vertex: HashSet<usize> = HashSet::new();
    let mut active_cycles_number: Vec<usize> = (0..cycles.len()).collect();

    if !hub_registry.hub_vertices.is_empty() {
        active_cycles_number.sort_by(|&a, &b| {
            cycle_hub_score(&cycles[b], hub_registry, g).cmp(&cycle_hub_score(&cycles[a], hub_registry, g))
        });
    }

    while merged {
        let (_new_block_clauses, new_merged, merged_numbers, new_cycle) =
            merge_cycles(&cycles, encoder, g, contractor, hub_registry, block_method, balanced, &mut cache_vertex, &active_cycles_number, opt);
        merged = new_merged;

        if merged {
            cycles.push(new_cycle);
            let mut remove_indices = [merged_numbers.0, merged_numbers.1];
            remove_indices.sort_unstable_by(|a, b| b.cmp(a));
            for &idx in &remove_indices {
                active_cycles_number.remove(idx);
            }
            active_cycles_number.push(cycles.len() - 1);
            if !hub_registry.hub_vertices.is_empty() {
                active_cycles_number.sort_by(|&a, &b| {
                    cycle_hub_score(&cycles[b], hub_registry, g).cmp(&cycle_hub_score(&cycles[a], hub_registry, g))
                });
            }
        }

        // Try candidate 3-opt merge when 2-opt cannot merge further
        if !merged && three_opt == 1 && active_cycles_number.len() >= 3 {
            let (_three_block_clauses, three_merged, three_indices, three_cycle) =
                merge_three_cycles(&cycles, encoder, g, contractor, hub_registry, block_method, balanced, &active_cycles_number);
            if three_merged {
                cycles.push(three_cycle);
                let (ia, ib, ic) = three_indices;
                let mut remove_indices = [ia, ib, ic];
                remove_indices.sort_unstable_by(|a, b| b.cmp(a));
                for &idx in &remove_indices {
                    active_cycles_number.remove(idx);
                }
                active_cycles_number.push(cycles.len() - 1);
                merged = true;
                cache_vertex.clear();
                if !hub_registry.hub_vertices.is_empty() {
                    active_cycles_number.sort_by(|&a, &b| {
                        cycle_hub_score(&cycles[b], hub_registry, g).cmp(&cycle_hub_score(&cycles[a], hub_registry, g))
                    });
                }
                continue;
            }
        }

        if active_cycles_number.len() == 1 || !merged {
            break;
        }
    }

    let mut active_cycles = get_active_cycles(&cycles, &active_cycles_number);
    if active_cycles.len() > 1 {
        let absorbed = CycleChainAbsorber::absorb_all(&active_cycles, g, contractor, hub_registry);
        if absorbed.len() < active_cycles.len() {
            println!("CycleChainAbsorber: merged from {} to {} subcycles (giant cycle len {})", active_cycles.len(), absorbed.len(), absorbed[0].len());
            active_cycles = absorbed;
        } else {
            let old_absorbed = SubcycleAbsorber::absorb_subcycles(&active_cycles, g, contractor, hub_registry);
            if old_absorbed.len() < active_cycles.len() {
                println!("SubcycleAbsorber: merged from {} to {} subcycles (giant cycle len {})", active_cycles.len(), old_absorbed.len(), old_absorbed[0].len());
                active_cycles = old_absorbed;
            }
        }
    }

    let block_clauses = if active_cycles.len() == 1 && active_cycles[0].len() == g.adjacency_list.len() {
        Vec::new()
    } else {
        match opt {
            3 => get_blocking_clauses(&active_cycles, encoder, g, block_method, balanced),
            2 => {
                let mut cl = get_blocking_clauses(sol_cycles, encoder, g, block_method, balanced);
                cl.extend(get_blocking_clauses(&active_cycles, encoder, g, block_method, balanced));
                cl
            }
            _ => get_blocking_clauses(sol_cycles, encoder, g, block_method, balanced),
        }
    };

    println!("number of connected cycles = {}", cycles.len() - sol_cycles.len());
    println!("number of merged cycles = {}", active_cycles.len());
    println!("merged cycle lengths map (length:number) = {:?}", map_cycle_lengths(&active_cycles));

    (block_clauses, active_cycles)
}

fn merge_cycles(
    cycles: &Vec<Vec<i32>>,
    encoder: &mut Encoder,
    g: &Graph,
    contractor: &Degree2Contractor,
    hub_registry: &HubRegistry,
    block_method: i32,
    balanced: i32,
    cache_vertex: &mut HashSet<usize>,
    active_cycles_number: &Vec<usize>,
    opt: i32,
) -> (Vec<Clause>, bool, (usize, usize), Vec<i32>) {
    //(block_clauses,merged,(merged_number1,merged_number2),new_cycle)
    
    for i in 0..active_cycles_number.len(){
        let left = active_cycles_number[i];
        if !cache_vertex.contains(&left){
            for j in i+1..active_cycles_number.len(){
                let right = active_cycles_number[j];

                match swap_node(&cycles[left],&cycles[right],&g, contractor, hub_registry){
                    Some(new_cycle) =>{
                    let new_block_clauses = get_blocking_clauses(&vec!(new_cycle.clone()), encoder, g, block_method, balanced);
                    return (new_block_clauses,true,(i,j),new_cycle)
                    }
                    None =>{
                        continue
                    }
                }
                
            }
            cache_vertex.insert(left);
        }
        if opt == 4 || opt == 5{
            return (vec!(),false,(0,0),vec!())
        }
    }
    
    (vec!(),false,(0,0),vec!())
}


fn swap_node(
    cycle1: &Vec<i32>,
    cycle2: &Vec<i32>,
    g: &Graph,
    contractor: &Degree2Contractor,
    _hub_registry: &HubRegistry,
) -> Option<Vec<i32>> {
    let n2 = cycle2.len();
    if cycle1.len() < 3 || n2 < 3 {
        return None;
    }

    let mut cycle2_pos: HashMap<i32, usize> = HashMap::with_capacity(n2);
    for (idx, &v) in cycle2.iter().enumerate() {
        cycle2_pos.insert(v, idx);
    }

    let is_protected = |u: i32, v: i32| -> bool {
        contractor.chain_map.contains_key(&(u, v)) || contractor.chain_map.contains_key(&(v, u))
    };

    let has_edge = |u: i32, v: i32| -> bool {
        g.adjacency_list.get(&u).map_or(false, |nbrs| nbrs.contains(&v))
    };

    for i in 0..cycle1.len() {
        let u1 = cycle1[i];
        let v1 = cycle1[(i + 1) % cycle1.len()];
        if is_protected(u1, v1) {
            continue;
        }

        let adjs_of_left_head = match g.adjacency_list.get(&u1) {
            Some(a) => a,
            None => continue,
        };

        for &u2 in adjs_of_left_head {
            if let Some(&j) = cycle2_pos.get(&u2) {
                let v2_fwd = cycle2[(j + 1) % n2];
                let v2_rev = cycle2[(j + n2 - 1) % n2];

                if has_edge(v1, v2_fwd) && !is_protected(u2, v2_fwd) {
                    return cycle_join(cycle1, cycle2, i, j, true);
                }

                if has_edge(v1, v2_rev) && !is_protected(u2, v2_rev) {
                    return cycle_join(cycle1, cycle2, i, j, false);
                }
            }
        }
    }
    None
}

fn cycle_join(cycle1:&Vec<i32>,cycle2:&Vec<i32>,i:usize,j:usize,reverse:bool) -> Option<Vec<i32>>{
    let mut new_cycle = Vec::new();

    if reverse{
        // cycle1のindex iまでを追加
        new_cycle.extend(&cycle1[0..=i]);

        // cycle2のindex jから逆順にindex 0までの要素を追加
        new_cycle.extend(cycle2[0..=j].iter().rev());
        if j != cycle2.len()-1{
        // cycle2のindexの最後から、j+1までの要素を逆順に追加
            new_cycle.extend(cycle2[j+1..].iter().rev());
        }
        
        if i != cycle1.len()-1{
        // cycle1のindex i+1から最後までをnew_cycleに追加
            new_cycle.extend(&cycle1[i+1..]);
        }
    }else{
        new_cycle.extend(&cycle1[0..=i]);
        new_cycle.extend(&cycle2[j..]);
        if j != 0{
            new_cycle.extend(&cycle2[0..=j-1]);
        }
        if i != cycle1.len()-1{
            new_cycle.extend(&cycle1[i+1..]);
        }
    }

    Some(new_cycle)
}

/// Try to merge three directed cycles by a 3-edge swap.
/// Config A (0): C1 -> C2 -> C3 -> C1  (u1->v2, u2->v3, u3->v1)
/// Try to merge three directed cycles by a 3-edge swap.
/// Config A (0): C1 -> C2 -> C3 -> C1  (u1->v2, u2->v3, u3->v1)
/// Config B (1): C1 -> C3 -> C2 -> C1  (u1->v3, u3->v2, u2->v1)
fn swap_three_nodes(
    c1: &Vec<i32>,
    c2: &Vec<i32>,
    c3: &Vec<i32>,
    g: &Graph,
    contractor: &Degree2Contractor,
    _hub_registry: &HubRegistry,
) -> Option<Vec<i32>> {
    let n2 = c2.len();
    let n3 = c3.len();
    if c1.len() < 3 || n2 < 3 || n3 < 3 {
        return None;
    }

    let mut c2_pos: HashMap<i32, usize> = HashMap::with_capacity(n2);
    for (idx, &v) in c2.iter().enumerate() {
        c2_pos.insert(v, idx);
    }

    let mut c3_pos: HashMap<i32, usize> = HashMap::with_capacity(n3);
    for (idx, &v) in c3.iter().enumerate() {
        c3_pos.insert(v, idx);
    }

    let is_protected = |u: i32, v: i32| -> bool {
        contractor.chain_map.contains_key(&(u, v)) || contractor.chain_map.contains_key(&(v, u))
    };

    let has_edge = |u: i32, v: i32| -> bool {
        g.adjacency_list.get(&u).map_or(false, |nbrs| nbrs.contains(&v))
    };

    for i in 0..c1.len() {
        let u1 = c1[i];
        let v1 = c1[(i + 1) % c1.len()];
        if is_protected(u1, v1) {
            continue;
        }

        let adjs_u1 = match g.adjacency_list.get(&u1) {
            Some(a) => a,
            None => continue,
        };

        // Check Config A: u1 -> v2, u2 -> v3, u3 -> v1
        for &v2 in adjs_u1 {
            if let Some(&v2_idx) = c2_pos.get(&v2) {
                let j = (v2_idx + n2 - 1) % n2;
                let u2 = c2[j];
                if is_protected(u2, v2) {
                    continue;
                }

                if let Some(adjs_u2) = g.adjacency_list.get(&u2) {
                    for &v3 in adjs_u2 {
                        if let Some(&v3_idx) = c3_pos.get(&v3) {
                            let k = (v3_idx + n3 - 1) % n3;
                            let u3 = c3[k];
                            if is_protected(u3, v3) {
                                continue;
                            }

                            if has_edge(u3, v1) {
                                return cycle_join_three(c1, c2, c3, i, j, k, 0);
                            }
                        }
                    }
                }
            }
        }

        // Check Config B: u1 -> v3, u3 -> v2, u2 -> v1
        for &v3 in adjs_u1 {
            if let Some(&v3_idx) = c3_pos.get(&v3) {
                let k = (v3_idx + n3 - 1) % n3;
                let u3 = c3[k];
                if is_protected(u3, v3) {
                    continue;
                }

                if let Some(adjs_u3) = g.adjacency_list.get(&u3) {
                    for &v2 in adjs_u3 {
                        if let Some(&v2_idx) = c2_pos.get(&v2) {
                            let j = (v2_idx + n2 - 1) % n2;
                            let u2 = c2[j];
                            if is_protected(u2, v2) {
                                continue;
                            }

                            if has_edge(u2, v1) {
                                return cycle_join_three(c1, c2, c3, i, j, k, 1);
                            }
                        }
                    }
                }
            }
        }
    }

    None
}

/// Reconstruct a single merged cycle from three directed cycles given cut positions.
fn cycle_join_three(
    c1: &Vec<i32>, c2: &Vec<i32>, c3: &Vec<i32>,
    i: usize, j: usize, k: usize,
    config: u8,
) -> Option<Vec<i32>> {
    let mut new_cycle = Vec::new();
    if config == 0 {
        new_cycle.extend(&c1[0..=i]);
        new_cycle.extend(&c2[j+1..]);
        new_cycle.extend(&c2[..=j]);
        new_cycle.extend(&c3[k+1..]);
        new_cycle.extend(&c3[..=k]);
        if i + 1 < c1.len() {
            new_cycle.extend(&c1[i+1..]);
        }
    } else {
        new_cycle.extend(&c1[0..=i]);
        new_cycle.extend(&c3[k+1..]);
        new_cycle.extend(&c3[..=k]);
        new_cycle.extend(&c2[j+1..]);
        new_cycle.extend(&c2[..=j]);
        if i + 1 < c1.len() {
            new_cycle.extend(&c1[i+1..]);
        }
    }
    Some(new_cycle)
}

/// Try to merge a triplet of active cycles using candidate graph filtering.
fn merge_three_cycles(
    cycles: &Vec<Vec<i32>>,
    encoder: &mut Encoder,
    g: &Graph,
    contractor: &Degree2Contractor,
    hub_registry: &HubRegistry,
    block_method: i32,
    balanced: i32,
    active_cycles_number: &Vec<usize>,
) -> (Vec<Clause>, bool, (usize, usize, usize), Vec<i32>) {
    let n = active_cycles_number.len();

    let mut vertex_to_active: HashMap<i32, usize> = HashMap::new();
    for (active_idx, &cycle_idx) in active_cycles_number.iter().enumerate() {
        for &v in &cycles[cycle_idx] {
            vertex_to_active.insert(v, active_idx);
        }
    }
    let mut cycle_neighbors: Vec<HashSet<usize>> = vec![HashSet::new(); n];
    for (active_idx, &cycle_idx) in active_cycles_number.iter().enumerate() {
        for &u in &cycles[cycle_idx] {
            if let Some(adjs) = g.adjacency_list.get(&u) {
                for &v in adjs {
                    if let Some(&neighbor_active) = vertex_to_active.get(&v) {
                        if neighbor_active != active_idx {
                            cycle_neighbors[active_idx].insert(neighbor_active);
                        }
                    }
                }
            }
        }
    }

    for a in 0..n {
        let neighbors_a: Vec<usize> = cycle_neighbors[a].iter()
            .filter(|&&b| b > a)
            .cloned()
            .collect();
        for &b in &neighbors_a {
            for &c in &cycle_neighbors[b] {
                if c <= b { continue; }
                if !cycle_neighbors[a].contains(&c) { continue; }
                let ci = active_cycles_number[a];
                let cj = active_cycles_number[b];
                let ck = active_cycles_number[c];
                if let Some(new_cycle) = swap_three_nodes(&cycles[ci], &cycles[cj], &cycles[ck], g, contractor, hub_registry) {
                    let new_block_clauses = get_blocking_clauses(&vec!(new_cycle.clone()), encoder, g, block_method, balanced);
                    return (new_block_clauses, true, (a, b, c), new_cycle);
                }
            }
        }
    }
    (vec![], false, (0, 0, 0), vec![])
}

/// Inject Partial MTZ ordering constraints for a set of vertices K.
/// Uses ORDER encoding: for each v in K, creates n-1 boolean variables
/// o[v][t] where o[v][t] = 1 iff position(v) >= t+1.
///
/// Prevents ANY subtour entirely within K\{source}.
fn inject_partial_mtz(
    k_vertices: &Vec<i32>,
    g: &Graph,
    encoder: &mut Encoder,
    _n: usize,
) -> Vec<Clause> {
    
    let mut clauses: Vec<Clause> = Vec::new();
    let source = k_vertices[0];  // Pick first vertex as source
    let k_set: HashSet<i32> = k_vertices.iter().cloned().collect();
    
    let mut order_vars: HashMap<i32, Vec<Lit>> = HashMap::new();
    let k_len = k_vertices.len();
    for &v in k_vertices.iter() {
        let mut vars = Vec::new();
        // Allocate order variables corresponding to positions within K (1..k_len-1)
        for _t in 0..k_len-1 {
            let lit = encoder.instance.new_lit();
            vars.push(lit);
        }
        order_vars.insert(v, vars);
    }
    
    // 1. Monotonicity: o[v][t+1] → o[v][t]
    //    Clause: ¬o[v][t+1] ∨ o[v][t]
    for &v in k_vertices.iter() {
        let vars = &order_vars[&v];
        for t in 0..vars.len()-1 {
            clauses.push(clause!(!vars[t+1], vars[t]));
        }
    }
    
    // 3. Source constraint: position(source) = 0 → o[source][t] = 0 for all t
    for lit in order_vars[&source].iter() {
        clauses.push(clause!(!*lit));
    }
    
    // 4. MTZ constraints: for each directed edge (u,v) where u,v ∈ K\{source}
    //    x_{u,v} = 1 → position(v) >= position(u) + 1
    //    Clause: ¬x_{u,v} ∨ ¬o[u][t] ∨ o[v][t+1]  for t = 0..n-3
    //    Also:   ¬x_{u,v} ∨ o[v][0]  (position(v) >= 1 when arc is selected)
    for &u in k_vertices.iter() {
        if let Some(adjs) = g.adjacency_list.get(&u) {
            for &v in adjs.iter() {
                if !k_set.contains(&v) { continue; }  // Only edges within K
                if v == source && u == source { continue; } // Skip self-loops
                
                let x_uv = match encoder.graph_lit_map.get(&(u, v)) {
                    Some(lit) => *lit,
                    None => continue,  // Edge doesn't exist in encoding
                };
                
                if u == source {
                    // Edge from source: x_{s,v} → position(v) >= 1
                    let o_v = &order_vars[&v];
                    clauses.push(clause!(!x_uv, o_v[0]));
                } else if v == source {
                    // Edge TO source: no MTZ constraint needed
                    // (source can be at position 0 and receive edge from position n-1)
                    continue;
                } else {
                    // Edge between non-source K vertices
                    let o_u = &order_vars[&u];
                    let o_v = &order_vars[&v];
                    
                    // x_{u,v} → position(v) >= 1
                    clauses.push(clause!(!x_uv, o_v[0]));
                    
                    // x_{u,v} ∧ o[u][t] → o[v][t+1] for t = 0..n-3
                    for t in 0..o_u.len()-1 {
                        clauses.push(clause!(!x_uv, !o_u[t], o_v[t+1]));
                    }
                }
            }
        }
    }
    
    clauses
}

fn get_blocking_clauses(
    cycles: &Vec<Vec<i32>>,
    encoder: &Encoder,
    g: &Graph,
    _block_method: i32,
    _balanced: i32,
) -> Vec<Clause> {
    let mut options = CutSelectorOptions::default();
    if cycles.len() <= 20 {
        options.base_max_cuts = 50;
        options.max_cycle_len_threshold = usize::MAX;
        options.tiny_cycle_boundary_len = usize::MAX;
        options.enable_reverse_blocking = true;
        options.enable_incoming_boundary_cuts = true;
    } else if cycles.len() <= 50 {
        let total_v = g.adjacency_list.len();
        options.base_max_cuts = 40;
        options.max_cycle_len_threshold = if total_v > 0 { total_v / 2 } else { usize::MAX };
        options.tiny_cycle_boundary_len = 32;
    }
    let (mut clauses, selected) = CutSelector::select_and_generate_cuts(cycles, g, encoder, &options);
    if !selected.is_empty() {
        println!("CutSelector: selected {}/{} subcycles (generated {} budgeted clauses)", selected.len(), cycles.len(), clauses.len());
    }

    if cycles.len() >= 2 && cycles.len() <= 16 {
        let hemi_cuts = HemisphereSplicer::generate_hemisphere_crossing_cuts(cycles, g, encoder);
        if !hemi_cuts.is_empty() {
            println!("HemisphereSplicer: generated {} bi-partition crossing cut clauses", hemi_cuts.len());
            clauses.extend(hemi_cuts);
        }
    }

    clauses
}

#[allow(dead_code)]
fn cegar_blocking_clauses(cycle:&Vec<i32>,lit_map:&BTreeMap<(i32,i32),Lit>)-> Vec<Clause>{
    let mut clauses =  Vec::new();
    // for cycle in sol_cycles.iter() {
    let len = cycle.len();
    // 順方向
    let mut clause = rustsat::types::Clause::new();
    for i in 0..len {
        let lit = lit_map.get(&(cycle[i], cycle[(i+1)%len])).unwrap();
        clause.add(!*lit);
    }
    clauses.push(clause);

    // 逆方向
    if len != 2{
        let mut clause = rustsat::types::Clause::new();
        for i in (0..len).rev() {
            let lit = lit_map.get(&(cycle[i], cycle[(i+len-1)%len])).unwrap();
            clause.add(!*lit);
        }
        clauses.push(clause);
    }
    // }
    clauses

}

#[allow(dead_code)]
fn asp_blocking_clauses(cycle:&Vec<i32>,encoder: &mut Encoder,g:&Graph, method: i32,balanced:i32) -> Vec<Clause>{
    let mut clauses = Vec::new();
    if method != 4{
        // for cycle in sol_cycles {
        // cycleごとに節を作る
        let mut clause1 = rustsat::types::Clause::new();
        let mut clause2 = rustsat::types::Clause::new();
        for u in cycle.iter() {
            for adjs in g.adjacency_list.get(u).iter(){
                // cycleの中の頂点と、その頂点に接続している頂点のなかで、cycleに入っていないものとの辺を節の中に加える
                for v in adjs.iter(){
                    if !cycle.contains(v){
                        //閉路から出ていく辺と閉路へと入っていく辺両方を同じ節へ追加する
                        if method == 1{
                            let lit1 = encoder.graph_lit_map.get(&(*u,*v)).unwrap();
                            let lit2 = encoder.graph_lit_map.get(&(*v,*u)).unwrap();
                            clause1.extend([*lit1,*lit2]);
                        //閉路から出ていく辺と閉路へと入っていく辺を別々の節へ追加する
                        }else if method == 2{
                            let lit1 = encoder.graph_lit_map.get(&(*u,*v)).unwrap();
                            let lit2 = encoder.graph_lit_map.get(&(*v,*u)).unwrap();
                            clause1.add(*lit1);
                            clause2.add(*lit2);
                        //閉路から出ていく辺のみを節へ追加する
                        }else if method == 3{
                            let lit = encoder.graph_lit_map.get(&(*u,*v)).unwrap();
                            clause1.add(*lit);
                        }
                    }

                }
            }
        }
        if balanced == 0 {
            clauses.push(clause1);
            if clause2.len() != 0{
                clauses.push(clause2);
            }
        }else if balanced == 1 {
            let lits1:Vec<Lit> = clause1.iter().cloned().collect();
            let lits2:Vec<Lit> = clause2.iter().cloned().collect();
            let (adder_clause1,s) = encoder.bailleux_tortalize(lits1.to_vec(),&vec!());
            let (adder_clause2,_) = encoder.bailleux_tortalize(lits2.to_vec(),&s);
            clauses.extend(adder_clause1);
            clauses.extend(adder_clause2);
            clauses.push(clause!(s[0]));

        }else{
            let lits1:Vec<Lit> = clause1.iter().cloned().collect();
            let lits2:Vec<Lit> = clause2.iter().cloned().collect();
            let (adder_clause1,s) = encoder.bailleux_tortalize(lits1.to_vec(),&vec!());
            let (adder_clause2,_) = encoder.bailleux_tortalize(lits2.to_vec(),&s);
            clauses.extend(adder_clause1);
            clauses.extend(adder_clause2);
            clauses.push(clause!(s[0]));
            clauses.push(clause1);
            clauses.push(clause2);
        }
        // }
    }else{
        let highest_v = g.get_highest_degree_vertex();
        // for cycle in sol_cycles {
        //次数が一番高い頂点が含まれてる閉路のみブロック節を追加する
        if cycle.contains(&highest_v){
            let mut clause1 = rustsat::types::Clause::new();
            let mut clause2 = rustsat::types::Clause::new();
            for u in cycle.iter() {
                for adjs in g.adjacency_list.get(u).iter(){
                    // cycleの中の頂点と、その頂点に接続している頂点のなかで、cycleに入っていないものとの辺を節の中に加える
                    for v in adjs.iter(){
                        if !cycle.contains(v){
                        //閉路から出ていく辺と閉路へと入っていく辺を別々の節へ追加する
                        let lit1 = encoder.graph_lit_map.get(&(*u,*v)).unwrap();
                        let lit2 = encoder.graph_lit_map.get(&(*v,*u)).unwrap();
                        clause1.add(*lit1);
                        clause2.add(*lit2);
                        }
                    }
                }
            }
            clauses.push(clause1);
            clauses.push(clause2);
        }
        // }
    }
    clauses
}

    // 要素を追加する関数
fn _add_to_set(set: &mut HashSet<(usize, usize)>, a: usize, b: usize) {
    let pair = if a < b { (a, b) } else { (b, a) };
    set.insert(pair);
}

fn _contains_in_set(set: &HashSet<(usize, usize)>, a: usize, b: usize) -> bool {
    let pair = if a < b { (a, b) } else { (b, a) };
    set.contains(&pair)
}

fn get_active_cycles(cycles: &Vec<Vec<i32>>, active_cycles_number: &Vec<usize>) -> Vec<Vec<i32>> {
    active_cycles_number.iter()
        .map(|&index| cycles[index].clone())
        .collect()
}

fn map_cycle_lengths(cycles: &Vec<Vec<i32>>) -> BTreeMap<usize, i32> {
    let mut length_map = BTreeMap::new();

    for cycle in cycles {
        let length = cycle.len();
        *length_map.entry(length).or_insert(0) += 1;
    }

    length_map
}

pub fn get_boundary_cut_clauses(
    cycle: &[i32],
    encoder: &mut Encoder,
    g: &Graph,
    total_vertices: usize,
    balanced: i32,
) -> Vec<Clause> {
    let mut clauses = Vec::new();
    let cycle_set: HashSet<i32> = cycle.iter().cloned().collect();

    // Technique A3: If |C| > |V| / 2, use complementary set S = V \ C
    let (target_set, is_complementary) = if cycle.len() > total_vertices / 2 && total_vertices > 0 {
        let all_v: HashSet<i32> = g.adjacency_list.keys().cloned().collect();
        let comp_set: HashSet<i32> = all_v.difference(&cycle_set).cloned().collect();
        (comp_set, true)
    } else {
        (cycle_set, false)
    };

    let mut clause_out = rustsat::types::Clause::new();
    let mut clause_in = rustsat::types::Clause::new();

    // Technique A1: Iterate over vertices in target_set and only collect boundary cut edges
    for &u in &target_set {
        if let Some(adjs) = g.adjacency_list.get(&u) {
            for &v in adjs {
                if !target_set.contains(&v) {
                    if let Some(lit_out) = encoder.graph_lit_map.get(&(u, v)) {
                        clause_out.add(*lit_out);
                    }
                    if let Some(lit_in) = encoder.graph_lit_map.get(&(v, u)) {
                        clause_in.add(*lit_in);
                    }
                }
            }
        }
    }

    let (c1, c2) = if is_complementary {
        // By duality: delta^+(V \ C) = delta^-(C) and delta^-(V \ C) = delta^+(C)
        (clause_in, clause_out)
    } else {
        (clause_out, clause_in)
    };

    if balanced == 0 {
        if !c1.is_empty() {
            clauses.push(c1);
        }
        if !c2.is_empty() {
            clauses.push(c2);
        }
    } else if balanced == 1 {
        let lits1: Vec<Lit> = c1.iter().cloned().collect();
        let lits2: Vec<Lit> = c2.iter().cloned().collect();
        let (adder_clause1, s) = encoder.bailleux_tortalize(lits1.to_vec(), &vec![]);
        let (adder_clause2, _) = encoder.bailleux_tortalize(lits2.to_vec(), &s);
        clauses.extend(adder_clause1);
        clauses.extend(adder_clause2);
        if !s.is_empty() {
            clauses.push(clause!(s[0]));
        }
    } else {
        let lits1: Vec<Lit> = c1.iter().cloned().collect();
        let lits2: Vec<Lit> = c2.iter().cloned().collect();
        let (adder_clause1, s) = encoder.bailleux_tortalize(lits1.to_vec(), &vec![]);
        let (adder_clause2, _) = encoder.bailleux_tortalize(lits2.to_vec(), &s);
        clauses.extend(adder_clause1);
        clauses.extend(adder_clause2);
        if !s.is_empty() {
            clauses.push(clause!(s[0]));
        }
        if !c1.is_empty() {
            clauses.push(c1);
        }
        if !c2.is_empty() {
            clauses.push(c2);
        }
    }

    clauses
}

pub fn get_induced_subgraph_sec_clauses(
    cycle: &[i32],
    encoder: &Encoder,
    g: &Graph,
) -> Vec<Clause> {
    let mut clauses = Vec::new();
    let len = cycle.len();
    if len < 3 || len > 6 {
        return clauses;
    }

    let cycle_set: HashSet<i32> = cycle.iter().cloned().collect();

    // 1. Standard forward & reverse exclusion for active cycle C
    let mut fwd_clause = rustsat::types::Clause::new();
    for i in 0..len {
        if let Some(lit) = encoder.graph_lit_map.get(&(cycle[i], cycle[(i + 1) % len])) {
            fwd_clause.add(!*lit);
        }
    }
    if !fwd_clause.is_empty() {
        clauses.push(fwd_clause);
    }

    if len != 2 {
        let mut rev_clause = rustsat::types::Clause::new();
        for i in (0..len).rev() {
            if let Some(lit) = encoder.graph_lit_map.get(&(cycle[i], cycle[(i + len - 1) % len])) {
                rev_clause.add(!*lit);
            }
        }
        if !rev_clause.is_empty() {
            clauses.push(rev_clause);
        }
    }

    // 2. Search for internal chords in G[C] to forbid chord subtours
    // For small |C| <= 6, enumerate chord paths
    for i in 0..len {
        let u = cycle[i];
        if let Some(adjs) = g.adjacency_list.get(&u) {
            for &v in adjs {
                if cycle_set.contains(&v) {
                    let next_u = cycle[(i + 1) % len];
                    let prev_u = cycle[(i + len - 1) % len];
                    // If (u, v) is a chord (not consecutive in C)
                    if v != next_u && v != prev_u {
                        // Forbid shortcut triangle (u, next_u, v) if (next_u, v) is an edge
                        if let Some(next_adjs) = g.adjacency_list.get(&next_u) {
                            if next_adjs.contains(&v) {
                                if let (Some(l1), Some(l2), Some(l3)) = (
                                    encoder.graph_lit_map.get(&(u, next_u)),
                                    encoder.graph_lit_map.get(&(next_u, v)),
                                    encoder.graph_lit_map.get(&(v, u)),
                                ) {
                                    clauses.push(clause!(!*l1, !*l2, !*l3));
                                }
                                if let (Some(l1), Some(l2), Some(l3)) = (
                                    encoder.graph_lit_map.get(&(u, v)),
                                    encoder.graph_lit_map.get(&(v, next_u)),
                                    encoder.graph_lit_map.get(&(next_u, u)),
                                ) {
                                    clauses.push(clause!(!*l1, !*l2, !*l3));
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    clauses
}



// enum MySolver<'a> {
//     Minisat(rustsat_minisat::core::Minisat),
//     Kissat(rustsat_kissat::Kissat<'a>), 
//     CaDiCaL(rustsat_cadical::CaDiCaL<'a, 'a>),
// }


// impl Solve for MySolver<'_> {
//     fn solve(&mut self) -> Result<SolverResult, SolverError>{
//         match self {
//             MySolver::Minisat(solver) => solver.solve(),
//             MySolver::Kissat(solver) => solver.solve(), 
//             MySolver::CaDiCaL(solver) => solver.solve(),
//         }
//     }

//     fn full_solution(&self) -> Result<Assignment, SolverError>{
//         match self {
//             MySolver::Minisat(solver) => solver.full_solution(),
//             MySolver::Kissat(solver) => solver.full_solution(), 
//             MySolver::CaDiCaL(solver) => solver.full_solution(),
//         }
//     }

//     fn add_cnf(&mut self, cnf: Cnf) -> SolveMightFail{
//         match self {
//             MySolver::Minisat(solver) => solver.add_cnf(cnf),
//             MySolver::Kissat(solver) => solver.add_cnf(cnf),
//             MySolver::CaDiCaL(solver) => solver.add_cnf(cnf),
//         }
//     }

//     fn signature(&self) -> &'static str {
//         match self {
//             MySolver::Minisat(solver) => solver.signature(),
//             MySolver::Kissat(solver) => solver.signature(),
//             MySolver::CaDiCaL(solver) => solver.signature(),
//         }
//     }

//     fn lit_val(&self, lit: rustsat::types::Lit) -> Result<rustsat::types::TernaryVal, rustsat::solvers::SolverError> { 
//         match self {
//             MySolver::Minisat(solver) => solver.lit_val(lit),
//             MySolver::Kissat(solver) => solver.lit_val(lit),
//             MySolver::CaDiCaL(solver) => solver.lit_val(lit),
//         }
//     }

//     fn add_clause(&mut self, clause: rustsat::types::Clause) -> Result<(), rustsat::solvers::SolverError>{
//         match self {
//             MySolver::Minisat(solver) => solver.add_clause(clause),
//             MySolver::Kissat(solver) => solver.add_clause(clause),
//             MySolver::CaDiCaL(solver) => solver.add_clause(clause),
//         }
//     }

//     // Implement other required methods in a similar way...
// }

// impl<'a> Extend<rustsat::types::Clause> for MySolver<'a> {
//     fn extend<T: IntoIterator<Item = rustsat::types::Clause>>(&mut self, iter: T) {
//         match self {
//             MySolver::Minisat(solver) => solver.extend(iter),
//             MySolver::Kissat(solver) => solver.extend(iter),
//             MySolver::CaDiCaL(solver) => solver.extend(iter),
//         }
//     }
// }

// impl SolveStats for MySolver<'_>{
//     fn stats(&self) -> SolverStats{
//         match self {
//             MySolver::Minisat(solver) => solver.stats(),
//             MySolver::Kissat(solver) => solver.stats(),
//             MySolver::CaDiCaL(solver) => solver.stats(),
//         }
//     }
// }

// impl<'a> MySolver<'a>{
//     fn set_configuration(&mut self, config:Config) -> SolveMightFail{
//         match self{
//             MySolver::CaDiCaL(solver) => solver.set_configuration(config),
//             _ => Err(SolverError::Api("このソルバーではset_configrationは使用できません。".to_string()))
//         }
//     }
// }

#[cfg(test)]
mod tests_blocking_enhancements {
    use super::*;

    #[test]
    fn test_boundary_cut_complementary_equivalence() {
        // Build a 6-vertex cycle graph 1-2-3-4-5-6-1 with a chord 1-4
        let mut g = Graph::new();
        g.add_edge(1, 2);
        g.add_edge(2, 3);
        g.add_edge(3, 4);
        g.add_edge(4, 5);
        g.add_edge(5, 6);
        g.add_edge(6, 1);
        g.add_edge(1, 4); // chord

        let mut encoder = Encoder::new();
        encoder.encode(&g, 1, 0, 0, 0, 0, 0);

        // Subcycle C = [1, 2, 3, 4] (|C| = 4 > 6/2 = 3) -> Complementary S = [5, 6]
        let c = vec![1, 2, 3, 4];
        let clauses = get_boundary_cut_clauses(&c, &mut encoder, &g, 6, 0);
        assert!(!clauses.is_empty());
        // Verify both out-cut and in-cut clauses exist
        assert_eq!(clauses.len(), 2);

        // Compare with direct calculation (total_vertices = 0 disables complementation)
        let clauses_direct = get_boundary_cut_clauses(&c, &mut encoder, &g, 0, 0);
        assert_eq!(clauses_direct.len(), 2);

        let set_comp_0: HashSet<Lit> = clauses[0].iter().cloned().collect();
        let set_dir_0: HashSet<Lit> = clauses_direct[0].iter().cloned().collect();
        assert_eq!(set_comp_0, set_dir_0);

        let set_comp_1: HashSet<Lit> = clauses[1].iter().cloned().collect();
        let set_dir_1: HashSet<Lit> = clauses_direct[1].iter().cloned().collect();
        assert_eq!(set_comp_1, set_dir_1);
    }

    #[test]
    fn test_induced_subgraph_chord_cycles() {
        // 4-vertex graph 1-2-3-4 with chord (1,3)
        let mut g = Graph::new();
        g.add_edge(1, 2);
        g.add_edge(2, 3);
        g.add_edge(3, 4);
        g.add_edge(4, 1);
        g.add_edge(1, 3); // chord (1, 3)

        let mut encoder = Encoder::new();
        encoder.encode(&g, 1, 0, 0, 0, 0, 0);

        let c = vec![1, 2, 3, 4];
        let clauses = get_induced_subgraph_sec_clauses(&c, &encoder, &g);
        // Must generate at least the main cycle exclusion + chord triangle exclusion
        assert!(!clauses.is_empty());
        // Forward (1) + Reverse (1) for 4-cycle, plus chord triangles (forward & reverse for each of 2 triangles = 4 clauses) -> 6 clauses
        assert_eq!(clauses.len(), 6);

        // Out of bounds cycles (|C| < 3 or |C| > 6) should return empty
        let c_small = vec![1, 2];
        assert!(get_induced_subgraph_sec_clauses(&c_small, &encoder, &g).is_empty());

        let c_large = vec![1, 2, 3, 4, 5, 6, 7];
        assert!(get_induced_subgraph_sec_clauses(&c_large, &encoder, &g).is_empty());
    }

    #[test]
    fn test_hub_aware_swap_node() {
        // Build graph with hub 1 connected to vertices 2..=31
        let mut g = Graph::new();
        for v in 2..=31 {
            g.add_edge(1, v);
            let next_v = if v == 31 { 2 } else { v + 1 };
            g.add_edge(v, next_v);
        }
        let (contracted_g, contractor) = Degree2Contractor::contract(&g);
        let hub_registry = HubRegistry::new(&contracted_g);
        assert!(hub_registry.is_hub_vertex(1));

        let c1 = vec![1, 2, 3];
        let c2 = vec![4, 5, 6];
        let merged = swap_node(&c1, &c2, &contracted_g, &contractor, &hub_registry);
        assert!(merged.is_some());
    }

    #[test]
    fn test_two_opt_hub_prioritization() {
        let mut g = Graph::new();
        for v in 2..=31 {
            g.add_edge(1, v);
            let next_v = if v == 31 { 2 } else { v + 1 };
            g.add_edge(v, next_v);
        }
        g.add_edge(6, 4);
        g.add_edge(9, 7);
        g.add_edge(5, 8);
        g.add_edge(5, 9);
        let (contracted_g, contractor) = Degree2Contractor::contract(&g);
        let hub_registry = HubRegistry::new(&contracted_g);

        let sol_cycles = vec![
            vec![7, 8, 9],
            vec![1, 2, 3],
            vec![4, 5, 6],
        ];
        let mut encoder = Encoder::new();
        encoder.encode(&contracted_g, 1, 0, 0, 0, 0, 0);

        let (_clauses, merged_cycles) = two_opt(
            &sol_cycles,
            &mut encoder,
            &contracted_g,
            &contractor,
            &hub_registry,
            3,
            0,
            1,
            0,
        );
        // All 9 vertices should be merged into a single cycle
        assert_eq!(merged_cycles.len(), 1);
        assert_eq!(merged_cycles[0].len(), 9);
    }

    #[test]
    fn test_short_cycle_pruning_triangles_and_quads() {
        // 5-vertex graph: 1-2-3-4-5-1 with chord 1-3
        // Triangles: {1, 2, 3} (2 clauses)
        // 4-cycles: {1, 3, 4, 5} with edges 1-3, 3-4, 4-5, 5-1 (2 clauses)
        let mut g = Graph::new();
        g.add_edge(1, 2);
        g.add_edge(2, 3);
        g.add_edge(3, 1);
        g.add_edge(3, 4);
        g.add_edge(4, 5);
        g.add_edge(5, 1);

        let mut encoder = Encoder::new();
        let mut cnf = encoder.encode(&g, 1, 0, 0, 0, 0, 0);
        let base_clause_count = cnf.len();

        let added = add_global_short_cycle_cuts(&g, &encoder, &mut cnf, 4);
        assert_eq!(added, 4);
        assert_eq!(cnf.len(), base_clause_count + 4);

        // Verify with max_cycle_len = 3 (only triangles)
        let mut cnf3 = encoder.encode(&g, 1, 0, 0, 0, 0, 0);
        let added3 = add_global_short_cycle_cuts(&g, &encoder, &mut cnf3, 3);
        assert_eq!(added3, 2);

        // Verify with max_cycle_len = 2 (none)
        let mut cnf2 = encoder.encode(&g, 1, 0, 0, 0, 0, 0);
        let added2 = add_global_short_cycle_cuts(&g, &encoder, &mut cnf2, 2);
        assert_eq!(added2, 0);

        // Verify safety on 3-vertex triangle graph (must NOT prune N=3)
        let mut g3 = Graph::new();
        g3.add_edge(1, 2);
        g3.add_edge(2, 3);
        g3.add_edge(3, 1);
        let mut encoder3 = Encoder::new();
        let mut cnf_g3 = encoder3.encode(&g3, 1, 0, 0, 0, 0, 0);
        let added_g3 = add_global_short_cycle_cuts(&g3, &encoder3, &mut cnf_g3, 4);
        assert_eq!(added_g3, 0);

        // Verify safety on 4-vertex 4-cycle graph (must NOT prune N=4)
        let mut g4 = Graph::new();
        g4.add_edge(1, 2);
        g4.add_edge(2, 3);
        g4.add_edge(3, 4);
        g4.add_edge(4, 1);
        let mut encoder4 = Encoder::new();
        let mut cnf_g4 = encoder4.encode(&g4, 1, 0, 0, 0, 0, 0);
        let added_g4 = add_global_short_cycle_cuts(&g4, &encoder4, &mut cnf_g4, 4);
        assert_eq!(added_g4, 0);
    }

    #[test]
    fn test_cluster_cut_constraints() {
        // Build a graph with 25 hubs (1..=25, complete clique so degree = 24 >= 20)
        // and 50 satellite vertices (101..=150)
        let mut g = Graph::new();
        // 25-clique for hubs
        for i in 1..=25 {
            for j in (i + 1)..=25 {
                g.add_edge(i, j);
            }
        }
        // Satellite cluster of size 50: 101..=150
        for v in 101..150 {
            g.add_edge(v, v + 1);
        }
        g.add_edge(150, 101);

        // Boundary cut edges from satellite cluster to hubs
        g.add_edge(101, 1);
        g.add_edge(110, 2);
        g.add_edge(120, 3);
        g.add_edge(130, 4);
        g.add_edge(140, 5);

        let mut encoder = Encoder::new();
        let mut cnf = encoder.encode(&g, 1, 0, 0, 0, 0, 0);
        let base_clause_count = cnf.len();

        let added = add_cluster_cut_constraints(&g, &encoder, &mut cnf);
        assert!(added >= 2, "Must add at least out-cut and in-cut clauses for cluster >= 50, got {}", added);
        assert_eq!(cnf.len(), base_clause_count + added);

        // Verify safety on small graph (< 50 satellite vertices)
        let mut g_small = Graph::new();
        for i in 1..=25 {
            for j in (i + 1)..=25 {
                g_small.add_edge(i, j);
            }
        }
        for v in 101..120 {
            g_small.add_edge(v, v + 1);
        }
        g_small.add_edge(120, 101);
        g_small.add_edge(101, 1);
        g_small.add_edge(110, 2);

        let mut encoder_small = Encoder::new();
        let mut cnf_small = encoder_small.encode(&g_small, 1, 0, 0, 0, 0, 0);
        let added_small = add_cluster_cut_constraints(&g_small, &encoder_small, &mut cnf_small);
        assert_eq!(added_small, 0, "Must not add cluster cuts when satellite clusters are < 50");
    }

    #[test]
    fn test_limit_conflicts_with_assumptions() {
        use rustsat_cadical::CaDiCaL;
        use rustsat::solvers::{LimitConflicts, Solve, SolveIncremental, SolverResult};
        let mut solver = CaDiCaL::default();
        let lit1 = Lit::positive(0);
        let lit2 = Lit::positive(1);

        // Add clause: lit1 or lit2
        solver.add_clause(clause![lit1, lit2]).unwrap();

        // Limit conflicts and solve with assumptions
        let _ = solver.limit_conflicts(Some(5000));
        let res = solver.solve_assumps(&[!lit1, !lit2]);
        let _ = solver.limit_conflicts(None);

        assert_eq!(res.unwrap(), SolverResult::Unsat);
    }
}
