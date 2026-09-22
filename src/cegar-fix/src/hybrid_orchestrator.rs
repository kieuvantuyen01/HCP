use crate::auto_classifier::{AutoTopologyClassifier, TargetTrack};
use crate::contraction::Degree2Contractor;
use crate::graph::Graph;
use crate::hcp_solver;
use crate::hub_registry::HubRegistry;
use crate::tour_verifier::TourVerifier;
use crate::two_tier_orchestrator::{TwoTierOptions, TwoTierOrchestrator};
use std::time::Instant;

#[derive(Debug, Clone)]
pub struct HybridOptions {
    pub auto_mode: bool,
    pub structural_decomposition: bool,
    pub degree2_contraction: bool,
    pub cycle_repair: bool,
    pub timeout_secs: f64,
    pub output_tour: Option<String>,
    pub macro_gadget: bool,
    pub bounded_freezer: bool,
    pub alternating_engine: bool,
}

impl Default for HybridOptions {
    fn default() -> Self {
        Self {
            auto_mode: true,
            structural_decomposition: true,
            degree2_contraction: true,
            cycle_repair: true,
            timeout_secs: 1800.0,
            output_tour: None,
            macro_gadget: false,
            bounded_freezer: false,
            alternating_engine: true,
        }
    }
}

impl HybridOptions {
    pub fn apply_ablation_mode(&mut self, mode: usize) {
        match mode {
            0 => {
                self.macro_gadget = false;
                self.bounded_freezer = false;
            }
            1 => {
                self.macro_gadget = false;
                self.bounded_freezer = true;
            }
            2 => {
                self.macro_gadget = true;
                self.bounded_freezer = false;
            }
            3 => {
                self.macro_gadget = true;
                self.bounded_freezer = true;
            }
            _ => {}
        }
    }
}

pub struct HybridOrchestrator;

impl HybridOrchestrator {
    pub fn solve(g: &Graph, options: &HybridOptions) -> Option<Vec<i32>> {
        let features = AutoTopologyClassifier::extract_features(g);
        let track = if options.macro_gadget || options.bounded_freezer {
            TargetTrack::GadgetInterfaceParity
        } else if options.auto_mode {
            AutoTopologyClassifier::classify(&features)
        } else {
            TargetTrack::GadgetInterfaceParity
        };

        println!(
            "AutoClassifier: N={}, M={}, Density={:.2}, Hubs={} -> Track: {:?}",
            features.n, features.m, features.density, features.hub_count, track
        );
        println!(
            "Structural decomposition: {}",
            if options.structural_decomposition { "enabled" } else { "disabled" }
        );
        println!(
            "Degree-2 contraction: {}",
            if options.degree2_contraction { "enabled" } else { "disabled" }
        );
        println!(
            "Cycle repair: {}",
            if options.cycle_repair { "enabled" } else { "disabled" }
        );

        let tour = match track {
            TargetTrack::B1LadderTwoTier if options.structural_decomposition => {
                let tt_opts = TwoTierOptions {
                    timeout_secs: options.timeout_secs,
                    output_tour: options.output_tour.clone(),
                    cycle_repair: options.cycle_repair,
                };
                let res = TwoTierOrchestrator::solve(g, &tt_opts);
                if res.is_some() {
                    res
                } else {
                    println!("TwoTier decomposition finished without tour. Falling back to Global CEGAR solver...");
                    Self::solve_global(g, options, TargetTrack::GeneralCaDiCaL)
                }
            }
            TargetTrack::B1LadderTwoTier => {
                println!("TwoTier decomposition skipped by experimental profile; using the matched global CEGAR path.");
                Self::solve_global(g, options, TargetTrack::GeneralCaDiCaL)
            }
            _ => Self::solve_global(g, options, track),
        };

        if let Some(ref t) = tour {
            if let Some(ref out_path) = options.output_tour {
                if let Err(e) = TourVerifier::write_tsplib_hcp(t, "tour", out_path) {
                    eprintln!("Warning: failed to write tour to {}: {}", out_path, e);
                } else {
                    println!("Wrote certified tour to {}", out_path);
                }
            }
        }

        tour
    }

    fn solve_global(g: &Graph, options: &HybridOptions, track: TargetTrack) -> Option<Vec<i32>> {
        // Pre-flight check: cut-vertices or disconnected components
        if g.has_articulation_points() {
            println!("Graph has cut-vertex or is disconnected.");
            println!("s UNSATISFIABLE");
            return None;
        }
        let mut uncontracted_g = g.clone();
        let pruned = uncontracted_g.prune_degree2_triangles();
        if pruned > 0 {
            println!("Pruned {} degree-2 triangle shortcut edges", pruned);
        }
        let (contracted_g, contractor) = if options.degree2_contraction {
            Degree2Contractor::contract(&uncontracted_g)
        } else {
            (
                uncontracted_g.clone(),
                Degree2Contractor::identity(&uncontracted_g),
            )
        };

        if let Some(cycle) = contractor.is_direct_cycle.clone() {
            println!("Graph is a single 2-regular Hamiltonian cycle.");
            let uncontracted = contractor.uncontract_cycle(&cycle);
            let line = uncontracted
                .iter()
                .map(|i| i.to_string())
                .collect::<Vec<String>>()
                .join(" ");
            println!();
            println!("solution: ");
            println!("{}\n", line);
            println!("s SATISFIABLE");
            if let Some(ref out_path) = options.output_tour {
                let _ = TourVerifier::write_tsplib_hcp(&uncontracted, "tour", out_path);
            }
            return Some(uncontracted);
        }

        if contractor.is_infeasible {
            println!("Infeasible degree-2 structure detected.");
            println!("s UNSATISFIABLE");
            return None;
        }

        if contractor.contracted_vertices_count < contractor.original_vertices_count {
            println!(
                "Degree-2 contraction: compressed graph from {} to {} vertices (reduced by {}%)",
                contractor.original_vertices_count,
                contractor.contracted_vertices_count,
                (contractor.original_vertices_count - contractor.contracted_vertices_count) * 100
                    / contractor.original_vertices_count
            );
        }

        let hub_reg = HubRegistry::new(&contracted_g);
        if !hub_reg.hub_vertices.is_empty() {
            println!(
                "Dense Hub optimization: detected {} hub vertices (sample: {:?})",
                hub_reg.hub_vertices.len(),
                &hub_reg.hub_vertices[..hub_reg.hub_vertices.len().min(5)]
            );
        }

        let start = Instant::now();
        let mg = if options.macro_gadget { 1 } else { 0 };
        let bf = if options.bounded_freezer { 1 } else { 0 };
        let alt = if options.alternating_engine { 1 } else { 0 };
        match track {
            TargetTrack::SnarkKeyBridge => {
                // CaDiCaL encoding: -e 0 -b 3 -l 1 --three-opt 1 --set-configration 1
                hcp_solver::solve_hamilton_with_decomposition(
                    contracted_g,
                    &contractor,
                    &hub_reg,
                    0, 0, 3, 0, 1, 1, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 60, 200,
                    options.timeout_secs,
                    start,
                    "",
                    mg,
                    bf,
                    alt,
                    options.structural_decomposition,
                    options.cycle_repair,
                )
            }
            TargetTrack::GadgetInterfaceParity | TargetTrack::B2SinzChainSMT => {
                // Sinz encoding: -e 1 -b 3 -y 2 -t 3 -l 1 --three-opt 1 --set-configration 1
                hcp_solver::solve_hamilton_with_decomposition(
                    contracted_g,
                    &contractor,
                    &hub_reg,
                    0, 1, 3, 2, 3, 1, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 60, 200,
                    options.timeout_secs,
                    start,
                    "",
                    mg,
                    bf,
                    alt,
                    options.structural_decomposition,
                    options.cycle_repair,
                )
            }
            _ => {
                // General CaDiCaL encoding: -e 0 -b 3 -l 1
                hcp_solver::solve_hamilton_with_decomposition(
                    contracted_g,
                    &contractor,
                    &hub_reg,
                    0, 0, 3, 0, 1, 1, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 60, 200,
                    options.timeout_secs,
                    start,
                    "",
                    mg,
                    bf,
                    alt,
                    options.structural_decomposition,
                    options.cycle_repair,
                )
            }
        }
    }
}
