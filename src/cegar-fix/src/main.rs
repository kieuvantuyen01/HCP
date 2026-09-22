mod chained_lk;
mod contraction;
mod encoder;
mod file_operations;
mod graph;
mod hcp_solver;
mod hub_registry;
mod hub_sub_hcp;
mod ils_patcher;
mod macro_solver;
mod matching_patcher;
pub mod modular_solver;
pub mod modular_tree;
mod options;
mod parallel_sub_hcp;
mod patching;
pub mod stem_cycle_patcher;
pub mod two_tier_decomposer;
pub mod pinpointed_strip_solver;
pub mod global_demand_coordinator;
pub mod macro_splicer;
pub mod two_tier_orchestrator;
pub mod staged_subcycle_filter;
pub mod dual_cut_generator;
pub mod staged_lazy_smt_solver;
pub mod subcycle_absorber;
pub mod bridge_cut_generator;
pub mod backbone_freezer;
pub mod auto_classifier;
pub mod cycle_chain_absorber;
pub mod tour_verifier;
pub mod hybrid_orchestrator;
pub mod snark_bridge;
pub mod gadget_parity;
pub mod component_meta_graph;
pub mod macro_mtz_encoder;
pub mod cut_selector;
pub mod solver_reseeder;
pub mod hemisphere_splicer;
pub mod static_cycle_cutter;
pub mod boundary_alternating_patcher;
pub mod metagraph_router;
pub mod parallel_sat_portfolio;
pub mod macro_cycle_stitcher;
pub mod giant_cycle_stitcher;
pub mod macro_crt_encoder;
pub mod macro_lfsr_encoder;
pub mod transitive_macro_splicer;
pub mod interface_port_synchronizer;
pub mod inverse_3sat_synthesizer;
pub mod hub_hierarchical_decomposer;
pub mod multi_opt_sat_splicer;
pub mod empirical_backbone_cutter;
pub mod cnf_subsumer;
pub mod twin_giant_splicer;
pub mod macro_component_splicer;
pub mod sat_macro_patcher;
pub mod gadget_path_absorber;
pub mod incremental_sat;
pub mod macro_crossover_splicer;
pub mod quotient_block_cutter;
pub mod bipartite_module_detector;
pub mod module_dual_path_extractor;
pub mod module_state_cnf_encoder;
pub mod balanced_pair_cutset;
pub mod localized_sat_repair;
pub mod alternating_port_engine;
pub mod port_corridor_lns;
pub mod modular_ring_dp_solver;


use contraction::Degree2Contractor;
use hub_registry::HubRegistry;
use std::time::Instant;
use log::info;

fn main() {
    env_logger::init();
    info!("プログラム開始");
    let instant = Instant::now();

    let matches = options::get_options();

    // solver,encodingのオプションをintで受け取る
    let solver = matches.value_of_t::<i32>("solver").unwrap_or(0);
    let encoding = matches.value_of_t::<i32>("encoding").unwrap_or(0);
    let blocking = matches.value_of_t::<i32>("blocking").unwrap_or(0);
    let symmetry = matches.value_of_t::<i32>("symmetry").unwrap_or(0);
    let two_opt = matches.value_of_t::<i32>("2-opt").unwrap_or(0);
    let three_opt = matches.value_of_t::<i32>("three-opt").unwrap_or(0);
    let loop_prohibition = matches.value_of_t::<i32>("loop-prohibition").unwrap_or(0);
    let cnf_normalize = matches.value_of_t::<i32>("cnf-normalize").unwrap_or(0);
    let balanced = matches.value_of_t::<i32>("balanced").unwrap_or(0);
    let de_arcify = matches.value_of_t::<i32>("de-arcify").unwrap_or(0);
    let config = matches.value_of_t::<i32>("set-configration").unwrap_or(0);
    let degree_order = matches.value_of_t::<i32>("degree-order").unwrap_or(0);
    let arcs_order = matches.value_of_t::<i32>("arcs-order").unwrap_or(0);
    let cegar_fallback = matches.value_of_t::<i32>("cegar-fallback").unwrap_or(0);
    let mtz_stall = matches.value_of_t::<i32>("mtz-stall").unwrap_or(0);
    let adaptive_escalation = matches.value_of_t::<i32>("adaptive-escalation").unwrap_or(1);
    let sub_hcp_timeout = matches.value_of_t::<u64>("sub-hcp-timeout").unwrap_or(60);
    let max_cluster_size = matches.value_of_t::<usize>("max-cluster-size").unwrap_or(500);
    let is_two_tier = matches.is_present("two-tier") && matches.value_of("two-tier").map_or(true, |v| v != "0");
    let is_staged_smt = matches.is_present("staged-smt") && matches.value_of("staged-smt").map_or(true, |v| v != "0");
    let auto_mode = matches.value_of("auto").map_or(true, |v| v != "0");
    let structural_decomposition = matches
        .value_of("structural-decomposition")
        .map_or(true, |v| v != "0");
    let degree2_contraction = matches
        .value_of("degree2-contraction")
        .map_or(true, |v| v != "0");
    let cycle_repair = matches
        .value_of("cycle-repair")
        .map_or(true, |v| v != "0");
    let macro_gadget = matches.is_present("macro-gadget") && matches.value_of("macro-gadget").map_or(true, |v| v != "0");
    let bounded_freezer = matches.is_present("bounded-freezer") && matches.value_of("bounded-freezer").map_or(true, |v| v != "0");
    let alternating_engine = if matches.is_present("no-alternating-engine") {
        false
    } else if matches.is_present("alternating-engine") {
        matches.value_of("alternating-engine").map_or(true, |v| v != "0")
    } else {
        true
    };
    let timeout_secs = matches.value_of_t::<f64>("timeout").unwrap_or(1800.0);
    let output_tour_path = matches.value_of("output-tour").map(|s| s.to_string());
    // solver,encodingのオプションを&strで受け取る
    let input_filename = matches
        .value_of("input")
        .or_else(|| matches.value_of("positional_input"))
        .unwrap_or("default");
    let output_foldername = matches.value_of("output").unwrap_or("default");

    println!("solve {}", input_filename);
    // let g = instance();
    let mut g = file_operations::input_to_graph(input_filename);
    if de_arcify != 0{
        g.remove_redundant_arcs();
    }
    if is_staged_smt {
        let opt = staged_lazy_smt_solver::StagedLazySmtOptions {
            max_batch_size: 500,
            timeout_secs,
            output_path: output_tour_path,
        };
        let res = staged_lazy_smt_solver::solve_staged_lazy_smt(&g, &opt);
        if let Some(tour) = res {
            println!("s SATISFIABLE");
            print!("solution: \n");
            for v in &tour {
                print!("{} ", v);
            }
            println!();
            println!("overall time = {:?}", instant.elapsed());
        } else {
            println!("s UNSATISFIABLE");
            println!("overall time = {:?}", instant.elapsed());
        }
        return;
    }
    if is_two_tier {
        let opt = two_tier_orchestrator::TwoTierSolverOptions {
            timeout_secs,
            max_iterations: 50_000,
            enable_patching: cycle_repair,
            output_path: output_tour_path,
        };
        let res = two_tier_orchestrator::solve_graph_two_tier(&g, &opt);
        if let Some(tour) = res {
            println!("s SATISFIABLE");
            print!("solution: \n");
            for v in &tour {
                print!("{} ", v);
            }
            println!();
            println!("overall time = {:?}", instant.elapsed());
        } else {
            println!("s UNSATISFIABLE");
            println!("overall time = {:?}", instant.elapsed());
        }
        return;
    }

    let has_manual_overrides = matches.is_present("encoding")
        || matches.is_present("blocking")
        || matches.is_present("symmetry")
        || matches.is_present("2-opt")
        || matches.is_present("three-opt")
        || matches.is_present("set-configration");

    let mut hybrid_opts = hybrid_orchestrator::HybridOptions {
        auto_mode: auto_mode && !macro_gadget && !bounded_freezer && !matches.is_present("ablation"),
        structural_decomposition,
        degree2_contraction,
        cycle_repair,
        timeout_secs,
        output_tour: output_tour_path.clone(),
        macro_gadget,
        bounded_freezer,
        alternating_engine,
    };
    if let Some(abl_str) = matches.value_of("ablation") {
        if let Ok(mode) = abl_str.parse::<usize>() {
            hybrid_opts.apply_ablation_mode(mode);
            println!("Ablation mode active: C{} (macro_gadget={}, bounded_freezer={})", mode, hybrid_opts.macro_gadget, hybrid_opts.bounded_freezer);
        }
    }

    if (auto_mode || matches.is_present("ablation") || hybrid_opts.macro_gadget || hybrid_opts.bounded_freezer) && !has_manual_overrides {
        let res = hybrid_orchestrator::HybridOrchestrator::solve(&g, &hybrid_opts);
        if res.is_none() {
            println!("s UNSATISFIABLE");
        }
        println!("overall time = {:?}", instant.elapsed());
        return;
    }
    if g.has_articulation_points() {
        println!("Graph has cut-vertex or is disconnected.");
        println!("s UNSATISFIABLE");
        return;
    }
    let pruned = g.prune_degree2_triangles();
    if pruned > 0 {
        println!("Pruned {} degree-2 triangle shortcut edges", pruned);
    }
    let (contracted_g, contractor) = if degree2_contraction {
        Degree2Contractor::contract(&g)
    } else {
        (g.clone(), Degree2Contractor::identity(&g))
    };

    if let Some(cycle) = contractor.is_direct_cycle {
        println!("Graph is a single 2-regular Hamiltonian cycle.");
        print!("solution: \n");
        for v in &cycle {
            print!("{} ", v);
        }
        println!();
        println!("s SATISFIABLE");
        println!("overall time = {:?}", instant.elapsed());
        return;
    }

    if contractor.is_infeasible {
        println!("Infeasible degree-2 structure detected.");
        println!("s UNSATISFIABLE");
        return;
    }

    if contractor.contracted_vertices_count < contractor.original_vertices_count {
        println!(
            "Degree-2 contraction: compressed graph from {} to {} vertices (reduced by {}%)",
            contractor.original_vertices_count,
            contractor.contracted_vertices_count,
            (contractor.original_vertices_count - contractor.contracted_vertices_count) * 100 / contractor.original_vertices_count
        );
    }
    let hub_registry = HubRegistry::new(&contracted_g);
    if !hub_registry.hub_vertices.is_empty() {
        println!(
            "Dense Hub optimization: detected {} hub vertices (sample: {:?})",
            hub_registry.hub_vertices.len(),
            &hub_registry.hub_vertices[..hub_registry.hub_vertices.len().min(5)]
        );
    }
    let time1 = instant.elapsed();
    // println!("encodhing time = {:?} sec",instant.elapsed().as_secs());
    // let instant2 = Instant::now();

    // println!("solver={},encoding={}",solver,encoding);
    // println!("{:?}",g);
    println!("file input time = {:?}", time1);
    let tour = hcp_solver::solve_hamilton_with_decomposition(contracted_g, &contractor, &hub_registry, solver, encoding, blocking, symmetry, two_opt, loop_prohibition, cnf_normalize, balanced, de_arcify,config,degree_order,arcs_order,three_opt,cegar_fallback,mtz_stall,adaptive_escalation,sub_hcp_timeout,max_cluster_size,timeout_secs,instant,output_foldername, if hybrid_opts.macro_gadget { 1 } else { 0 }, if hybrid_opts.bounded_freezer { 1 } else { 0 }, if hybrid_opts.alternating_engine { 1 } else { 0 }, structural_decomposition, cycle_repair);
    if let Some(ref t) = tour {
        if let Some(ref out_path) = output_tour_path {
            if let Err(e) = tour_verifier::TourVerifier::write_tsplib_hcp(t, "tour", out_path) {
                eprintln!("Warning: failed to write tour to {}: {}", out_path, e);
            } else {
                println!("Wrote certified tour to {}", out_path);
            }
        }
    }
    let time2 = instant.elapsed() - time1;

    // println!("solving time = {:?} sec",instant2.elapsed().as_secs());
    println!("solving time = {:?}", time2);
    println!("overall time = {:?}", instant.elapsed());
    info!("プログラム終了");
}

fn _instance() -> graph::Graph {
    let mut g = graph::Graph::new();
    g.add_edge(1, 2);
    g.add_edge(1, 8);
    g.add_edge(2, 3);
    g.add_edge(3, 4);
    g.add_edge(3, 7);
    g.add_edge(4, 5);
    g.add_edge(4, 6);
    g.add_edge(4, 8);
    g.add_edge(5, 6);
    g.add_edge(6, 7);
    g.add_edge(7, 8);
    g
}

// fn add_cl(lit_names: BTreeMap<Lit, String>, instance: &mut SatInstance) {
//     let mut map_i: BTreeMap<String, Vec<Lit>> = BTreeMap::new();
//     let mut map_j: BTreeMap<String, Vec<Lit>> = BTreeMap::new();

//     for (lit, name) in lit_names.iter() {
//         let parts: Vec<&str> = name.split("_").collect();
//         map_i.entry(parts[0].to_string()).or_insert(Vec::new()).push(lit.clone());
//         map_j.entry(parts[1].to_string()).or_insert(Vec::new()).push(lit.clone());
//     }

//     for (_, lits) in map_i.iter() {
//         for i in 0..lits.len() {
//             for j in i+1..lits.len() {
//                 instance.add_binary(!lits[i], !lits[j]);
//             }
//         }
//         instance.add_clause(lits.as_slice().into());
//     }

//     for (_, lits) in map_j.iter() {
//         for i in 0..lits.len() {
//             for j in i+1..lits.len() {
//                 instance.add_binary(!lits[i], !lits[j]);
//             }
//         }
//         instance.add_clause(lits.as_slice().into());
//     }
// }

// fn check_hamilton(edges: Vec<&str>) -> bool {
//     let mut graph: BTreeMap<&str,&str> = std::collections::BTreeMap::new();
//     for edge in edges {
//         let nodes: Vec<&str> = edge.split('_').collect();
//         let from = nodes[0];
//         let to = nodes[1];
//         graph.insert(from, to);
//     }

//     let start_node = graph.keys().next().unwrap();
//     let mut current_node = start_node;
//     let mut visited = std::collections::BTreeSet::new();
//     loop {
//         visited.insert(current_node);
//         current_node = match graph.get(current_node) {
//             Some(node) => node,
//             None => break,
//         };
//         if visited.contains(current_node) {
//             break;
//         }
//     }

//     visited.len() == graph.len() && current_node == start_node
// }

// fn solve_hamilton(lit_names: BTreeMap<Lit, String>,solver: &mut rustsat_minisat::core::Minisat){

//     let res = solver.solve().unwrap();
//     if res == SolverResult::Sat{
//         let sol = solver.full_solution().unwrap();
//         let true_lits: Vec<Lit> = lit_names.iter().filter_map(|(lit, _)| if sol[lit.var()] == TernaryVal::True { Some(lit.clone()) } else { None }).collect();
//         let true_lits_names: Vec<&str> = true_lits.iter().map(|lit| lit_names.get(lit).unwrap().as_str()).collect();
//         if check_hamilton(true_lits_names.clone()){
//             println!("{:?}",true_lits_names);
//         }else{
//             // let block_clause:Vec<Lit> = true_lits.iter().map(|&lit|!lit).collect();
//             let mut block_clause = rustsat::types::Clause::new();
//             for lit in true_lits.iter(){
//                 block_clause.add(!*lit);
//             }
//             // let block_clause:rustsat::types::Clause = true_lits.iter().collect();
//             let _ = solver.add_clause(block_clause);
//             solve_hamilton(lit_names,solver)
//         }
//     }else{
//         println!("UNSAT");
//     }

// }
