use crate::core::graph::Graph;
use crate::core::solver_utils::create_solver_with_deadline;
use rustsat::clause;
use rustsat::instances::{BasicVarManager, ManageVars};
use rustsat::solvers::{Solve, SolverResult};
use rustsat::types::{Clause, Lit, TernaryVal};
use std::collections::{HashMap, HashSet};
use std::time::{Duration, Instant};

/// Checks if graph `g` is 3-regular (cubic): all vertices have degree exactly 3.
pub fn can_solve_directed_cubic(g: &Graph) -> bool {
    let n = g.adjacency_list.len();
    if n < 4 {
        return false;
    }
    g.adjacency_list.values().all(|nbrs| nbrs.len() == 3)
}

/// Solves Hamiltonian cycle on 3-regular graph `g` using Directed 2-Loop SAT-CEGAR.
pub fn solve_directed_cubic(g: &Graph, timeout_secs: f64) -> Result<Vec<i32>, String> {
    solve_directed_cubic_with_seed(g, timeout_secs, 1)
}

pub fn solve_directed_cubic_with_seed(
    g: &Graph,
    timeout_secs: f64,
    seed: i32,
) -> Result<Vec<i32>, String> {
    let deadline = Instant::now() + Duration::from_secs_f64(timeout_secs);
    let mut nodes: Vec<i32> = g.adjacency_list.keys().copied().collect();
    nodes.sort_unstable();
    let n = nodes.len();
    if n < 3 {
        return Err("Graph has fewer than 3 vertices".to_string());
    }

    // Map each directed edge (u, v) where {u, v} in E to a SAT literal
    let mut var_mgr = BasicVarManager::default();
    let mut edge_vars: HashMap<(i32, i32), Lit> = HashMap::new();

    for &u in &nodes {
        if let Some(nbrs) = g.adjacency_list.get(&u) {
            for &v in nbrs {
                let lit = var_mgr.new_var().pos_lit();
                edge_vars.insert((u, v), lit);
            }
        }
    }

    let mut solver = create_solver_with_deadline(deadline);
    let _ = solver.set_option("seed", seed.max(1));

    // 1. Out-degree == 1 for every vertex u:
    //    \sum_{v in N(u)} x_{(u, v)} = 1
    for &u in &nodes {
        let nbrs = match g.adjacency_list.get(&u) {
            Some(n) => n,
            None => return Err("UNSAT".to_string()),
        };
        let out_lits: Vec<Lit> = nbrs.iter().map(|&v| edge_vars[&(u, v)]).collect();
        if out_lits.is_empty() {
            return Err("UNSAT".to_string());
        }
        // At-least-1
        let _ = solver.add_clause(Clause::from_iter(out_lits.iter().copied()));
        // At-most-1
        for i in 0..out_lits.len() {
            for j in (i + 1)..out_lits.len() {
                let _ = solver.add_clause(clause![!out_lits[i], !out_lits[j]]);
            }
        }
    }

    // 2. In-degree == 1 for every vertex v:
    //    \sum_{u in N(v)} x_{(u, v)} = 1
    for &v in &nodes {
        let nbrs = match g.adjacency_list.get(&v) {
            Some(n) => n,
            None => return Err("UNSAT".to_string()),
        };
        let in_lits: Vec<Lit> = nbrs.iter().map(|&u| edge_vars[&(u, v)]).collect();
        if in_lits.is_empty() {
            return Err("UNSAT".to_string());
        }
        // At-least-1
        let _ = solver.add_clause(Clause::from_iter(in_lits.iter().copied()));
        // At-most-1
        for i in 0..in_lits.len() {
            for j in (i + 1)..in_lits.len() {
                let _ = solver.add_clause(clause![!in_lits[i], !in_lits[j]]);
            }
        }
    }

    // 3. 2-loop prohibition:
    //    For every undirected edge {u, v}, \neg x_{(u, v)} \lor \neg x_{(v, u)}
    for &u in &nodes {
        if let Some(nbrs) = g.adjacency_list.get(&u) {
            for &v in nbrs {
                if u < v {
                    let lit_uv = edge_vars[&(u, v)];
                    let lit_vu = edge_vars[&(v, u)];
                    let _ = solver.add_clause(clause![!lit_uv, !lit_vu]);
                }
            }
        }
    }

    // CEGAR loop
    let mut iter_count = 0;
    loop {
        if Instant::now() >= deadline {
            return Err("TIMEOUT".to_string());
        }

        let res = match solver.solve() {
            Ok(r) => r,
            Err(e) => return Err(format!("CaDiCaL error: {:?}", e)),
        };

        match res {
            SolverResult::Unsat => return Err("UNSAT".to_string()),
            SolverResult::Interrupted => return Err("TIMEOUT".to_string()),
            SolverResult::Sat => {
                iter_count += 1;
                let sol = solver
                    .full_solution()
                    .map_err(|e| format!("Failed to get full solution: {:?}", e))?;

                // Map each vertex to its directed successor
                let mut succ: HashMap<i32, i32> = HashMap::with_capacity(n);
                for &u in &nodes {
                    if let Some(nbrs) = g.adjacency_list.get(&u) {
                        for &v in nbrs {
                            let lit = edge_vars[&(u, v)];
                            if sol.lit_value(lit) == TernaryVal::True {
                                succ.insert(u, v);
                                break;
                            }
                        }
                    }
                }

                // Extract directed cycles
                let mut visited = HashSet::with_capacity(n);
                let mut cycles: Vec<Vec<i32>> = Vec::new();

                for &start_u in &nodes {
                    if visited.contains(&start_u) {
                        continue;
                    }
                    let mut cyc = Vec::new();
                    let mut curr = start_u;
                    while !visited.contains(&curr) {
                        visited.insert(curr);
                        cyc.push(curr);
                        match succ.get(&curr) {
                            Some(&nxt) => curr = nxt,
                            None => break,
                        }
                    }
                    if !cyc.is_empty() {
                        cycles.push(cyc);
                    }
                }

                // Termination: single cycle containing all vertices
                if cycles.len() == 1 && cycles[0].len() == n {
                    println!(
                        "[directed_cegar] Solved 3-regular cubic graph in {} iterations!",
                        iter_count
                    );
                    return Ok(cycles.pop().unwrap());
                }

                // Add directed cut clauses: for each cycle C, at least one edge must leave C
                // \bigvee_{u in C, v not in C, (u, v) in E} x_{(u, v)}
                for cyc in &cycles {
                    let cyc_set: HashSet<i32> = cyc.iter().copied().collect();
                    let mut cut_lits = Vec::new();

                    for &u in cyc {
                        if let Some(nbrs) = g.adjacency_list.get(&u) {
                            for &v in nbrs {
                                if !cyc_set.contains(&v) {
                                    if let Some(&lit) = edge_vars.get(&(u, v)) {
                                        cut_lits.push(lit);
                                    }
                                }
                            }
                        }
                    }

                    if cut_lits.is_empty() {
                        return Err("UNSAT".to_string());
                    }
                    let _ = solver.add_clause(Clause::from_iter(cut_lits));
                }
            }
        }
    }
}
