use crate::core::graph::Graph;
use crate::core::solver_utils::create_solver_with_deadline;
use crate::core::tour_verifier::TourVerifier;
use rustsat::clause;
use rustsat::instances::{BasicVarManager, ManageVars};
use rustsat::solvers::{Solve, SolverResult};
use rustsat::types::{Clause, Lit, TernaryVal};
use rustsat_cadical::CaDiCaL;
use std::collections::{HashMap, HashSet};
use std::time::{Duration, Instant};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct CegarConfig {
    pub enable_repair: bool,
    pub seed: i32,
}

impl Default for CegarConfig {
    fn default() -> Self {
        Self {
            enable_repair: true,
            seed: 1,
        }
    }
}

/// Extracts disjoint Eulerian cycles from an active adjacency map where every vertex has degree 2.
fn extract_subcycles(nodes: &[i32], active_adj: &HashMap<i32, Vec<i32>>) -> Vec<Vec<i32>> {
    let mut visited = HashSet::new();
    let mut cycles = Vec::new();

    for &start_u in nodes {
        if visited.contains(&start_u) {
            continue;
        }
        let mut cyc = Vec::new();
        let mut curr = start_u;
        let mut prev: Option<i32> = None;

        while !visited.contains(&curr) {
            visited.insert(curr);
            cyc.push(curr);
            let nbrs = match active_adj.get(&curr) {
                Some(n) => n,
                None => break,
            };
            let next = if let Some(p) = prev {
                if nbrs.len() > 1 && nbrs[0] == p {
                    nbrs[1]
                } else if !nbrs.is_empty() {
                    nbrs[0]
                } else {
                    break;
                }
            } else if !nbrs.is_empty() {
                nbrs[0]
            } else {
                break;
            };
            prev = Some(curr);
            curr = next;
        }
        if !cyc.is_empty() {
            cycles.push(cyc);
        }
    }

    cycles
}

/// Solves Hamiltonian cycle on `g` using deterministic SAT-CEGAR with 2-opt merging and DFJ cuts.
pub fn solve_cycle(g: &Graph, timeout_secs: f64) -> Result<Vec<i32>, String> {
    solve_cycle_with_config(g, timeout_secs, CegarConfig::default())
}

pub fn solve_cycle_with_config(
    g: &Graph,
    timeout_secs: f64,
    config: CegarConfig,
) -> Result<Vec<i32>, String> {
    solve_cycle_with_forced_edges_and_config(g, timeout_secs, &HashSet::new(), config)
}

/// Solves Hamiltonian cycle on `g` requiring all edges in `forced_edges` to be selected.
pub fn solve_cycle_with_forced_edges(
    g: &Graph,
    timeout_secs: f64,
    forced_edges: &HashSet<(i32, i32)>,
) -> Result<Vec<i32>, String> {
    solve_cycle_with_forced_edges_and_config(g, timeout_secs, forced_edges, CegarConfig::default())
}

pub fn solve_cycle_with_forced_edges_and_config(
    g: &Graph,
    timeout_secs: f64,
    forced_edges: &HashSet<(i32, i32)>,
    config: CegarConfig,
) -> Result<Vec<i32>, String> {
    let deadline = Instant::now() + Duration::from_secs_f64(timeout_secs);
    let mut nodes: Vec<i32> = g.adjacency_list.keys().copied().collect();
    nodes.sort_unstable();
    let n = nodes.len();
    if n < 3 {
        return Err("Graph has fewer than 3 vertices".to_string());
    }

    // Collect deduplicated undirected edges and sorted adjacency sets
    let mut edges: Vec<(i32, i32)> = Vec::new();
    let mut adj_sorted: HashMap<i32, Vec<i32>> = HashMap::new();
    let mut adj_sets: HashMap<i32, HashSet<i32>> = HashMap::new();

    for &u in &nodes {
        let mut neighbors: Vec<i32> = g.adjacency_list[&u].clone();
        neighbors.sort_unstable();
        neighbors.dedup();
        for &w in &neighbors {
            if u < w {
                edges.push((u, w));
            }
        }
        adj_sets.insert(u, neighbors.iter().copied().collect());
        adj_sorted.insert(u, neighbors);
    }
    edges.sort_unstable();

    // Map each undirected edge to a SAT literal
    let mut var_mgr = BasicVarManager::default();
    let mut edge_vars: HashMap<(i32, i32), Lit> = HashMap::new();

    for &(u, v) in &edges {
        let lit = var_mgr.new_var().pos_lit();
        edge_vars.insert((u, v), lit);
        edge_vars.insert((v, u), lit);
    }

    let mut solver = create_solver_with_deadline(deadline);
    let _ = solver.set_option("seed", config.seed.max(1));
    let mut base_clauses: Vec<Clause> = Vec::new();

    // Add exact degree-2 constraints for each vertex
    for &u in &nodes {
        let mut inc: Vec<Lit> = Vec::new();
        if let Some(nbrs) = adj_sorted.get(&u) {
            for &v in nbrs {
                if let Some(&lit) = edge_vars.get(&(u, v)) {
                    inc.push(lit);
                }
            }
        }

        let deg = inc.len();
        if deg < 2 {
            return Err("UNSAT".to_string());
        } else if deg == 2 {
            base_clauses.push(clause![inc[0]]);
            base_clauses.push(clause![inc[1]]);
        } else {
            // At-least-1
            base_clauses.push(Clause::from_iter(inc.iter().copied()));

            // At-least-2: for each edge i, clause containing all other incident edges
            for i in 0..deg {
                let mut cl = Vec::with_capacity(deg - 1);
                for j in 0..deg {
                    if i != j {
                        cl.push(inc[j]);
                    }
                }
                base_clauses.push(Clause::from_iter(cl));
            }

            // At-most-2
            if deg <= 32 {
                for i in 0..deg {
                    for j in (i + 1)..deg {
                        for k in (j + 1)..deg {
                            let mut cl = Clause::new();
                            cl.add(!inc[i]);
                            cl.add(!inc[j]);
                            cl.add(!inc[k]);
                            base_clauses.push(cl);
                        }
                    }
                }
            } else {
                add_at_most_2(&mut solver, &mut var_mgr, &inc);
            }
        }
    }

    // Force each specified edge as a unit clause
    for &fe in forced_edges {
        let e = (fe.0.min(fe.1), fe.0.max(fe.1));
        if let Some(&lit) = edge_vars.get(&e) {
            base_clauses.push(clause![lit]);
        }
    }

    for cl in &base_clauses {
        let _ = solver.add_clause_ref(cl);
    }

    let mut forbidden_edges: HashSet<(i32, i32)> = HashSet::new();
    for &fe in forced_edges {
        forbidden_edges.insert((fe.0.min(fe.1), fe.0.max(fe.1)));
    }

    let t_start = Instant::now();
    let mut iter_count = 0;

    // CEGAR loop with DFJ subcycle cuts and 2-opt merge acceleration
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

                let mut active_adj: HashMap<i32, Vec<i32>> = HashMap::new();
                for &v in &nodes {
                    active_adj.insert(v, Vec::new());
                }

                for &(u, v) in &edges {
                    if let Some(&lit) = edge_vars.get(&(u, v)) {
                        if sol.lit_value(lit) == TernaryVal::True {
                            active_adj.entry(u).or_default().push(v);
                            active_adj.entry(v).or_default().push(u);
                        }
                    }
                }

                let cycles = extract_subcycles(&nodes, &active_adj);

                if iter_count % 10 == 0 || cycles.len() <= 8 {
                    let largest = cycles.iter().map(|c| c.len()).max().unwrap_or(0);
                    println!(
                        "  [cegar_engine] Iter {:3} ({:.2}s): {} cycles (largest={})",
                        iter_count,
                        t_start.elapsed().as_secs_f64(),
                        cycles.len(),
                        largest
                    );
                }

                if cycles.len() == 1 && cycles[0].len() == n {
                    let (ok, err) = TourVerifier::verify(g, &cycles[0]);
                    if ok {
                        return Ok(cycles[0].clone());
                    } else {
                        return Err(format!("Tour verification failed: {}", err));
                    }
                }

                // Heuristic 2-opt/3-opt merge acceleration for 2..=32 cycles
                if config.enable_repair && cycles.len() >= 2 && cycles.len() <= 32 {
                    if let Some(merged) = safe_2opt_merge(&cycles, &adj_sets, &forbidden_edges) {
                        if merged.len() == n {
                            let (ok, _) = TourVerifier::verify(g, &merged);
                            if ok {
                                println!(
                                    "  [cegar_engine] 2-opt merged {} cycles into 1 valid tour at iter {} ({:.2}s)!",
                                    cycles.len(),
                                    iter_count,
                                    t_start.elapsed().as_secs_f64()
                                );
                                return Ok(merged);
                            }
                        }
                    }
                }

                // Add DFJ cuts and cycle blocking clauses
                for cyc in &cycles {
                    if cyc.len() < n {
                        let cyc_set: HashSet<i32> = cyc.iter().copied().collect();
                        let mut cut_lits = Vec::new();
                        for &u in cyc {
                            if let Some(nbrs) = adj_sorted.get(&u) {
                                for &v in nbrs {
                                    if !cyc_set.contains(&v) {
                                        if let Some(&lit) = edge_vars.get(&(u, v)) {
                                            cut_lits.push(lit);
                                        }
                                    }
                                }
                            }
                        }
                        cut_lits.sort_unstable();
                        cut_lits.dedup();
                        if cut_lits.is_empty() {
                            return Err("UNSAT".to_string());
                        }
                        let cl = Clause::from_iter(cut_lits.clone());
                        let _ = solver.add_clause(cl);

                        if cut_lits.len() <= 10 {
                            for i in 0..cut_lits.len() {
                                let mut cl = Vec::with_capacity(cut_lits.len());
                                cl.push(!cut_lits[i]);
                                for j in 0..cut_lits.len() {
                                    if i != j {
                                        cl.push(cut_lits[j]);
                                    }
                                }
                                let _ = solver.add_clause(Clause::from_iter(cl));
                            }
                        }

                        // Cycle edge blocking clause: \bigvee_{e in C} \neg e
                        let mut block_lits = Vec::new();
                        for i in 0..cyc.len() {
                            let u = cyc[i];
                            let v = cyc[(i + 1) % cyc.len()];
                            if let Some(&lit) = edge_vars.get(&(u, v)) {
                                block_lits.push(!lit);
                            }
                        }
                        if !block_lits.is_empty() {
                            let cl = Clause::from_iter(block_lits);
                            let _ = solver.add_clause(cl);
                        }
                    }
                }
            }
        }
    }
}

/// Solves Hamiltonian path between `port_u` and `port_v` on `g` with forced edges.
pub fn solve_path_with_forced_edges(
    g: &Graph,
    port_u: i32,
    port_v: i32,
    timeout_secs: f64,
    forced_edges: &HashSet<(i32, i32)>,
) -> Result<Vec<i32>, String> {
    solve_path_with_forced_edges_and_config(
        g,
        port_u,
        port_v,
        timeout_secs,
        forced_edges,
        CegarConfig::default(),
    )
}

pub fn solve_path_with_forced_edges_and_config(
    g: &Graph,
    port_u: i32,
    port_v: i32,
    timeout_secs: f64,
    forced_edges: &HashSet<(i32, i32)>,
    config: CegarConfig,
) -> Result<Vec<i32>, String> {
    if port_u == port_v {
        return Err("Path ports must be distinct".to_string());
    }

    let max_id = g.adjacency_list.keys().copied().max().unwrap_or(0);
    let dummy_w = max_id + 1;

    let mut augmented_g = g.clone();
    augmented_g.add_edge(dummy_w, port_u);
    augmented_g.add_edge(dummy_w, port_v);

    let cycle =
        solve_cycle_with_forced_edges_and_config(&augmented_g, timeout_secs, forced_edges, config)?;

    let dummy_pos = cycle
        .iter()
        .position(|&x| x == dummy_w)
        .ok_or_else(|| "Dummy vertex missing from cycle".to_string())?;

    let cycle_len = cycle.len();
    let mut raw_path = Vec::with_capacity(cycle_len - 1);
    for i in 1..cycle_len {
        raw_path.push(cycle[(dummy_pos + i) % cycle_len]);
    }

    if raw_path.first() == Some(&port_u) && raw_path.last() == Some(&port_v) {
        Ok(raw_path)
    } else if raw_path.first() == Some(&port_v) && raw_path.last() == Some(&port_u) {
        raw_path.reverse();
        Ok(raw_path)
    } else {
        Err("Cycle orientation mismatch with specified ports".to_string())
    }
}

/// Solves Hamiltonian path between `port_u` and `port_v` on `g`.
/// Reduces the path problem to Hamiltonian cycle by adding a dummy vertex `w`
/// connected strictly to `port_u` and `port_v`.
pub fn solve_path(
    g: &Graph,
    port_u: i32,
    port_v: i32,
    timeout_secs: f64,
) -> Result<Vec<i32>, String> {
    solve_path_with_forced_edges(g, port_u, port_v, timeout_secs, &HashSet::new())
}

/// Merges two cycles if there exist vertices u1, u2 on c1 and v1, v2 on c2
/// such that deleting (u1, u2) and (v1, v2) and reconnecting forms a single cycle.
/// Never deletes edges present in `forbidden_delete`.
fn merge_two_cycles(
    c1: &[i32],
    c2: &[i32],
    adj: &HashMap<i32, HashSet<i32>>,
    forbidden_delete: &HashSet<(i32, i32)>,
) -> Option<Vec<i32>> {
    let n1 = c1.len();
    let n2 = c2.len();
    if n1 == 0 || n2 == 0 {
        return None;
    }

    let c2_pos: HashMap<i32, usize> = c2.iter().enumerate().map(|(idx, &v)| (v, idx)).collect();

    for i in 0..n1 {
        let u1 = c1[i];
        let u2 = c1[(i + 1) % n1];
        let e1 = (u1.min(u2), u1.max(u2));
        if forbidden_delete.contains(&e1) {
            continue;
        }

        let u1_nbrs = match adj.get(&u1) {
            Some(s) => s,
            None => continue,
        };
        let u2_nbrs = match adj.get(&u2) {
            Some(s) => s,
            None => continue,
        };

        for &v1 in u1_nbrs {
            if let Some(&j) = c2_pos.get(&v1) {
                // Forward orientation: cut (v1, v2) where v2 = c2[(j + 1) % n2]
                let v2 = c2[(j + 1) % n2];
                let e2 = (v1.min(v2), v1.max(v2));
                if !forbidden_delete.contains(&e2) && u2_nbrs.contains(&v2) {
                    let mut tour = Vec::with_capacity(n1 + n2);
                    for k in 1..=n1 {
                        tour.push(c1[(i + k) % n1]);
                    }
                    for k in 0..n2 {
                        let idx = (j + n2 - (k % n2)) % n2;
                        tour.push(c2[idx]);
                    }
                    return Some(tour);
                }

                // Backward orientation: cut (v0, v1) where v0 = c2[(j + n2 - 1) % n2]
                let v0 = c2[(j + n2 - 1) % n2];
                let e0 = (v0.min(v1), v0.max(v1));
                if !forbidden_delete.contains(&e0) && u2_nbrs.contains(&v0) {
                    let j_prev = (j + n2 - 1) % n2;
                    let mut tour = Vec::with_capacity(n1 + n2);
                    for k in 1..=n1 {
                        tour.push(c1[(i + k) % n1]);
                    }
                    for k in 0..n2 {
                        let idx = (j_prev + 1 + k) % n2;
                        tour.push(c2[idx]);
                    }
                    return Some(tour);
                }
            }
        }
    }

    None
}

/// 3-opt cycle merge:
/// Splits c1 at (B, C) and (D, A), with chord (C, A) in c1,
/// and splits c2 at (x, y) with cross-edges (B, x) and (y, D).
/// Inserts c2 into c1 and reverses one of c1's segments using chord (C, A).
fn merge_two_cycles_3opt(
    c1: &[i32],
    c2: &[i32],
    adj: &HashMap<i32, HashSet<i32>>,
    forbidden_delete: &HashSet<(i32, i32)>,
) -> Option<Vec<i32>> {
    let n1 = c1.len();
    let n2 = c2.len();
    if n1 < 4 || n2 == 0 {
        return None;
    }
    let pos1: HashMap<i32, usize> = c1.iter().enumerate().map(|(i, &u)| (u, i)).collect();

    for idx2 in 0..n2 {
        let x = c2[idx2];
        for &(y, rev) in &[
            (c2[(idx2 + 1) % n2], false),
            (c2[(idx2 + n2 - 1) % n2], true),
        ] {
            let e_xy = (x.min(y), x.max(y));
            if forbidden_delete.contains(&e_xy) {
                continue;
            }
            let nbrs_x: Vec<i32> = adj
                .get(&x)
                .map(|s| s.iter().filter(|v| pos1.contains_key(v)).copied().collect())
                .unwrap_or_default();
            let nbrs_y: Vec<i32> = adj
                .get(&y)
                .map(|s| s.iter().filter(|v| pos1.contains_key(v)).copied().collect())
                .unwrap_or_default();
            if nbrs_x.is_empty() || nbrs_y.is_empty() {
                continue;
            }

            for &b in &nbrs_x {
                let idx_b = pos1[&b];
                let c = c1[(idx_b + 1) % n1];
                let e_bc = (b.min(c), b.max(c));
                if forbidden_delete.contains(&e_bc) {
                    continue;
                }

                for &d in &nbrs_y {
                    if d == b || d == c {
                        continue;
                    }
                    let idx_d = pos1[&d];
                    let a = c1[(idx_d + 1) % n1];
                    if a == b || a == c || a == d {
                        continue;
                    }
                    let e_da = (d.min(a), d.max(a));
                    if forbidden_delete.contains(&e_da) {
                        continue;
                    }

                    if adj.get(&c).map_or(false, |s| s.contains(&a)) {
                        let mut p1 = Vec::new();
                        let mut curr = (idx_d + 1) % n1;
                        loop {
                            p1.push(c1[curr]);
                            if c1[curr] == b {
                                break;
                            }
                            curr = (curr + 1) % n1;
                        }

                        let mut p2 = Vec::new();
                        let mut curr2 = (idx_b + 1) % n1;
                        loop {
                            p2.push(c1[curr2]);
                            if c1[curr2] == d {
                                break;
                            }
                            curr2 = (curr2 + 1) % n1;
                        }

                        if p1.len() + p2.len() == n1 {
                            let mut q_path = Vec::with_capacity(n2);
                            if !rev {
                                for k in 0..n2 {
                                    q_path.push(c2[(idx2 + n2 - (k % n2)) % n2]);
                                }
                            } else {
                                for k in 0..n2 {
                                    q_path.push(c2[(idx2 + k) % n2]);
                                }
                            }
                            p2.reverse();
                            let mut res = p1;
                            res.extend(q_path);
                            res.extend(p2);
                            return Some(res);
                        }
                    }
                }
            }
        }
    }
    None
}

/// Attempts greedy 2-opt pairwise merge on disjoint cycles.
/// Never deletes any edge in `forbidden_delete`.
fn safe_2opt_merge(
    cycles: &[Vec<i32>],
    adj: &HashMap<i32, HashSet<i32>>,
    forbidden_delete: &HashSet<(i32, i32)>,
) -> Option<Vec<i32>> {
    let mut curr = cycles.to_vec();
    let mut merged_any = true;

    while merged_any && curr.len() > 1 {
        merged_any = false;
        curr.sort_by(|a, b| b.len().cmp(&a.len()));

        let num_cycles = curr.len();
        let mut merge_step = None;

        'search: for i in 0..num_cycles {
            for j in (i + 1)..num_cycles {
                if let Some(res) = merge_two_cycles(&curr[i], &curr[j], adj, forbidden_delete) {
                    merge_step = Some((i, j, res));
                    break 'search;
                }
            }
        }

        if let Some((i, j, res)) = merge_step {
            curr.remove(j);
            curr[i] = res;
            merged_any = true;
        }
    }

    let mut merged_3opt = true;
    while merged_3opt && curr.len() > 1 && curr.len() <= 6 {
        merged_3opt = false;
        curr.sort_by(|a, b| b.len().cmp(&a.len()));
        let num_cycles = curr.len();
        let mut merge_step = None;
        'search_3opt: for i in 0..num_cycles {
            for j in (i + 1)..num_cycles {
                if let Some(res) = merge_two_cycles_3opt(&curr[i], &curr[j], adj, forbidden_delete)
                {
                    merge_step = Some((i, j, res));
                    break 'search_3opt;
                }
                if let Some(res) = merge_two_cycles_3opt(&curr[j], &curr[i], adj, forbidden_delete)
                {
                    merge_step = Some((j, i, res));
                    break 'search_3opt;
                }
            }
        }
        if let Some((i, j, res)) = merge_step {
            let max_idx = i.max(j);
            let min_idx = i.min(j);
            curr.remove(max_idx);
            curr[min_idx] = res;
            merged_3opt = true;

            let mut merged_2opt = true;
            while merged_2opt && curr.len() > 1 {
                merged_2opt = false;
                curr.sort_by(|a, b| b.len().cmp(&a.len()));
                let n_c = curr.len();
                let mut step_2opt = None;
                'search_inner: for a in 0..n_c {
                    for b in (a + 1)..n_c {
                        if let Some(r) = merge_two_cycles(&curr[a], &curr[b], adj, forbidden_delete)
                        {
                            step_2opt = Some((a, b, r));
                            break 'search_inner;
                        }
                    }
                }
                if let Some((a, b, r)) = step_2opt {
                    curr.remove(b);
                    curr[a] = r;
                    merged_2opt = true;
                }
            }
        }
    }

    if curr.len() == 1 {
        Some(curr.into_iter().next().unwrap())
    } else {
        None
    }
}

/// Sequential counter encoding for at-most-2 cardinality constraint.
pub fn add_at_most_2(solver: &mut CaDiCaL, var_mgr: &mut BasicVarManager, edge_lits: &[Lit]) {
    let n = edge_lits.len();
    if n <= 2 {
        return;
    }
    if n <= 32 {
        for i in 0..n {
            for j in (i + 1)..n {
                for k in (j + 1)..n {
                    let _ = solver.add_clause(clause![!edge_lits[i], !edge_lits[j], !edge_lits[k]]);
                }
            }
        }
        return;
    }

    let mut s: Vec<Vec<Lit>> = Vec::with_capacity(n - 1);
    for _ in 0..(n - 1) {
        let s0 = var_mgr.new_var().pos_lit();
        let s1 = var_mgr.new_var().pos_lit();
        s.push(vec![s0, s1]);
    }

    // Base clauses (i = 0)
    let _ = solver.add_clause(clause![!edge_lits[0], s[0][0]]);
    let _ = solver.add_clause(clause![!s[0][1]]);

    // Step clauses (i = 1..n-1)
    for i in 1..(n - 1) {
        let _ = solver.add_clause(clause![!s[i - 1][0], s[i][0]]);
        let _ = solver.add_clause(clause![!s[i - 1][1], s[i][1]]);

        let _ = solver.add_clause(clause![!edge_lits[i], s[i][0]]);
        let _ = solver.add_clause(clause![!edge_lits[i], !s[i - 1][0], s[i][1]]);
        let _ = solver.add_clause(clause![!edge_lits[i], !s[i - 1][1]]);
    }

    // Final clause for last literal
    let _ = solver.add_clause(clause![!edge_lits[n - 1], !s[n - 2][1]]);
}
