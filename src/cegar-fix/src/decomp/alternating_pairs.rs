use crate::core::graph::Graph;
use crate::core::tour_verifier::TourVerifier;
use crate::decomp::spqr_series::{contract_series_chains, expand_series_tour};
use rustsat::clause;
use rustsat::instances::{BasicVarManager, Cnf, ManageVars};
use rustsat::solvers::{ControlSignal, Solve, SolverResult, Terminate};
use rustsat::types::{Clause, Lit};
use rustsat_cadical::CaDiCaL;
use std::collections::{HashMap, HashSet};
use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::{mpsc, Arc};
use std::thread;
use std::time::{Duration, Instant};

enum WorkerCmd {
    Solve,
    AddCnf(Cnf),
    Stop,
}

enum WorkerMsg {
    Solution(usize, Vec<Lit>),
    #[allow(dead_code)]
    Unsat,
    Cancelled,
}

fn is_valid_dir_cycle(cyc: &[usize], dir_adj: &[Vec<usize>]) -> bool {
    let n = cyc.len();
    if n < 3 {
        return false;
    }
    for i in 0..n {
        let u = cyc[i];
        let v = cyc[(i + 1) % n];
        if !dir_adj[u].contains(&v) {
            return false;
        }
    }
    true
}

fn merge_two_dir_cycles(c1: &[usize], c2: &[usize], dir_adj: &[Vec<usize>]) -> Option<Vec<usize>> {
    let n1 = c1.len();
    let n2 = c2.len();
    if n1 == 0 || n2 == 0 {
        return None;
    }
    let c2_pos: HashMap<usize, usize> = c2.iter().enumerate().map(|(idx, &v)| (v, idx)).collect();

    for i in 0..n1 {
        let u1 = c1[i];
        let u2 = c1[(i + 1) % n1];
        for &v1 in &dir_adj[u1] {
            if let Some(&j) = c2_pos.get(&v1) {
                let j_prev = (j + n2 - 1) % n2;
                let v_prev = c2[j_prev];
                if dir_adj[v_prev].contains(&u2) {
                    let mut merged = Vec::with_capacity(n1 + n2);
                    for k in 0..=i {
                        merged.push(c1[k]);
                    }
                    for k in 0..n2 {
                        merged.push(c2[(j + k) % n2]);
                    }
                    for k in (i + 1)..n1 {
                        merged.push(c1[k]);
                    }
                    if is_valid_dir_cycle(&merged, dir_adj) {
                        return Some(merged);
                    }
                }
            }
        }
    }
    None
}

fn pairwise_dir_merge(cycles: &[Vec<usize>], dir_adj: &[Vec<usize>]) -> Vec<Vec<usize>> {
    let mut curr = cycles.to_vec();
    let mut merged_any = true;
    while merged_any && curr.len() > 1 {
        merged_any = false;
        curr.sort_by(|a, b| b.len().cmp(&a.len()));
        let n_cycles = curr.len();
        let mut merge_step = None;
        'search: for i in 0..n_cycles {
            for j in (i + 1)..n_cycles {
                if let Some(res) = merge_two_dir_cycles(&curr[i], &curr[j], dir_adj) {
                    merge_step = Some((i, j, res));
                    break 'search;
                }
                if let Some(res) = merge_two_dir_cycles(&curr[j], &curr[i], dir_adj) {
                    merge_step = Some((j, i, res));
                    break 'search;
                }
            }
        }
        if let Some((i, j, res)) = merge_step {
            let max_idx = i.max(j);
            let min_idx = i.min(j);
            curr.remove(max_idx);
            curr[min_idx] = res;
            merged_any = true;
        }
    }
    curr
}

struct AlternatingPairGraph {
    contracted_g: Graph,
    chain_map: HashMap<(i32, i32), Vec<i32>>,
    pairs: Vec<(i32, i32)>,
    v_partner: HashMap<i32, i32>,
    color: HashMap<i32, u8>,
}

const MIN_DEG2_FRACTION: f64 = 0.15;

fn extract_alternating_pairs(raw_g: &Graph) -> Option<AlternatingPairGraph> {
    let deg2_count = raw_g
        .adjacency_list
        .values()
        .filter(|nbrs| nbrs.len() == 2)
        .count();
    let n = raw_g.adjacency_list.len();
    let deg2_fraction = deg2_count as f64 / n as f64;
    if deg2_fraction < MIN_DEG2_FRACTION {
        return None;
    }

    let series_decomp = contract_series_chains(raw_g);
    let g = series_decomp.contracted_g;
    let chain_map = series_decomp.chain_map;

    let mut v_partner: HashMap<i32, i32> = HashMap::new();
    let mut pairs: Vec<(i32, i32)> = Vec::new();
    let mut sorted_chains: Vec<(i32, i32)> = chain_map.keys().copied().collect();
    sorted_chains.sort_unstable();

    for (u, w) in sorted_chains {
        if u < w {
            pairs.push((u, w));
            v_partner.insert(u, w);
            v_partner.insert(w, u);
        }
    }

    if pairs.is_empty() {
        return None;
    }

    let mut color: HashMap<i32, u8> = HashMap::new();
    let (u0, w0) = pairs[0];
    color.insert(u0, 0);
    color.insert(w0, 1);
    let mut q = vec![u0, w0];

    while let Some(curr) = q.pop() {
        let curr_c = color[&curr];

        if let Some(&vp) = v_partner.get(&curr) {
            if let Some(&existing_color) = color.get(&vp) {
                if existing_color == curr_c {
                    return None;
                }
            } else {
                color.insert(vp, 1 - curr_c);
                q.push(vp);
            }
        }

        if let Some(nbrs) = g.adjacency_list.get(&curr) {
            let vp = v_partner.get(&curr).copied().unwrap_or(-1);
            for &nxt in nbrs {
                if nxt != vp {
                    if let Some(&existing_color) = color.get(&nxt) {
                        if existing_color == curr_c {
                            return None;
                        }
                    } else {
                        color.insert(nxt, 1 - curr_c);
                        q.push(nxt);
                    }
                }
            }
        }
    }

    if g.adjacency_list.len() != pairs.len() * 2 || color.len() != g.adjacency_list.len() {
        return None;
    }

    Some(AlternatingPairGraph {
        contracted_g: g,
        chain_map,
        pairs,
        v_partner,
        color,
    })
}

pub fn can_solve_alternating_pairs(raw_g: &Graph) -> bool {
    extract_alternating_pairs(raw_g).is_some()
}

pub fn solve_alternating_pairs(raw_g: &Graph, timeout_secs: f64) -> Option<Vec<i32>> {
    solve_alternating_pairs_with_config(raw_g, timeout_secs, 3, true, 1)
}

pub fn solve_alternating_pairs_with_config(
    raw_g: &Graph,
    timeout_secs: f64,
    num_workers: usize,
    enable_repair: bool,
    seed: i32,
) -> Option<Vec<i32>> {
    let t_start = Instant::now();
    let deadline = t_start + Duration::from_secs_f64(timeout_secs);

    let apg = extract_alternating_pairs(raw_g)?;

    let AlternatingPairGraph {
        contracted_g: g,
        chain_map,
        pairs,
        v_partner,
        color,
    } = apg;

    let n_dir = pairs.len();
    let mut node_to_id: HashMap<i32, usize> = HashMap::new();
    let mut id_to_pair: HashMap<usize, (i32, i32)> = HashMap::new();

    for (idx, &(u, w)) in pairs.iter().enumerate() {
        node_to_id.insert(u, idx);
        node_to_id.insert(w, idx);
        let u_c = color.get(&u).copied().unwrap_or(0);
        let in_v = if u_c == 0 { u } else { w };
        let out_v = if u_c == 1 { u } else { w };
        id_to_pair.insert(idx, (in_v, out_v));
    }

    let mut dir_adj: Vec<Vec<usize>> = vec![Vec::new(); n_dir];
    let mut arc_set: HashSet<(usize, usize)> = HashSet::new();

    for (&u, &c) in &color {
        if c == 1 {
            let u_id = node_to_id[&u];
            if let Some(nbrs) = g.adjacency_list.get(&u) {
                let mut sorted_nbrs = nbrs.clone();
                sorted_nbrs.sort_unstable();
                let vp = v_partner.get(&u).copied().unwrap_or(0);
                for v in sorted_nbrs {
                    if v != vp {
                        if let Some(&v_id) = node_to_id.get(&v) {
                            if arc_set.insert((u_id, v_id)) {
                                dir_adj[u_id].push(v_id);
                            }
                        }
                    }
                }
            }
        }
    }

    let mut var_mgr = BasicVarManager::default();
    let mut arc_lit_map: HashMap<(usize, usize), Lit> = HashMap::new();
    for &arc in &arc_set {
        arc_lit_map.insert(arc, var_mgr.new_var().pos_lit());
    }

    let mut in_arcs: Vec<Vec<usize>> = vec![Vec::new(); n_dir];
    for &(u, v) in &arc_set {
        in_arcs[v].push(u);
    }

    let mut cnf = Cnf::new();
    for u in 0..n_dir {
        let out_lits: Vec<Lit> = dir_adj[u].iter().map(|&v| arc_lit_map[&(u, v)]).collect();
        cnf.add_clause(Clause::from_iter(out_lits.iter().copied()));
        for i in 0..out_lits.len() {
            for j in (i + 1)..out_lits.len() {
                cnf.add_clause(clause!(!out_lits[i], !out_lits[j]));
            }
        }
    }
    for v in 0..n_dir {
        let in_lits: Vec<Lit> = in_arcs[v].iter().map(|&u| arc_lit_map[&(u, v)]).collect();
        cnf.add_clause(Clause::from_iter(in_lits.iter().copied()));
        for i in 0..in_lits.len() {
            for j in (i + 1)..in_lits.len() {
                cnf.add_clause(clause!(!in_lits[i], !in_lits[j]));
            }
        }
    }
    for &(u, v) in &arc_set {
        if u < v && arc_set.contains(&(v, u)) {
            cnf.add_clause(clause!(!arc_lit_map[&(u, v)], !arc_lit_map[&(v, u)]));
        }
    }

    // Static short cycle cuts
    let mut static_3c = Vec::new();
    for u in 0..n_dir {
        for &v in &dir_adj[u] {
            for &w in &dir_adj[v] {
                if w != u && dir_adj[w].contains(&u) {
                    if u < v && u < w {
                        static_3c.push(vec![u, v, w]);
                    }
                }
            }
        }
    }
    for c in &static_3c {
        let lits = vec![
            !arc_lit_map[&(c[0], c[1])],
            !arc_lit_map[&(c[1], c[2])],
            !arc_lit_map[&(c[2], c[0])],
        ];
        cnf.add_clause(Clause::from_iter(lits));
    }

    let num_workers = num_workers.clamp(1, 3);
    let (tx_res, rx_res) = mpsc::channel::<WorkerMsg>();
    let mut worker_senders: Vec<mpsc::Sender<WorkerCmd>> = Vec::new();
    let mut cancel_flags: Vec<Arc<AtomicBool>> = Vec::new();
    let mut worker_handles = Vec::new();

    for worker_id in 0..num_workers {
        let (tx_cmd, rx_cmd) = mpsc::channel::<WorkerCmd>();
        worker_senders.push(tx_cmd);
        let cancel_flag = Arc::new(AtomicBool::new(false));
        cancel_flags.push(cancel_flag.clone());

        let tx_res_clone = tx_res.clone();
        let cnf_clone = cnf.clone();
        let worker_seed = seed
            .max(1)
            .saturating_add((worker_id as i32).saturating_mul(104_729));

        let handle = thread::spawn(move || {
            let mut solver = CaDiCaL::default();
            let _ = solver.set_option("chrono", 1);
            let _ = solver.set_option("seed", worker_seed);
            let restart_interval = [50, 100, 200][worker_id];
            let _ = solver.set_option("restartint", restart_interval);

            let _ = solver.add_cnf(cnf_clone.clone());

            let cancel_ref = cancel_flag.clone();
            solver.attach_terminator(move || {
                if cancel_ref.load(Ordering::Relaxed) || Instant::now() >= deadline {
                    ControlSignal::Terminate
                } else {
                    ControlSignal::Continue
                }
            });

            while let Ok(cmd) = rx_cmd.recv() {
                match cmd {
                    WorkerCmd::Solve => {
                        cancel_flag.store(false, Ordering::SeqCst);
                        let res = solver.solve();
                        if cancel_flag.load(Ordering::Relaxed) {
                            let _ = tx_res_clone.send(WorkerMsg::Cancelled);
                        } else if matches!(res, Ok(SolverResult::Sat)) {
                            let sol = solver.full_solution().unwrap();
                            let _ = tx_res_clone
                                .send(WorkerMsg::Solution(worker_id, sol.into_iter().collect()));
                        } else if matches!(res, Ok(SolverResult::Unsat)) {
                            let _ = tx_res_clone.send(WorkerMsg::Unsat);
                        } else {
                            let _ = tx_res_clone.send(WorkerMsg::Cancelled);
                        }
                    }
                    WorkerCmd::AddCnf(cuts) => {
                        let _ = solver.add_cnf(cuts);
                    }
                    WorkerCmd::Stop => break,
                }
            }
        });
        worker_handles.push(handle);
    }

    let mut accumulated_cuts = Cnf::new();
    let mut result_tour = None;

    for round in 1..=500 {
        let t_round = Instant::now();

        if t_start.elapsed().as_secs_f64() > timeout_secs {
            break;
        }

        for flag in &cancel_flags {
            flag.store(false, Ordering::SeqCst);
        }
        for tx in &worker_senders {
            let _ = tx.send(WorkerCmd::Solve);
        }

        let mut winning_sol = None;
        let mut finished_count = 0;

        while finished_count < num_workers {
            match rx_res.recv() {
                Ok(WorkerMsg::Solution(w_id, sol)) => {
                    if winning_sol.is_none() {
                        winning_sol = Some(sol);
                        for (idx, flag) in cancel_flags.iter().enumerate() {
                            if idx != w_id {
                                flag.store(true, Ordering::SeqCst);
                            }
                        }
                    }
                    finished_count += 1;
                }
                Ok(WorkerMsg::Unsat) | Ok(WorkerMsg::Cancelled) => {
                    finished_count += 1;
                }
                Err(_) => break,
            }
        }

        if winning_sol.is_none() {
            break;
        }

        let sol = winning_sol.unwrap();
        let sol_set: HashSet<Lit> = sol.into_iter().collect();

        let mut successor: HashMap<usize, usize> = HashMap::new();
        for (&(u, v), &lit) in &arc_lit_map {
            if sol_set.contains(&lit) {
                successor.insert(u, v);
            }
        }

        let mut visited = HashSet::new();
        let mut cycles = Vec::new();
        for start in 0..n_dir {
            if !visited.contains(&start) {
                let mut c = Vec::new();
                let mut curr = start;
                while !visited.contains(&curr) {
                    visited.insert(curr);
                    c.push(curr);
                    if let Some(&nxt) = successor.get(&curr) {
                        curr = nxt;
                    } else {
                        break;
                    }
                }
                if !c.is_empty() {
                    cycles.push(c);
                }
            }
        }

        cycles.sort_by_key(|c| std::cmp::Reverse(c.len()));
        let raw_cycles = cycles.clone();

        if enable_repair {
            cycles = pairwise_dir_merge(&cycles, &dir_adj);
        }

        let lens: Vec<usize> = cycles.iter().map(|c| c.len()).collect();
        let top5: Vec<usize> = lens.iter().take(5).copied().collect();
        println!(
            "[alternating_pairs] Round {:2} ({:?}): {} cycles (largest={}, top5={:?})",
            round,
            t_round.elapsed(),
            cycles.len(),
            cycles[0].len(),
            top5
        );

        if cycles.len() == 1 || cycles[0].len() == n_dir {
            let final_dir_tour = &cycles[0];
            let mut contracted_tour = Vec::with_capacity(n_dir * 2);
            for &dir_v in final_dir_tour {
                let (in_v, out_v) = id_to_pair[&dir_v];
                contracted_tour.push(in_v);
                contracted_tour.push(out_v);
            }
            let full_tour = expand_series_tour(&contracted_tour, &chain_map);
            let (valid, _err) = TourVerifier::verify(raw_g, &full_tour);
            if valid {
                result_tour = Some(full_tour);
                break;
            }
        }

        let mut cuts = Cnf::new();
        for c in &raw_cycles {
            if c.len() < n_dir {
                let mut lits = Vec::new();
                for i in 0..c.len() {
                    let u = c[i];
                    let v = c[(i + 1) % c.len()];
                    if let Some(&lit) = arc_lit_map.get(&(u, v)) {
                        lits.push(!lit);
                    }
                }
                if !lits.is_empty() {
                    cuts.add_clause(Clause::from_iter(lits));
                }

                let c_set: HashSet<usize> = c.iter().copied().collect();
                let mut out_lits = Vec::new();
                for &u in c {
                    for &v in &dir_adj[u] {
                        if !c_set.contains(&v) {
                            if let Some(&lit) = arc_lit_map.get(&(u, v)) {
                                out_lits.push(lit);
                            }
                        }
                    }
                }
                out_lits.sort_unstable();
                out_lits.dedup();
                if !out_lits.is_empty() {
                    cuts.add_clause(Clause::from_iter(out_lits.iter().copied()));
                }

                let mut in_lits = Vec::new();
                for &v in c {
                    for &u in &in_arcs[v] {
                        if !c_set.contains(&u) {
                            if let Some(&lit) = arc_lit_map.get(&(u, v)) {
                                in_lits.push(lit);
                            }
                        }
                    }
                }
                in_lits.sort_unstable();
                in_lits.dedup();
                if !in_lits.is_empty() {
                    cuts.add_clause(Clause::from_iter(in_lits.iter().copied()));
                }
            }
        }

        for cl in cuts.iter() {
            accumulated_cuts.add_clause(cl.clone());
        }

        for tx in &worker_senders {
            let _ = tx.send(WorkerCmd::AddCnf(cuts.clone()));
        }
    }

    for tx in &worker_senders {
        let _ = tx.send(WorkerCmd::Stop);
    }
    for handle in worker_handles {
        let _ = handle.join();
    }

    result_tour
}
