use cegar_fix::pipeline::options::Options;
use cegar_fix::pipeline::solver_pipeline::{run_batch, solve_single_graph_with_config};
use std::time::Instant;

fn main() {
    env_logger::init();
    let instant = Instant::now();
    let opts = Options::parse_from_args();

    if opts.batch_mode {
        println!(
            "Running batch solver for graphs {}..={} with {} workers (timeout: {:.1}s, ablation: {}, seed: {})...",
            opts.batch_start,
            opts.batch_end,
            opts.workers,
            opts.timeout,
            opts.ablation,
            opts.seed,
        );
        run_batch(&opts);
        return;
    }

    println!("solve {}", opts.graph_file);
    match solve_single_graph_with_config(
        &opts.graph_file,
        opts.timeout,
        opts.output_tour_file.as_deref(),
        opts.ablation,
        opts.seed,
    ) {
        Ok((tour, elapsed, _vertices)) => {
            println!("s SATISFIABLE");
            println!("solution: ");
            for v in &tour {
                print!("{} ", v);
            }
            println!();
            println!("overall time = {:.6}s", elapsed);
        }
        Err(e) => {
            if e.contains("UNSAT") || e.contains("Infeasible") || e.contains("cut-vertex") {
                println!("s UNSATISFIABLE");
            } else {
                println!("s UNKNOWN");
            }
            eprintln!("Solver output/error: {}", e);
            println!("overall time = {:?}", instant.elapsed());
        }
    }
}
