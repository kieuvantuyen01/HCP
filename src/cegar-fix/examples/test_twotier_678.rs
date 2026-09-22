use std::time::Instant;
use cegar_fix::file_operations;
use cegar_fix::two_tier_orchestrator::{TwoTierOrchestrator, TwoTierOptions};

fn main() {
    let raw_g = file_operations::input_to_graph("FHCPCS-col/graph678.col");
    let options = TwoTierOptions {
        timeout_secs: 10.0,
        output_tour: None,
        cycle_repair: true,
    };
    let t0 = Instant::now();
    let tour_opt = TwoTierOrchestrator::solve(&raw_g, &options);
    println!("TwoTier solve result in {:?}: is_some={}", t0.elapsed(), tour_opt.is_some());
}
