use crate::pipeline::ablation::AblationConfig;
use clap::{App, Arg};

#[derive(Debug, Clone, PartialEq)]
pub struct Options {
    pub batch_mode: bool,
    pub batch_start: usize,
    pub batch_end: usize,
    pub workers: usize,
    pub checkpoint: String,
    pub graph_file: String,
    pub timeout: f64,
    pub output_tour_file: Option<String>,
    pub ablation: AblationConfig,
    pub seed: i32,
}

impl Default for Options {
    fn default() -> Self {
        Self {
            batch_mode: false,
            batch_start: 0,
            batch_end: 0,
            workers: 1,
            checkpoint: "scratch/batch_1001_results.json".to_string(),
            graph_file: String::new(),
            timeout: 1800.0,
            output_tour_file: None,
            ablation: AblationConfig::Full,
            seed: 1,
        }
    }
}

impl Options {
    pub fn parse_from_args() -> Self {
        Self::try_parse_from_args().unwrap_or_else(|e| {
            eprintln!("Error: {}", e);
            std::process::exit(1);
        })
    }

    pub fn try_parse_from_args() -> Result<Self, String> {
        let matches = get_options();
        Self::try_from_matches(&matches)
    }

    pub fn from_matches(matches: &clap::ArgMatches) -> Self {
        Self::try_from_matches(matches).unwrap_or_else(|e| {
            eprintln!("Error: {}", e);
            std::process::exit(1);
        })
    }

    pub fn try_from_matches(matches: &clap::ArgMatches) -> Result<Self, String> {
        let timeout = matches.value_of_t::<f64>("timeout").unwrap_or(1800.0);
        let output_tour_file = matches.value_of("output-tour").map(|s| s.to_string());
        let ablation = matches
            .value_of("ablation")
            .unwrap_or("full")
            .parse::<AblationConfig>()?;
        let seed = matches
            .value_of("seed")
            .unwrap_or("1")
            .parse::<i32>()
            .map_err(|_| "Invalid --seed; expected a positive 32-bit integer".to_string())?;
        if seed <= 0 {
            return Err("Invalid --seed; expected a positive 32-bit integer".to_string());
        }

        // Check if batch mode is requested
        if let Some(mut batch_vals) = matches.values_of("batch") {
            let start_str = batch_vals
                .next()
                .ok_or_else(|| "Missing START for --batch".to_string())?;
            let end_str = batch_vals.next().unwrap_or(start_str);
            let start: usize = start_str
                .parse()
                .map_err(|_| "Invalid START graph ID for --batch".to_string())?;
            let end: usize = end_str
                .parse()
                .map_err(|_| "Invalid END graph ID for --batch".to_string())?;
            let workers: usize = matches
                .value_of("workers")
                .and_then(|w| w.parse().ok())
                .unwrap_or(2);
            let checkpoint = matches
                .value_of("checkpoint")
                .unwrap_or("scratch/batch_1001_results.json")
                .to_string();

            return Ok(Options {
                batch_mode: true,
                batch_start: start,
                batch_end: end,
                workers,
                checkpoint,
                graph_file: String::new(),
                timeout,
                output_tour_file,
                ablation,
                seed,
            });
        }

        let input_filename = matches
            .value_of("input")
            .or_else(|| matches.value_of("positional_input"));
        let graph_file = match input_filename {
            Some(f) => f.to_string(),
            None => {
                return Err(
                    "No input graph specified. Provide -i <FILE> or use --batch <START> <END>."
                        .to_string(),
                );
            }
        };

        Ok(Options {
            batch_mode: false,
            batch_start: 0,
            batch_end: 0,
            workers: 1,
            checkpoint: String::new(),
            graph_file,
            timeout,
            output_tour_file,
            ablation,
            seed,
        })
    }
}

pub fn get_options() -> clap::ArgMatches {
    get_options_app().get_matches()
}

pub fn get_options_app() -> clap::App<'static> {
    App::new("HCP Solver")
        .version("1.0")
        .author("Me <me@example.com>")
        .about("Solves Hamiltonian cycles in graphs")
        .arg(
            Arg::with_name("solver")
                .short('s')
                .long("solver")
                .value_name("n")
                .help("Solver: 0: minisat (default), 1: kissat, 2: cadical")
                .takes_value(true),
        )
        .arg(
            Arg::with_name("encoding")
                .short('e')
                .long("encoding")
                .value_name("n")
                .help("Encoding method: 0: binominal (default), 1: sinz, 2: adder, 3: advanced sinz, 4: product + binominal, 5: product recursive, 6: ladder")
                .takes_value(true),
        )
        .arg(
            Arg::with_name("input")
                .short('i')
                .long("input")
                .value_name("FILE NAME")
                .help("Input file")
                .takes_value(true),
        )
        .arg(
            Arg::with_name("positional_input")
                .value_name("INPUT FILE")
                .help("Input file path (positional)")
                .takes_value(true)
                .index(1),
        )
        .arg(
            Arg::with_name("output-tour")
                .short('o')
                .long("output-tour")
                .alias("output")
                .value_name("FILE")
                .help("Output certified HCP tour file path (TSPLIB format)")
                .takes_value(true),
        )
        .arg(
            Arg::with_name("blocking")
                .short('b')
                .long("block")
                .value_name("n")
                .help("Blocking method")
                .takes_value(true),
        )
        .arg(
            Arg::with_name("symmetry")
                .short('y')
                .long("symmetry")
                .value_name("n")
                .help("Symmetry blocking method")
                .takes_value(true),
        )
        .arg(
            Arg::with_name("2-opt")
                .long("two-opt")
                .value_name("n")
                .help("2-opt method")
                .takes_value(true),
        )
        .arg(
            Arg::with_name("three-opt")
                .short('x')
                .long("three-opt")
                .value_name("n")
                .help("Restricted 3-opt method:\n 0: Disabled (default)\n 1: Enabled")
                .takes_value(true),
        )
        .arg(
            Arg::with_name("loop-prohibition")
                .short('l')
                .long("loop")
                .value_name("n")
                .help("Loop prohibition")
                .takes_value(true),
        )
        .arg(
            Arg::with_name("cnf-normalize")
                .short('n')
                .long("normalize")
                .value_name("n")
                .help("CNF normalization")
                .takes_value(true),
        )
        .arg(
            Arg::with_name("balanced")
                .short('c')
                .long("balanced")
                .value_name("n")
                .help("Block clauses balanced")
                .takes_value(true),
        )
        .arg(
            Arg::with_name("de-arcify")
                .short('d')
                .long("de-arcify")
                .value_name("n")
                .help("Remove redundant arcs")
                .takes_value(true),
        )
        .arg(
            Arg::with_name("set-configration")
                .long("set-configration")
                .value_name("n")
                .help("cadical set configuration")
                .takes_value(true),
        )
        .arg(
            Arg::with_name("degree-order")
                .short('r')
                .long("degree_order")
                .value_name("n")
                .help("clauses order")
                .takes_value(true),
        )
        .arg(
            Arg::with_name("arcs-order")
                .short('a')
                .long("arc_order")
                .value_name("n")
                .help("literal number order")
                .takes_value(true),
        )
        .arg(
            Arg::with_name("cegar-fallback")
                .short('f')
                .long("cegar-fallback")
                .value_name("n")
                .help("CEGAR hard blocking fallback option")
                .takes_value(true),
        )
        .arg(
            Arg::with_name("mtz-stall")
                .long("mtz-stall")
                .value_name("n")
                .help("Partial MTZ injection stall threshold")
                .takes_value(true),
        )
        .arg(
            Arg::with_name("adaptive-escalation")
                .short('A')
                .long("adaptive-escalation")
                .value_name("n")
                .help("Adaptive stall-based escalation strategy")
                .takes_value(true),
        )
        .arg(
            Arg::with_name("sub-hcp-timeout")
                .long("sub-hcp-timeout")
                .value_name("n")
                .help("Sub-HCP solver timeout per cluster in seconds")
                .takes_value(true),
        )
        .arg(
            Arg::with_name("max-cluster-size")
                .long("max-cluster-size")
                .value_name("n")
                .help("Maximum vertices per cluster for sub-HCP solving")
                .takes_value(true),
        )
        .arg(
            Arg::with_name("two-tier")
                .long("two-tier")
                .value_name("n")
                .help("Two-tier demand-coordinated solver")
                .takes_value(true)
                .min_values(0),
        )
        .arg(
            Arg::with_name("staged-smt")
                .long("staged-smt")
                .value_name("n")
                .help("Staged-length lazy SMT solver")
                .takes_value(true)
                .min_values(0),
        )
        .arg(
            Arg::with_name("timeout")
                .short('t')
                .long("timeout")
                .value_name("SECONDS")
                .help("Timeout in seconds (default: 1800.0)")
                .takes_value(true),
        )
        .arg(
            Arg::with_name("batch")
                .long("batch")
                .value_names(&["START", "END"])
                .number_of_values(2)
                .help("Run batch solver on graph ID range [START, END]")
                .takes_value(true),
        )
        .arg(
            Arg::with_name("workers")
                .long("workers")
                .value_name("N")
                .help("Number of worker threads for batch mode (default: 2 or available CPUs)")
                .takes_value(true),
        )
        .arg(
            Arg::with_name("checkpoint")
                .long("checkpoint")
                .value_name("PATH")
                .help("JSON checkpoint file path (default: scratch/batch_1001_results.json)")
                .takes_value(true)
                .default_value("scratch/batch_1001_results.json"),
        )
        .arg(
            Arg::with_name("auto")
                .long("auto")
                .value_name("n")
                .help("Auto topology classification and hybrid solver routing:\n 0: Disabled\n 1: Enabled (default)")
                .takes_value(true)
                .default_value("1"),
        )
        .arg(
            Arg::with_name("macro-gadget")
                .long("macro-gadget")
                .value_name("n")
                .help("Macro-gadget state encoding (H1)")
                .takes_value(true)
                .min_values(0),
        )
        .arg(
            Arg::with_name("bounded-freezer")
                .long("bounded-freezer")
                .value_name("n")
                .help("Topologically bounded backbone freezing (H2-Refined)")
                .takes_value(true)
                .min_values(0),
        )
        .arg(
            Arg::with_name("ablation")
                .long("ablation")
                .value_name("NAME")
                .help("Ablation: full, no-decomposition, no-dispatch, no-series, no-two-cut, no-repair, undirected-cubic, or one-alternating-worker")
                .takes_value(true)
                .default_value("full"),
        )
        .arg(
            Arg::with_name("seed")
                .long("seed")
                .value_name("N")
                .help("Positive CaDiCaL random seed (default: 1)")
                .takes_value(true)
                .default_value("1"),
        )
        .arg(
            Arg::with_name("alternating-engine")
                .long("alternating-engine")
                .value_name("n")
                .help("Alternating Port Engine for cycle compression")
                .takes_value(true)
                .min_values(0),
        )
        .arg(
            Arg::with_name("no-alternating-engine")
                .long("no-alternating-engine")
                .help("Disable Alternating Port Engine")
                .takes_value(false),
        )
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_options_missing_input_returns_err() {
        let app = get_options_app();
        let matches = app.try_get_matches_from(vec!["cegar-fix"]).unwrap();
        let res = Options::try_from_matches(&matches);
        assert!(res.is_err());
        assert!(res.unwrap_err().contains("No input graph specified"));
    }

    #[test]
    fn test_options_single_graph_success() {
        let app = get_options_app();
        let matches = app
            .try_get_matches_from(vec!["cegar-fix", "-i", "test.col", "-t", "15"])
            .unwrap();
        let res = Options::try_from_matches(&matches);
        assert!(res.is_ok());
        let opts = res.unwrap();
        assert_eq!(opts.graph_file, "test.col");
        assert_eq!(opts.timeout, 15.0);
        assert!(!opts.batch_mode);
        assert_eq!(opts.ablation, AblationConfig::Full);
        assert_eq!(opts.seed, 1);
    }

    #[test]
    fn test_options_batch_mode_success() {
        let app = get_options_app();
        let matches = app
            .try_get_matches_from(vec!["cegar-fix", "--batch", "1", "10", "--workers", "4"])
            .unwrap();
        let res = Options::try_from_matches(&matches);
        assert!(res.is_ok());
        let opts = res.unwrap();
        assert!(opts.batch_mode);
        assert_eq!(opts.batch_start, 1);
        assert_eq!(opts.batch_end, 10);
        assert_eq!(opts.workers, 4);
    }

    #[test]
    fn test_options_named_ablation_and_seed() {
        let app = get_options_app();
        let matches = app
            .try_get_matches_from(vec![
                "cegar-fix",
                "-i",
                "test.col",
                "--ablation",
                "Proposed-NoRepair",
                "--seed",
                "42",
            ])
            .unwrap();
        let opts = Options::try_from_matches(&matches).unwrap();
        assert_eq!(opts.ablation, AblationConfig::NoRepair);
        assert_eq!(opts.seed, 42);
    }
}
