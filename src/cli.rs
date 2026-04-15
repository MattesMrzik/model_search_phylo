use std::fmt::{self, Display};
use std::path::PathBuf;

use clap::Parser;
use ftail::Ftail;
use log::LevelFilter;
use ntimestamp::Timestamp;

use phylo::evolutionary_models::FrequencyOptimisation;
use phylo::optimisers::StopCondition;

use anyhow::{Context, Result};

#[derive(Clone, clap::ValueEnum, Copy, Debug)]
#[allow(clippy::upper_case_acronyms)]
pub enum GapHandling {
    PIP,
    TKF91,
    TKF92,
    Missing,
}

#[derive(Clone, clap::ValueEnum, Copy, Debug)]
#[allow(clippy::upper_case_acronyms)]
pub enum SubstModelId {
    WAG,
    HIVB,
    BLOSUM,
    JC69,
    K80,
    HKY,
    TN93,
    GTR,
    NONE,
}

impl Display for SubstModelId {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{:?}", self)
    }
}

#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
pub struct Cli {
    /// Output folder
    #[arg(short = 'd', long, value_name = "OUTPUT_FOLDER", default_value = ".")]
    pub out_folder: PathBuf,

    /// Run identifier for output
    #[arg(short, long, value_name = "RUN_NAME")]
    pub run_name: Option<String>,

    /// Sequence file in fasta format
    #[arg(short, long, value_name = "SEQ_FILE")]
    pub seq_file: PathBuf,

    /// Tree file in newick format (required for model optimization)
    #[arg(short, long, value_name = "TREE_FILE")]
    pub tree_file: PathBuf,

    /// Sequence evolution model
    #[arg(short, long, value_name = "MODEL", ignore_case = true)]
    pub model: SubstModelId,

    /// Sequence evolution model parameters.
    #[arg(short = 'p', long, value_name = "MODEL_PARAMS", num_args = 0..)]
    pub params: Vec<f64>,

    /// Sequence evolution model stationary frequencies
    #[arg(short = 'f', long, value_name = "FREQUENCIES", num_args = 0..)]
    pub freqs: Vec<f64>,

    /// Frequency optimisation method: fixed, empirical, estimated
    #[arg(
        short = 'o',
        long,
        value_name = "FREQ_OPTIMISATION",
        default_value = "empirical"
    )]
    pub freq_opt: FrequencyOptimisation,

    /// Gap handling method
    #[arg(
        short,
        long,
        value_name = "GAP_HANDLING",
        ignore_case = true,
        default_value = "pip"
    )]
    pub gap_handling: GapHandling,

    /// Epsilon value for numerical optimisation
    #[arg(short, long, value_name = "EPSILON", default_value = "1e-5")]
    pub epsilon: f64,

    /// PRNG seed
    #[arg(long = "seed", value_name = "PRNG_SEED")]
    pub prng_seed: Option<u64>,

    /// Do not add a timestamp to the output folder and files
    #[arg(long, default_value_t = false)]
    pub no_timestamp: bool,

    /// Log level for the run
    #[arg(
        short = 'l',
        long,
        value_name = "LOG_LEVEL",
        default_value = "info",
        ignore_case = true
    )]
    pub log_level: LevelFilter,
}

pub struct Config {
    pub timestamp: Timestamp,
    #[allow(dead_code)]
    pub out_fldr: PathBuf,
    pub out_logl: PathBuf,
    pub out_params: PathBuf,
    #[allow(dead_code)]
    pub run_id: String,
    pub seq_file: PathBuf,
    pub input_tree: PathBuf,
    pub model: SubstModelId,
    pub params: Vec<f64>,
    pub freqs: Vec<f64>,
    pub freq_opt: FrequencyOptimisation,
    pub gap_handling: GapHandling,
    pub stop_condition: StopCondition,
    pub prng_seed: Option<u64>,
    #[allow(dead_code)]
    pub no_timestamp: bool,
    #[allow(dead_code)]
    pub log_level: LevelFilter,
}

impl Config {
    pub fn setup(cli: Cli) -> Result<Self> {
        let timestamp = Timestamp::now();
        let run_id = match (cli.run_name, cli.no_timestamp) {
            (Some(name), false) => format!("{}_{}", name, timestamp.as_u64()),
            (Some(name), true) => name,
            (None, false) => format!("{}", timestamp.as_u64()),
            (None, true) => "model_search".to_string(),
        };

        let out_fldr = cli.out_folder.join(format!("{run_id}_out"));
        std::fs::create_dir_all(&out_fldr).context("Failed to create output directory")?;

        Ftail::new()
            .datetime_format("%H:%M:%S")
            .console(cli.log_level)
            .single_file(
                out_fldr.join(format!("{run_id}.log")).to_str().unwrap(),
                true,
                LevelFilter::Debug,
            )
            .init()
            .context("Failed to initialize logger")?;

        let out_logl = out_fldr.join(format!("{run_id}_logl.out"));
        let out_params = out_fldr.join(format!("{run_id}_params.json"));

        Ok(Config {
            timestamp,
            out_fldr,
            out_logl,
            out_params,
            run_id,
            seq_file: cli.seq_file,
            input_tree: cli.tree_file,
            model: cli.model,
            params: cli.params,
            freqs: cli.freqs,
            freq_opt: cli.freq_opt,
            gap_handling: cli.gap_handling,
            stop_condition: StopCondition::epsilon(cli.epsilon),
            prng_seed: cli.prng_seed,
            no_timestamp: cli.no_timestamp,
            log_level: cli.log_level,
        })
    }
}
