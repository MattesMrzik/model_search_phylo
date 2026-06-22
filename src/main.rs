use std::fs::File;
use std::io::Write;

use anyhow::{Context, Result};
use clap::Parser;
use log::info;
use serde::Serialize;

use phylo::alphabets::{Alphabet, DNA_ALPHABET, PROTEIN_ALPHABET};
use phylo::likelihood::{ModelSearchCost, TreeSearchCost};
use phylo::optimisers::ModelOptimiser;
use phylo::phylo_info::PhyloInfoBuilder;
use phylo::pip_model::{PIPCostBuilder, PIPModel};
use phylo::random::DefaultGenerator;
use phylo::substitution_models::{
    dna_models::*, protein_models::*, SubstModel, SubstitutionCostBuilder,
};
use phylo::tkf_model::{
    TKF91CostBuilder, TKF91IndelCostBuilder, TKF92CostBuilder, TKF92IndelCostBuilder,
};

mod cli;
use crate::cli::{Cli, Config, GapHandling as Gap, SubstModelId};

macro_rules! extract_result {
    ($res:expr) => {{
        let res = $res;
        let cost = &res.cost;
        let params: Vec<f64> = (0..cost.param_count()).map(|i| cost.param(i)).collect();
        let freqs = cost.freqs().clone();
        (
            ModelSearchCost::cost(cost),
            res.final_cost,
            params,
            freqs,
            cost.tree().clone(),
        )
    }};
}

macro_rules! run_model_search {
    ($macro:ident, $cfg:expr, $info:expr) => {
        match $cfg.model {
            SubstModelId::JC69 => extract_result!($macro!(JC69, $cfg, $info)),
            SubstModelId::K80 => extract_result!($macro!(K80, $cfg, $info)),
            SubstModelId::HKY => extract_result!($macro!(HKY, $cfg, $info)),
            SubstModelId::TN93 => extract_result!($macro!(TN93, $cfg, $info)),
            SubstModelId::GTR => extract_result!($macro!(GTR, $cfg, $info)),
            SubstModelId::WAG => extract_result!($macro!(WAG, $cfg, $info)),
            SubstModelId::HIVB => extract_result!($macro!(HIVB, $cfg, $info)),
            SubstModelId::BLOSUM => extract_result!($macro!(BLOSUM, $cfg, $info)),
            SubstModelId::NONE => unreachable!("Model::NONE is only valid for TKF91/TKF92"),
        }
    };
}

macro_rules! pip_optimisation {
    ($model:ty, $cfg:expr, $info:expr) => {
        ModelOptimiser::with_stop_condition(
            PIPCostBuilder::new(PIPModel::<$model>::new(&$cfg.freqs, &$cfg.params), $info)
                .build()?,
            $cfg.freq_opt,
            $cfg.stop_condition,
        )
        .run()?
    };
}

macro_rules! tkf91_optimisation {
    ($model:ty, $cfg:expr, $info:expr) => {{
        let split = 2.min($cfg.params.len());

        let core = &$cfg.params[..split];
        let rest = &$cfg.params[split..];

        ModelOptimiser::with_stop_condition(
            TKF91CostBuilder::new(core, SubstModel::<$model>::new(&$cfg.freqs, rest), $info)
                .build()?,
            $cfg.freq_opt,
            $cfg.stop_condition,
        )
        .run()?
    }};
}

macro_rules! tkf92_optimisation {
    ($model:ty, $cfg:expr, $info:expr) => {{
        let split = 3.min($cfg.params.len());

        let core = &$cfg.params[..split];
        let rest = &$cfg.params[split..];

        ModelOptimiser::with_stop_condition(
            TKF92CostBuilder::new(core, SubstModel::<$model>::new(&$cfg.freqs, rest), $info)
                .build()?,
            $cfg.freq_opt,
            $cfg.stop_condition,
        )
        .run()?
    }};
}

macro_rules! subst_optimisation {
    ($model:ty, $cfg:expr, $info:expr) => {
        ModelOptimiser::with_stop_condition(
            SubstitutionCostBuilder::new(
                SubstModel::<$model>::new(&$cfg.freqs, &$cfg.params),
                $info,
            )
            .build()?,
            $cfg.freq_opt,
            $cfg.stop_condition,
        )
        .run()?
    };
}

macro_rules! tkf91_indel_optimisation {
    ($cfg:expr, $info:expr) => {
        ModelOptimiser::with_stop_condition(
            TKF91IndelCostBuilder::new(&$cfg.params, $info).build()?,
            $cfg.freq_opt,
            $cfg.stop_condition,
        )
        .run()?
    };
}

macro_rules! tkf92_indel_optimisation {
    ($cfg:expr, $info:expr) => {
        ModelOptimiser::with_stop_condition(
            TKF92IndelCostBuilder::new(&$cfg.params, $info).build()?,
            $cfg.freq_opt,
            $cfg.stop_condition,
        )
        .run()?
    };
}

fn main() -> Result<()> {
    let cli = Cli::parse();
    let cfg = Config::setup(cli)?;

    info!("Starting model parameter search on a fixed tree.");
    info!("Input sequences: {}", cfg.seq_file.display());
    info!("Input tree: {}", cfg.input_tree.display());

    let alphabet: &'static Alphabet = match cfg.model {
        SubstModelId::JC69
        | SubstModelId::K80
        | SubstModelId::HKY
        | SubstModelId::TN93
        | SubstModelId::GTR
        | SubstModelId::NONE => {
            info!("Assuming DNA sequences");
            &DNA_ALPHABET
        }
        SubstModelId::WAG | SubstModelId::HIVB | SubstModelId::BLOSUM => {
            info!("Assuming protein sequences");
            &PROTEIN_ALPHABET
        }
    };

    let mut rng = DefaultGenerator::new(cfg.prng_seed.unwrap_or_else(|| cfg.timestamp.as_u64()));

    let (log_likelihood, final_cost, params, freqs, _final_tree) = match cfg.gap_handling {
        Gap::TKF91 | Gap::TKF92 => {
            let info = PhyloInfoBuilder::new(&cfg.seq_file)
                .tree_file(Some(&cfg.input_tree))
                .alphabet(Some(alphabet))
                .build_with_ancestors_w_rng(&mut rng)
                .context("Failed to build PhyloInfo with ancestors for TKF model")?;

            match cfg.gap_handling {
                Gap::TKF91 => {
                    if matches!(cfg.model, SubstModelId::NONE) {
                        extract_result!(tkf91_indel_optimisation!(&cfg, info.clone()))
                    } else {
                        run_model_search!(tkf91_optimisation, &cfg, info.clone())
                    }
                }
                Gap::TKF92 => {
                    if matches!(cfg.model, SubstModelId::NONE) {
                        extract_result!(tkf92_indel_optimisation!(&cfg, info.clone()))
                    } else {
                        run_model_search!(tkf92_optimisation, &cfg, info.clone())
                    }
                }
                _ => unreachable!(),
            }
        }
        Gap::PIP | Gap::Missing => {
            if matches!(cfg.model, SubstModelId::NONE) {
                anyhow::bail!("Model::NONE is only compatible with TKF91 and TKF92");
            }
            let info = PhyloInfoBuilder::new(&cfg.seq_file)
                .tree_file(Some(&cfg.input_tree))
                .alphabet(Some(alphabet))
                .build_w_rng(&mut rng)
                .context("Failed to build PhyloInfo for PIP/Substitution model")?;

            match cfg.gap_handling {
                Gap::PIP => run_model_search!(pip_optimisation, &cfg, info.clone()),
                Gap::Missing => run_model_search!(subst_optimisation, &cfg, info.clone()),
                _ => unreachable!(),
            }
        }
    };
    assert_eq!(log_likelihood, final_cost);

    info!("Optimization complete.");
    info!("Final log-likelihood: {}", log_likelihood);
    // info!("Saving optimized tree to {}", cfg.out_tree.display());

    // write_newick_to_file(std::slice::from_ref(&final_tree), &cfg.out_tree)
    //     .context("Failed to write optimized tree to file")?;

    let mut out_logl =
        File::create(&cfg.out_logl).context("Failed to create log-likelihood file")?;
    writeln!(out_logl, "{}", final_cost).context("Failed to write log-likelihood to file")?;

    #[derive(Serialize)]
    struct ParamsOutput {
        params: Vec<f64>,
        freqs: Vec<f64>,
    }
    let params_output = ParamsOutput {
        params,
        freqs: freqs.iter().cloned().collect(),
    };
    let mut out_params = File::create(&cfg.out_params).context("Failed to create params file")?;
    serde_json::to_writer_pretty(&mut out_params, &params_output)
        .context("Failed to write params to file")?;

    Ok(())
}
