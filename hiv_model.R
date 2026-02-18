#' Create HIV model and interventions for Kenya
#'
#' R translation of hiv_model.py, using reticulate to call stisim/starsim.
#' Requires: install.packages("reticulate")
#'           devtools::install_github("starsimhub/rstarsim")
#'
#' Set working directory to the repo root before sourcing this file.

# %% Imports and settings
library(reticulate)
library(starsim)

load_starsim()
sti <- import("stisim")
sc  <- import("sciris")
np  <- import("numpy")
pd  <- import("pandas")


# Kenya-specific parameters
sim_pars <- list(
  start         = 1985L,
  stop          = 2030L,
  n_agents      = 10000,
  use_migration = TRUE
)

sti_pars <- list(
  hiv = list(
    beta_m2f      = 0.012,
    eff_condom    = 0.95,
    rel_init_prev = 0.5
  )
)

nw_pars <- list(
  prop_f0     = 0.8,
  prop_m0     = 0.75,
  f1_conc     = 0.15,
  m1_conc     = 0.15,
  p_pair_form = 0.5
)


make_custom_interventions <- function(test_years = NULL) {
  #' Create custom interventions for Kenya: HIV testing, ART, and PrEP.
  #'
  #' ART is created here (rather than auto-loaded from art_coverage.csv)
  #' to allow setting future_coverage.
  #'
  #' @param test_years Integer vector of years for testing coverage. Default: 1990:2050.
  #' @return List of intervention instances.

  if (is.null(test_years)) {
    test_years <- 1990L:2050L
  }

  scaleup_end   <- min(2020L, tail(test_years, 1))
  scaleup_years <- test_years[1]:scaleup_end
  n_scaleup     <- length(scaleup_years)
  n_future      <- length(test_years) - n_scaleup

  fsw_prob     <- c(seq(0, 0.75, length.out = n_scaleup), seq(0.75, 0.85, length.out = n_future))
  gp_prob      <- c(seq(0, 0.1,  length.out = n_scaleup), seq(0.1,  0.1,  length.out = n_future))
  low_cd4_prob <- c(seq(0, 0.85, length.out = n_scaleup), seq(0.85, 0.95, length.out = n_future))

  fsw_eligibility <- function(sim) {
    sim$networks$structuredsexual$fsw & !sim$diseases$hiv$diagnosed & !sim$diseases$hiv$on_art
  }

  other_eligibility <- function(sim) {
    !sim$networks$structuredsexual$fsw & !sim$diseases$hiv$diagnosed & !sim$diseases$hiv$on_art
  }

  low_cd4_eligibility <- function(sim) {
    (sim$diseases$hiv$cd4 < 200L) & !sim$diseases$hiv$diagnosed
  }

  # Testing
  testing <- list(
    sti$HIVTest(years = test_years, test_prob_data = fsw_prob,
                name = "fsw_testing", eligibility = fsw_eligibility, label = "fsw_testing"),
    sti$HIVTest(years = test_years, test_prob_data = gp_prob,
                name = "other_testing", eligibility = other_eligibility, label = "other_testing"),
    sti$HIVTest(years = test_years, test_prob_data = low_cd4_prob,
                name = "low_cd4_testing", eligibility = low_cd4_eligibility, label = "low_cd4_testing")
  )

  # ART
  data_path <- file.path(getwd(), "data")
  n_art <- pd$read_csv(file.path(data_path, "n_art.csv"))$set_index("year")
  art   <- sti$ART(coverage_data = n_art, future_coverage = list(year = 2024L, prop = 0.97))

  # PrEP
  prep <- sti$Prep(
    coverage = c(0, 0.01, 0.5, 0.8),
    years    = c(2004L, 2005L, 2015L, 2025L),
    eff_prep = 0.8
  )

  c(testing, list(art, prep))
}


make_sim_pars <- function(sim, calib_pars) {
  #' Apply calibration parameters to a simulation.
  #'
  #' @param sim An sti.Sim instance (may be uninitialized).
  #' @param calib_pars Named list of calibration parameter values.
  #' @return The modified sim.

  if (!sim$initialized) sim$init()
  hiv <- sim$diseases$hiv
  nw  <- sim$networks$structuredsexual

  for (k in names(calib_pars)) {
    pars_val <- calib_pars[[k]]

    if (k == "rand_seed") {
      sim$pars$rand_seed <- pars_val
      next
    }
    if (k %in% c("index", "mismatch")) next

    if (is.list(pars_val)) {
      v <- pars_val$value
    } else if (is.numeric(pars_val)) {
      v <- pars_val
    } else {
      stop(paste("Parameter", k, "not recognized"))
    }

    if (grepl("hiv_", k)) {
      k <- sub("hiv_", "", k)
      hiv$pars[[k]] <- v
    } else if (grepl("nw_", k)) {
      k <- sub("nw_", "", k)
      if (grepl("pair_form", k)) {
        nw$pars[[k]]$set(v)
      } else {
        nw$pars[[k]] <- v
      }
    } else {
      stop(paste("Parameter", k, "not recognized"))
    }
  }

  sim
}


make_sim <- function(...) {
  #' Create a Kenya HIV simulation.
  #'
  #' Uses data_path to auto-load init_prev and condom_use data via DataLoader.
  #' Custom interventions (testing, ART with future_coverage, PrEP) are created

  #' separately and merged with any user-provided interventions.
  #'
  #' @param ... Override any Sim parameters (e.g. verbose, rand_seed, stop, analyzers, interventions).
  #' @return An sti.Sim instance.

  kwargs <- list(...)

  intvs      <- make_custom_interventions()
  user_intvs <- if (!is.null(kwargs$interventions)) as.list(kwargs$interventions) else list()
  kwargs$interventions <- NULL

  # Default analyzers
  user_analyzers <- if (!is.null(kwargs$analyzers)) as.list(kwargs$analyzers) else list()
  kwargs$analyzers <- NULL
  analyzers <- c(user_analyzers, list(sti$sw_stats(diseases = list("hiv"))))

  sim_args <- c(
    list(
      location      = "kenya",
      diseases      = "hiv",
      data_path     = file.path(getwd(), "data"),
      sim_pars      = sim_pars,
      nw_pars       = nw_pars,
      sti_pars      = sti_pars,
      interventions = c(intvs, user_intvs),
      analyzers     = analyzers
    ),
    kwargs
  )

  do.call(sti$Sim, sim_args)
}


run_msim <- function(use_calib = TRUE, n_pars = 1L, do_save = TRUE) {
  #' Run multiple simulations, optionally applying calibration parameters.
  #'
  #' @param use_calib Whether to apply saved calibration parameters.
  #' @param n_pars Number of parameter sets to run.
  #' @param do_save Whether to save results.
  #' @return List of completed sims.

  calib <- if (use_calib) sc$loadobj("results/kenya_hiv_calib.obj") else NULL

  sims <- list()
  for (par_idx in seq_len(n_pars) - 1L) {  # 0-indexed to match Python calibration df
    sim <- make_sim(verbose = -1L)
    if (use_calib) {
      calib_pars <- calib$df$iloc[par_idx]$to_dict()
      sim$init()
      sim <- make_sim_pars(sim, calib_pars)
      message(sprintf("Using calibration parameters for index %d", par_idx))
    }
    sim$par_idx <- par_idx
    sims <- c(sims, list(sim))
  }
  sims <- ss$parallel(sims)$sims

  if (do_save) {
    dfs <- list()
    for (sim in iterate(sims)) {
      par_idx <- sim$par_idx
      df <- sim$to_df(resample = "year", use_years = TRUE, sep = ".")
      df["res_no"] <- par_idx
      dfs <- c(dfs, list(df))
    }
    df <- pd$concat(dfs)
    sc$saveobj("results/msim.df", df)
  }

  sims
}


save_stats <- function(sims, resfolder = "results") {
  #' Save age/sex stratified epi stats and SW stats.

  dfs <- list()
  for (sim in iterate(sims)) {
    par_idx  <- sim$par_idx
    age_bins <- sim$diseases$hiv$age_bins
    n_bins   <- length(age_bins)

    for (sex in c("f", "m")) {
      sex_label <- if (sex == "f") "Female" else "Male"
      for (i in seq_len(n_bins - 1)) {
        ab1 <- age_bins[i]
        ab2 <- age_bins[i + 1]
        age <- if (ab1 == 65) "65+" else paste0(ab1, "-", ab2)

        prev_key <- sprintf("prevalence_%s_%s_%s", sex, ab1, ab2)
        inf_key  <- sprintf("new_infections_%s_%s_%s", sex, ab1, ab2)

        dd <- data.frame(
          age            = age,
          sex            = sex_label,
          prevalence     = as.numeric(tail(sim$results[["hiv"]][[prev_key]], 1)),
          new_infections = as.numeric(mean(tail(sim$results[["hiv"]][[inf_key]], 120))),
          par_idx        = par_idx
        )
        dfs <- c(dfs, list(dd))
      }
    }
  }
  epi_df <- do.call(rbind, dfs)
  sc$saveobj(file.path(resfolder, "epi_df.df"), epi_df)

  # Save SW stats (from first sim)
  first_sim <- NULL
  for (sim in iterate(sims)) {
    if (sim$par_idx == 0) { first_sim <- sim; break }
  }
  sw_res <- first_sim$results[["sw_stats"]]
  sw_df  <- sw_res$to_df(resample = "year", use_years = TRUE, sep = ".")
  sc$saveobj(file.path(resfolder, "sw_df.df"), sw_df)
}


# %% Run as a script (when sourced interactively or via Rscript)
run_main <- function() {

  # SETTINGS
  debug     <- FALSE
  seed      <- 1L
  do_save   <- TRUE
  do_run    <- TRUE
  do_plot   <- TRUE
  use_calib <- TRUE

  to_run <- c(
    # "run_sim",
    "run_msim"
  )

  if ("run_sim" %in% to_run) {

    if (do_run) {
      sim <- make_sim(rand_seed = seed, verbose = 1/12)
      if (use_calib) {
        calib      <- sc$loadobj("results/kenya_hiv_calib.obj")
        calib_pars <- calib$df$iloc[0L]$to_dict()
        sim$init()
        sim <- make_sim_pars(sim, calib_pars)
        message("Using calibration parameters")
      }
      sim$run()
      df       <- sim$to_df(resample = "year", use_years = TRUE, sep = ".")
      df$index <- df[["timevec"]]
      if (do_save) {
        sc$saveobj("results/kenya_sim.df", df)
        sc$saveobj("results/kenya.sim", sim)
      }
    } else {
      df <- sc$loadobj("results/kenya_sim.df")
    }

    if (do_plot) {
      # Use Python plotting (or implement ggplot2 equivalent)
      source("plot_sims.R", local = TRUE)
      plot_hiv_sims(df, start_year = 1985L, title = "hiv_plots")
    }
  }

  if ("run_msim" %in% to_run) {
    n_pars <- if (!debug) 50L else 2L
    if (do_run) {
      sims <- run_msim(use_calib = use_calib, n_pars = n_pars, do_save = do_save)
    } else {
      sims <- NULL
    }

    if (do_save && !is.null(sims)) {
      save_stats(sims, resfolder = "results")
    }
  }
}

# Uncomment to run:
# run_main()