#' Test that the HIV Kenya model runs successfully.
#'
#' Run from repo root:
#'   Rscript tests/test_model.R
#' Or with testthat:
#'   Rscript -e "testthat::test_file('tests/test_model.R')"

library(testthat)

# Source the model; if running from tests/, step up to repo root first
if (!file.exists("hiv_model.R")) setwd("..")
source("hiv_model.R")


test_that("HIV model runs and produces valid results", {

  # Create and run the sim
  sim <- make_sim(n_agents = 1000L)
  sim$run()

  # Extract results
  res  <- sim$results$hiv
  prev <- as.numeric(res$prevalence$values)
  art  <- as.numeric(res$n_on_art$values)

  # Check prevalence
  expect_true(all(prev > 0), label = "Expect nonzero prevalence at all timepoints")
  expect_gt(prev[length(prev)], prev[1], label = "Expect prevalence to increase during the sim")

  # Check ART
  expect_equal(art[1], 0, label = "Expect no one on ART at simulation start")
  expect_gt(art[length(art)], 0, label = "Expect people on ART at simulation end")
})
