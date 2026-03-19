#!/bin/bash
# Install R dependencies for the HIV Kenya model.
# Run from the repo root:  bash install_R.sh
#
# This installs the R packages (reticulate, devtools, testthat, rstarsim)
# and sets up a Python environment managed by R with starsim and stisim.
#
# Usage:
#   bash install_R.sh                 # Full install (including Python env)
#   bash install_R.sh --skip-python   # R packages only (Python already installed)

set -e

echo "Installing R packages (reticulate, devtools, testthat)..."
Rscript -e 'install.packages(c("reticulate", "devtools", "testthat"), repos = "https://cloud.r-project.org")'

echo "Installing rstarsim from GitHub..."
Rscript -e 'devtools::install_github("starsimhub/rstarsim")'

if [ "$1" != "--skip-python" ]; then
  echo "Initializing Starsim (creates a Python environment if needed)..."
  Rscript -e 'library(starsim); init_starsim()'

  echo "Installing STIsim into the R-managed Python environment..."
  Rscript -e 'library(reticulate); reticulate::py_install("stisim", pip = TRUE)'
fi

echo "Done. You can now run: Rscript tests/test_model.R"
