#' Plot HIV simulation results for Kenya
#'
#' R translation of plot_sims.py. Produces a 2x3 panel of key HIV epi outputs
#' with validation data overlay. Uses base R graphics.
#'
#' Expects a data frame (from sim$to_df()) with columns like n_alive,
#' hiv.n_infected, hiv.prevalence_15_49, hiv.new_infections, hiv.new_deaths,
#' hiv.n_diagnosed, hiv.n_on_art, indexed by year.

library(reticulate)

# Load validation data -- use the calibration CSV which has the same columns
.load_hiv_data <- function(location = "kenya") {
  fpath <- file.path(getwd(), "data", paste0(location, "_hiv_data.csv"))
  if (!file.exists(fpath)) {
    fpath <- file.path(getwd(), "data", paste0(location, "_hiv_calib.csv"))
  }
  read.csv(fpath, stringsAsFactors = FALSE)
}


#' Extract y-values for plotting depending on single vs multi sim
#'
#' @param df Data frame of simulation results
#' @param which "single" or "multi"
#' @param resname Column name (e.g. "hiv.n_infected")
#' @return Numeric vector
get_y <- function(df, which, resname) {
  if (which == "single") {
    return(df[[resname]])
  } else if (which == "multi") {
    return(df[[paste0(resname, ".50%")]])
  }
}


#' Format large numbers with SI suffixes (K, M) for axis labels
si_format <- function(x) {
  ifelse(is.na(x), "",
    ifelse(abs(x) >= 1e6, paste0(round(x / 1e6, 1), "M"),
      ifelse(abs(x) >= 1e3, paste0(round(x / 1e3, 0), "K"),
        as.character(round(x, 1)))))
}


#' Plot HIV simulation results
#'
#' @param df Data frame from sim$to_df(), with years as index or in a "timevec" column.
#'   For single sims, columns are like "n_alive", "hiv.n_infected", etc.
#'   For multi sims, columns include quantile suffixes like "hiv.n_infected.50%".
#' @param start_year Start year for x-axis (default 2000)
#' @param end_year End year for x-axis (default 2025)
#' @param which "single" or "multi"
#' @param title File name prefix for saved figure
#' @param location Location name for loading data (default "kenya")
plot_hiv_sims <- function(df, start_year = 2000, end_year = 2025, which = "single",
                          title = "hiv_plots", location = "kenya") {

  # Handle reticulate pandas DataFrame -- convert to R data.frame
  if (inherits(df, "pandas.core.frame.DataFrame")) {
    years <- as.numeric(df$index$values)
    df_r <- as.data.frame(py_to_r(df))
    df_r$year <- years
  } else if (is.data.frame(df)) {
    if ("timevec" %in% names(df)) {
      df_r <- df
      df_r$year <- df_r$timevec
    } else if ("index" %in% names(df)) {
      df_r <- df
      df_r$year <- df_r$index
    } else {
      df_r <- df
      df_r$year <- as.numeric(rownames(df))
    }
  } else {
    stop("df must be a pandas DataFrame or R data.frame")
  }

  # Load validation data
  hiv_data <- .load_hiv_data(location)
  # The calib CSV uses "time" not "year"
  if ("time" %in% names(hiv_data) && !"year" %in% names(hiv_data)) {
    names(hiv_data)[names(hiv_data) == "time"] <- "year"
  }
  hiv_data <- hiv_data[hiv_data$year >= start_year & hiv_data$year <= end_year, ]

  # Subset model output to year range
  dfplot <- df_r[df_r$year >= start_year & df_r$year <= end_year, ]
  x <- dfplot$year

  # Create output directory
  dir.create("figures", showWarnings = FALSE)

  # Set up 2x3 layout
  outfile <- file.path("figures", paste0(title, start_year, "_", which, ".png"))
  png(outfile, width = 1800, height = 700, res = 100)
  par(mfrow = c(2, 3), mar = c(4, 4.5, 2.5, 1), oma = c(0, 0, 0, 0), cex = 1.1)

  col_model <- "#1f77b4"
  col_data  <- "black"

  # Helper to make SI-formatted y-axis
  si_axis <- function(ymax) {
    ticks <- pretty(c(0, ymax))
    axis(2, at = ticks, labels = si_format(ticks))
  }

  # --- Panel 1: Population size ---
  resname <- "n_alive"
  y <- get_y(dfplot, which, resname)
  ymax <- max(c(y, hiv_data[[resname]]), na.rm = TRUE) * 1.1
  plot(x, y, type = "l", col = col_model, lwd = 2,
       xlab = "", ylab = "", main = "Population size",
       ylim = c(0, ymax), yaxt = "n")
  si_axis(ymax)
  points(hiv_data$year, hiv_data[[resname]], pch = 16, col = col_data)
  legend("topleft", legend = c("Data", "Modeled"), col = c(col_data, col_model),
         pch = c(16, NA), lty = c(NA, 1), lwd = c(NA, 2), bty = "n")

  # --- Panel 2: PLHIV ---
  resname <- "hiv.n_infected"
  y <- get_y(dfplot, which, resname)
  ymax <- max(c(y, hiv_data[[resname]]), na.rm = TRUE) * 1.1
  plot(x, y, type = "l", col = col_model, lwd = 2,
       xlab = "", ylab = "", main = "PLHIV",
       ylim = c(0, ymax), yaxt = "n")
  si_axis(ymax)
  points(hiv_data$year, hiv_data[[resname]], pch = 16, col = col_data)

  # --- Panel 3: HIV prevalence ---
  resname <- "hiv.prevalence_15_49"
  y <- get_y(dfplot, which, resname) * 100
  data_prev <- hiv_data[[resname]] * 100
  ymax <- max(c(y, data_prev), na.rm = TRUE) * 1.1
  plot(x, y, type = "l", col = col_model, lwd = 2,
       xlab = "", ylab = "", main = "HIV prevalence (%)",
       ylim = c(0, ymax))
  points(hiv_data$year, data_prev, pch = 16, col = col_data)

  # --- Panel 4: HIV infections ---
  resname <- "hiv.new_infections"
  y <- get_y(dfplot, which, resname)
  ymax <- max(c(y, hiv_data[[resname]]), na.rm = TRUE) * 1.1
  plot(x, y, type = "l", col = col_model, lwd = 2,
       xlab = "", ylab = "", main = "HIV infections",
       ylim = c(0, ymax), yaxt = "n")
  si_axis(ymax)
  points(hiv_data$year, hiv_data[[resname]], pch = 16, col = col_data)

  # --- Panel 5: HIV deaths ---
  resname <- "hiv.new_deaths"
  y <- get_y(dfplot, which, resname)
  ymax <- max(c(y, hiv_data[[resname]]), na.rm = TRUE) * 1.1
  plot(x, y, type = "l", col = col_model, lwd = 2,
       xlab = "", ylab = "", main = "HIV-related deaths",
       ylim = c(0, ymax), yaxt = "n")
  si_axis(ymax)
  points(hiv_data$year, hiv_data[[resname]], pch = 16, col = col_data)

  # --- Panel 6: Diagnosed and treated (90-90-90) ---
  resnames <- c("hiv.n_infected" = "PLHIV", "hiv.n_diagnosed" = "Dx", "hiv.n_on_art" = "Treated")
  cols_909090 <- c("#1f77b4", "#ff7f0e", "#2ca02c")

  # Find y-range across all three series
  ymax <- 0
  for (rname in names(resnames)) {
    yy <- get_y(dfplot, which, rname)
    ymax <- max(ymax, max(yy, na.rm = TRUE))
  }
  ymax <- ymax * 1.1

  first <- TRUE
  for (i in seq_along(resnames)) {
    rname <- names(resnames)[i]
    yy <- get_y(dfplot, which, rname)
    if (first) {
      plot(x, yy, type = "l", col = cols_909090[i], lwd = 2,
           xlab = "", ylab = "", main = "Diagnosed and treated",
           ylim = c(0, ymax), yaxt = "n")
      si_axis(ymax)
      first <- FALSE
    } else {
      lines(x, yy, col = cols_909090[i], lwd = 2)
    }
  }
  points(hiv_data$year, hiv_data[["hiv.n_infected"]], pch = 16, col = col_data)
  legend("topleft", legend = unname(resnames), col = cols_909090,
         lty = 1, lwd = 2, bty = "n")

  dev.off()
  message(sprintf("Figure saved to %s", outfile))

  invisible(NULL)
}
