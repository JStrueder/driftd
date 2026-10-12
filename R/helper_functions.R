
#' @title Extract Individual Participant Parameter Estimates
#'
#' @description
#' A getter function to safely extract the table of individual-level parameter
#' estimates from a \code{driftd_results} object.
#'
#' @param x An object of class \code{driftd_results}, typically the output
#'   from the \code{\link{driftd}} function.
#'
#' @return A \code{data.frame} or \code{tibble} where each row represents a
#'   participant's fit and columns contain the estimated diffusion model
#'   parameters, along with convergence information and fit statistics (e.g., LLE, AIC, BIC).
#'
#' @details
#' This function performs a class check to ensure the input object is valid,
#' preventing errors that might occur if a user tries to access the internal
#' structure of an object that wasn't created by \code{\link{driftd}}.
#'
#' @seealso
#' \code{\link{driftd}} for the main fitting function,
#' \code{\link{summary.driftd_results}} for group-level summary statistics.
#'
#' @export
individual_estimates <- function(x) {
  if (!inherits(x, "driftd_results")) {
    stop("Object must be of class 'driftd_results'")
  }
  return(x$individual_estimates)
}




#' @title Internal: Retrieve Diffusion Model Parameters for a Specific Condition
#'
#' @description
#' A highly specific internal helper function used to map the estimated
#' parameter columns from the \code{individual_estimates} table back to a simple
#' list of parameter values needed for a simulation in a given experimental condition.
#'
#' @param participant_params A single row (\code{data.frame} or \code{tibble})
#'   containing all estimated parameters for one participant.
#' @param condition A \code{character} string specifying the experimental
#'   condition cell, as stored in the \code{cond} column of the fitted data
#'   (e.g., \code{"neg_0.45"}).
#' @param settings A \code{list} containing the model fitting settings,
#'   particularly \code{settings$cell_parameters}, which records the name of the
#'   parameter that applies to each condition cell.
#'
#' @return A \code{list} containing the seven primary DDM parameter values
#'   (\code{a}, \code{t0}, \code{z}, \code{v}, \code{st0}, \code{sz}, \code{sv})
#'   for the specified condition.
#'
#' @details
#' For each parameter, the column to use is looked up in
#' \code{settings$cell_parameters}, which \code{\link{driftd}} builds by matching
#' trials to parameters on the actual condition values. Fixed parameters (recorded
#' as \code{NA}) are read from the column named after the parameter itself
#' (e.g., \code{st0}).
#'
#' **Note:** This function is primarily intended for internal use by the
#' \code{\link{predict.driftd_results}} and \code{\link{simulate.driftd_results}}, and is not typically called by end-users.
#'
#' @keywords internal
get_parameters_for_condition <- function(participant_params, condition, settings) {

  if (is.null(settings$cell_parameters)) {
    stop("This driftd_results object was created by an older version of driftd. ",
         "Please refit the model.")
  }

  params <- list(a = NA_real_, t0 = NA_real_, z = NA_real_, v = NA_real_,
                 st0 = NA_real_, sz = NA_real_, sv = NA_real_)

  cell <- settings$cell_parameters[settings$cell_parameters$.cell == condition, , drop = FALSE]

  for (param_name in names(params)) {
    column <- if (nrow(cell) == 1) cell[[param_name]] else NA_character_
    if (is.na(column)) column <- param_name

    value <- participant_params[[column]]
    if (!(is.numeric(value) && length(value) == 1 && is.finite(value))) {
      stop("No valid estimate of ", param_name, " for condition '", condition, "'.")
    }
    params[[param_name]] <- value
  }

  return(params)
}




# Internal: EZ-diffusion estimates (Wagenmakers, van der Maas, & Grasman, 2007,
# Psychonomic Bulletin & Review, 14, 3-22) of drift rate, boundary separation, and
# non-decision time from the proportion of upper-boundary responses and the mean and
# variance of the RTs (in seconds), with the diffusion coefficient s = 1 used by rtdists.
# Vectorised over condition cells. Used to derive starting values for the full fit.
#
# EZ assumes an unbiased starting point and no between-trial variability, so RTs at both
# boundaries have the same distribution; RTs from both responses are therefore pooled.
# Proportions of 0, .5 and 1 are moved by 1/(2n), following the paper's edge correction.
# Cells where EZ is undefined (fewer than 2 trials, zero RT variance) return NA.
#' @noRd
ez_diffusion <- function(prop_upper, mean_rt, var_rt, n_trials) {
  edge <- 1 / (2 * n_trials)
  p <- pmin(pmax(prop_upper, edge), 1 - edge)
  p <- ifelse(abs(p - 0.5) < 1e-12, 0.5 + edge, p)
  L <- stats::qlogis(p)
  x <- L * (L * p^2 - L * p + p - 0.5) / var_rt
  v <- sign(p - 0.5) * x^(1 / 4)
  a <- L / v
  y <- -v * a
  mdt <- (a / (2 * v)) * (1 - exp(y)) / (1 + exp(y))
  out <- data.frame(v = v, a = a, t0 = mean_rt - mdt)
  out[!(n_trials >= 2 & is.finite(var_rt) & var_rt > 0), ] <- NA
  out[!is.finite(as.matrix(out))] <- NA
  out
}




# Internal: element k (k >= 1) of the van der Corput sequence in the given prime base.
# Using a different prime per dimension gives the Halton sequence, a deterministic
# quasi-random sequence that covers [0, 1)^d more evenly than independent random draws.
#' @noRd
halton <- function(k, base) {
  f <- 1
  r <- 0
  while (k > 0) {
    f <- f / base
    r <- r + f * (k %% base)
    k <- k %/% base
  }
  r
}
