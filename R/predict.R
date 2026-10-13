#' @title Model Predictions for a Diffusion Model Fit
#'
#' @description
#' Compares the observed data with the predictions of the fitted diffusion model, for
#' every participant and condition cell: the probability of each response and the
#' quantiles of its RT distribution. Predictions are computed exactly from the fitted
#' parameters. Optionally, prediction intervals show how much the observed statistics are
#' expected to vary under the model, given the number of trials in each cell.
#'
#' @param object An object of class \code{driftd_results}, the output of \code{\link{driftd}}.
#' @param probs Numeric vector of RT quantiles to compare. Default is the conventional
#'   \code{c(0.1, 0.3, 0.5, 0.7, 0.9)}.
#' @param intervals Logical: compute prediction intervals by simulation? Default TRUE.
#' @param n_sim Integer: number of simulated replicate datasets for the prediction intervals.
#'   Default 200.
#' @param level Coverage of the prediction intervals. Default 0.95.
#' @param min_n Integer: minimum number of responses of a given type in a cell for its RT
#'   quantiles to be computed. Below this, observed (and simulated) quantiles are \code{NA}
#'   and only the response probability is compared. Default 5.
#' @param parameter_uncertainty Logical: if TRUE, each simulated replicate uses parameters
#'   drawn from a multivariate normal distribution with the estimates as means and
#'   \code{object$vcov} as covariance (a parametric bootstrap), so the prediction intervals
#'   also reflect estimation uncertainty. Parameters without a standard error are kept at
#'   their estimates, and draws are truncated to the parameter bounds used in fitting.
#'   Default FALSE.
#' @param seed Optional integer seed for the simulations, for reproducible intervals. The
#'   random number generator state is restored afterwards.
#' @param ... Not used.
#'
#' @return A data frame of class \code{driftd_prediction} with one row per participant,
#'   condition cell, response (\code{"upper"}, \code{"lower"}) and statistic, and columns:
#'   \itemize{
#'     \item the participant ID column, \code{.cell} (a label of the condition cell) and the condition columns;
#'     \item \code{response}, \code{statistic} (\code{"probability"} or \code{"quantile"}) and
#'       \code{quantile} (the quantile level, \code{NA} for probabilities);
#'     \item \code{n_cell} (trials in the cell) and \code{n_response} (trials with this response);
#'     \item \code{observed}, \code{predicted}, and the prediction interval \code{lower} and
#'       \code{upper} (\code{NA} if \code{intervals = FALSE}). Probabilities are between 0
#'       and 1; RT quantiles are in seconds (plots show milliseconds by default);
#'     \item \code{.row}, a row identifier that links each row to its simulated replicates
#'       (used by \code{aggregate()} and \code{plot()}; keep it when subsetting).
#'   }
#'   The result is a data frame and can be filtered, reordered and summarized with base R or
#'   dplyr, and \code{aggregate()}, \code{summary()} and \code{plot()} still work on subsets.
#'   Use \code{\link{aggregate.driftd_prediction}} for group-level (averaged) values,
#'   \code{\link{summary.driftd_prediction}} for a per-participant misfit summary, and
#'   \code{\link{plot.driftd_prediction}} for plots.
#'
#' @details
#' \strong{Predictions.} The predicted probability of each response and the predicted RT
#' quantiles of each response (the quantiles of its defective distribution, scaled by its
#' probability) are computed exactly with \code{\link[rtdists]{pdiffusion}} and
#' \code{\link[rtdists]{qdiffusion}}, so they contain no simulation noise. Observed quantiles
#' use \code{\link[stats]{quantile}} (type 7).
#'
#' \strong{Prediction intervals.} For each participant and cell, \code{n_sim} datasets with
#' the observed number of trials are simulated from the fitted model
#' (\code{\link[rtdists]{rdiffusion}}), and the same statistics are computed for each. The
#' interval is the central \code{level} range of the simulated statistics: under the model,
#' an observed statistic falls outside it with probability of about \code{1 - level}. Cells
#' with few trials therefore get wide intervals.
#'
#' \strong{Group level.} Group-level values are averages of the participant-level statistics
#' (Vincentizing; Ratcliff, 1979), computed for observed and predicted values over the same
#' participants, rather than statistics of the pooled trials.
#'
#' @references
#' Ratcliff, R. (1979). Group reaction time distributions and an analysis of distribution
#' statistics. \emph{Psychological Bulletin, 86}, 446-461.
#'
#' Ratcliff, R., & McKoon, G. (2008). The diffusion decision model: Theory and data for
#' two-choice decision tasks. \emph{Neural Computation, 20}, 873-922.
#'
#' @seealso \code{\link{plot.driftd_prediction}}, \code{\link{simulate.driftd_results}}
#'
#' @importFrom stats predict
#' @export
predict.driftd_results <- function(object, probs = c(0.1, 0.3, 0.5, 0.7, 0.9), intervals = TRUE,
                                   n_sim = 200, level = 0.95, min_n = 5,
                                   parameter_uncertainty = FALSE, seed = NULL, ...) {
  cols <- fit_columns(object)
  if (!is.numeric(probs) || any(probs <= 0 | probs >= 1)) stop("probs must be between 0 and 1")
  if (level <= 0 || level >= 1) stop("level must be between 0 and 1")
  if (intervals && n_sim < 20) stop("n_sim should be at least 20 for prediction intervals")
  probs <- sort(probs)
  if (!is.null(seed)) local_seed(seed)

  d <- object$data
  id <- cols$sub_id
  est <- object$individual_estimates
  precision <- if (is.null(object$settings$precision)) 3 else object$settings$precision
  upper_resp <- d[[cols$response]] == "upper"

  stat_info <- data.frame(response = rep(c("upper", "lower"), each = 1 + length(probs)),
                          statistic = rep(c("probability", rep("quantile", length(probs))), 2),
                          quantile = rep(c(NA, probs), 2),
                          stringsAsFactors = FALSE)
  n_stat <- nrow(stat_info)

  cell_rows <- split(seq_len(nrow(d)), list(d[[id]], d$.cell), drop = TRUE, sep = "\r")
  out <- vector("list", length(cell_rows))
  reps <- vector("list", length(cell_rows))

  for (k in seq_along(cell_rows)) {
    idx <- cell_rows[[k]]
    pid <- d[[id]][idx[1]]
    cd <- d$.cell[idx[1]]
    rt <- d[[cols$rt]][idx]
    up <- upper_resp[idx]
    est_row <- est[est[[id]] == pid, , drop = FALSE]
    pars <- get_parameters_for_condition(est_row, cd, object$settings)

    obs <- cell_statistics(rt, up, probs, min_n)
    pred <- predicted_statistics(pars, probs, precision)
    n_resp <- rep(c(sum(up), sum(!up)), each = 1 + length(probs))

    out[[k]] <- data.frame(ID = pid, .cell = cd, stat_info, n_cell = length(idx), n_response = n_resp,
                           observed = obs, predicted = pred, stringsAsFactors = FALSE)

    if (intervals) {
      par_list <- if (parameter_uncertainty) {
        lapply(draw_parameters(est_row, object$vcov[[as.character(pid)]], n_sim, min(d[[cols$rt]][d[[id]] == pid])),
               get_parameters_for_condition, condition = cd, settings = object$settings)
      } else {
        list(pars)
      }
      reps[[k]] <- replicate_statistics(par_list, length(idx), n_sim, probs, min_n)
    }
  }

  res <- do.call(rbind, out)
  if (intervals) {
    rep_mat <- do.call(rbind, reps)
    bounds <- interval_bounds(rep_mat, level)
    res$lower <- bounds[, 1]
    res$upper <- bounds[, 2]
  } else {
    rep_mat <- NULL
    res$lower <- NA_real_
    res$upper <- NA_real_
  }

  # add the condition columns and order the rows by participant, then by cell
  cells <- object$settings$cells
  res <- cbind(res[, c("ID", ".cell")], cells[match(res$.cell, cells$.cell), setdiff(names(cells), ".cell"), drop = FALSE],
               res[, setdiff(names(res), c("ID", ".cell"))])
  ord <- order(match(res$ID, unique(est[[id]])), match(res$.cell, cells$.cell), seq_len(nrow(res)))
  res <- res[ord, ]
  if (!is.null(rep_mat)) rep_mat <- rep_mat[ord, , drop = FALSE]
  rownames(res) <- NULL
  names(res)[1] <- id
  # .row links each row to its simulated replicates, so they stay matched after any
  # subsetting or reordering of the rows (base R or dplyr)
  res$.row <- seq_len(nrow(res))
  if (!is.null(rep_mat)) rownames(rep_mat) <- res$.row

  structure(res, class = c("driftd_prediction", "data.frame"), replicates = rep_mat,
            driftd = list(sub_id = id, condition_columns = object$settings$condition_columns,
                          probs = probs, level = level, min_n = min_n, n_sim = n_sim,
                          intervals = intervals, parameter_uncertainty = parameter_uncertainty))
}




#' @title Group-Level Averages of Model Predictions
#'
#' @description
#' Averages participant-level observed and predicted statistics from
#' \code{\link{predict.driftd_results}} across participants (Vincentizing), and optionally
#' across the levels of condition columns not listed in \code{by}.
#'
#' @param x An object of class \code{driftd_prediction}.
#' @param by Character vector of condition columns to keep. Statistics are averaged over
#'   participants and over all other condition columns. Default: all condition columns,
#'   i.e. one row per condition cell.
#' @param ... Not used.
#'
#' @return A data frame with the \code{by} columns, \code{response}, \code{statistic},
#'   \code{quantile}, \code{n_participants} (number of participants contributing),
#'   \code{observed}, \code{predicted} (averages over participants), and \code{lower} and
#'   \code{upper}, the prediction interval for individual participants (see Details).
#'
#' @details
#' Only participant-cells with an observed value contribute to each average, and observed
#' and predicted values are averaged over the same participant-cells. RT quantiles are
#' therefore averaged only over participants with at least \code{min_n} responses of that
#' type (see \code{\link{predict.driftd_results}}).
#'
#' The prediction interval shows where the statistics of \emph{individual participants} are
#' expected to lie under the fitted model, not the uncertainty of the average. The simulated
#' replicates of all participants in a cell (\code{n_sim} per participant, averaged over the
#' condition columns not in \code{by} like the observed values) are pooled, and the interval
#' is their central \code{level} range. It therefore reflects the differences between
#' participants as well as trial-to-trial variability (and estimation uncertainty if
#' \code{parameter_uncertainty = TRUE} in \code{predict()}), and about \code{level} of the
#' individual participants' observed values should fall inside it. Because each participant is
#' fitted separately, this interval mainly shows the expected spread; misfit of individual
#' participants is better judged with their own intervals (\code{plot(..., participant = )},
#' \code{\link{summary.driftd_prediction}}).
#'
#' @importFrom stats aggregate
#' @export
aggregate.driftd_prediction <- function(x, by = NULL, ...) {
  info <- attr(x, "driftd")
  if (is.null(info)) stop("x must be the output of predict() on a driftd fit")
  if (is.null(by)) by <- info$condition_columns
  bad <- setdiff(by, info$condition_columns)
  if (length(bad) > 0) {
    stop("by must name condition columns of the fit (", paste(info$condition_columns, collapse = ", "),
         "); not: ", paste(bad, collapse = ", "))
  }

  use <- !is.na(x$observed)
  x_use <- as.data.frame(x)[use, , drop = FALSE]
  keys <- x_use[, c(by, "response", "statistic", "quantile"), drop = FALSE]
  group <- do.call(paste, c(lapply(keys, as.character), sep = "\r"))
  group <- factor(group, levels = unique(group))

  n_part <- tapply(x_use[[info$sub_id]], group, function(v) length(unique(v)))
  res <- keys[!duplicated(group), , drop = FALSE]
  res$n_participants <- as.integer(n_part)
  res$observed <- as.vector(tapply(x_use$observed, group, mean))
  res$predicted <- as.vector(tapply(x_use$predicted, group, mean))

  reps <- replicates_for(x, use)
  if (!is.null(reps)) {
    # Prediction intervals for individual participants: each participant's simulated
    # replicates (averaged over the condition columns not in `by`, like the observed values)
    # are pooled over participants, and the interval is the central `level` range of the
    # pooled values. It shows where individual participants' statistics are expected to lie.
    pgroup <- paste(x_use[[info$sub_id]], group, sep = "\r")
    pgroup <- factor(pgroup, levels = unique(pgroup))
    sums <- rowsum(reps, pgroup, na.rm = TRUE, reorder = FALSE)
    counts <- rowsum((!is.na(reps)) * 1, pgroup, reorder = FALSE)
    pmeans <- sums / counts
    pmeans[counts == 0] <- NA
    pmeans_group <- group[match(levels(pgroup), as.character(pgroup))]
    bounds <- t(vapply(split(seq_len(nrow(pmeans)), pmeans_group), function(rows) {
      v <- pmeans[rows, , drop = FALSE]
      ok <- !is.na(v)
      if (sum(ok) < length(v) / 2) return(c(NA_real_, NA_real_))
      alpha <- (1 - info$level) / 2
      stats::quantile(v[ok], c(alpha, 1 - alpha), names = FALSE)
    }, numeric(2)))
    res$lower <- bounds[, 1]
    res$upper <- bounds[, 2]
  } else {
    res$lower <- NA_real_
    res$upper <- NA_real_
  }
  rownames(res) <- NULL
  res
}




#' @title Summarize Model Misfit per Participant
#'
#' @description
#' For each participant, counts how many observed statistics (response probabilities and RT
#' quantiles across all cells) fall outside their prediction intervals.
#'
#' @param object An object of class \code{driftd_prediction} with prediction intervals.
#' @param ... Not used.
#'
#' @return A data frame with one row per participant, sorted from the worst to the best
#'   fitting: \code{n_statistics} (statistics with an observed value and an interval),
#'   \code{n_outside} and \code{prop_outside} (number and proportion outside the interval),
#'   and \code{mean_abs_dev} (mean absolute deviation between observed and predicted RT
#'   quantiles, in seconds).
#'
#' @details
#' If the model describes a participant's data well, the proportion of statistics outside
#' the intervals should be close to \code{1 - level} (e.g. 5%). Statistics within a
#' participant are correlated, so this is a descriptive screening measure, not a test.
#'
#' @export
summary.driftd_prediction <- function(object, ...) {
  info <- attr(object, "driftd")
  if (!info$intervals) stop("summary() needs prediction intervals; use predict(..., intervals = TRUE)")
  ok <- !is.na(object$observed) & !is.na(object$lower)
  outside <- ok & (object$observed < object$lower | object$observed > object$upper)
  id <- object[[info$sub_id]]
  ids <- unique(id)
  is_q <- object$statistic == "quantile" & ok
  res <- data.frame(id = ids,
                    n_statistics = as.vector(tapply(ok, factor(id, levels = ids), sum)),
                    n_outside = as.vector(tapply(outside, factor(id, levels = ids), sum)))
  res$prop_outside <- res$n_outside / res$n_statistics
  res$mean_abs_dev <- as.vector(tapply(ifelse(is_q, abs(object$observed - object$predicted), NA),
                                       factor(id, levels = ids), mean, na.rm = TRUE))
  names(res)[1] <- info$sub_id
  res <- res[order(-res$prop_outside, -res$mean_abs_dev), ]
  rownames(res) <- NULL
  res
}




#' @title Simulate Data from a Diffusion Model Fit
#'
#' @description
#' Simulates replicate datasets from the fitted model, with the same participants, condition
#' cells and number of trials per cell as the observed data.
#'
#' @param object An object of class \code{driftd_results}, the output of \code{\link{driftd}}.
#' @param nsim Integer: number of replicate datasets. Default 1.
#' @param seed Optional integer seed. The random number generator state is restored afterwards.
#' @param parameter_uncertainty Logical: draw the parameters of each replicate from their
#'   estimated sampling distribution (see \code{\link{predict.driftd_results}}). Default FALSE.
#' @param ... Not used.
#'
#' @return A data frame with columns \code{sim} (replicate number), the participant ID,
#'   \code{.cell} (a label of the condition cell), the condition columns, the response (a factor with levels
#'   \code{"lower"}, \code{"upper"}) and the RT (in seconds), using the column names of the
#'   original data.
#'
#' @importFrom stats simulate
#' @export
simulate.driftd_results <- function(object, nsim = 1, seed = NULL, parameter_uncertainty = FALSE, ...) {
  cols <- fit_columns(object)
  if (!is.null(seed)) local_seed(seed)
  d <- object$data
  id <- cols$sub_id
  est <- object$individual_estimates
  cells <- object$settings$cells
  cell_rows <- split(seq_len(nrow(d)), list(d[[id]], d$.cell), drop = TRUE, sep = "\r")

  sims <- lapply(cell_rows, function(idx) {
    pid <- d[[id]][idx[1]]
    cd <- d$.cell[idx[1]]
    est_row <- est[est[[id]] == pid, , drop = FALSE]
    n <- length(idx)
    if (parameter_uncertainty) {
      draws <- draw_parameters(est_row, object$vcov[[as.character(pid)]], nsim, min(d[[cols$rt]][d[[id]] == pid]))
      s <- do.call(rbind, lapply(draws, function(row) {
        simulate_cell(get_parameters_for_condition(row, cd, object$settings), n)
      }))
    } else {
      s <- simulate_cell(get_parameters_for_condition(est_row, cd, object$settings), n * nsim)
    }
    data.frame(sim = rep(seq_len(nsim), each = n), ID = pid, .cell = cd,
               response = s$response, rt = s$rt, stringsAsFactors = FALSE)
  })
  res <- do.call(rbind, sims)
  res <- res[order(res$sim, match(res$ID, unique(est[[id]])), match(res$.cell, cells$.cell)), ]
  res <- cbind(res[, c("sim", "ID", ".cell")],
               cells[match(res$.cell, cells$.cell), setdiff(names(cells), ".cell"), drop = FALSE],
               res[, c("response", "rt")])
  res$response <- factor(as.character(res$response), levels = c("lower", "upper"))
  names(res)[names(res) == "ID"] <- id
  names(res)[names(res) == "response"] <- cols$response
  names(res)[names(res) == "rt"] <- cols$rt
  rownames(res) <- NULL
  res
}




#' @title Subset Model Predictions
#'
#' @description
#' Subsetting a \code{driftd_prediction} object with \code{[} keeps the information needed by
#' \code{aggregate()}, \code{summary()} and \code{plot()}, for row and column subsets alike.
#'
#' @param x An object of class \code{driftd_prediction}.
#' @param ... Indices passed to the data frame method.
#'
#' @return The subset, as a \code{driftd_prediction} object if it is a data frame.
#' @keywords internal
#' @export
`[.driftd_prediction` <- function(x, ...) {
  out <- NextMethod()
  if (is.data.frame(out)) {
    attr(out, "driftd") <- attr(x, "driftd")
    attr(out, "replicates") <- attr(x, "replicates")
    class(out) <- class(x)
  }
  out
}




# ---- Internal helpers ---------------------------------------------------------------------

# Simulated replicate statistics for the rows `use` of a prediction table, looked up by the
# .row column (so they match after subsetting or reordering). NULL if unavailable.
#' @noRd
replicates_for <- function(x, use) {
  reps <- attr(x, "replicates")
  if (is.null(reps)) return(NULL)
  if (!".row" %in% names(x)) {
    warning("The .row column is missing, so prediction intervals for averages cannot be computed.",
            call. = FALSE)
    return(NULL)
  }
  reps[match(x$.row[use], as.integer(rownames(reps))), , drop = FALSE]
}

# Column names stored by driftd(); fails clearly for objects from older versions.
#' @noRd
fit_columns <- function(object) {
  if (!inherits(object, "driftd_results")) stop("object must be of class 'driftd_results' from driftd()")
  if (is.null(object$settings$columns) || is.null(object$settings$cells)) {
    stop("This driftd_results object was created by an older version of driftd. Please refit the model.")
  }
  object$settings$columns
}

# Set the RNG seed in the calling function and restore the previous RNG state when it exits.
#' @noRd
local_seed <- function(seed, envir = parent.frame()) {
  had_seed <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
  old_seed <- if (had_seed) get(".Random.seed", envir = globalenv(), inherits = FALSE) else NULL
  do.call(on.exit, list(substitute({
    if (HAD) assign(".Random.seed", OLD, envir = globalenv())
    else if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) rm(".Random.seed", envir = globalenv())
  }, list(HAD = had_seed, OLD = old_seed)), add = TRUE), envir = envir)
  set.seed(seed)
}

# Response probabilities and RT quantiles of one condition cell, in the order: upper
# probability, upper quantiles, lower probability, lower quantiles. Quantiles of a response
# need at least min_n trials with that response, otherwise they are NA.
#' @noRd
cell_statistics <- function(rt, upper, probs, min_n) {
  n <- length(rt)
  stats_for <- function(sel) {
    q <- if (sum(sel) >= min_n) stats::quantile(rt[sel], probs, names = FALSE) else rep(NA_real_, length(probs))
    c(sum(sel) / n, q)
  }
  c(stats_for(upper), stats_for(!upper))
}

# Exact predicted response probabilities and RT quantiles (same order as cell_statistics).
#' @noRd
predicted_statistics <- function(p, probs, precision = 3) {
  args <- list(a = p$a, v = p$v, t0 = p$t0, z = p$z * p$a, sz = p$sz, sv = p$sv, st0 = p$st0,
               precision = precision)
  p_upper <- do.call(rtdists::pdiffusion, c(list(rt = Inf, response = "upper"), args))
  p_resp <- c(upper = p_upper, lower = 1 - p_upper)
  q_for <- function(resp) {
    pr <- p_resp[[resp]]
    if (!is.finite(pr) || pr < 1e-6) return(rep(NA_real_, length(probs)))
    tryCatch(suppressWarnings(do.call(rtdists::qdiffusion, c(list(p = probs * pr, response = resp), args))),
             error = function(e) rep(NA_real_, length(probs)))
  }
  c(p_resp[["upper"]], q_for("upper"), p_resp[["lower"]], q_for("lower"))
}

# Simulate n trials of one cell; returns a data frame with rt and response.
#' @noRd
simulate_cell <- function(p, n) {
  rtdists::rdiffusion(n, a = p$a, v = p$v, t0 = p$t0, z = p$z * p$a, sz = p$sz, sv = p$sv, st0 = p$st0)
}

# Statistics of n_sim simulated replicates of one cell (columns). par_list holds either one
# parameter set (used for all replicates) or one per replicate.
#' @noRd
replicate_statistics <- function(par_list, n_cell, n_sim, probs, min_n) {
  n_stat <- 2 * (1 + length(probs))
  if (length(par_list) == 1) {
    s <- simulate_cell(par_list[[1]], n_cell * n_sim)
    rt <- matrix(s$rt, nrow = n_cell)
    up <- matrix(s$response == "upper", nrow = n_cell)
    vapply(seq_len(n_sim), function(r) cell_statistics(rt[, r], up[, r], probs, min_n), numeric(n_stat))
  } else {
    vapply(par_list, function(p) {
      s <- simulate_cell(p, n_cell)
      cell_statistics(s$rt, s$response == "upper", probs, min_n)
    }, numeric(n_stat))
  }
}

# Central `level` interval of each row of a replicate matrix; NA if fewer than half of the
# replicates have a value.
#' @noRd
interval_bounds <- function(rep_mat, level) {
  alpha <- (1 - level) / 2
  t(apply(rep_mat, 1, function(r) {
    r <- r[!is.na(r)]
    if (length(r) < ncol(rep_mat) / 2) return(c(NA_real_, NA_real_))
    stats::quantile(r, c(alpha, 1 - alpha), names = FALSE)
  }))
}

# n_sim parameter sets for one participant, drawn from a multivariate normal distribution
# with the estimates as means and vcov as covariance, truncated to the fitting bounds.
# Returns a list of one-row data frames in the layout of individual_estimates.
#' @noRd
draw_parameters <- function(est_row, vcov, n_sim, min_rt) {
  if (is.null(vcov)) {
    warning("No variance-covariance matrix for participant ", est_row[[1]],
            "; using the point estimates for all replicates.", call. = FALSE)
    return(rep(list(est_row), n_sim))
  }
  free <- colnames(vcov)
  L <- tryCatch(chol(vcov), error = function(e) NULL)
  if (is.null(L)) {
    warning("Variance-covariance matrix of participant ", est_row[[1]],
            " is not positive definite; using the point estimates for all replicates.", call. = FALSE)
    return(rep(list(est_row), n_sim))
  }
  draws <- matrix(stats::rnorm(n_sim * length(free)), n_sim) %*% L
  draws <- sweep(draws, 2, unlist(est_row[free]), "+")
  type <- sub("_.*", "", free)
  lower <- c(a = 0.01, t0 = 0, z = 0, v = -Inf, st0 = 0, sz = 0, sv = 0)[type]
  upper <- c(a = Inf, t0 = 0.999 * min_rt, z = 1, v = Inf, st0 = Inf, sz = 1, sv = Inf)[type]
  draws <- sweep(sweep(draws, 2, lower, pmax), 2, upper, pmin)
  lapply(seq_len(n_sim), function(r) {
    row <- est_row
    row[free] <- as.list(draws[r, ])
    row
  })
}
