#' @title Parameter Recovery for a Diffusion Model Fit
#'
#' @description
#' Checks whether the parameters of a fitted model can be recovered with the design and
#' number of trials of the data: data are simulated from each participant's estimates
#' (the "true" parameters), the same model is fitted to the simulated data, and the
#' recovered parameters are compared with the true ones.
#'
#' @param fit An object of class \code{driftd_results}, the output of \code{\link{driftd}}.
#' @param n_rep Integer: number of simulated datasets (each with all participants). Default 1.
#' @param seed Optional integer seed for the simulations. The random number generator state
#'   is restored afterwards.
#' @param compute_se Logical: compute standard errors in the refits, which adds the coverage
#'   of 95% confidence intervals to the summary. Default FALSE (faster).
#' @param verbose Logical: report progress. Default TRUE.
#'
#' @return An object of class \code{driftd_recovery}, a list with:
#'   \itemize{
#'     \item \code{data}: one row per replication, participant and estimated parameter with
#'       the parameter type, its condition columns, \code{name}, \code{true}, \code{recovered},
#'       \code{se} (if \code{compute_se}), and the convergence code of the refit.
#'     \item \code{summary}: one row per parameter (\code{name}) with the correlation between
#'       true and recovered values across participants (\code{r}), the mean difference
#'       (\code{bias}), the root mean squared error (\code{rmse}) and, if \code{compute_se},
#'       the proportion of 95% confidence intervals containing the true value
#'       (\code{coverage}).
#'     \item \code{summary_by_type}: the same, pooled over the conditions of each parameter type.
#'     \item \code{error_correlation}: correlations between the recovery errors (recovered
#'       minus true) of all pairs of parameters, across participants and replications. Large
#'       values indicate trade-offs: when one parameter is over- or underestimated, the other
#'       tends to be as well.
#'   }
#'
#' @details
#' Parameter recovery shows how precisely the parameters can be estimated given the model,
#' the design and the number of trials (Wilson & Collins, 2019). Using each participant's
#' estimates as the true parameters keeps them in a realistic range. Recovery is good if
#' recovered values correlate highly with the true values, with little bias, and if the
#' recovered value of one parameter does not depend on the true values of other parameters.
#' Correlations need several participants; with few participants, increase \code{n_rep}.
#'
#' Each replication refits the model for every participant, so this takes about as long as
#' the original fit, times \code{n_rep}.
#'
#' @references
#' Wilson, R. C., & Collins, A. G. E. (2019). Ten simple rules for the computational modeling
#' of behavioral data. \emph{eLife, 8}, e49547.
#'
#' @seealso \code{\link{recover_models}}, \code{\link{simulate.driftd_results}}
#' @export
recover_parameters <- function(fit, n_rep = 1, seed = NULL, compute_se = FALSE, verbose = TRUE) {
  fit_columns(fit)
  check_refittable(fit)
  if (!is.null(seed)) local_seed(seed)
  id <- fit$settings$columns$sub_id
  truth <- estimates(fit, include_fixed = FALSE)

  reps <- lapply(seq_len(n_rep), function(r) {
    if (verbose) message("Replication ", r, " of ", n_rep, ": simulating and refitting ", nrow(fit$individual_estimates),
                         " participants")
    sim <- simulate(fit, nsim = 1)
    refit <- refit_model(fit, sim, compute_se = compute_se)
    rec <- estimates(refit, include_fixed = FALSE)
    out <- truth
    out$converge <- NULL
    out$fixed <- NULL
    names(out)[names(out) == "estimate"] <- "true"
    m <- match(paste(out[[id]], out$name), paste(rec[[id]], rec$name))
    out$recovered <- rec$estimate[m]
    out$se <- if (compute_se) rec$se[m] else NULL
    out$converge <- rec$converge[m]
    cbind(rep = r, out)
  })
  data <- do.call(rbind, reps)
  rownames(data) <- NULL

  structure(list(data = data,
                 summary = recovery_summary(data, "name", compute_se),
                 summary_by_type = recovery_summary(data, "parameter", compute_se),
                 error_correlation = error_correlation(data, id),
                 n_rep = n_rep, sub_id = id),
            class = "driftd_recovery")
}




#' @title Model Recovery for a Set of Diffusion Models
#'
#' @description
#' Checks whether the models being compared can be distinguished with the design and number
#' of trials of the data: data are simulated from each fitted model, every model is fitted to
#' every simulated dataset, and the model selected by AIC and BIC is recorded.
#'
#' @param ... Two or more objects of class \code{driftd_results} fitted to the same data,
#'   as for \code{\link{compare_models}}. Arguments can be named.
#' @param names Optional character vector of model names.
#' @param n_rep Integer: number of simulated datasets per generating model. Default 1.
#' @param seed Optional integer seed for the simulations. The random number generator state
#'   is restored afterwards.
#' @param verbose Logical: report progress. Default TRUE.
#'
#' @return An object of class \code{driftd_model_recovery}, a list with:
#'   \itemize{
#'     \item \code{confusion}: for \code{"AIC"} and \code{"BIC"}, a matrix with the
#'       generating models in rows and the selected models in columns, giving the proportion
#'       of simulated participants for whom each model was selected. Rows sum to 1; a good
#'       design has values close to 1 on the diagonal.
#'     \item \code{confusion_group}: the same at the group level (summed criteria), as the
#'       proportion of the \code{n_rep} simulated datasets.
#'     \item \code{comparisons}: the \code{\link{compare_models}} results for every generating
#'       model and replication.
#'   }
#'
#' @details
#' Model recovery shows whether the data could, in principle, tell the models apart, and
#' whether the selection criterion is biased towards simpler or more complex models (Wilson &
#' Collins, 2019). Simulated data use each participant's estimates under the generating model,
#' with the observed design and trial counts.
#'
#' This fits every model \code{n_rep} times to data from every model, i.e. (number of
#' models)^2 x \code{n_rep} fits of all participants, which can take a long time.
#'
#' @references
#' Wilson, R. C., & Collins, A. G. E. (2019). Ten simple rules for the computational modeling
#' of behavioral data. \emph{eLife, 8}, e49547.
#'
#' @seealso \code{\link{compare_models}}, \code{\link{recover_parameters}}
#' @export
recover_models <- function(..., names = NULL, n_rep = 1, seed = NULL, verbose = TRUE) {
  fits <- list(...)
  if (length(fits) < 2) stop("recover_models() needs at least two fits")
  if (is.null(names)) {
    names <- names(fits)
    exprs <- vapply(as.list(substitute(list(...)))[-1], function(e) paste(deparse(e), collapse = ""), character(1))
    if (is.null(names)) names <- exprs else names[names == ""] <- exprs[names == ""]
  }
  # validates the fits and checks that they share the data
  invisible(do.call(compare_models, c(unname(fits), list(names = names))))
  invisible(lapply(fits, check_refittable))
  if (!is.null(seed)) local_seed(seed)

  n_models <- length(fits)
  design <- recovery_design(fits)
  if (verbose) message("Model recovery: ", n_models^2 * n_rep, " fits of all participants")
  comparisons <- list()
  counts <- list(AIC = matrix(0, n_models, n_models, dimnames = list(generating = names, selected = names)),
                 BIC = matrix(0, n_models, n_models, dimnames = list(generating = names, selected = names)))
  group <- counts
  for (g in seq_len(n_models)) {
    for (r in seq_len(n_rep)) {
      if (verbose) message("Generating model '", names[g], "', replication ", r, " of ", n_rep)
      sim <- simulate_trials(fits[[g]], design)
      refits <- lapply(fits, refit_model, data = sim, compute_se = FALSE)
      cmp <- do.call(compare_models, c(unname(refits), list(names = names)))
      comparisons[[paste(names[g], r, sep = "_rep")]] <- cmp
      for (crit in c("AIC", "BIC")) {
        best <- cmp$table[[paste0("n_best_", crit)]]
        counts[[crit]][g, ] <- counts[[crit]][g, ] + best
        winner <- which.min(cmp$table[[crit]])
        group[[crit]][g, winner] <- group[[crit]][g, winner] + 1
      }
    }
  }
  to_prop <- function(m) sweep(m, 1, rowSums(m), "/")
  structure(list(confusion = lapply(counts, to_prop), confusion_group = lapply(group, to_prop),
                 comparisons = comparisons, n_rep = n_rep, n_participants = nrow(fits[[1]]$individual_estimates)),
            class = "driftd_model_recovery")
}




#' @title Print Recovery Results
#' @description Prints the summary of a parameter or model recovery.
#' @param x An object of class \code{driftd_recovery} or \code{driftd_model_recovery}.
#' @param digits Number of decimal places. Default 2.
#' @param ... Not used.
#' @return \code{x}, invisibly.
#' @export
print.driftd_recovery <- function(x, digits = 2, ...) {
  n_id <- length(unique(x$data[[x$sub_id]]))
  cat("Parameter recovery (", n_id, " participants, ", x$n_rep, " simulated dataset",
      if (x$n_rep > 1) "s", ")\n\n", sep = "")
  round_df <- function(df) {
    num <- vapply(df, is.numeric, logical(1))
    df[num] <- lapply(df[num], function(v) round(v, digits))
    df
  }
  cat("By parameter type (pooled over conditions):\n")
  print(round_df(x$summary_by_type), row.names = FALSE)
  cat("\nBy parameter:\n")
  print(round_df(x$summary), row.names = FALSE)
  n_bad <- sum(x$data$converge != 0, na.rm = TRUE)
  if (n_bad > 0) cat("\nNote:", n_bad, "recovered values come from refits that did not converge.\n")
  invisible(x)
}

#' @rdname print.driftd_recovery
#' @export
print.driftd_model_recovery <- function(x, digits = 2, ...) {
  cat("Model recovery (", x$n_participants, " participants, ", x$n_rep, " simulated dataset",
      if (x$n_rep > 1) "s", " per model)\n", sep = "")
  for (crit in c("BIC", "AIC")) {
    cat("\nSelected by ", crit, ", per simulated participant (rows: generating model):\n", sep = "")
    print(round(x$confusion[[crit]], digits))
  }
  cat("\nSelected by summed BIC, per simulated dataset:\n")
  print(round(x$confusion_group$BIC, digits))
  invisible(x)
}




#' @title Plot Recovery Results
#'
#' @description
#' For parameter recovery: recovered against true values, one panel per parameter type
#' (or per parameter), with the identity line and the correlation in the panel title.
#' For model recovery: the confusion matrix as a heatmap.
#'
#' @param x An object of class \code{driftd_recovery} or \code{driftd_model_recovery}.
#' @param by For parameter recovery: \code{"parameter"} (default, one panel per parameter
#'   type, conditions pooled) or \code{"name"} (one panel per parameter).
#' @param criterion For model recovery: \code{"BIC"} (default) or \code{"AIC"}.
#' @param level For model recovery: \code{"participant"} (default) or \code{"group"}.
#' @param ... Not used.
#' @return A \code{ggplot} object.
#' @export
plot.driftd_recovery <- function(x, by = c("parameter", "name"), ...) {
  by <- match.arg(by)
  d <- x$data[!is.na(x$data$recovered), ]
  r <- tapply(seq_len(nrow(d)), d[[by]], function(k) suppressWarnings(stats::cor(d$true[k], d$recovered[k])))
  levels_order <- unique(d[[by]])
  d$panel <- factor(paste0(d[[by]], "  (r = ", formatC(r[d[[by]]], format = "f", digits = 2), ")"),
                    levels = paste0(levels_order, "  (r = ", formatC(r[levels_order], format = "f", digits = 2), ")"))
  # invisible points on the diagonal so that every panel has equal x and y ranges
  rng <- do.call(rbind, lapply(split(d, d$panel, drop = TRUE), function(s) {
    data.frame(panel = s$panel[1], v = range(c(s$true, s$recovered)))
  }))
  ggplot2::ggplot(d, ggplot2::aes(x = .data$true, y = .data$recovered)) +
    ggplot2::geom_blank(data = rng, ggplot2::aes(x = .data$v, y = .data$v)) +
    ggplot2::geom_abline(intercept = 0, slope = 1, colour = "grey50", linetype = "dashed") +
    ggplot2::geom_point(colour = level_colors[1], alpha = 0.7, size = 1.8) +
    ggplot2::facet_wrap("panel", scales = "free") +
    ggplot2::labs(x = "True value (estimate from the data)", y = "Recovered value",
                  title = "Parameter recovery",
                  caption = "Each point is one simulated participant. Dashed line: perfect recovery.") +
    ggplot2::theme_bw() +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank())
}

#' @rdname plot.driftd_recovery
#' @export
plot.driftd_model_recovery <- function(x, criterion = c("BIC", "AIC"), level = c("participant", "group"), ...) {
  criterion <- match.arg(criterion)
  level <- match.arg(level)
  m <- if (level == "participant") x$confusion[[criterion]] else x$confusion_group[[criterion]]
  d <- as.data.frame(as.table(m), stringsAsFactors = TRUE)
  names(d) <- c("generating", "selected", "p")
  d$generating <- factor(d$generating, levels = rev(rownames(m)))
  ggplot2::ggplot(d, ggplot2::aes(x = .data$selected, y = .data$generating, fill = .data$p)) +
    ggplot2::geom_tile(colour = "white") +
    ggplot2::geom_text(ggplot2::aes(label = formatC(.data$p, format = "f", digits = 2),
                                    colour = .data$p > 0.5), show.legend = FALSE) +
    ggplot2::scale_fill_gradient(low = "#F2F2F2", high = level_colors[1], limits = c(0, 1), name = "Proportion") +
    ggplot2::scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = "black")) +
    ggplot2::labs(x = paste0("Selected model (", criterion, ")"), y = "Generating model",
                  title = "Model recovery",
                  caption = wrap_caption(paste("Proportion of simulated",
                                               if (level == "participant") "participants" else "datasets",
                                               "for which each model was selected. Rows sum to 1."), width = 60)) +
    ggplot2::theme_minimal()
}




# ---- Internal helpers ---------------------------------------------------------------------

# Stop unless the fit stores the model specification needed to refit it.
#' @noRd
check_refittable <- function(fit) {
  if (is.null(fit$settings$model) || is.null(fit$settings$parameters)) {
    stop("This driftd_results object was created by an older version of driftd. Please refit the model.",
         call. = FALSE)
  }
  invisible(TRUE)
}

# Fit the model specification of `fit` to new data in the layout returned by simulate().
#' @noRd
refit_model <- function(fit, data, compute_se = FALSE) {
  cols <- fit$settings$columns
  args <- c(list(data = data, sub_id = cols$sub_id, response = cols$response, rt = cols$rt),
            fit$settings$model, list(compute_se = compute_se, ask = FALSE))
  out <- NULL
  utils::capture.output(out <- suppressMessages(do.call(driftd, args)))
  out
}

# The observed trials (in their original order) with the participant, response and RT columns
# and every condition column used by any of the fits, so that data simulated from one model
# can be fitted with all of them.
#' @noRd
recovery_design <- function(fits) {
  design <- data.frame(row = seq_len(nrow(fits[[1]]$data)))
  for (f in fits) {
    cols <- f$settings$columns
    for (col in c(cols$sub_id, cols$response, cols$rt, f$settings$condition_columns)) {
      if (!col %in% names(design)) design[[col]] <- f$data[[col]]
    }
  }
  design$row <- NULL
  design
}

# Simulate one response and RT for every trial of `design` from the model `fit` (fitted to
# the same trials in the same order), keeping all other columns.
#' @noRd
simulate_trials <- function(fit, design) {
  cols <- fit$settings$columns
  d <- fit$data
  id <- cols$sub_id
  est <- fit$individual_estimates
  rt <- numeric(nrow(d))
  resp <- character(nrow(d))
  for (idx in split(seq_len(nrow(d)), list(d[[id]], d$.cell), drop = TRUE, sep = "\r")) {
    pars <- get_parameters_for_condition(est[est[[id]] == d[[id]][idx[1]], , drop = FALSE], d$.cell[idx[1]],
                                         fit$settings)
    s <- simulate_cell(pars, length(idx))
    rt[idx] <- s$rt
    resp[idx] <- as.character(s$response)
  }
  design[[cols$rt]] <- rt
  design[[cols$response]] <- factor(resp, levels = c("lower", "upper"))
  design
}

# Recovery statistics per group of rows (parameter name or type).
#' @noRd
recovery_summary <- function(data, by, compute_se) {
  groups <- unique(data[[by]])
  res <- do.call(rbind, lapply(groups, function(g) {
    d <- data[data[[by]] == g & !is.na(data$recovered), ]
    err <- d$recovered - d$true
    row <- data.frame(group = g, parameter = d$parameter[1], n = nrow(d),
                      r = if (nrow(d) > 2 && stats::sd(d$true) > 0) stats::cor(d$true, d$recovered) else NA_real_,
                      bias = mean(err), rmse = sqrt(mean(err^2)), stringsAsFactors = FALSE)
    if (compute_se) {
      ok <- !is.na(d$se)
      row$coverage <- if (any(ok)) mean(abs(err[ok]) <= 1.96 * d$se[ok]) else NA_real_
    }
    row
  }))
  if (by == "parameter") {
    res$group <- NULL
  } else {
    names(res)[names(res) == "group"] <- "name"
  }
  res
}

# Correlations between the recovery errors (recovered - true) of all pairs of parameters,
# across participants and replications.
#' @noRd
error_correlation <- function(data, id) {
  key <- paste(data$rep, data[[id]], sep = "\r")
  params <- unique(data$name)
  err <- tapply(data$recovered - data$true, list(key, factor(data$name, levels = params)), function(v) v[1])
  out <- suppressWarnings(stats::cor(err[, params, drop = FALSE], use = "pairwise.complete.obs"))
  dimnames(out) <- list(params, params)
  out
}
