#' @title Compare Diffusion Models Fitted to the Same Data
#'
#' @description
#' Compares two or more \code{\link{driftd}} fits of the same data, e.g. a model in which a
#' manipulation affects only the drift rate against one in which it also affects the
#' boundary separation. Reports information criteria summed across participants, how many
#' participants each model describes best, and likelihood-ratio tests for nested models.
#'
#' @param ... Two or more objects of class \code{driftd_results}, fitted to the same data.
#'   Arguments can be named (e.g. \code{compare_models(drift = fit1, boundary = fit2)}).
#' @param names Optional character vector of model names. Default: the argument names, or
#'   the expressions passed.
#'
#' @return An object of class \code{driftd_comparison}, a list with:
#'   \itemize{
#'     \item \code{table}: one row per model with the number of free parameters per
#'       participant (\code{k}), and summed over participants: \code{logLik}, \code{AIC},
#'       \code{BIC}; their differences from the best model (\code{dAIC}, \code{dBIC}); AIC and
#'       BIC weights (\code{w_AIC}, \code{w_BIC}); the number of participants for whom the model
#'       has the lowest AIC and BIC (\code{n_best_AIC}, \code{n_best_BIC}); and the number of
#'       participants whose fit did not converge (\code{n_not_converged}).
#'     \item \code{participants}: one row per participant and model with \code{k},
#'       \code{n_trials}, \code{logLik}, \code{AIC}, \code{BIC}, their differences from the
#'       participant's best model, and the convergence code.
#'     \item \code{lrt}: one row per pair of nested models (\code{restricted} nested in
#'       \code{full}) with the summed likelihood-ratio statistic \code{chisq}, its degrees of
#'       freedom \code{df} (difference in free parameters times number of participants),
#'       \code{p}, the number of participants with a significant individual test at .05
#'       (\code{n_significant}), the number of participants for whom the full model fits worse
#'       than the restricted one (\code{n_negative}, which indicates an optimization problem),
#'       and \code{boundary} (TRUE if the restricted model fixes a variability parameter at 0).
#'       \code{NULL} if no models are nested.
#'   }
#'
#' @details
#' Each participant is fitted separately, so the log-likelihood of a model for the whole data
#' set is the sum of the participants' log-likelihoods, and AIC and BIC are summed likewise
#' (each participant's BIC uses that participant's number of trials). Lower values are better.
#' AIC tends to favor more complex models than BIC, especially with many trials; reporting
#' both, together with the number of participants best described by each model, shows how
#' robust a conclusion is. Weights are computed from the summed criteria.
#'
#' A model is nested in another if every parameter it estimates is also estimated by the
#' other, varying across the same or a subset of the condition columns, and every parameter
#' it fixes is either estimated by the other or fixed at the same value. For nested models,
#' twice the summed log-likelihood difference is compared with a chi-square distribution with
#' degrees of freedom equal to the summed difference in free parameters. If the restricted
#' model fixes \code{sv}, \code{sz} or \code{st0} at 0 (the boundary of their range), the
#' chi-square reference distribution is not exact and the test is conservative (Self & Liang,
#' 1987); this is flagged in the \code{boundary} column.
#'
#' Models must be fitted to identical data (same participants, trials and RTs). Participants
#' whose fits did not converge are included, but counted in \code{n_not_converged}; consider
#' inspecting them (\code{fit$fit_info}) before interpreting the comparison.
#'
#' @references
#' Self, S. G., & Liang, K.-Y. (1987). Asymptotic properties of maximum likelihood estimators
#' and likelihood ratio tests under nonstandard conditions. \emph{Journal of the American
#' Statistical Association, 82}, 605-610.
#'
#' @examples
#' \dontrun{
#' fit_a   <- driftd(speed_acc, sub_id = "id", vary_a = "condition", vary_v = "stimulus",
#'                   fix_st0 = FALSE)
#' fit_at0 <- driftd(speed_acc, sub_id = "id", vary_a = "condition", vary_t0 = "condition",
#'                   vary_v = "stimulus", fix_st0 = FALSE)
#' cmp <- compare_models(a = fit_a, a_t0 = fit_at0)
#' cmp
#' cmp$participants
#' }
#'
#' @seealso \code{\link{driftd}}
#' @importFrom stats pchisq
#' @export
compare_models <- function(..., names = NULL) {
  fits <- list(...)
  if (length(fits) < 2) stop("compare_models() needs at least two fits")
  if (!all(vapply(fits, inherits, logical(1), "driftd_results"))) {
    stop("All arguments must be driftd_results objects (outputs of driftd())")
  }
  invisible(lapply(fits, fit_columns))
  if (any(vapply(fits, function(f) is.null(f$settings$parameters), logical(1)))) {
    stop("A driftd_results object was created by an older version of driftd. Please refit the model.")
  }
  if (is.null(names)) {
    names <- names(fits)
    exprs <- vapply(as.list(substitute(list(...)))[-1], function(e) paste(deparse(e), collapse = ""), character(1))
    if (is.null(names)) names <- exprs else names[names == ""] <- exprs[names == ""]
  }
  if (length(names) != length(fits) || anyDuplicated(names)) stop("names must be unique, one per model")
  check_same_data(fits, names)

  id <- fits[[1]]$settings$columns$sub_id
  d <- fits[[1]]$data
  ids <- fits[[1]]$individual_estimates[[id]]
  n_trials <- as.vector(table(factor(as.character(d[[id]]), levels = as.character(ids))))

  # per participant and model
  per <- do.call(rbind, lapply(seq_along(fits), function(m) {
    f <- fits[[m]]
    est <- f$individual_estimates[match(as.character(ids), as.character(f$individual_estimates[[id]])), ]
    # parameters estimated for each participant (a parameter of a condition a participant
    # has no trials in is not estimated for them)
    estimated <- f$settings$parameters$name[!f$settings$parameters$fixed]
    k <- rowSums(!is.na(as.matrix(est[estimated])))
    data.frame(ID = ids, model = names[m], k = k, n_trials = n_trials, logLik = est$LLE,
               AIC = 2 * k - 2 * est$LLE, BIC = k * log(n_trials) - 2 * est$LLE,
               converge = est$converge, stringsAsFactors = FALSE)
  }))
  per$model <- factor(per$model, levels = names)
  best_aic <- tapply(per$AIC, as.character(per$ID), min)[as.character(per$ID)]
  best_bic <- tapply(per$BIC, as.character(per$ID), min)[as.character(per$ID)]
  per$dAIC <- per$AIC - best_aic
  per$dBIC <- per$BIC - best_bic
  per <- per[order(match(per$ID, ids), per$model), ]
  rownames(per) <- NULL

  # group level
  tab <- data.frame(model = names,
                    k = vapply(fits, function(f) sum(!f$settings$parameters$fixed), numeric(1)),
                    logLik = as.vector(tapply(per$logLik, per$model, sum)),
                    AIC = as.vector(tapply(per$AIC, per$model, sum)),
                    BIC = as.vector(tapply(per$BIC, per$model, sum)),
                    stringsAsFactors = FALSE)
  tab$dAIC <- tab$AIC - min(tab$AIC)
  tab$dBIC <- tab$BIC - min(tab$BIC)
  tab$w_AIC <- exp(-tab$dAIC / 2) / sum(exp(-tab$dAIC / 2))
  tab$w_BIC <- exp(-tab$dBIC / 2) / sum(exp(-tab$dBIC / 2))
  tab$n_best_AIC <- as.vector(tapply(per$dAIC == 0, per$model, sum))
  tab$n_best_BIC <- as.vector(tapply(per$dBIC == 0, per$model, sum))
  tab$n_not_converged <- as.vector(tapply(per$converge != 0, per$model, sum))
  names(per)[1] <- id

  structure(list(table = tab, participants = per, lrt = nested_tests(fits, names, per),
                 n_participants = length(ids)),
            class = "driftd_comparison")
}




#' @title Print a Model Comparison
#' @description Prints the group-level comparison, the likelihood-ratio tests of nested models,
#'   and notes on convergence and possible optimization problems.
#' @param x An object of class \code{driftd_comparison}.
#' @param digits Number of decimal places. Default 2.
#' @param ... Not used.
#' @return \code{x}, invisibly.
#' @export
print.driftd_comparison <- function(x, digits = 2, ...) {
  tab <- x$table
  cat("Model comparison (", x$n_participants, " participants; criteria summed over participants)\n\n", sep = "")
  shown <- tab
  num <- vapply(shown, is.numeric, logical(1))
  shown[num] <- lapply(shown[num], function(v) round(v, digits))
  print(shown, row.names = FALSE)
  cat("\nLowest AIC: ", tab$model[which.min(tab$AIC)], "   Lowest BIC: ", tab$model[which.min(tab$BIC)], "\n", sep = "")

  if (!is.null(x$lrt)) {
    cat("\nLikelihood-ratio tests of nested models:\n")
    lrt <- x$lrt
    lrt$chisq <- round(lrt$chisq, digits)
    lrt$p <- format.pval(lrt$p, digits = 3, eps = 1e-4)
    print(lrt, row.names = FALSE)
    if (any(lrt$boundary)) {
      cat("Note: boundary = TRUE means the restricted model fixes a variability parameter at 0;",
          "the test is then conservative.\n")
    }
    if (any(lrt$n_negative > 0)) {
      cat("Note: for some participants the full model fits worse than the restricted model (n_negative),",
          "which indicates that an optimization did not find the maximum.\n")
    }
  }
  if (any(tab$n_not_converged > 0)) {
    cat("\nNote: some fits did not converge (n_not_converged); see fit$fit_info.\n")
  }
  invisible(x)
}




#' @title Plot a Model Comparison per Participant
#'
#' @description
#' Shows, for each participant, how much worse each model is than that participant's best
#' model (difference in AIC or BIC; 0 = best). Reveals whether a group-level conclusion holds
#' for most participants or is driven by a few.
#'
#' @param x An object of class \code{driftd_comparison}.
#' @param criterion \code{"BIC"} (default) or \code{"AIC"}.
#' @param ... Not used.
#' @return A \code{ggplot} object.
#' @export
plot.driftd_comparison <- function(x, criterion = c("BIC", "AIC"), ...) {
  criterion <- match.arg(criterion)
  per <- x$participants
  id <- names(per)[1]
  per$delta <- per[[paste0("d", criterion)]]
  # participants sorted by how clearly the overall best model wins for them
  best <- x$table$model[which.min(x$table[[criterion]])]
  ord <- per[per$model == best, ]
  ord <- ord[order(-ord$delta), ]
  per$label <- factor(as.character(per[[id]]), levels = as.character(ord[[id]]))
  ggplot2::ggplot(per, ggplot2::aes(x = .data$delta, y = .data$label, colour = .data$model, shape = .data$model)) +
    ggplot2::geom_point(size = 2.2) +
    ggplot2::scale_x_continuous(trans = "pseudo_log", breaks = c(0, 2, 6, 10, 30, 100, 300, 1000, 3000, 10000)) +
    ggplot2::scale_colour_manual(values = rep(level_colors, length.out = nlevels(per$model)), name = "Model") +
    ggplot2::scale_shape_manual(values = rep(c(16, 17, 15, 18, 4, 8, 3, 1), length.out = nlevels(per$model)),
                                name = "Model") +
    ggplot2::labs(x = paste0("\u0394", criterion, " relative to the participant's best model (0 = best)"), y = id,
                  title = paste("Model comparison per participant:", criterion),
                  caption = wrap_caption(paste0("Lower is better. Differences above 10 are commonly read as strong ",
                                                "evidence against a model for that participant."))) +
    ggplot2::theme_bw() +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank())
}




# ---- Internal helpers ---------------------------------------------------------------------

# Stop unless all fits used the same participants, trials and RTs.
#' @noRd
check_same_data <- function(fits, names) {
  key <- function(f) {
    cols <- f$settings$columns
    list(id = as.character(f$data[[cols$sub_id]]), rt = round(f$data[[cols$rt]], 10),
         response = as.character(f$data[[cols$response]]))
  }
  ref <- key(fits[[1]])
  for (m in seq_along(fits)[-1]) {
    if (!identical(key(fits[[m]]), ref)) {
      stop("Models '", names[1], "' and '", names[m], "' were not fitted to the same data ",
           "(participants, trials, responses and RTs must be identical).", call. = FALSE)
    }
  }
}

# Is model a nested in model b? Returns NA if not nested, otherwise whether the test is at a
# boundary (a fixes sv, sz or st0 at 0 where b estimates it).
#' @noRd
nesting <- function(a, b) {
  boundary <- FALSE
  for (type in c("a", "t0", "z", "v", "st0", "sz", "sv")) {
    a_fixed <- a$settings$fixed_parameters[[type]]
    b_fixed <- b$settings$fixed_parameters[[type]]
    if (!a_fixed) {
      if (b_fixed) return(NA)
      if (!all(a$settings$varying_parameters[[type]] %in% b$settings$varying_parameters[[type]])) return(NA)
    } else if (b_fixed) {
      if (!isTRUE(all.equal(a$individual_estimates[[type]][1], b$individual_estimates[[type]][1]))) return(NA)
    } else if (type %in% c("st0", "sz", "sv") && a$individual_estimates[[type]][1] == 0) {
      boundary <- TRUE
    }
  }
  if (sum(!a$settings$parameters$fixed) >= sum(!b$settings$parameters$fixed)) return(NA)
  boundary
}

# Likelihood-ratio tests for all pairs of nested models (NULL if there are none).
#' @noRd
nested_tests <- function(fits, names, per) {
  rows <- list()
  for (i in seq_along(fits)) for (j in seq_along(fits)) {
    if (i == j) next
    boundary <- nesting(fits[[i]], fits[[j]])
    if (is.na(boundary)) next
    restricted <- per[per$model == names[i], ]
    full <- per[per$model == names[j], ]
    lr <- 2 * (full$logLik - restricted$logLik)
    df_each <- full$k - restricted$k
    p_each <- ifelse(df_each > 0, stats::pchisq(pmax(lr, 0), pmax(df_each, 1), lower.tail = FALSE), 1)
    rows[[length(rows) + 1]] <- data.frame(
      restricted = names[i], full = names[j], chisq = sum(lr), df = sum(df_each),
      p = stats::pchisq(sum(lr), sum(df_each), lower.tail = FALSE),
      n_significant = sum(p_each < 0.05),
      n_negative = sum(lr < -1e-6), boundary = boundary, stringsAsFactors = FALSE)
  }
  if (length(rows) == 0) NULL else do.call(rbind, rows)
}
