#' @title Parameter Estimates in Analysis-Ready Format
#'
#' @description
#' Returns the individual parameter estimates of a \code{\link{driftd}} fit as a data frame
#' for further analyses (e.g., t-tests, ANOVAs or mixed models on the parameters), reports,
#' or custom plots. Condition values are taken from the original condition columns, so the
#' result works for any design and any level names.
#'
#' @param fit An object of class \code{driftd_results}, the output of \code{\link{driftd}}.
#' @param format The layout of the result:
#'   \itemize{
#'     \item \code{"long"} (default): one row per participant and parameter (e.g. each
#'       condition-specific drift rate is one row), with the parameter type, the condition
#'       values it applies to, its estimate and standard error.
#'     \item \code{"cells"}: one row per participant and condition cell, with the condition
#'       columns and one column per parameter type (\code{a}, \code{t0}, \code{z}, \code{v},
#'       \code{st0}, \code{sz}, \code{sv}) holding the value that applies in that cell. Ready
#'       for e.g. \code{aov(v ~ prob * val + Error(ID / (prob * val)), data = ...)}.
#'   }
#' @param include_fixed Logical: include parameters that were fixed rather than estimated
#'   (\code{"long"} format only). Default TRUE; they are marked in the \code{fixed} column.
#'
#' @return A data frame.
#'
#'   For \code{format = "long"}: the participant ID column, \code{parameter} (the parameter
#'   type, e.g. \code{"v"}), one column per condition column (the condition level the
#'   parameter applies to, or \code{NA} if the parameter does not vary across that column),
#'   \code{name} (the column name in \code{fit$individual_estimates}, e.g.
#'   \code{"v_0.4_neg"}), \code{estimate}, \code{se} (\code{NA} if not available, see
#'   \code{\link{driftd}}), \code{fixed}, and \code{converge} (the participant's
#'   convergence code).
#'
#'   For \code{format = "cells"}: the participant ID column, the condition columns,
#'   \code{n_trials} (the participant's number of trials in the cell), the seven parameters,
#'   and \code{converge}. Parameters that do not vary across conditions repeat the same
#'   value in every cell of a participant.
#'
#' @details
#' \code{z} is the relative starting point (between 0 and 1), and \code{t0}, \code{st0} are
#' in seconds, as in \code{\link{driftd}}.
#'
#' In a two-stage analysis, the estimates are treated as data in a second analysis across
#' participants. Participants whose fit did not converge (\code{converge != 0}) should
#' usually be inspected or excluded first.
#'
#' @examples
#' \dontrun{
#' fit <- driftd(speed_acc, sub_id = "id", vary_a = "condition", vary_v = "stimulus",
#'               fix_st0 = FALSE)
#' est <- estimates(fit)                     # one row per participant and parameter
#' cells <- estimates(fit, format = "cells") # one row per participant and condition
#' aggregate(cbind(a, v) ~ condition + stimulus, data = cells, FUN = mean)
#' }
#'
#' @seealso \code{\link{driftd}}
#' @export
estimates <- function(fit, format = c("long", "cells"), include_fixed = TRUE) {
  format <- match.arg(format)
  cols <- fit_columns(fit)
  if (is.null(fit$settings$parameters)) {
    stop("This driftd_results object was created by an older version of driftd. Please refit the model.")
  }
  id <- cols$sub_id
  est <- fit$individual_estimates
  conditions <- fit$settings$condition_columns
  if (format == "long") estimates_long(fit, est, id, conditions, include_fixed) else
    estimates_cells(fit, est, id, conditions)
}




# ---- Internal helpers ---------------------------------------------------------------------

#' @noRd
estimates_long <- function(fit, est, id, conditions, include_fixed) {
  info <- fit$settings$parameters
  if (!include_fixed) info <- info[!info$fixed, , drop = FALSE]
  se_tab <- fit$standard_errors
  n_par <- nrow(info)
  n_id <- nrow(est)

  rows <- rep(seq_len(n_par), times = n_id)
  res <- data.frame(ID = rep(est[[id]], each = n_par), parameter = info$parameter[rows], stringsAsFactors = FALSE)
  for (column in conditions) res[[column]] <- info[[column]][rows]
  res$name <- info$name[rows]
  res$estimate <- as.vector(t(as.matrix(est[, info$name, drop = FALSE])))
  res$se <- vapply(seq_len(nrow(res)), function(k) {
    name <- res$name[k]
    if (is.null(se_tab) || !(name %in% names(se_tab))) return(NA_real_)
    se_tab[[name]][match(res$ID[k], se_tab[[id]])]
  }, numeric(1))
  res$fixed <- info$fixed[rows]
  res$converge <- rep(est$converge, each = n_par)
  names(res)[1] <- id
  rownames(res) <- NULL
  res
}

#' @noRd
estimates_cells <- function(fit, est, id, conditions) {
  d <- fit$data
  cells <- fit$settings$cells
  present <- unique(d[, c(id, ".cell")])
  counts <- table(factor(paste(d[[id]], d$.cell, sep = "\r"), levels = paste(present[[id]], present$.cell, sep = "\r")))
  types <- c("a", "t0", "z", "v", "st0", "sz", "sv")

  values <- t(vapply(seq_len(nrow(present)), function(k) {
    row <- est[est[[id]] == present[[id]][k], , drop = FALSE]
    unlist(get_parameters_for_condition(row, present$.cell[k], fit$settings))[types]
  }, numeric(length(types))))

  res <- data.frame(ID = present[[id]], stringsAsFactors = FALSE)
  cell_rows <- match(present$.cell, cells$.cell)
  for (column in conditions) res[[column]] <- cells[[column]][cell_rows]
  res$n_trials <- as.integer(counts)
  res <- cbind(res, as.data.frame(values))
  res$converge <- est$converge[match(present[[id]], est[[id]])]

  # order by participant (as in the estimates) and then by the condition values
  res <- res[do.call(order, c(list(match(res$ID, est[[id]])), lapply(conditions, function(col) res[[col]]))), ]
  names(res)[1] <- id
  rownames(res) <- NULL
  res
}
