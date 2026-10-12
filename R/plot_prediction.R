#' @title Plot Observed Data Against Model Predictions
#'
#' @description
#' Graphical checks of how well the fitted diffusion model describes the data, based on the
#' output of \code{\link{predict.driftd_results}}. All plots work for any number of condition
#' columns, levels and cells, and return \code{ggplot} objects that can be modified further.
#'
#' @param x An object of class \code{driftd_prediction}.
#' @param type The plot type:
#'   \itemize{
#'     \item \code{"qp"} (default): quantile-probability plot (Ratcliff & McKoon, 2008). Each
#'       condition cell contributes its RT quantiles of upper responses at x = P(upper) and of
#'       lower responses at x = P(lower), so all cells and both responses appear in one panel.
#'       Observed values are points, predictions are crosses connected across cells for each
#'       quantile, and bars show prediction intervals.
#'     \item \code{"scatter"}: observed against predicted values, one panel per statistic
#'       (P(upper) and the RT quantiles in \code{quantiles} for each response), with the
#'       identity line. Points outside their prediction interval are highlighted.
#'     \item \code{"cdf"}: defective cumulative distributions per cell: the RT quantiles of each
#'       response plotted against cumulative probability, scaled by the response probability,
#'       so both how often and how fast each response is given are shown.
#'     \item \code{"participants"}: proportion of each participant's statistics outside their
#'       prediction intervals (see \code{\link{summary.driftd_prediction}}), to spot poorly
#'       fitted participants.
#'   }
#' @param participant Optional participant ID: plot only this participant's data and
#'   predictions instead of group averages.
#' @param by Character vector of condition columns to keep when averaging over participants
#'   (see \code{\link{aggregate.driftd_prediction}}). Default: all condition columns.
#' @param facet Optional character vector of condition columns (a subset of \code{by}) to split
#'   the plot into panels (\code{"qp"} and \code{"scatter"}).
#' @param color Optional condition column used to color points (\code{"scatter"}).
#' @param select Optional named list to restrict the plot to some conditions, e.g.
#'   \code{list(val = "pos", prob = c(0.4, 0.6))}.
#' @param level For \code{"scatter"}: \code{"group"} (default, one point per condition cell,
#'   averaged over participants) or \code{"participant"} (one point per participant and cell).
#' @param quantiles RT quantiles shown in \code{"scatter"}. Default: 0.1, 0.5 and 0.9 (those
#'   of them that were predicted).
#' @param max_cells Maximum number of cells (panels) shown by \code{"cdf"}. Default 16; use
#'   \code{select} or \code{by} to choose which cells to show.
#' @param intervals Logical: show prediction intervals (if they were computed). Default TRUE.
#' @param ... Not used.
#'
#' @return A \code{ggplot} object.
#'
#' @details
#' Group-level plots show averages of participant-level statistics (Vincentizing); see
#' \code{\link{aggregate.driftd_prediction}}. RT quantiles of a response are only shown where
#' participants had at least \code{min_n} such responses (see \code{\link{predict.driftd_results}}).
#'
#' @references
#' Ratcliff, R., & McKoon, G. (2008). The diffusion decision model: Theory and data for
#' two-choice decision tasks. \emph{Neural Computation, 20}, 873-922.
#'
#' @seealso \code{\link{predict.driftd_results}}
#'
#' @importFrom rlang .data
#' @export
plot.driftd_prediction <- function(x, type = c("qp", "scatter", "cdf", "participants"), participant = NULL,
                                   by = NULL, facet = NULL, color = NULL, select = NULL,
                                   level = c("group", "participant"), quantiles = c(0.1, 0.5, 0.9),
                                   max_cells = 16, intervals = TRUE, ...) {
  type <- match.arg(type)
  level <- match.arg(level)
  info <- attr(x, "driftd")
  intervals <- intervals && isTRUE(info$intervals)
  if (is.null(by)) by <- info$condition_columns
  check_columns(c(facet, color), by, "facet and color must be among the condition columns kept in `by`")

  x <- filter_prediction(x, participant, select)

  switch(type,
         qp = plot_qp(aggregate(x, by = by), facet, intervals, participant),
         scatter = plot_scatter(if (level == "group") aggregate(x, by = by) else x, info, facet, color,
                                quantiles, intervals, level, participant),
         cdf = plot_cdf(aggregate(x, by = by), by, max_cells, participant),
         participants = plot_participants(x, info))
}




# ---- Internal plotting helpers --------------------------------------------------------------

# Colors for upper and lower responses (Okabe-Ito, colorblind-safe)
response_colors <- c(upper = "#0072B2", lower = "#D55E00")

#' @noRd
check_columns <- function(cols, allowed, msg) {
  bad <- setdiff(cols, allowed)
  if (length(bad) > 0) stop(msg, " (not: ", paste(bad, collapse = ", "), ")", call. = FALSE)
}

# Restrict predictions to one participant and/or to the conditions in `select`, keeping the
# attributes (and the matching rows of the replicate matrix).
#' @noRd
filter_prediction <- function(x, participant, select) {
  info <- attr(x, "driftd")
  keep <- rep(TRUE, nrow(x))
  if (!is.null(participant)) {
    if (!any(x[[info$sub_id]] == participant)) stop("No participant with ID ", participant, call. = FALSE)
    keep <- keep & x[[info$sub_id]] == participant
  }
  if (!is.null(select)) {
    check_columns(names(select), info$condition_columns, "select must name condition columns of the fit")
    for (col in names(select)) keep <- keep & x[[col]] %in% select[[col]]
    if (!any(keep)) stop("No data left after `select`", call. = FALSE)
  }
  reps <- attr(x, "replicates")
  out <- x[keep, , drop = FALSE]
  attr(out, "replicates") <- if (is.null(reps)) NULL else reps[keep, , drop = FALSE]
  attr(out, "driftd") <- info
  class(out) <- class(x)
  out
}

# Readable cell labels such as "val = pos, prob = 0.4"
#' @noRd
cell_labels <- function(df, by) {
  if (length(by) == 0) return(rep("all conditions", nrow(df)))
  do.call(paste, c(lapply(by, function(col) paste(col, "=", df[[col]])), sep = ", "))
}

#' @noRd
title_suffix <- function(participant) {
  if (is.null(participant)) "group averages" else paste("participant", participant)
}

#' @noRd
plot_qp <- function(g, facet, intervals, participant) {
  key <- c(setdiff(names(g), c("statistic", "quantile", "n_participants", "observed", "predicted", "lower", "upper")))
  probs <- g[g$statistic == "probability", c(key, "observed", "predicted")]
  names(probs)[names(probs) %in% c("observed", "predicted")] <- c("p_obs", "p_pred")
  qs <- merge(g[g$statistic == "quantile", ], probs, by = key)
  qs$quantile_label <- factor(format(qs$quantile))

  p <- ggplot2::ggplot(qs) +
    ggplot2::geom_line(ggplot2::aes(x = .data$p_pred, y = .data$predicted, group = .data$quantile_label),
                       colour = "grey60", linewidth = 0.4)
  if (intervals) {
    p <- p + ggplot2::geom_errorbar(ggplot2::aes(x = .data$p_pred, ymin = .data$lower, ymax = .data$upper),
                                    colour = "grey60", width = 0.015, linewidth = 0.3, na.rm = TRUE)
  }
  p <- p +
    ggplot2::geom_point(ggplot2::aes(x = .data$p_pred, y = .data$predicted), shape = 4, colour = "grey30", size = 1.8) +
    ggplot2::geom_point(ggplot2::aes(x = .data$p_obs, y = .data$observed, colour = .data$response), size = 2) +
    ggplot2::scale_colour_manual(values = response_colors, name = "Response") +
    ggplot2::scale_x_continuous(limits = c(0, 1)) +
    ggplot2::labs(x = "Response probability", y = "RT quantile (s)",
                  title = paste0("Quantile-probability plot (", title_suffix(participant), ")"),
                  caption = paste0("Points: observed. Crosses and lines: predicted (quantiles ",
                                   paste(levels(qs$quantile_label), collapse = ", "), ")",
                                   if (intervals) ". Bars: prediction intervals." else ".")) +
    ggplot2::theme_bw()
  if (length(facet) > 0) p <- p + ggplot2::facet_wrap(facet, labeller = "label_both")
  p
}

#' @noRd
plot_scatter <- function(df, info, facet, color, quantiles, intervals, level, participant) {
  df <- as.data.frame(df)
  quantiles <- intersect(round(quantiles, 6), round(info$probs, 6))
  keep <- (df$statistic == "probability" & df$response == "upper") |
    (df$statistic == "quantile" & round(df$quantile, 6) %in% quantiles)
  df <- df[keep & !is.na(df$observed), ]
  df$panel <- ifelse(df$statistic == "probability", "P(upper)",
                     paste0(df$response, ": RT quantile ", format(df$quantile)))
  df$panel <- factor(df$panel, levels = c("P(upper)", paste0(rep(c("upper", "lower"), each = length(quantiles)),
                                                              ": RT quantile ", format(rep(quantiles, 2)))))
  if (intervals) {
    df$fit <- ifelse(is.na(df$lower), "no interval",
                     ifelse(df$observed < df$lower | df$observed > df$upper, "outside interval", "inside interval"))
  }

  # invisible points on the diagonal so that every panel has equal x and y ranges
  rng <- do.call(rbind, lapply(split(df, df$panel, drop = TRUE), function(s) {
    r <- range(c(s$observed, s$predicted, if (intervals) c(s$lower, s$upper)), na.rm = TRUE)
    data.frame(panel = s$panel[1], v = r)
  }))

  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$predicted, y = .data$observed)) +
    ggplot2::geom_blank(data = rng, ggplot2::aes(x = .data$v, y = .data$v)) +
    ggplot2::geom_abline(intercept = 0, slope = 1, colour = "grey50", linetype = "dashed")
  if (intervals) {
    p <- p + ggplot2::geom_errorbar(ggplot2::aes(ymin = .data$lower, ymax = .data$upper),
                                    colour = "grey70", width = 0, na.rm = TRUE)
  }
  if (!is.null(color)) {
    df[[color]] <- factor(df[[color]])
    if (intervals) {
      p <- p + ggplot2::geom_point(data = df, ggplot2::aes(colour = .data[[color]], shape = .data$fit), size = 2) +
        ggplot2::scale_shape_manual(values = c("inside interval" = 16, "outside interval" = 4, "no interval" = 1),
                                    name = NULL)
    } else {
      p <- p + ggplot2::geom_point(data = df, ggplot2::aes(colour = .data[[color]]), size = 2)
    }
  } else if (intervals) {
    p <- p + ggplot2::geom_point(ggplot2::aes(colour = .data$fit), size = 2) +
      ggplot2::scale_colour_manual(values = c("inside interval" = "grey30", "outside interval" = "#D55E00",
                                              "no interval" = "grey70"), name = NULL)
  } else {
    p <- p + ggplot2::geom_point(colour = "grey30", size = 2)
  }
  wrap <- if (length(facet) > 0) ggplot2::facet_wrap(c("panel", facet), scales = "free") else
    ggplot2::facet_wrap("panel", scales = "free")
  what <- if (!is.null(participant)) paste("participant", participant) else
    if (level == "group") "group averages per condition cell" else "one point per participant and condition cell"
  p + wrap +
    ggplot2::labs(x = "Predicted", y = "Observed", title = paste0("Observed vs predicted (", what, ")"),
                  caption = paste0("Dashed line: observed = predicted.",
                                   if (intervals) " Bars: prediction intervals." else "")) +
    ggplot2::theme_bw()
}

#' @noRd
plot_cdf <- function(g, by, max_cells, participant) {
  g <- as.data.frame(g)
  g$cell <- cell_labels(g, by)
  cells <- unique(g$cell)
  if (length(cells) > max_cells) {
    message("Showing the first ", max_cells, " of ", length(cells), " cells. Use `select` or `by` ",
            "to choose cells, or increase `max_cells`.")
    cells <- cells[seq_len(max_cells)]
    g <- g[g$cell %in% cells, ]
  }
  g$cell <- factor(g$cell, levels = cells)
  probs <- g[g$statistic == "probability", c("cell", "response", "observed", "predicted")]
  names(probs)[3:4] <- c("p_obs", "p_pred")
  qs <- merge(g[g$statistic == "quantile", ], probs, by = c("cell", "response"))
  qs$cum_obs <- qs$quantile * qs$p_obs
  qs$cum_pred <- qs$quantile * qs$p_pred
  qs <- qs[order(qs$cell, qs$response, qs$quantile), ]

  ggplot2::ggplot(qs, ggplot2::aes(colour = .data$response)) +
    ggplot2::geom_line(ggplot2::aes(x = .data$predicted, y = .data$cum_pred, group = .data$response), linewidth = 0.6) +
    ggplot2::geom_point(ggplot2::aes(x = .data$observed, y = .data$cum_obs), size = 1.8) +
    ggplot2::facet_wrap("cell") +
    ggplot2::scale_colour_manual(values = response_colors, name = "Response") +
    ggplot2::labs(x = "RT (s)", y = "Cumulative probability",
                  title = paste0("Defective cumulative distributions (", title_suffix(participant), ")"),
                  caption = paste("Points: observed RT quantiles. Lines: predicted.",
                                  "Each curve ends at the probability of that response.")) +
    ggplot2::theme_bw()
}

#' @noRd
plot_participants <- function(x, info) {
  s <- summary(x)
  id <- info$sub_id
  s$label <- factor(as.character(s[[id]]), levels = rev(as.character(s[[id]])))
  ggplot2::ggplot(s, ggplot2::aes(x = .data$prop_outside, y = .data$label)) +
    ggplot2::geom_vline(xintercept = 1 - info$level, colour = "grey50", linetype = "dashed") +
    ggplot2::geom_point(size = 2) +
    ggplot2::expand_limits(x = 0) +
    ggplot2::labs(x = "Proportion of statistics outside the prediction interval", y = id,
                  title = "Model misfit per participant",
                  caption = paste0("Dashed line: proportion expected by chance (", 1 - info$level,
                                   ") if the model describes the data.")) +
    ggplot2::theme_bw()
}
