#' @title Plot Observed Data Against Model Predictions
#'
#' @description
#' Graphical checks of how well the fitted diffusion model describes the data, based on the
#' output of \code{\link{predict.driftd_results}}. All plots work for any number of condition
#' columns and levels, and return \code{ggplot} objects that can be modified further.
#'
#' @param x An object of class \code{driftd_prediction}.
#' @param type The plot type:
#'   \itemize{
#'     \item \code{"rt"} (default): RT quantiles. RT on the x-axis and quantile level on the
#'       y-axis, one column of panels per response. Crosses are observed quantiles averaged
#'       over participants, large points are model predictions, horizontal bars are
#'       prediction intervals for individual participants, and small faint points are
#'       individual participants.
#'     \item \code{"choice"}: probability of one response (see \code{response}) across the
#'       levels of a condition column (\code{x}). Crosses are observed probabilities
#'       averaged over participants, large points and dashed lines are predictions, vertical
#'       bars are prediction intervals for individual participants, and small faint points are
#'       individual participants.
#'     \item \code{"scatter"}: observed against predicted values, one panel per statistic
#'       (P(upper) and the RT quantiles in \code{quantiles} for each response), with the
#'       identity line and prediction intervals. A compact overview of all cells, e.g. for
#'       supplementary material.
#'     \item \code{"participants"}: proportion of each participant's statistics outside their
#'       prediction intervals (see \code{\link{summary.driftd_prediction}}), to spot poorly
#'       fitted participants.
#'     \item \code{"qp"}: quantile-probability plot (Ratcliff & McKoon, 2008): RT quantiles of
#'       each response plotted against its response probability, all cells in one panel. Most
#'       readable for designs with a single condition column that orders difficulty.
#'   }
#'   As in Ratcliff's quantile-probability plots, data are drawn as crosses and model
#'   predictions as circles (joined by lines where there is an order), with prediction
#'   intervals around the predictions.
#' @param participant Optional participant ID: plot only this participant's data and
#'   predictions instead of group averages.
#' @param by Character vector of condition columns to keep when averaging over participants
#'   (see \code{\link{aggregate.driftd_prediction}}). Default: all condition columns.
#'   Statistics are averaged over the condition columns not listed.
#' @param color Condition column shown by color (\code{"rt"}, \code{"choice"}, and as point
#'   color in \code{"scatter"}). By default, \code{"rt"} and \code{"choice"} use the condition
#'   column with the fewest levels (for \code{"choice"}, among those not on the x-axis).
#' @param panels Condition columns that split the plot into panels. By default, all condition
#'   columns in \code{by} not used for \code{x_var} or \code{color}. In \code{"rt"} plots, panels
#'   are arranged in rows, with one column per response.
#' @param x_var For \code{"choice"}: the condition column on the x-axis. Default: the condition
#'   column with the most levels.
#' @param response For \code{"choice"}: the response whose probability is shown,
#'   \code{"upper"} (default) or \code{"lower"}.
#' @param response_labels Optional named character vector of labels for the responses, e.g.
#'   \code{c(upper = "correct", lower = "error")}.
#' @param individuals Logical: show individual participants as faint points in \code{"rt"}
#'   and \code{"choice"} plots. Default TRUE.
#' @param rt_unit Unit for RTs in \code{"rt"} plots: \code{"ms"} (default) or \code{"s"}.
#' @param select Optional named list to restrict the plot to some conditions, e.g.
#'   \code{list(val = "pos", prob = c(0.4, 0.6))}.
#' @param level For \code{"scatter"}: \code{"group"} (default, one point per condition cell,
#'   averaged over participants) or \code{"participant"} (one point per participant and cell).
#' @param quantiles RT quantiles shown in \code{"scatter"}. Default: 0.1, 0.5 and 0.9 (those
#'   of them that were predicted).
#' @param intervals Logical: show prediction intervals (if they were computed). Default TRUE.
#' @param ... Not used.
#'
#' @return A \code{ggplot} object.
#'
#' @details
#' Group-level values are averages of participant-level statistics (Vincentizing). Prediction
#' intervals in group-level plots are intervals for individual participants: about 95% (with
#' the default \code{level}) of the participants' observed values should fall inside them.
#' See \code{\link{aggregate.driftd_prediction}}. RT quantiles of a response are only shown where
#' participants had at least \code{min_n} such responses (see \code{\link{predict.driftd_results}}).
#'
#' Reporting response probabilities and RT quantiles in separate figures keeps each figure
#' simple; \code{"choice"} and \code{"rt"} plots are designed to be used together.
#'
#' @references
#' Ratcliff, R., & McKoon, G. (2008). The diffusion decision model: Theory and data for
#' two-choice decision tasks. \emph{Neural Computation, 20}, 873-922.
#'
#' @seealso \code{\link{predict.driftd_results}}
#'
#' @importFrom rlang .data
#' @export
plot.driftd_prediction <- function(x, type = c("rt", "choice", "scatter", "participants", "qp"),
                                   participant = NULL, by = NULL, color = NULL, panels = NULL, x_var = NULL,
                                   response = c("upper", "lower"), response_labels = NULL,
                                   individuals = TRUE, rt_unit = c("ms", "s"), select = NULL,
                                   level = c("group", "participant"), quantiles = c(0.1, 0.5, 0.9),
                                   intervals = TRUE, ...) {
  type <- match.arg(type)
  level <- match.arg(level)
  response <- match.arg(response)
  rt_unit <- match.arg(rt_unit)
  info <- attr(x, "driftd")
  if (is.null(info)) stop("x must be the output of predict() on a driftd fit")
  intervals <- intervals && isTRUE(info$intervals)
  if (is.null(by)) by <- info$condition_columns
  check_columns(by, info$condition_columns, "by must name condition columns of the fit")
  check_columns(c(color, panels, x_var), by, "x_var, color and panels must be among the condition columns kept in `by`")
  labels <- c(upper = "upper", lower = "lower")
  if (!is.null(response_labels)) labels[names(response_labels)] <- response_labels

  x <- filter_prediction(x, participant, select)
  individuals <- individuals && is.null(participant)

  switch(type,
         rt = plot_rt(x, aggregate(x, by = by), info, by, color, panels, individuals, intervals,
                      rt_unit, labels, participant),
         choice = plot_choice(x, aggregate(x, by = by), info, by, x_var, color, panels, response,
                              individuals, intervals, labels, participant),
         scatter = plot_scatter(if (level == "group") aggregate(x, by = by) else x, info, panels, color,
                                quantiles, intervals, level, participant),
         participants = plot_participants(x, info),
         qp = plot_qp(aggregate(x, by = by), panels, intervals, participant))
}




# ---- Internal plotting helpers --------------------------------------------------------------

# Colors for upper and lower responses (Okabe-Ito, colorblind-safe)
response_colors <- c(upper = "#0072B2", lower = "#D55E00")

#' @noRd
check_columns <- function(cols, allowed, msg) {
  bad <- setdiff(cols, allowed)
  if (length(bad) > 0) stop(msg, " (not: ", paste(bad, collapse = ", "), ")", call. = FALSE)
}

# Restrict predictions to one participant and/or to the conditions in `select` (replicates
# stay matched through the .row column).
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
  x[keep, , drop = FALSE]
}

#' @noRd
title_suffix <- function(participant) {
  if (is.null(participant)) "group averages" else paste("participant", participant)
}

# Okabe-Ito colors (colorblind-safe) for condition levels
level_colors <- c("#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9", "#F0E442", "#000000")

# Assign the condition columns in `by` to the x-axis (choice plots), color and panels.
# Defaults: x = the column with the most levels, color = the remaining column with the
# fewest levels, panels = all others. Every column in `by` gets a role, so no cells overlap.
#' @noRd
assign_roles <- function(g, by, x_var = NULL, color = NULL, panels = NULL, need_x = FALSE) {
  n_levels <- vapply(by, function(col) length(unique(g[[col]])), integer(1))
  if (need_x && is.null(x_var) && length(by) > 0) {
    free <- setdiff(by, c(color, panels))
    if (length(free) > 0) x_var <- free[which.max(n_levels[free])]
  }
  if (is.null(color)) {
    free <- setdiff(by, c(x_var, panels))
    if (length(free) > 0) color <- free[which.min(n_levels[free])]
  }
  if (length(intersect(c(x_var, color), panels)) > 0 || (!is.null(x_var) && identical(x_var, color))) {
    stop("x_var, color and panels must be different condition columns", call. = FALSE)
  }
  panels <- c(panels, setdiff(by, c(x_var, color, panels)))
  list(x = x_var, color = color, panels = panels)
}

# Color scale and legend title for the color role (a single color if there is none).
#' @noRd
color_scale <- function(color, levels) {
  if (is.null(color)) return(ggplot2::scale_colour_manual(values = level_colors[1], guide = "none"))
  values <- if (length(levels) <= length(level_colors)) level_colors[seq_along(levels)] else
    grDevices::hcl.colors(length(levels), "Dark 3")
  ggplot2::scale_colour_manual(values = stats::setNames(values, levels), name = color)
}

#' @noRd
color_column <- function(df, color) {
  if (is.null(color)) factor(rep("all", nrow(df))) else factor(df[[color]], levels = sort(unique(df[[color]])))
}

# Wrap long captions so they fit the plot width
# Caption text for the prediction intervals in group and single-participant plots
#' @noRd
interval_caption <- function(participant, level) {
  paste0(" Bars: ", round(100 * level), "% prediction intervals",
         if (is.null(participant)) " for individual participants." else ".")
}

#' @noRd
wrap_caption <- function(text, width = 90) paste(strwrap(text, width = width), collapse = "\n")

#' @noRd
panel_labeller <- function(df, cols) {
  if (length(cols) == 0) return(rep("", nrow(df)))
  labels <- do.call(paste, c(lapply(cols, function(col) paste(col, "=", df[[col]])), sep = " | "))
  # panels in the order of the condition levels (factor levels, or sorted values)
  ord <- do.call(order, unname(lapply(cols, function(col) df[[col]])))
  factor(labels, levels = unique(labels[ord]))
}

# Observed participant-level values, averaged within each participant over the condition
# columns not in `by` (so there is one point per participant per plotted cell).
#' @noRd
individual_values <- function(x, by, info) {
  df <- as.data.frame(x)
  df <- df[!is.na(df$observed), ]
  keys <- c(info$sub_id, by, "response", "statistic", "quantile")
  if (setequal(by, info$condition_columns)) return(df[, c(keys, "observed")])
  group <- do.call(paste, c(lapply(df[keys], as.character), sep = "\r"))
  res <- df[!duplicated(group), keys, drop = FALSE]
  res$observed <- as.vector(tapply(df$observed, factor(group, levels = unique(group)), mean))
  res
}

#' @noRd
plot_rt <- function(x, g, info, by, color, panels, individuals, intervals, rt_unit, labels, participant) {
  roles <- assign_roles(g, by, color = color, panels = panels)
  scale <- if (rt_unit == "ms") 1000 else 1
  prep <- function(df) {
    df <- as.data.frame(df)
    df <- df[df$statistic == "quantile" & !is.na(df$observed), ]
    df$col <- color_column(df, roles$color)
    df$panel_row <- panel_labeller(df, roles$panels)
    df$response_label <- factor(labels[df$response], levels = labels[c("upper", "lower")])
    # dodge the colors vertically within each quantile row
    k <- nlevels(df$col)
    df$y <- df$quantile + (as.integer(df$col) - (k + 1) / 2) * min(0.05, 0.12 / k)
    df
  }
  gq <- prep(g)
  levels_col <- levels(gq$col)

  p <- ggplot2::ggplot(gq, ggplot2::aes(y = .data$y, colour = .data$col))
  if (individuals) {
    iq <- prep(individual_values(x, by, info))
    iq$col <- factor(as.character(iq$col), levels = levels_col)
    p <- p + ggplot2::geom_point(data = iq, ggplot2::aes(x = .data$observed * scale), alpha = 0.3, size = 1)
  }
  if (intervals) {
    p <- p + ggplot2::geom_segment(ggplot2::aes(x = .data$lower * scale, xend = .data$upper * scale, yend = .data$y),
                                   linewidth = 0.5, na.rm = TRUE)
  }
  p <- p +
    ggplot2::geom_point(ggplot2::aes(x = .data$predicted * scale), size = 2.8) +
    ggplot2::geom_point(ggplot2::aes(x = .data$observed * scale), shape = 4, size = 2.4, stroke = 1.1,
                        colour = "black") +
    color_scale(roles$color, levels_col) +
    ggplot2::scale_y_continuous(breaks = info$probs) +
    ggplot2::labs(x = paste0("RT (", rt_unit, ")"), y = "Quantile",
                  title = paste0("RT quantiles (", title_suffix(participant), ")"),
                  caption = wrap_caption(paste0("Crosses: observed", if (is.null(participant)) " (average over participants)",
                                   ". Large points: model predictions.",
                                   if (intervals) interval_caption(participant, info$level) else "",
                                   if (individuals) " Small points: individual participants." else ""))) +
    ggplot2::theme_bw() +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank())
  if (length(roles$panels) > 0) {
    p + ggplot2::facet_grid(panel_row ~ response_label)
  } else {
    p + ggplot2::facet_grid(. ~ response_label)
  }
}

#' @noRd
plot_choice <- function(x, g, info, by, x_var, color, panels, response, individuals, intervals, labels,
                        participant) {
  roles <- assign_roles(g, by, x_var = x_var, color = color, panels = panels, need_x = TRUE)
  prep <- function(df) {
    df <- as.data.frame(df)
    df <- df[df$statistic == "probability" & df$response == response, ]
    df$col <- color_column(df, roles$color)
    df$xval <- if (is.null(roles$x)) factor(rep("all conditions", nrow(df))) else
      if (is.numeric(df[[roles$x]])) df[[roles$x]] else factor(df[[roles$x]], levels = sort(unique(df[[roles$x]])))
    df$panel <- panel_labeller(df, roles$panels)
    df
  }
  gp <- prep(g)
  levels_col <- levels(gp$col)
  numeric_x <- is.numeric(gp$xval)
  dodge_width <- if (numeric_x && length(unique(gp$xval)) > 1) 0.3 * min(diff(sort(unique(gp$xval)))) else 0.3
  dodge <- ggplot2::position_dodge(width = dodge_width)

  p <- ggplot2::ggplot(gp, ggplot2::aes(x = .data$xval, colour = .data$col, group = .data$col))
  if (individuals) {
    ip <- prep(individual_values(x, by, info))
    ip$col <- factor(as.character(ip$col), levels = levels_col)
    p <- p + ggplot2::geom_point(data = ip, ggplot2::aes(y = .data$observed), position = dodge, alpha = 0.3, size = 1)
  }
  if (intervals) {
    p <- p + ggplot2::geom_errorbar(ggplot2::aes(ymin = .data$lower, ymax = .data$upper), position = dodge,
                                    width = dodge_width / 2, linewidth = 0.5, na.rm = TRUE)
  }
  if (length(unique(gp$xval)) > 1) {
    p <- p + ggplot2::geom_line(ggplot2::aes(y = .data$predicted), position = dodge, linewidth = 0.5, linetype = "dashed")
  }
  p <- p +
    ggplot2::geom_point(ggplot2::aes(y = .data$predicted), position = dodge, size = 2.8) +
    ggplot2::geom_point(ggplot2::aes(y = .data$observed), position = dodge, shape = 4, size = 2.4, stroke = 1.1,
                        colour = "black") +
    color_scale(roles$color, levels_col) +
    ggplot2::scale_y_continuous(limits = c(0, 1)) +
    ggplot2::labs(x = if (is.null(roles$x)) NULL else roles$x, y = paste0("P(", labels[[response]], ")"),
                  title = paste0("Response probabilities (", title_suffix(participant), ")"),
                  caption = wrap_caption(paste0("Crosses: observed", if (is.null(participant)) " (average over participants)",
                                   ". Large points and dashed lines: model predictions.",
                                   if (intervals) interval_caption(participant, info$level) else "",
                                   if (individuals) " Small points: individual participants." else ""))) +
    ggplot2::theme_bw() +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank())
  if (numeric_x) p <- p + ggplot2::scale_x_continuous(breaks = sort(unique(gp$xval)))
  if (length(roles$panels) > 0) p <- p + ggplot2::facet_wrap("panel")
  p
}

#' @noRd
plot_qp <- function(g, panels, intervals, participant) {
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
    ggplot2::geom_point(ggplot2::aes(x = .data$p_pred, y = .data$predicted), shape = 1, colour = "grey30", size = 1.8) +
    ggplot2::geom_point(ggplot2::aes(x = .data$p_obs, y = .data$observed, colour = .data$response), shape = 4,
                        size = 2, stroke = 1.1) +
    ggplot2::scale_colour_manual(values = response_colors, name = "Response") +
    ggplot2::scale_x_continuous(limits = c(0, 1)) +
    ggplot2::labs(x = "Response probability", y = "RT quantile (s)",
                  title = paste0("Quantile-probability plot (", title_suffix(participant), ")"),
                  caption = wrap_caption(paste0("Crosses: observed. Circles and lines: predicted (quantiles ",
                                   paste(levels(qs$quantile_label), collapse = ", "), ")",
                                   if (intervals) ". Bars: prediction intervals." else "."))) +
    ggplot2::theme_bw()
  if (length(panels) > 0) p <- p + ggplot2::facet_wrap(panels, labeller = "label_both")
  p
}

#' @noRd
plot_scatter <- function(df, info, panels, color, quantiles, intervals, level, participant) {
  df <- as.data.frame(df)
  quantiles <- intersect(round(quantiles, 6), round(info$probs, 6))
  keep <- (df$statistic == "probability" & df$response == "upper") |
    (df$statistic == "quantile" & round(df$quantile, 6) %in% quantiles)
  df <- df[keep & !is.na(df$observed), ]
  df$panel <- ifelse(df$statistic == "probability", "P(upper)",
                     paste0(df$response, ": RT quantile ", format(df$quantile)))
  df$panel <- factor(df$panel, levels = c("P(upper)", paste0(rep(c("upper", "lower"), each = length(quantiles)),
                                                              ": RT quantile ", format(rep(quantiles, 2)))))

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
    p <- p + ggplot2::geom_point(data = df, ggplot2::aes(colour = .data[[color]]), size = 2)
  } else {
    p <- p + ggplot2::geom_point(colour = "grey30", size = 2)
  }
  wrap <- if (length(panels) > 0) ggplot2::facet_wrap(c("panel", panels), scales = "free") else
    ggplot2::facet_wrap("panel", scales = "free")
  what <- if (!is.null(participant)) paste("participant", participant) else
    if (level == "group") "group averages per condition cell" else "one point per participant and condition cell"
  p + wrap +
    ggplot2::labs(x = "Predicted", y = "Observed", title = paste0("Observed vs predicted (", what, ")"),
                  caption = paste0("Dashed line: observed = predicted.",
                                   if (intervals) interval_caption(if (level == "group") participant else "", info$level)
                                   else "")) +
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
