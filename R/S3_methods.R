
#' @title Print Diffusion Model Fitting Results
#'
#' @description
#' Provides a brief summary when displaying an object of class \code{driftd_results}.
#'
#' @param x An object of class \code{driftd_results}, typically the output
#'   from the \code{\link{driftd}} function.
#' @param ... Further arguments passed to or from other methods (e.g., passed
#'   to \code{print.data.frame} for the head of parameter estimates).
#'
#' @return Invisibly returns the input object \code{x}. Called for its side
#'   effect of printing a summary to the console.
#'
#' @details
#' This method prints key summary information, including the number of
#' participants, number of conditions, and the head (first 6 rows) of the
#' individual parameter estimates table.
#'
#' @seealso
#' \code{\link{driftd}} for the main fitting function,
#' \code{\link{summary.driftd_results}} for a more detailed output, and
#' \code{\link{plot.driftd_results}} for visualization.
#' @importFrom utils head
#'
#' @export
print.driftd_results <- function(x, ...) {
  cat("Diffusion Model Results\n")
  cat("=======================\n")
  cat("Participants:", x$settings$n_participants, "\n")
  cat("Conditions:", x$settings$n_conditions, "\n")
  cat("\nParameter Estimates (first 6 participants):\n")
  print(head(x$individual_estimates), ...)
  cat("\nUse summary() for detailed results or plot() for visualization.\n")
  invisible(x)
}




#' @title Summarize Diffusion Model Fitting Results
#'
#' @description
#' Generates a detailed summary of the fitting process and parameter estimates
#' for an object of class \code{driftd_results}.
#'
#' @param object An object of class \code{driftd_results}, typically the output
#'   from the \code{\link{driftd}} function.
#' @param ... Further arguments passed to or from other methods.
#'
#' @return A \code{data.frame} or \code{tibble} containing the aggregated
#'   **Mean** and **Standard Deviation (SD)** for all estimated and fixed
#'   diffusion model parameters across participants. The function also prints
#'   convergence information to the console.
#'
#' @details
#' This method is the primary way to inspect the overall performance of the
#' model fit. It reports the total number of participants and conditions,
#' a table of optimizer convergence codes, and the summary table of group-level
#' parameter statistics (Mean and SD).
#'
#' @seealso
#' \code{\link{driftd}} for the main fitting function,
#' \code{\link{print.driftd_results}} for a brief console output.
#'
#' @export
summary.driftd_results <- function(object, ...) {
  cat("Diffusion Model Summary\n")
  cat("=======================\n")
  cat("Participants:", object$settings$n_participants, "\n")
  cat("Conditions:", object$settings$n_conditions, "\n")
  cat("\nConvergence:\n")
  print(object$convergence)
  cat("\nParameter Summary:\n")
  return(object$summary)  # Your existing summary table
}




#' @title Visualize Diffusion Model Fitting Results
#'
#' @description
#' Plots various diagnostic and summary visualizations for an object of class
#' \code{driftd_results}.
#'
#' @param x An object of class \code{driftd_results}, typically the output
#'   from the \code{\link{driftd}} function.
#' @param type A \code{character} string specifying the type of plot to generate.
#'   Must be one of:
#'   \itemize{
#'     \item \code{"parameters"} (default): Generates a 2x2 panel plot showing
#'       the distribution of the individual participant estimates for the main
#'       DDM parameters (\code{a}, \code{t0}, \code{z}, \code{v}).
#'     \item \code{"convergence"}: Generates a bar plot summarizing the optimizer
#'       convergence codes across all participants.
#'   }
#' @param boxplot_width A numeric value (default is \code{0.8}) used to control
#'   the width of the boxplots in the parameter visualization.
#' @param ... Further graphical parameters passed to the underlying plotting
#'   functions (e.g., \code{\link{boxplot}} or \code{\link{barplot}}).
#'
#' @return The function is called for its side effect of creating a plot. It
#'   invisibly returns \code{NULL}.
#'
#' @details
#' The default plot (\code{type = "parameters"}) displays **Boxplots** and
#' **Jittered Data Points** for the boundary separation (\code{a}), non-decision
#' time (\code{t0}), starting point bias (\code{z}), and drift rate (\code{v}).
#' The plots are adjusted to display parameter estimates across experimental
#' conditions if the parameter was fit as condition-dependent. For the starting
#' point (\code{z}) and drift rate (\code{v}) plots, a reference line is added
#' at \code{0.5} and \code{0}, respectively, for context.
#'
#' The \code{type = "convergence"} plot helps assess the quality of the fit by
#' showing how many participants reached successful convergence (code 0).
#'
#' @seealso
#' \code{\link{driftd}} for the main fitting function,
#' \code{\link{print.driftd_results}} and \code{\link{summary.driftd_results}}
#' for textual output.
#'
#' @importFrom graphics par title boxplot stripchart abline barplot text plot
#' @export
plot.driftd_results <- function(x, type = "parameters", boxplot_width = 0.8, ...) {
  if (type == "parameters") {
    # Extract parameter estimates and settings
    param_data <- x$individual_estimates
    settings <- x$settings

    # Set up 2x2 layout
    par(mfrow = c(2, 2), mar = c(5, 4, 3, 2), oma = c(0, 0, 2, 0))

    # Get parameter columns (excluding ID, converge, LLE, AIC, BIC)
    param_cols <- names(param_data)[!names(param_data) %in% c("ID", "converge", "LLE", "AIC", "BIC")]

    # Function to extract condition from parameter name
    extract_condition <- function(param_name) {
      if (grepl("_", param_name)) {
        return(sub("^[^_]*_", "", param_name))
      } else {
        return("Overall")
      }
    }

    # Function to extract parameter type
    extract_param_type <- function(param_name) {
      return(sub("_.*", "", param_name))
    }

    # Prepare data for plotting
    plot_data <- data.frame()
    for (col in param_cols) {
      param_type <- extract_param_type(col)
      condition <- extract_condition(col)

      # Skip if not one of our target parameters
      if (!param_type %in% c("a", "t0", "z", "v")) next

      values <- param_data[[col]]
      plot_data <- rbind(plot_data,
                         data.frame(Parameter = param_type,
                                    Condition = condition,
                                    Value = values))
    }

    # Calculate actual ranges to ensure we don't cut off data
    a_range <- range(plot_data$Value[plot_data$Parameter == "a"], na.rm = TRUE)
    t0_range <- range(plot_data$Value[plot_data$Parameter == "t0"], na.rm = TRUE)
    z_range <- range(plot_data$Value[plot_data$Parameter == "z"], na.rm = TRUE)
    v_range <- range(plot_data$Value[plot_data$Parameter == "v"], na.rm = TRUE)

    # Set y-axis limits with safety margins
    a_ylim <- c(0, max(3, a_range[2] * 1.1))
    t0_ylim <- c(0, max(1, t0_range[2] * 1.1))
    z_ylim <- c(0, 1)
    v_ylim <- c(min(-3, v_range[1] * 1.1), max(3, v_range[2] * 1.1))

    # Adjust boxplot width based on number of conditions
    # Fewer conditions = narrower boxes, more conditions = wider boxes
    get_boxplot_width <- function(data, base_width = boxplot_width) {
      n_conditions <- length(unique(data$Condition))
      if (n_conditions <= 2) {
        return(base_width * 0.4)  # Narrower for 1-2 conditions
      } else if (n_conditions <= 4) {
        return(base_width * 0.6)  # Slightly narrower for 3-4 conditions
      } else {
        return(base_width)        # Default for 5+ conditions
      }
    }

    # Panel 1: Parameter a
    a_data <- plot_data[plot_data$Parameter == "a", ]
    if (nrow(a_data) > 0) {
      # Check if a varies by conditions
      vary_a <- settings$varying_parameters$a
      boxwex_a <- get_boxplot_width(a_data)

      if (length(vary_a) == 0) {
        # a doesn't vary - use parameter name on x-axis
        boxplot(Value ~ Parameter, data = a_data,
                main = "Boundary Separation (a)",
                xlab = "", ylab = "Value",
                ylim = a_ylim, boxwex = boxwex_a,
                col = "lightblue", ...)
      } else {
        # a varies - use condition names on x-axis
        boxplot(Value ~ Condition, data = a_data,
                main = "Boundary Separation (a)",
                xlab = "Condition", ylab = "Value",
                ylim = a_ylim, boxwex = boxwex_a,
                col = "lightblue", ...)
      }
      # Add points for individual participants
      if (length(vary_a) == 0) {
        stripchart(Value ~ Parameter, data = a_data,
                   vertical = TRUE, method = "jitter",
                   pch = 16, col = rgb(0, 0, 0, 0.3), add = TRUE)
      } else {
        stripchart(Value ~ Condition, data = a_data,
                   vertical = TRUE, method = "jitter",
                   pch = 16, col = rgb(0, 0, 0, 0.3), add = TRUE)
      }
    } else {
      plot(0, type = 'n', xaxt = 'n', yaxt = 'n', ann = FALSE)
      text(1, 0, "No a parameter estimates", cex = 1.2)
    }

    # Panel 2: Parameter t0
    t0_data <- plot_data[plot_data$Parameter == "t0", ]
    if (nrow(t0_data) > 0) {
      # Check if t0 varies by conditions
      vary_t0 <- settings$varying_parameters$t0
      boxwex_t0 <- get_boxplot_width(t0_data)

      if (length(vary_t0) == 0) {
        # t0 doesn't vary - use parameter name on x-axis
        boxplot(Value ~ Parameter, data = t0_data,
                main = "Non-decision Time (t0)",
                xlab = "", ylab = "Value",
                ylim = t0_ylim, boxwex = boxwex_t0,
                col = "lightgreen", ...)
      } else {
        # t0 varies - use condition names on x-axis
        boxplot(Value ~ Condition, data = t0_data,
                main = "Non-decision Time (t0)",
                xlab = "Condition", ylab = "Value",
                ylim = t0_ylim, boxwex = boxwex_t0,
                col = "lightgreen", ...)
      }
      if (length(vary_t0) == 0) {
        stripchart(Value ~ Parameter, data = t0_data,
                   vertical = TRUE, method = "jitter",
                   pch = 16, col = rgb(0, 0, 0, 0.3), add = TRUE)
      } else {
        stripchart(Value ~ Condition, data = t0_data,
                   vertical = TRUE, method = "jitter",
                   pch = 16, col = rgb(0, 0, 0, 0.3), add = TRUE)
      }
    } else {
      plot(0, type = 'n', xaxt = 'n', yaxt = 'n', ann = FALSE)
      text(1, 0, "No t0 parameter estimates", cex = 1.2)
    }

    # Panel 3: Parameter z
    z_data <- plot_data[plot_data$Parameter == "z", ]
    if (nrow(z_data) > 0) {
      # Check if z varies by conditions
      vary_z <- settings$varying_parameters$z
      boxwex_z <- get_boxplot_width(z_data)

      if (length(vary_z) == 0) {
        # z doesn't vary - use parameter name on x-axis
        boxplot(Value ~ Parameter, data = z_data,
                main = "Starting Point (z)",
                xlab = "", ylab = "Value",
                ylim = z_ylim, boxwex = boxwex_z,
                col = "lightyellow", ...)
      } else {
        # z varies - use condition names on x-axis
        boxplot(Value ~ Condition, data = z_data,
                main = "Starting Point (z)",
                xlab = "Condition", ylab = "Value",
                ylim = z_ylim, boxwex = boxwex_z,
                col = "lightyellow", ...)
      }
      abline(h = 0.5, lty = 2, col = "red", lwd = 2)
      if (length(vary_z) == 0) {
        stripchart(Value ~ Parameter, data = z_data,
                   vertical = TRUE, method = "jitter",
                   pch = 16, col = rgb(0, 0, 0, 0.3), add = TRUE)
      } else {
        stripchart(Value ~ Condition, data = z_data,
                   vertical = TRUE, method = "jitter",
                   pch = 16, col = rgb(0, 0, 0, 0.3), add = TRUE)
      }
    } else {
      plot(0, type = 'n', xaxt = 'n', yaxt = 'n', ann = FALSE)
      text(1, 0, "No z parameter estimates", cex = 1.2)
    }

    # Panel 4: Parameter v
    v_data <- plot_data[plot_data$Parameter == "v", ]
    if (nrow(v_data) > 0) {
      # Check if v varies by conditions
      vary_v <- settings$varying_parameters$v
      boxwex_v <- get_boxplot_width(v_data)

      if (length(vary_v) == 0) {
        # v doesn't vary - use parameter name on x-axis
        boxplot(Value ~ Parameter, data = v_data,
                main = "Drift Rate (v)",
                xlab = "", ylab = "Value",
                ylim = v_ylim, boxwex = boxwex_v,
                col = "lightcoral", ...)
      } else {
        # v varies - use condition names on x-axis
        boxplot(Value ~ Condition, data = v_data,
                main = "Drift Rate (v)",
                xlab = "Condition", ylab = "Value",
                ylim = v_ylim, boxwex = boxwex_v,
                col = "lightcoral", ...)
      }
      abline(h = 0, lty = 2, col = "red", lwd = 2)
      if (length(vary_v) == 0) {
        stripchart(Value ~ Parameter, data = v_data,
                   vertical = TRUE, method = "jitter",
                   pch = 16, col = rgb(0, 0, 0, 0.3), add = TRUE)
      } else {
        stripchart(Value ~ Condition, data = v_data,
                   vertical = TRUE, method = "jitter",
                   pch = 16, col = rgb(0, 0, 0, 0.3), add = TRUE)
      }
    } else {
      plot(0, type = 'n', xaxt = 'n', yaxt = 'n', ann = FALSE)
      text(1, 0, "No v parameter estimates", cex = 1.2)
    }

    # Add overall title
    title(main = "Diffusion Model Parameter Estimates", outer = TRUE)

    # Reset layout
    par(mfrow = c(1, 1))

  } else if (type == "convergence") {
    # Existing convergence plot code
    barplot(x$convergence, main = "Model Convergence",
            xlab = "Convergence Code", ylab = "Number of Participants")
  } else {
    stop("Unknown plot type. Use 'parameters' or 'convergence'")
  }
}




#' @title Plot Posterior Predictive Checks (PPC)
#'
#' @description
#' Visualizes the discrepancy between observed data and model-simulated data
#' (posterior predictive checks) contained within a \code{driftd_ppc} object.
#'
#' @param x An object of class \code{driftd_ppc}, typically the output
#'   from the \code{\link{posterior_predict}} function.
#' @param type A \code{character} string specifying the type of PPC visualization.
#'   Must be one of:
#'   \itemize{
#'     \item \code{"summary"} (default): Plots aggregated fit statistics (e.g., mean RT, accuracy)
#'       comparing observed vs. simulated data across all participants, split by experimental cells.
#'     \item \code{"participant"}: Plots the observed vs. simulated fit for each
#'       individual participant (requires the helper function \code{plot_ppc_participant}).
#'     \item \code{"density"}: Plots the full reaction time (RT) density distributions
#'       for observed and simulated data, possibly filtered by condition.
#'   }
#' @param condition A \code{character} string or \code{NULL} (default). If
#'   \code{type = "density"} is selected, this can specify a particular
#'   experimental condition to plot. If \code{NULL}, the density is plotted
#'   for all aggregated data.
#' @param ... Further graphical parameters passed to the underlying plotting
#'   functions.
#'
#' @return The function is called for its side effect of creating a plot or plots.
#'   It invisibly returns \code{NULL}.
#'
#' @details
#' The primary goal of this function is to assess the **goodness-of-fit** of the
#' diffusion model by visually comparing the raw data to data generated using the
#' individual parameter estimates.
#'
#' The three plot types offer different levels of detail:
#' \itemize{
#'   \item \strong{Summary:} Provides a high-level view (e.g., scatter plots of observed vs. simulated means).
#'   \item \strong{Participant:} Provides participant-level diagnostics (e.g., Q-Q plots).
#'   \item \strong{Density:} Compares the full shape of the reaction time distributions (e.g., via density plots).
#' }
#'
#' This method relies on three internal helper functions (\code{plot_ppc_summary},
#' \code{plot_ppc_participant}, and \code{plot_ppc_density}) which are
#' defined elsewhere in the package.
#'
#' @seealso
#' \code{\link{posterior_predict}} for generating the input object.
#'
#' @export
plot.driftd_ppc <- function(x, type = "summary", condition = NULL, ...) {

  if (!inherits(x, "driftd_ppc")) {
    stop("x must be of class 'driftd_ppc' from posterior_predict() function")
  }

  # Extract data and column names
  observed_data <- x$observed_data
  simulated_data <- x$simulated_data
  column_info <- x$column_info

  # Use the stored column names
  sub_id <- column_info$sub_id
  response <- column_info$response
  rt <- column_info$rt

  if (type == "summary") {
    plot_ppc_summary(observed_data, simulated_data, sub_id, response, rt)
  } else if (type == "participant") {
    plot_ppc_participant(observed_data, simulated_data, sub_id, response, rt)
  } else if (type == "density") {
    plot_ppc_density(observed_data, simulated_data, sub_id, response, rt, condition)
  } else {
    stop("type must be 'summary', 'participant', or 'density'")
  }
}
