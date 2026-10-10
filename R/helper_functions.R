
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




#' @title Perform Posterior Predictive Checks for Diffusion Model Fits
#'
#' @description
#' Uses the individual parameter estimates from a \code{driftd} fit to simulate
#' new data and generate a posterior predictive check (PPC) object.
#'
#' @param object An object of class \code{driftd_results}, typically the output
#'   from the \code{\link{driftd}} function.
#' @param n_sim An \code{integer} specifying the number of simulated datasets
#'   to generate for each participant (or each condition cell). Defaults to \code{1}.
#'
#' @return An object of class \code{driftd_ppc}, which is a \code{list}
#'   containing:
#'   \itemize{
#'     \item \code{simulated_data}: A \code{data.frame} combining all simulated
#'       reaction times (RTs) and responses, with a column \code{n_sim}
#'       indicating the simulation run.
#'     \item \code{observed_data}: The original data used for fitting.
#'     \item \code{column_info}: Information about the original column names (ID, RT, Response).
#'     \item \code{n_sim}: The number of simulations performed.
#'     \item \code{participants}: A list of unique participant IDs.
#'   }
#'
#' @details
#' This function assesses the **goodness-of-fit** of the diffusion model by
#' generating simulated data using the estimated parameters. The simulated data
#' can then be compared to the observed data via plotting functions.
#'
#' The simulation is performed cell-by-cell (participant by condition) using the
#' **rtdists::rdiffusion** function. The total number of simulated trials for a
#' cell is equal to the observed number of trials for that cell multiplied by
#' \code{n_sim}.
#'
#' The function relies on the original data being stored in the \code{driftd_results}
#' object, and will stop if \code{object$data} is \code{NULL}.
#'
#' @seealso
#' \code{\link{driftd}} for the main fitting function,
#' \code{\link[rtdists]{rdiffusion}} for the underlying simulation function.
#'
#' @importFrom utils txtProgressBar setTxtProgressBar close
#' @importFrom rtdists rdiffusion
#' @export
posterior_predict <- function(object, n_sim = 1) {

  if (!inherits(object, "driftd_results")) {
    stop("object must be of class 'driftd_results' from driftd() function")
  }

  if (is.null(object$data)) {
    stop("driftd_results object does not contain data. Please refit model with updated driftd function.")
  }

  # Get data and parameter estimates
  data <- object$data
  param_estimates <- object$individual_estimates
  settings <- object$settings

  # Extract column names from the stored data
  sub_id <- names(data)[1]   # participant ID column
  cond_col <- "cond"         # condition column (always named cond)
  response_col <- names(data)[3]  # user's response column name
  rt_col <- names(data)[4]        # user's RT column name

  # Initialize results
  all_simulated <- data.frame()

  # Get unique participants
  participants <- unique(data[[sub_id]])

  cat("Simulating data for", length(participants), "participants...\n")
  pb <- txtProgressBar(min = 0, max = length(participants), style = 3)

  for (i in seq_along(participants)) {
    participant <- participants[i]

    # Participant-specific observed data
    participant_data <- data[data[[sub_id]] == participant, ]

    # Participant's parameter estimates
    participant_params <- param_estimates[param_estimates[[sub_id]] == participant, ]

    if (nrow(participant_params) == 0) {
      warning("No parameter estimates found for participant: ", participant)
      next
    }

    # For each condition cell
    for (cond in unique(participant_data[[cond_col]])) {
      n_trials <- sum(participant_data[[cond_col]] == cond)
      if (n_trials == 0) next

      for (sim in seq_len(n_sim)) {
        # Extract parameters for this condition
        params <- get_parameters_for_condition(participant_params, cond, settings)

        if (!all(sapply(params, function(x) is.numeric(x) && is.finite(x)))) {
          warning("Invalid parameters for participant ", participant, " condition '", cond, "'. Skipping.")
          next
        }

        # Simulate data
        sim_data <- rtdists::rdiffusion(
          n = n_trials,
          a = params$a,
          v = params$v,
          t0 = params$t0,
          z = params$z * params$a,  # Convert to absolute starting point
          st0 = params$st0,
          sz = params$sz,
          sv = params$sv
        )

        # rdiffusion returns columns: "rt", "response"
        # Rename them to match the observed data's column names
        names(sim_data)[names(sim_data) == "rt"] <- rt_col
        names(sim_data)[names(sim_data) == "response"] <- response_col

        # Add metadata
        sim_data$n_sim <- sim
        sim_data[[sub_id]] <- participant
        sim_data[[cond_col]] <- cond

        # Reorder columns to match the observed data order
        # Observed data order is: ID, cond, response, RT
        sim_data <- sim_data[, c("n_sim", sub_id, cond_col, response_col, rt_col)]

        all_simulated <- rbind(all_simulated, sim_data)
      }
    }

    setTxtProgressBar(pb, i)
  }

  close(pb)

  if (nrow(all_simulated) == 0) {
    warning("No data was simulated. Check for parameter estimation issues.")
  }

  # Create the PPC object with column information
  ppc_object <- list(
    simulated_data = all_simulated,
    observed_data = data,
    column_info = list(
      sub_id = sub_id,
      response = response_col,
      rt = rt_col
    ),
    n_sim = n_sim,
    participants = participants
  )

  class(ppc_object) <- "driftd_ppc"
  return(ppc_object)
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
#' \code{\link{posterior_predict}} function and is not typically called by end-users.
#'
#' @keywords internal
get_parameters_for_condition <- function(participant_params, condition, settings) {

  if (is.null(settings$cell_parameters)) {
    stop("This driftd_results object was created by an older version of driftd. ",
         "Please refit the model.")
  }

  # Defaults, used only if a parameter value cannot be found
  params <- list(
    a = 1.0, t0 = 0.2, z = 0.5, v = 0.0,
    st0 = 0.0, sz = 0.0, sv = 0.0
  )

  cell <- settings$cell_parameters[settings$cell_parameters$cond == condition, , drop = FALSE]

  for (param_name in names(params)) {
    column <- if (nrow(cell) == 1) cell[[param_name]] else NA_character_
    if (is.na(column)) column <- param_name

    value <- participant_params[[column]]
    if (is.numeric(value) && length(value) == 1 && is.finite(value)) {
      params[[param_name]] <- value
    } else {
      cat("  Issue finding parameter value. Using default parameter:", param_name, "=", params[[param_name]], "\n")
    }
  }

  return(params)
}




#' @title Internal: Plot Summary Statistics for PPC
#'
#' @description
#' Generates a set of scatter plots comparing key condition-level summary
#' statistics (proportion of upper responses, median RTs for both boundaries)
#' between the observed data and the model-simulated data.
#'
#' @param observed_data A \code{data.frame} of the original data.
#' @param simulated_data A \code{data.frame} of the simulated data generated
#'   by \code{\link{posterior_predict}}.
#' @param sub_id A \code{character} string specifying the participant ID column name.
#' @param response A \code{character} string specifying the response column name.
#' @param rt A \code{character} string specifying the reaction time column name.
#'
#' @return The function is called for its side effect of creating a multi-panel plot.
#'   It invisibly returns \code{NULL}.
#'
#' @details
#' This function computes three key summary statistics—**Proportion of Upper Responses**,
#' **Median RT for Upper Responses**, and **Median RT for Lower Responses**—for
#' each condition in both the observed and simulated datasets.
#'
#' It then plots the observed vs. simulated values for these statistics on three
#' separate panels. Points that fall close to the diagonal $y=x$ line (which is included
#' for reference) indicate a good fit for that specific summary statistic at the
#' condition level.
#'
#' **Note:** This is an internal helper called by \code{\link{plot.driftd_ppc}}
#' when \code{type = "summary"} is specified, and is not meant for direct user interaction.
#'
#' @importFrom dplyr group_by summarise left_join .data
#' @importFrom graphics par layout plot abline title legend
#' @importFrom stats median setNames
#' @importFrom grDevices rainbow rgb
#' @keywords internal
#' Plot condition-level summary comparison
plot_ppc_summary <- function(observed_data, simulated_data, sub_id, response, rt) {

  # Calculate summary statistics for observed data by condition
  observed_summary <- observed_data %>%
    dplyr::group_by(.data$cond) %>%
    dplyr::summarise(
      prop_upper = mean(.data[[response]] == "upper", na.rm = TRUE),
      median_rt_upper = median(.data[[rt]][.data[[response]] == "upper"], na.rm = TRUE),
      median_rt_lower = median(.data[[rt]][.data[[response]] == "lower"], na.rm = TRUE),
      .groups = 'drop'
    )

  # Calculate summary statistics for simulated data by condition
  simulated_summary <- simulated_data %>%
    dplyr::group_by(.data$cond, .data$n_sim) %>%
    dplyr::summarise(
      prop_upper = mean(.data[[response]] == "upper", na.rm = TRUE),
      median_rt_upper = median(.data[[rt]][.data[[response]] == "upper"], na.rm = TRUE),
      median_rt_lower = median(.data[[rt]][.data[[response]] == "lower"], na.rm = TRUE),
      .groups = 'drop'
    ) %>%
    dplyr::group_by(.data$cond) %>%
    dplyr::summarise(
      prop_upper = mean(.data$prop_upper, na.rm = TRUE),
      median_rt_upper = mean(.data$median_rt_upper, na.rm = TRUE),
      median_rt_lower = mean(.data$median_rt_lower, na.rm = TRUE),
      .groups = 'drop'
    )

  # Merge observed and simulated summaries
  comparison <- dplyr::left_join(observed_summary, simulated_summary,
                                 by = "cond", suffix = c("_obs", "_sim"))

  # Colors for each condition
  unique_conds <- unique(comparison$cond)
  colors <- setNames(rainbow(length(unique_conds)), unique_conds)

  # Range for RTs
  rt_range <- range(c(
    comparison$median_rt_upper_obs, comparison$median_rt_upper_sim,
    comparison$median_rt_lower_obs, comparison$median_rt_lower_sim
  ), na.rm = TRUE)
  rt_lim <- c(0, max(1, rt_range[2] * 1.1))

  # Set up layout with space for legend at bottom
  layout(matrix(c(1, 2, 3, 4, 4, 4), nrow = 2, ncol = 3, byrow = TRUE),
         heights = c(0.85, 0.15))
  par(mar = c(5, 5, 4, 2), oma = c(0, 0, 3, 0))

  # Plot 1
  plot(comparison$prop_upper_obs, comparison$prop_upper_sim,
       col = colors[comparison$cond], pch = 16, cex = 1.8,
       xlim = c(0, 1), ylim = c(0, 1),
       xlab = "Observed Proportion", ylab = "Simulated Proportion",
       main = "Proportion of Upper Boundary\nResponses",
       cex.lab = 1.3, cex.axis = 1.2, cex.main = 1.4)
  abline(0, 1, lty = 2, col = "gray", lwd = 2)

  # Plot 2
  plot(comparison$median_rt_upper_obs, comparison$median_rt_upper_sim,
       col = colors[comparison$cond], pch = 16, cex = 1.8,
       xlim = rt_lim, ylim = rt_lim,
       xlab = "Observed Median RT (s)", ylab = "Simulated Median RT (s)",
       main = "Median RT for Upper Boundary\nResponses",
       cex.lab = 1.3, cex.axis = 1.2, cex.main = 1.4)
  abline(0, 1, lty = 2, col = "gray", lwd = 2)

  # Plot 3
  plot(comparison$median_rt_lower_obs, comparison$median_rt_lower_sim,
       col = colors[comparison$cond], pch = 16, cex = 1.8,
       xlim = rt_lim, ylim = rt_lim,
       xlab = "Observed Median RT (s)", ylab = "Simulated Median RT (s)",
       main = "Median RT for Lower Boundary\nResponses",
       cex.lab = 1.3, cex.axis = 1.2, cex.main = 1.4)
  abline(0, 1, lty = 2, col = "gray", lwd = 2)

  # Legend panel (empty plot for legend)
  par(mar = c(0, 0, 0, 0))
  plot(1, type = "n", axes = FALSE, xlab = "", ylab = "")

  # Calculate optimal number of columns based on number of conditions
  n_conds <- length(unique_conds)
  ncol_legend <- ifelse(n_conds > 8, 4,
                        ifelse(n_conds > 4, 3, 1))

  legend("center", legend = unique_conds, col = colors, pch = 16,
         pt.cex = 1.4, cex = 1, ncol = ncol_legend)

  title(main = "Posterior Predictive Check: Observed vs Simulated Data",
        outer = TRUE, cex.main = 1.8)
  par(mfrow = c(1, 1), mar = c(5, 4, 4, 2) + 0.1)  # Reset to default
}




#' @title Internal: Plot Participant-Level Summary Statistics for PPC
#'
#' @description
#' Generates scatter plots comparing key summary statistics (accuracy and median RTs)
#' between the observed and simulated data, with a separate point representing
#' **each participant-by-condition cell**. This allows for assessment
#' of individual-level fit.
#'
#' @param observed_data A \code{data.frame} of the original data.
#' @param simulated_data A \code{data.frame} of the simulated data generated
#'   by \code{\link{posterior_predict}}.
#' @param sub_id A \code{character} string specifying the participant ID column name.
#' @param response A \code{character} string specifying the response column name.
#' @param rt A \code{character} string specifying the reaction time column name.
#'
#' @return Invisibly returns a \code{list} of **Pearson correlation coefficients**
#'   calculated between the observed and simulated values for the three metrics.
#'   The function is primarily called for its side effect of creating a multi-panel
#'   plot and printing the correlation summary to the console.
#'
#' @details
#' This function computes three key summary statistics—**Proportion of Upper Responses**,
#' **Median RT for Upper Responses**, and **Median RT for Lower Responses**—for
#' every unique combination of participant and condition.
#'
#' The resulting multi-panel plot shows participant level data (one data point for each
#' condition cell per participant). Points close to the $y=x$ line indicate that
#' the model successfully captured the behavioral statistics at the individual level.
#' The function also calculates and prints the correlation between observed and
#' simulated values for each metric, which serves as a quantitative measure of fit quality.
#'
#' **Note:** This is an internal helper called by \code{\link{plot.driftd_ppc}}
#' when \code{type = "participant"} is specified, and is not meant for direct user interaction.
#'
#' @importFrom dplyr group_by summarise left_join .data
#' @importFrom graphics par layout plot abline title legend
#' @importFrom stats median setNames cor
#' @importFrom grDevices rainbow rgb
#' @keywords internal
#' Plot participant-level comparison
plot_ppc_participant <- function(observed_data, simulated_data, sub_id, response, rt) {

  # Calculate summary statistics for observed data by participant and condition
  observed_summary <- observed_data %>%
    dplyr::group_by(.data[[sub_id]], .data$cond) %>%
    dplyr::summarise(
      prop_upper = mean(.data[[response]] == "upper", na.rm = TRUE),
      median_rt_upper = median(.data[[rt]][.data[[response]] == "upper"], na.rm = TRUE),
      median_rt_lower = median(.data[[rt]][.data[[response]] == "lower"], na.rm = TRUE),
      .groups = 'drop'
    )

  # Calculate summary statistics for simulated data by participant and condition
  simulated_summary <- simulated_data %>%
    dplyr::group_by(.data[[sub_id]], .data$cond, .data$n_sim) %>%
    dplyr::summarise(
      prop_upper = mean(.data[[response]] == "upper", na.rm = TRUE),
      median_rt_upper = median(.data[[rt]][.data[[response]] == "upper"], na.rm = TRUE),
      median_rt_lower = median(.data[[rt]][.data[[response]] == "lower"], na.rm = TRUE),
      .groups = 'drop'
    ) %>%
    dplyr::group_by(.data[[sub_id]], .data$cond) %>%
    dplyr::summarise(
      prop_upper = mean(.data$prop_upper, na.rm = TRUE),
      median_rt_upper = mean(.data$median_rt_upper, na.rm = TRUE),
      median_rt_lower = mean(.data$median_rt_lower, na.rm = TRUE),
      .groups = 'drop'
    )

  # Merge observed and simulated summaries
  comparison <- dplyr::left_join(observed_summary, simulated_summary,
                                 by = c(sub_id, "cond"), suffix = c("_obs", "_sim"))

  # Get unique conditions for coloring
  unique_conds <- unique(comparison$cond)
  colors <- rainbow(length(unique_conds))
  names(colors) <- unique_conds

  # Determine RT ranges
  rt_range_upper <- range(c(comparison$median_rt_upper_obs, comparison$median_rt_upper_sim), na.rm = TRUE)
  rt_range_lower <- range(c(comparison$median_rt_lower_obs, comparison$median_rt_lower_sim), na.rm = TRUE)
  rt_range <- range(c(rt_range_upper, rt_range_lower), na.rm = TRUE)
  rt_lim <- c(0, max(1, rt_range[2] * 1.1))

  # Set up layout with space for legend at bottom
  layout(matrix(c(1, 2, 3, 4, 4, 4), nrow = 2, ncol = 3, byrow = TRUE),
         heights = c(0.85, 0.15))
  par(mar = c(5, 5, 4, 2), oma = c(0, 0, 3, 0))

  # Plot 1: Proportion of upper boundary responses
  plot(comparison$prop_upper_obs, comparison$prop_upper_sim,
       col = colors[comparison$cond], pch = 16, cex = 1.2,
       xlim = c(0, 1), ylim = c(0, 1),
       xlab = "Observed Proportion", ylab = "Simulated Proportion",
       main = "Proportion of Upper Boundary\nResponses (Participant Level)",
       cex.lab = 1.3, cex.axis = 1.2, cex.main = 1.3)
  abline(a = 0, b = 1, lty = 2, col = "gray", lwd = 2)

  # Plot 2: Median RT for upper boundary responses
  plot(comparison$median_rt_upper_obs, comparison$median_rt_upper_sim,
       col = colors[comparison$cond], pch = 16, cex = 1.2,
       xlim = rt_lim, ylim = rt_lim,
       xlab = "Observed Median RT (s)", ylab = "Simulated Median RT (s)",
       main = "Median RT for Upper Boundary\nResponses (Participant Level)",
       cex.lab = 1.3, cex.axis = 1.2, cex.main = 1.3)
  abline(a = 0, b = 1, lty = 2, col = "gray", lwd = 2)

  # Plot 3: Median RT for lower boundary responses
  plot(comparison$median_rt_lower_obs, comparison$median_rt_lower_sim,
       col = colors[comparison$cond], pch = 16, cex = 1.2,
       xlim = rt_lim, ylim = rt_lim,
       xlab = "Observed Median RT (s)", ylab = "Simulated Median RT (s)",
       main = "Median RT for Lower Boundary\nResponses (Participant Level)",
       cex.lab = 1.3, cex.axis = 1.2, cex.main = 1.3)
  abline(a = 0, b = 1, lty = 2, col = "gray", lwd = 2)

  # Legend panel (empty plot for legend)
  par(mar = c(0, 0, 0, 0))
  plot(1, type = "n", axes = FALSE, xlab = "", ylab = "")

  # Calculate optimal number of columns based on number of conditions
  n_conds <- length(unique_conds)
  ncol_legend <- ifelse(n_conds > 8, 4,
                        ifelse(n_conds > 4, 3, 1))

  legend("center", legend = unique_conds, col = colors, pch = 16,
         pt.cex = 1.4, cex = 1, ncol = ncol_legend)

  # Add overall title
  title(main = "Posterior Predictive Check: Participant Level",
        outer = TRUE, cex.main = 1.6)
  par(mfrow = c(1, 1), mar = c(5, 4, 4, 2) + 0.1)  # Reset to default

  # Calculate and print correlations
  correlations <- list(
    prop_upper = cor(comparison$prop_upper_obs, comparison$prop_upper_sim, use = "complete.obs"),
    median_rt_upper = cor(comparison$median_rt_upper_obs, comparison$median_rt_upper_sim, use = "complete.obs"),
    median_rt_lower = cor(comparison$median_rt_lower_obs, comparison$median_rt_lower_sim, use = "complete.obs")
  )

  cat("Participant-level correlations between observed and simulated data:\n")
  cat("Proportion upper boundary:", round(correlations$prop_upper, 3), "\n")
  cat("Median RT upper boundary:", round(correlations$median_rt_upper, 3), "\n")
  cat("Median RT lower boundary:", round(correlations$median_rt_lower, 3), "\n")

  return(invisible(correlations))
}




#' @title Internal: Plot RT Density Distributions for PPC
#'
#' @description
#' Compares the full reaction time (RT) distributions between the observed data
#' and the model-simulated data by plotting smoothed density curves for both
#' the upper and lower response boundaries.
#'
#' @param observed_data A \code{data.frame} of the original data.
#' @param simulated_data A \code{data.frame} of the simulated data generated
#'   by \code{\link{posterior_predict}}.
#' @param sub_id A \code{character} string specifying the participant ID column name.
#' @param response A \code{character} string specifying the response column name.
#' @param rt A \code{character} string specifying the reaction time column name.
#' @param condition A \code{character} string or \code{NULL} (default). If
#'   specified, the density plot is filtered and limited to the data from
#'   this particular condition.
#'
#' @return The function is called for its side effect of creating a multi-panel plot.
#'   It invisibly returns \code{NULL}.
#'
#' @details
#' This is often the most visually informative PPC plot, as it directly compares
#' the **shape and spread** of the observed RT distribution (blue line) against
#' the distribution generated by the fitted model (red line).
#'
#' The function generates two plots:
#' \itemize{
#'   \item **Upper Boundary:** Density comparison for correct/designated responses.
#'   \item **Lower Boundary:** Density comparison for error/alternative responses.
#' }
#' Good model fit is indicated by the red simulated line closely overlapping the
#' blue observed line across the entire distribution, including the fast and slow tails.
#' The function can aggregate data across all participants and conditions, or focus
#' on a single condition if the \code{condition} argument is provided.
#'
#' **Note:** This is an internal helper called by \code{\link{plot.driftd_ppc}}
#' when \code{type = "density"} is specified, and is not meant for direct user interaction.
#'
#' @importFrom graphics par layout plot lines title legend
#' @importFrom stats density
#' @importFrom grDevices rainbow rgb
#' @keywords internal
#' Plot RT density distributions
plot_ppc_density <- function(observed_data, simulated_data, sub_id, response, rt, condition = NULL) {

  # Filter by condition if specified
  if (!is.null(condition)) {
    observed_data <- observed_data[observed_data$cond == condition, ]
    simulated_data <- simulated_data[simulated_data$cond == condition, ]
    plot_title <- paste("RT Distribution Comparison - Condition:", condition)
  } else {
    plot_title <- "RT Distribution Comparison (All Conditions)"
  }

  # Set up layout with space for legend at bottom
  layout(matrix(c(1, 2, 3, 3), nrow = 2, ncol = 2, byrow = TRUE),
         heights = c(0.85, 0.15))
  par(mar = c(5, 5, 4, 2), oma = c(0, 0, 3, 0))

  # Determine RT range for consistent x-axis
  upper_obs <- observed_data[[rt]][observed_data[[response]] == "upper"]
  upper_sim <- simulated_data[[rt]][simulated_data[[response]] == "upper"]
  lower_obs <- observed_data[[rt]][observed_data[[response]] == "lower"]
  lower_sim <- simulated_data[[rt]][simulated_data[[response]] == "lower"]

  all_rts <- c(upper_obs, upper_sim, lower_obs, lower_sim)
  rt_range <- range(all_rts, na.rm = TRUE)
  x_lim <- c(0, max(1, rt_range[2] * 1.1))

  # Plot 1: Upper boundary RT distributions
  dens_obs <- density(upper_obs, na.rm = TRUE)
  dens_sim <- density(upper_sim, na.rm = TRUE)

  y_max <- max(c(dens_obs$y, dens_sim$y)) * 1.1

  plot(dens_obs, col = "blue", lwd = 3,
       xlim = x_lim, ylim = c(0, y_max),
       xlab = "Response Time (s)", ylab = "Density",
       main = "Upper Boundary RT Distributions",
       cex.lab = 1.3, cex.axis = 1.2, cex.main = 1.4)
  lines(dens_sim, col = "red", lwd = 3)

  # Plot 2: Lower boundary RT distributions
  dens_obs_lower <- density(lower_obs, na.rm = TRUE)
  dens_sim_lower <- density(lower_sim, na.rm = TRUE)

  y_max_lower <- max(c(dens_obs_lower$y, dens_sim_lower$y)) * 1.1

  plot(dens_obs_lower, col = "blue", lwd = 3,
       xlim = x_lim, ylim = c(0, y_max_lower),
       xlab = "Response Time (s)", ylab = "Density",
       main = "Lower Boundary RT Distributions",
       cex.lab = 1.3, cex.axis = 1.2, cex.main = 1.4)
  lines(dens_sim_lower, col = "red", lwd = 3)

  # Legend panel (empty plot for legend)
  par(mar = c(0, 0, 0, 0))
  plot(1, type = "n", axes = FALSE, xlab = "", ylab = "")
  legend("center", legend = c("Observed", "Simulated"),
         col = c("blue", "red"), lwd = 3, cex = 1.0, ncol = 2)

  # Add overall title
  title(main = plot_title, outer = TRUE, cex.main = 1.6)
  par(mfrow = c(1, 1), mar = c(5, 4, 4, 2) + 0.1)  # Reset to default
}
