
#' @title Drift-Diffusion Model Fitting
#' @description Fits a drift-diffusion model (DDM) to response time (RT) and response data for individual participants using maximum likelihood estimation via \code{nlminb}.
#'
#' @param data A data frame containing the experimental data.
#' @param sub_id Character string specifying the column name for subject/participant ID. Default is "sub_id".
#' @param response Character string specifying the column name for responses (coded as 0/1 or 'lower'/'upper'). Default is "response".
#' @param rt Character string specifying the column name for reaction times. Default is "rt". RTs should be in seconds; if the mean RT is above 10, RTs are assumed to be in milliseconds and converted to seconds (with a message).
#' @param vary_a Character vector of column names indicating factors that modulate boundary separation (a). The model will estimate distinct values of (\code{a}) for each level combination of the specified variables.
#' @param vary_t0 Character vector of column names indicating factors that modulate non-decision time (t0). The model will estimate distinct values of (\code{t0}) for each level combination of the specified variables.
#' @param vary_z Character vector of column names indicating factors that modulate starting point (z). The model will estimate distinct values of (\code{z}) for each level combination of the specified variables.
#' @param vary_v Character vector of column names indicating factors that modulate drift rate (v). The model will estimate distinct values of (\code{v}) for each level combination of the specified variables.
#' @param vary_st0 Character vector of column names indicating factors that modulate trial-by-trial variability in non-decision time (st0). The model will estimate distinct values of (\code{st0}) for each level combination of the specified variables.
#' @param vary_sz Character vector of column names indicating factors that modulate trial-by-trial variability in starting point (sz). The model will estimate distinct values of (\code{sz}) for each level combination of the specified variables.
#' @param vary_sv Character vector of column names indicating factors that modulate trial-by-trial variability in drift rate (sv). The model will estimate distinct values of (\code{sv}) for each level combination of the specified variables.
#' @param fix_a Logical: If TRUE, \code{a} is fixed to \code{a_fixed} and will not be estimated.
#' @param fix_t0 Logical: If TRUE, \code{t0} is fixed to \code{t0_fixed} and will not be estimated.
#' @param fix_z Logical: If TRUE, \code{z} is fixed to \code{z_fixed} and will not be estimated.
#' @param fix_v Logical: If TRUE, \code{v} is fixed to \code{v_fixed} and will not be estimated.
#' @param fix_st0 Logical: If TRUE, \code{st0} is fixed to \code{st0_fixed} and will not be estimated. Default is TRUE.
#' @param fix_sz Logical: If TRUE, \code{sz} is fixed to \code{sz_fixed} and will not be estimated. Default is TRUE.
#' @param fix_sv Logical: If TRUE, \code{sv} is fixed to \code{sv_fixed} and will not be estimated. Default is TRUE.
#' @param a_fixed Numeric: Fixed value for \code{a} if \code{fix_a=TRUE}. Default is 1.
#' @param t0_fixed Numeric: Fixed value for \code{t0} if \code{fix_t0=TRUE}. Default is 0.2.
#' @param z_fixed Numeric: Fixed value for \code{z} (relative to \code{a}) if \code{fix_z=TRUE}. Default is 0.5.
#' @param v_fixed Numeric: Fixed value for \code{v} if \code{fix_v=TRUE}. Default is 0.5.
#' @param st0_fixed Numeric: Fixed value for \code{st0} if \code{fix_st0=TRUE}. Default is 0.
#' @param sz_fixed Numeric: Fixed value for \code{sz} if \code{fix_sz=TRUE}. Default is 0. As in 'rtdists', \code{sz} is the range of starting points on the scale of \code{a} (not relative to \code{a} like \code{z}), so it must be smaller than \code{2 * min(z, 1 - z) * a}.
#' @param sv_fixed Numeric: Fixed value for \code{sv} if \code{fix_sv=TRUE}. Default is 0.
#' @param precision Numeric: Precision parameter for \code{rtdists::ddiffusion}. Default is 3.
#' @param stop_on_error Logical: Passed to \code{rtdists::ddiffusion}. Default is FALSE.
#' @param print_starting_vals Logical: If TRUE, prints the derived starting parameter estimates for each participant. Default is FALSE.
#' @param compute_se Logical: If TRUE (default), computes standard errors of the estimated parameters from the numerical Hessian of the log-likelihood (see Details).
#' @param ask Logical: if TRUE (default) and R is running interactively, ask for confirmation before fitting when the pre-fit checks issue warnings (e.g. few trials per condition). Set to FALSE to always continue.
#' @param p_outlier Numeric between 0 and 1: the assumed proportion of outlier (contaminant) trials, as in HDDM (Wiecki, Sofer, & Frank, 2013). If above 0, each trial's likelihood is a mixture of the diffusion model density and a uniform outlier density, \code{(1 - p_outlier) * ddiffusion + p_outlier * w_outlier}, so a few trials the diffusion model cannot explain (e.g. fast guesses or lapses) cannot distort the estimates. \code{p_outlier} is fixed, not estimated. Default is 0 (no outlier component); HDDM recommends 0.05. See Details.
#' @param w_outlier Numeric: the density of the outlier distribution, as in HDDM. The default of 0.1 is HDDM's, a uniform distribution of RTs between 0 and 5 s with both responses equally likely (\code{1 / (2 * 5)}). Only used if \code{p_outlier > 0}.
#' @param max_retries Integer: Maximum number of additional optimizer runs for a participant whose fit does not converge, or whose estimate of \code{sv} or \code{sz} is at 0 (below 0.001, the same tolerance used for parameters at a bound when computing standard errors). Retries continue until the best converged solution has been reached from two independent starting points. Set to 0 to disable retries. Default is 5.
#'
#' @details
#' The function estimates DDM parameters (\code{a, t0, z, v, st0, sz, sv}) for each participant individually, allowing parameters to vary across experimental conditions specified in the \code{vary_} arguments. It uses the maximum likelihood method with the \code{rtdists::ddiffusion} function. Note that \code{z} is estimated as a proportion (0 to 1) and internally converted to \code{z * a} for \code{ddiffusion}. A parameter cannot be both fixed (\code{fix_} = TRUE) and varied across conditions (\code{vary_}); this combination is an error. \code{t0} is bounded above by the fastest RT it applies to (with \code{p_outlier > 0}, see Outliers below). It includes pre-checks for data quality.
#'
#' Starting values for \code{v}, \code{a} and \code{t0} are EZ-diffusion estimates (Wagenmakers, van der Maas, & Grasman, 2007) computed per participant and condition cell. If a fit does not converge, the optimizer is first restarted from where it stopped. If it still does not converge, or if an estimate of \code{sv} or \code{sz} is at 0 (below 0.001; where the likelihood is flat in these parameters, so the optimizer can stop there even when a positive value fits better), the model is refitted from further starting values: first with \code{sv}/\code{sz} moved to positive values, then alternative starting values spread over the plausible parameter space with a deterministic Halton sequence. This continues until the best converged solution has been found from two independent starting points or \code{max_retries} is reached. The fit with the highest likelihood is kept.
#'
#' \strong{Outliers.} With \code{p_outlier > 0}, driftd follows HDDM's robust likelihood (Ratcliff & Tuerlinckx, 2002; Wiecki et al., 2013). Every trial has a likelihood of at least \code{p_outlier * w_outlier}, so a trial the diffusion model cannot explain, such as a fast guess below \code{t0}, adds a bounded cost instead of pulling the parameters towards it. For the same reason, \code{t0} is then no longer bounded above by the fastest RT, but by the 10% quantile of the RTs it applies to, so that a few fast outliers cannot pull it down. The outlier component is part of the likelihood only: as in HDDM, \code{predict()}, \code{simulate()} and the recovery functions generate data from the diffusion model alone, so trials that the model treats as outliers show up as misfit in the model checks. The log-likelihood, AIC and BIC are those of the mixture, so only compare models fitted with the same \code{p_outlier}.
#'
#' Standard errors are the square roots of the diagonal of the inverse of the Hessian of the negative log-likelihood at the maximum likelihood estimates (the observed Fisher information), computed numerically with Richardson extrapolation (\code{numDeriv::hessian}) at a \code{precision} of at least 4. They are asymptotic (Wald) standard errors; approximate 95% confidence intervals are estimate +/- 1.96 SE. Parameters estimated at (within 0.001 of) a bound, such as \code{sv = 0} or \code{t0} at the fastest RT, get no standard error, because the usual asymptotic theory does not hold there; the standard errors of the other parameters are then computed with these held at their estimates. If the Hessian is not positive definite (the fit is not a proper maximum, or a parameter is not identified), no standard errors are reported for that participant.
#'
#' @return A list of class \code{driftd_results} containing:
#' \item{individual_estimates}{A data frame with estimated parameters, convergence status, Log-Likelihood (LLE), AIC, and BIC for each participant.}
#' \item{summary}{A data frame with mean and standard deviation of all estimated parameters across participants.}
#' \item{settings}{A list summarizing the model structure (varying and fixed parameters), including \code{cell_parameters}, a data frame recording which parameter applies in each condition cell.}
#' \item{data}{The processed data frame used for fitting.}
#' \item{convergence}{A table summarizing the convergence codes from \code{nlminb}.}
#' \item{standard_errors}{A data frame with the same layout as \code{individual_estimates} (one column per estimated parameter) containing the standard errors. \code{NA} for parameters at a bound, and for all parameters of a participant if the Hessian is not positive definite (see \code{fit_info$se_status}).}
#' \item{vcov}{A list with, per participant (named by ID), the estimated variance-covariance matrix of the parameters that have standard errors, or \code{NULL}. Needed for standard errors of functions of several parameters (e.g. differences between conditions).}
#' \item{fit_info}{A data frame with, per participant, the number of optimizer runs, the final convergence code, whether the best solution was reached from two independent starting points (\code{NA} if no retries were needed), the \code{nlminb} message, and the status of the standard error computation.}
#' @export
#'
#' @importFrom rtdists ddiffusion
#' @importFrom dplyr %>% group_by summarise n filter select mutate case_when left_join bind_cols across everything all_of
#' @importFrom rlang sym
#' @importFrom tibble tibble
#' @importFrom gt gt tab_spanner fmt_number
#' @importFrom stats nlminb sd quantile relevel setNames
#' @importFrom utils txtProgressBar setTxtProgressBar head
driftd <- function(data, sub_id = "sub_id", response = "response", rt = "rt",

                   vary_a = c(),
                   vary_t0 = c(),
                   vary_z = c(),
                   vary_v = c(),
                   vary_st0 = c(),
                   vary_sz = c(),
                   vary_sv = c(),

                   fix_a = FALSE,
                   fix_t0 = FALSE,
                   fix_z = FALSE,
                   fix_v = FALSE,
                   fix_st0 = TRUE,
                   fix_sz = TRUE,
                   fix_sv = TRUE,

                   a_fixed = 1,
                   t0_fixed = 0.2,
                   z_fixed = 0.5,
                   v_fixed = 0.5,
                   st0_fixed = 0,
                   sz_fixed = 0,
                   sv_fixed = 0,

                   precision = 3,
                   stop_on_error = FALSE,

                   print_starting_vals = FALSE,

                   p_outlier = 0,
                   w_outlier = 0.1,

                   max_retries = 5,

                   compute_se = TRUE,

                   ask = TRUE

)
{
  #SECTION 0: PRECHECKS AND WARNINGS

  # Check that data is dataframe
  if (!is.data.frame(data)) {
    stop("'data' must be a data frame")
  }

  # check that required columns are there
  required_cols <- c(sub_id, response, rt)
  missing_cols <- setdiff(required_cols, names(data))
  if (length(missing_cols) > 0) {
    stop("The following specified columns were not found in data: ",
         paste(missing_cols, collapse = ", "))
  }

  # Check for missing values in critical columns
  if (any(is.na(data[[sub_id]]))) {
    stop("NA values found in subject ID column.")
  }

  if (any(is.na(data[[response]]))) {
    stop("NA values found in response column.")
  }

  if (any(is.na(data[[rt]]))) {
    stop("NA values found in RT column.")
  }

  # Validate RT values
  if (any(data[[rt]] <= 0)) {
    stop("All RT values must be positive")
  }

  # rtdists works in seconds. RTs averaging over 10 are assumed to be milliseconds.
  if (mean(data[[rt]]) > 10) {
    message("RTs appear to be in milliseconds (mean RT = ", round(mean(data[[rt]]), 1),
            "). Converting to seconds; all estimates and simulated RTs are in seconds.")
    data[[rt]] <- data[[rt]] / 1000
  }

  # A parameter can either be fixed to one value or vary across conditions, not both
  vary_args <- list(a = vary_a, t0 = vary_t0, z = vary_z, v = vary_v,
                    st0 = vary_st0, sz = vary_sz, sv = vary_sv)
  fix_args <- c(a = fix_a, t0 = fix_t0, z = fix_z, v = fix_v,
                st0 = fix_st0, sz = fix_sz, sv = fix_sv)
  conflicting <- names(fix_args)[fix_args & lengths(vary_args) > 0]
  if (length(conflicting) > 0) {
    stop("These parameters are both fixed and set to vary across conditions: ",
         paste(conflicting, collapse = ", "), ". Set fix_", conflicting[1],
         " = FALSE to estimate it per condition, or remove vary_", conflicting[1], ".")
  }

  # Check that the condition columns exist and have no missing values
  vary_cols <- unique(unlist(vary_args))
  missing_vary <- setdiff(vary_cols, names(data))
  if (length(missing_vary) > 0) {
    stop("The following vary_ columns were not found in data: ",
         paste(missing_vary, collapse = ", "))
  }
  na_vary <- vary_cols[vapply(vary_cols, function(x) any(is.na(data[[x]])), logical(1))]
  if (length(na_vary) > 0) {
    stop("NA values found in condition column(s): ", paste(na_vary, collapse = ", "))
  }

  # The model specification, stored so that the same model can be refitted (e.g. to
  # simulated data for parameter and model recovery)
  model_spec <- list(vary_a = vary_a, vary_t0 = vary_t0, vary_z = vary_z, vary_v = vary_v,
                     vary_st0 = vary_st0, vary_sz = vary_sz, vary_sv = vary_sv,
                     fix_a = fix_a, fix_t0 = fix_t0, fix_z = fix_z, fix_v = fix_v,
                     fix_st0 = fix_st0, fix_sz = fix_sz, fix_sv = fix_sv,
                     a_fixed = a_fixed, t0_fixed = t0_fixed, z_fixed = z_fixed, v_fixed = v_fixed,
                     st0_fixed = st0_fixed, sz_fixed = sz_fixed, sv_fixed = sv_fixed,
                     precision = precision, stop_on_error = stop_on_error, max_retries = max_retries,
                     p_outlier = p_outlier, w_outlier = w_outlier)

  # Validate parameter bounds
  if (a_fixed <= 0) stop("a_fixed must be positive")
  if (t0_fixed < 0) stop("t0_fixed cannot be negative")
  # z = 0 or 1 puts the starting point on a boundary, where every trial has likelihood 0
  if (z_fixed <= 0 || z_fixed >= 1) stop("z_fixed must be strictly between 0 and 1")
  if (st0_fixed < 0) stop("st0_fixed cannot be negative")
  if (sz_fixed < 0 || sz_fixed > 1) stop("sz_fixed must be between 0 and 1")
  if (sv_fixed < 0) stop("sv_fixed cannot be negative")
  if (!is.numeric(p_outlier) || length(p_outlier) != 1 || is.na(p_outlier) || p_outlier < 0 || p_outlier >= 1) {
    stop("p_outlier must be a single number between 0 and 1 (0 = no outlier component)")
  }
  if (!is.numeric(w_outlier) || length(w_outlier) != 1 || is.na(w_outlier) || w_outlier <= 0) {
    stop("w_outlier must be a single positive number")
  }
  # sz is the range of starting points on the scale of a (as in rtdists), so all starting
  # points lie between the boundaries only if sz < 2 * min(z, 1 - z) * a. Otherwise every
  # trial has likelihood 0. This can be checked in advance when a, z and sz are all fixed.
  if (fix_sz && fix_a && fix_z && sz_fixed >= 2 * min(z_fixed, 1 - z_fixed) * a_fixed) {
    stop("sz_fixed (", sz_fixed, ") must be smaller than 2 * min(z_fixed, 1 - z_fixed) * a_fixed (",
         2 * min(z_fixed, 1 - z_fixed) * a_fixed, "), or starting points fall outside the boundaries")
  }

  run_pre_checks <- function(data, sub_id, response, rt, vary_a, vary_t0, vary_z, vary_v,
                             vary_st0, vary_sz, vary_sv, fix_st0, st0_fixed,
                             fix_sz, sz_fixed) {

    warnings_issued <- FALSE

    # Parameter checks
    if (fix_st0 && st0_fixed > 0.5) {
      message("WARNING: st0_fixed seems large (", st0_fixed, "). This may cause estimation issues.")
      warnings_issued <- TRUE
    }

    if (fix_sz && sz_fixed > 0.5) {
      message("WARNING: sz_fixed seems large (", sz_fixed, "). This may cause estimation issues.")
      warnings_issued <- TRUE
    }

    # FIXED: Check condition counts WITHIN participants
    check_conditions <- unique(unlist(c(vary_a, vary_t0, vary_z, vary_v, vary_st0, vary_sz, vary_sv)))
    check_conditions <- check_conditions[!is.na(check_conditions) & check_conditions != ""]

    if (length(check_conditions) > 0) {
      # Create condition column
      data$temp_cond <- do.call(paste, c(data[check_conditions], sep = "_"))

      # Count trials per participant per condition
      participant_condition_counts <- data %>%
        group_by(!!sym(sub_id), .data$temp_cond) %>%
        summarise(trial_count = n(), .groups = 'drop')

      # Find participant-condition combinations with low trial counts
      low_trial_combos <- participant_condition_counts %>%
        filter(.data$trial_count < 25)

      if (nrow(low_trial_combos) > 0) {
        message("WARNING: Some participant-condition combinations have very few trials (<25):")

        # Group by condition and list affected participants
        by_condition <- low_trial_combos %>%
          group_by(.data$temp_cond) %>%
          summarise(
            participants = paste(!!sym(sub_id), " (n=", .data$trial_count, ")", sep = "", collapse = ", "),
            .groups = 'drop'
          )

        for (i in seq_len(nrow(by_condition))) {
          message("  Condition '", by_condition$temp_cond[i], "': ", by_condition$participants[i])
        }

        message("This may cause estimation problems.")
        warnings_issued <- TRUE
      }

      # Remove temporary column
      data$temp_cond <- NULL
    }

    # Participant checks (total trials per participant). Counted on the IDs as text so
    # that unused levels of a factor ID column are not reported as participants with 0 trials.
    participant_counts <- table(as.character(data[[sub_id]]))
    small_participants <- participant_counts[participant_counts < 20]
    if (length(small_participants) > 0) {
      message("WARNING: Some participants have very few total trials (<20): ",
              paste(names(small_participants), " (n=", small_participants, ")", sep = "", collapse = ", "),
              ". This may cause estimation problems.")
      warnings_issued <- TRUE
    }

    # Participants who only ever responded at one boundary: without any responses at the
    # other boundary, boundary separation, drift rate and starting point are not well
    # identified (the likelihood can keep increasing as a and v grow).
    one_sided <- tapply(as.character(data[[response]]), as.character(data[[sub_id]]),
                        function(x) length(unique(x)) < 2)
    if (any(one_sided)) {
      message("WARNING: Some participants gave all their responses at one boundary: ",
              paste(names(one_sided)[one_sided], collapse = ", "),
              ". Their estimates of a, v and z are unreliable.")
      warnings_issued <- TRUE
    }

    # Very fast responses are usually anticipations. t0 cannot exceed the fastest RT it
    # applies to, so a single such response can pull a participant's t0 estimate down.
    fast <- table(as.character(data[[sub_id]])[data[[rt]] < 0.15])
    if (length(fast) > 0) {
      message("WARNING: Some participants have RTs below 0.15 s: ",
              paste(names(fast), " (n=", fast, ")", sep = "", collapse = ", "),
              ". t0 must be below a participant's fastest RT, so such responses can pull ",
              "the t0 estimate down. Consider excluding them before fitting.")
      warnings_issued <- TRUE
    }

    return(warnings_issued)
  }

  # Run the pre-checks
  cat("Running pre-fit checks...\n")
  warnings_issued <- run_pre_checks(data, sub_id, response, rt, vary_a, vary_t0, vary_z, vary_v,
                                    vary_st0, vary_sz, vary_sv, fix_st0, st0_fixed,
                                    fix_sz, sz_fixed)

  # Ask for confirmation
  if (warnings_issued && ask && interactive()) {
    cat("\n")
    msg <- readline("Do you want to continue with model fitting? (y/n): ")
    if (tolower(msg) != "y") {
      message("Model fitting cancelled by user.")
      return(invisible(NULL))
    }
    cat("Continuing with model fitting...\n\n")
  }

  # SECTION A: SETTING UP THE DATAFRAME

  # discern how many conditions there are that need to be considered when varying parameters
  conditions <- unlist(c(vary_a, vary_t0, vary_z, vary_v, vary_st0, vary_sz, vary_sv))
  conditions <- unique(conditions[!is.na(conditions) & conditions != ""])

  #create a list with sublists that hold the level names for each condition
  levels_list <- lapply(conditions, function(x) sort(unique(data[[x]])))
  names(levels_list) <- conditions

  names(a_fixed) <- "a"
  names(t0_fixed) <- "t0"
  names(z_fixed) <- "z"
  names(v_fixed) <- "v"
  names(st0_fixed) <- "st0"
  names(sz_fixed) <- "sz"
  names(sv_fixed) <- "sv"

  # create cond column that indicates which experimental cell a given trial came from
  # then paste this to the dataset and rename dataset df.rtdists for all further operations
  # rename and relevel responses to prep for ddiffusion later on
  if (length(conditions) < 1) {
    cond <- rep("NA", nrow(data))
  } else {
    cond <- do.call(paste, c(data[conditions], sep = "_"))
    if (length(unique(cond)) != nrow(unique(data[conditions]))) {
      stop("Some condition levels contain '_' in a way that gives two different conditions ",
           "the same label. Please recode these levels without '_'.")
    }
  }

  #before next steps, save original column names for processed_data frame
  S <- sub_id
  R <- response
  React <- rt

  ID <- data[[sub_id]]
  response <- data[[response]]
  RT <- data[[rt]]

  # Check if df.rtdists$response contains only 0 and 1 or upper lower
  if (all(response %in% c(0, 1))) {
    # Convert 0 to "lower" and 1 to "upper"
    response <- ifelse(response == 1, "upper", "lower")
  } else if (all(response %in% c("upper", "lower"))) {
    response <- response
  } else {
    stop("Response variable must be coded as 0/1 or 'lower'/'upper'")
  }

  df.rtdists = data.frame(ID, cond, response, RT)
  # rtdists reads a factor response by its level number (1 = lower, 2 = upper), so both
  # levels must always exist, even when every response is at the same boundary
  df.rtdists$response <- factor(df.rtdists$response, levels = c("lower", "upper"))

  # RP List. Will be used to loop over RPs later on
  RPlist <- unique(df.rtdists$ID)
  #RPlist <- RPlist[RPlist != 77]

  # SECTION B: GENERATING THE CONDITION CELLS EACH PARAMETER IS ESTIMATED FOR

  #order vary lists alphabetically
  vary_a <- sort(vary_a)
  vary_t0 <- sort(vary_t0)
  vary_z <- sort(vary_z)
  vary_v <- sort(vary_v)
  vary_st0 <- sort(vary_st0)
  vary_sz <- sort(vary_sz)
  vary_sv <- sort(vary_sv)
  vary_list <- list(a = vary_a, t0 = vary_t0, z = vary_z, v = vary_v,
                    st0 = vary_st0, sz = vary_sz, sv = vary_sv)

  # Function to generate all level combinations of the given condition columns, one row per
  # combination. The first column varies slowest (neg_0.4, neg_0.45, ..., pos_0.4, ...).
  generate_cells <- function(vary) {
    grid <- expand.grid(rev(levels_list[vary]), KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
    grid[vary]
  }

  # SECTION C: GENERATING STARTING PARAMETER LIST THAT WILL BE FED TO NLMINB

  # For each parameter, generate one name per condition cell it is estimated for, and keep
  # the cell's level values (param_cells) so trials can be matched to parameters by value.
  # A parameter that doesn't vary across any conditions is simply given its name.
  param_cells <- list()
  param_type <- c()
  for (par in names(vary_list)) {
    if (length(vary_list[[par]]) == 0) {
      param_cells[par] <- list(list())
      param_type[par] <- par
    } else {
      grid <- generate_cells(vary_list[[par]])
      cell_names <- paste0(par, "_", do.call(paste, c(grid, sep = "_")))
      for (r in seq_len(nrow(grid))) {
        param_cells[[cell_names[r]]] <- as.list(grid[r, , drop = FALSE])
      }
      param_type[cell_names] <- par
    }
  }
  if (length(param_cells) != length(param_type)) {
    stop("Some condition levels contain '_' in a way that gives two different conditions ",
         "the same parameter name. Please recode these levels without '_'.")
  }

  names_a <- names(param_type)[param_type == "a"]
  names_t0 <- names(param_type)[param_type == "t0"]
  names_z <- names(param_type)[param_type == "z"]
  names_v <- names(param_type)[param_type == "v"]
  names_st0 <- names(param_type)[param_type == "st0"]
  names_sz <- names(param_type)[param_type == "sz"]
  names_sv <- names(param_type)[param_type == "sv"]

  # determine how many times a parameter will be estimated
  length_a <- length(names_a)
  length_t0 <- length(names_t0)
  length_z <- length(names_z)
  length_v <- length(names_v)
  length_st0 <- length(names_st0)
  length_sz <- length(names_sz)
  length_sv <- length(names_sv)

  # create lists for each parameter, which contain the parameters starting value n times—n being the number of times it will be estimated
  start_a <- rep(a_fixed, length_a)
  start_t0 <- rep(t0_fixed, length_t0)
  start_z <- rep(z_fixed, length_z)
  start_v <- rep(v_fixed, length_v)
  start_st0 <- rep(st0_fixed, length_st0)
  start_sz <- rep(sz_fixed, length_sz)
  start_sv <- rep(sv_fixed, length_sv)

  # Initialize starting parameter list that will be fed into ddiffusion function
  start <- c(start_a, start_t0, start_z, start_v, start_st0, start_sz, start_sv)
  names(start) <- c(names_a, names_t0, names_z, names_v, names_st0, names_sz, names_sv)
  param_names <- c(names_a, names_t0, names_z, names_v, names_st0, names_sz, names_sv)
  col_names <- param_names

  # SECTION D: DEALING WITH FIXED PARAMETERS

  # check whether a parameter is fixed. If it is, remove it from the starting parameter list
  # Also, create new set of variables "n_fixed", which are either equal to the fixed parameter value or NULL, if parameter is not fixed
  if (fix_a) {
    start <- start[!names(start) %in% "a"]
    param_names <- param_names[!param_names %in% "a"]
  } else {
    a_fixed <- NULL
  }

  if (fix_t0 == TRUE) {
    start <- start[!names(start) %in% "t0"]
    param_names <- param_names[!param_names %in% "t0"]
  } else {
    t0_fixed <- NULL
  }

  if (fix_z == TRUE) {
    start <- start[!names(start) %in% "z"]
    param_names <- param_names[!param_names %in% "z"]
  } else {
    z_fixed <- NULL
  }

  if (fix_v == TRUE) {
    start <- start[!names(start) %in% "v"]
    param_names <- param_names[!param_names %in% "v"]
  } else {
    v_fixed <- NULL
  }

  if (fix_st0 == TRUE) {
    start <- start[!names(start) %in% "st0"]
    param_names <- param_names[!param_names %in% "st0"]
  } else {
    st0_fixed <- NULL
  }

  if (fix_sz == TRUE) {
    start <- start[!names(start) %in% "sz"]
    param_names <- param_names[!param_names %in% "sz"]
  } else {
    sz_fixed <- NULL
  }

  if (fix_sv == TRUE) {
    start <- start[!names(start) %in% "sv"]
    param_names <- param_names[!param_names %in% "sv"]
  } else {
    sv_fixed <- NULL
  }

  # list with fixed parameter values that will be used for naming fixed parameter columns later on
  fixed_par_vals <- c(a_fixed, t0_fixed, z_fixed, v_fixed, st0_fixed, sz_fixed, sv_fixed)

  # SECTION E: GENERATING TRACKING COLUMNS FOR ALL THE STARTING PARAMETERS
  # EACH STARTING PARAMETER GETS A COLUMN IN THE DATASET THAT EQUALS "TRUE", IF THE STARTING PARAMETER APPLIES FOR THAT TRIAL

  split_pars <- function(string) unlist(strsplit(string, "_"))
  list_with_split_pars <- lapply(param_names, split_pars)

  # A trial belongs to a parameter if its value in every condition column the parameter
  # varies across equals the parameter's level (always TRUE for non-varying parameters)
  tracking_cols <- lapply(param_names, function(par_name) {
    cell <- param_cells[[par_name]]
    contains_condition <- rep(TRUE, nrow(data))
    for (column in names(cell)) {
      contains_condition <- contains_condition & data[[column]] == cell[[column]]
    }
    name_i <- paste0("contains_condition_", par_name)
    tibble::as_tibble(setNames(list(contains_condition), name_i))
  })

  df.rtdists <- bind_cols(df.rtdists, tracking_cols)

  # Record which estimated parameter applies in each condition cell (NA if the parameter
  # is fixed). predict() and simulate() use this for each cell.
  cells <- unique(cond)
  first_trial <- match(cells, cond)
  cell_parameters <- data.frame(.cell = cells)
  for (par in names(vary_list)) {
    candidates <- param_names[param_type[param_names] == par]
    cell_parameters[[par]] <- vapply(first_trial, function(trial) {
      applies <- vapply(candidates, function(name) {
        df.rtdists[[paste0("contains_condition_", name)]][trial]
      }, logical(1))
      if (any(applies)) candidates[applies][1] else NA_character_
    }, character(1))
  }

  # SECTION F: CREATING THE WRAPPER FUNCTION FOR RTDISTS

  # Define wrapper calling rtdists. Wrapper will be given to optim to minimize
  # SECTION F: CREATING THE WRAPPER FUNCTION FOR RTDISTS
  ll_fun <- function(pars, RT, response,
                     a_start_val, t0_start_val, z_start_val, v_start_val, st0_start_val, sz_start_val, sv_start_val,
                     precision, stop_on_error,
                     param_names, contains_conditions_list, list_with_split_pars,
                     groups = NULL, p_outlier = 0, w_outlier = 0.1) {

    # Initialize base parameter vectors
    a   <- a_start_val
    t0  <- t0_start_val
    z   <- z_start_val
    v   <- v_start_val
    st0 <- st0_start_val
    sz  <- sz_start_val
    sv  <- sv_start_val

    # Apply condition-specific parameter updates
    for (p in seq_along(list_with_split_pars)) {
      parameter <- list_with_split_pars[[p]]
      condition <- contains_conditions_list[[p]]
      par_name  <- param_names[p]

      if (parameter[1] == "a") {
        a[condition] <- pars[par_name]
      } else if (parameter[1] == "t0") {
        t0[condition] <- pars[par_name]
      } else if (parameter[1] == "z") {
        z[condition] <- pars[par_name]
      } else if (parameter[1] == "v") {
        v[condition] <- pars[par_name]
      } else if (parameter[1] == "st0") {
        st0[condition] <- pars[par_name]
      } else if (parameter[1] == "sz") {
        sz[condition] <- pars[par_name]
      } else if (parameter[1] == "sv") {
        sv[condition] <- pars[par_name]
      }
    }

    # Compute likelihoods. Trials in the same group share all parameter values, so
    # ddiffusion is called once per group with scalar parameters. Passing trial-level
    # parameter vectors instead gives the same densities, but rtdists then has to find
    # the unique parameter combinations on every call, which dominates the fitting time.
    if (is.null(groups)) groups <- list(seq_along(RT))
    dens <- numeric(length(RT))
    for (g in groups) {
      k <- g[1]
      at <- function(x) if (length(x) == 1) x else x[k]
      dens[g] <- ddiffusion(
        rt = RT[g], response = response[g],
        a = at(a), v = at(v), t0 = at(t0), z = at(z) * at(a),
        st0 = at(st0), sz = at(sz), sv = at(sv),
        precision = precision, stop_on_error = stop_on_error
      )
    }

    # Outlier mixture as in HDDM: each trial comes from the diffusion model with probability
    # 1 - p_outlier, and otherwise from a uniform outlier distribution with density w_outlier
    if (p_outlier > 0) {
      dens[is.na(dens) | dens < 0] <- 0
      dens <- (1 - p_outlier) * dens + p_outlier * w_outlier
    }

    # Handle invalid densities. A density of 0 usually means the true density is too
    # small to represent (e.g. an RT just above t0), so it is replaced by the smallest
    # positive double. A larger floor (this used to be 1e-10) caps the cost of a poorly
    # fitted trial, and the optimizer can then improve the fit to the remaining trials
    # by effectively discarding one.
    if (any(is.na(dens) | dens <= 0)) {
      dens[dens <= 0 | is.na(dens)] <- .Machine$double.xmin
    }

    # Return negative log-likelihood
    return(-sum(log(dens)))
  }

  # SECTION G: SETTING UP THE OUTPUT TABLE FOR THE MODEL

  # Initialize estimated parameters tables
  Nparams <- length(start)
  Nfixed <- length(fixed_par_vals)
  paramEsts <- data.frame(matrix(NA, ncol=Nparams+Nfixed+5))
  names(paramEsts) <- c("ID", names(start), names(fixed_par_vals), "converge", "LLE", "AIC", "BIC")
  fit_info <- data.frame(ID = RPlist, optimizer_runs = NA_integer_,
                         converge = NA_integer_, best_replicated = NA, message = NA_character_,
                         se_status = NA_character_)
  seEsts <- data.frame(matrix(NA_real_, nrow = length(RPlist), ncol = Nparams + 1))
  names(seEsts) <- c("ID", names(start))
  vcov_list <- list()

  # SECTION H: FIT THE MODEL FOR EACH PARTICIPANT AND DUMP RESULTS INTO OUTPUT TABLE

  # Create bounds that match the actual parameters being optimized
  lower_bounds <- c()
  upper_bounds <- c()

  if (!fix_a && length_a > 0) {
    lower_bounds <- c(lower_bounds, rep(0.01, length_a))
    upper_bounds <- c(upper_bounds, rep(Inf, length_a))
  }
  if (!fix_t0 && length_t0 > 0) {
    lower_bounds <- c(lower_bounds, rep(0, length_t0))
    upper_bounds <- c(upper_bounds, rep(Inf, length_t0))
  }
  if (!fix_z && length_z > 0) {
    lower_bounds <- c(lower_bounds, rep(0, length_z))
    upper_bounds <- c(upper_bounds, rep(1, length_z))
  }
  if (!fix_v && length_v > 0) {
    lower_bounds <- c(lower_bounds, rep(-Inf, length_v))
    upper_bounds <- c(upper_bounds, rep(Inf, length_v))
  }
  if (!fix_st0 && length_st0 > 0) {
    lower_bounds <- c(lower_bounds, rep(0, length_st0))
    upper_bounds <- c(upper_bounds, rep(Inf, length_st0))
  }
  if (!fix_sz && length_sz > 0) {
    lower_bounds <- c(lower_bounds, rep(0, length_sz))
    upper_bounds <- c(upper_bounds, rep(1, length_sz))
  }
  if (!fix_sv && length_sv > 0) {
    lower_bounds <- c(lower_bounds, rep(0, length_sv))
    upper_bounds <- c(upper_bounds, rep(Inf, length_sv))
  }

  # starting message
  cat("Fitting diffusion model for", length(RPlist), "participants...\n")
  pb <- txtProgressBar(min = 0, max = length(RPlist), style = 3)

  ## Run model across all RPs
  loopCount <-  1
  # Run model
  for (i in RPlist) {

    # read in one participant
    data.to.model <- df.rtdists %>% filter(.data$ID==i)

    # turn the tracking columns into vectors and put them in a list that can be given to the rtdists wrapper
    target_columns <- grep("^contains_condition_", names(data.to.model), value = TRUE)
    contains_conditions_list <- data.to.model %>%
      dplyr::select(all_of(target_columns)) %>%
      dplyr::mutate(across(everything(), as.vector)) %>%
      dplyr::summarise(across(everything(), list)) %>%
      unlist(recursive = FALSE)

    # prepare starting parameters for rtdists function
    ### --- Compute condition-specific starting values (using "upper"/"lower" responses) --- ###

    # ---------- 1) RT units: already converted to seconds during the pre-checks ----------
    data.to.model$RT_sec <- data.to.model$RT

    # ---------- 2) Recode responses: "upper" -> +1, "lower" -> -1 ----------
    data.to.model$resp_num <- dplyr::case_when(
      data.to.model$response == "upper" ~  1,
      data.to.model$response == "lower" ~ -1,
      TRUE ~ NA_real_
    )

    # ---------- 3) Per-condition summaries ----------
    cond_summaries <- data.to.model %>%
      dplyr::group_by(.data$cond) %>%
      dplyr::summarise(
        n_trials = n(),
        meanRT = mean(.data$RT_sec, na.rm = TRUE),
        varRT  = stats::var(.data$RT_sec, na.rm = TRUE),
        minRT  = t0_cap(.data$RT_sec, p_outlier),
        q10RT  = quantile(.data$RT_sec, probs = 0.1, na.rm = TRUE, names = FALSE),
        prop_upper = mean(.data$response == "upper", na.rm = TRUE),
        .groups = "drop"
      )

    # ---------- 4) Starting guesses from EZ-diffusion ----------
    # EZ (Wagenmakers et al., 2007) gives closed-form estimates of v, a and t0 per condition
    # cell. t0 is kept below the fastest RT (the 10% quantile with an outlier component) so
    # every trial has a non-zero likelihood at the start. z starts unbiased, and the
    # variability parameters at small values.
    ez <- ez_diffusion(cond_summaries$prop_upper, cond_summaries$meanRT,
                       cond_summaries$varRT, cond_summaries$n_trials)
    cond_summaries <- cond_summaries %>%
      dplyr::mutate(
        v_guess = pmin(pmax(ez$v, -5), 5),
        a_guess = pmin(pmax(ez$a, 0.3), 5),
        t0_guess = pmin(ez$t0, 0.9 * .data$minRT),
        z_guess = 0.5,
        st0_guess = 0.1,
        sz_guess = 0.1,
        sv_guess = 0.1
      )

    # ---------- 5) Fallback starts where EZ is undefined or implausible ----------
    cond_summaries <- cond_summaries %>%
      dplyr::mutate(
        v_guess = ifelse(is.na(.data$v_guess), 0.5 * sign(.data$prop_upper - 0.5 + 1e-6), .data$v_guess),
        a_guess = ifelse(is.na(.data$a_guess), 1.0, .data$a_guess),
        # a non-positive EZ t0 (e.g. strongly skewed RTs) falls back to the fast-RT floor
        t0_guess = ifelse(is.na(.data$t0_guess) | .data$t0_guess <= 0,
                          pmin(.data$q10RT, 0.9 * .data$minRT), .data$t0_guess)
      )

    # ---------- 6) Join back to trials and create trial-level start vectors ----------
    data.to.model <- data.to.model %>%
      dplyr::left_join(cond_summaries, by = "cond")

    # --- Replicate across trials and respect any fixed flags ---
    a_start_val  <- if (fix_a)  a_fixed else data.to.model$a_guess
    v_start_val  <- if (fix_v)  v_fixed else data.to.model$v_guess
    t0_start_val <- if (fix_t0) t0_fixed else data.to.model$t0_guess
    z_start_val  <- if (fix_z)  z_fixed else data.to.model$z_guess
    st0_start_val <- if (fix_st0) st0_fixed else data.to.model$st0_guess
    sz_start_val  <- if (fix_sz)  sz_fixed else data.to.model$sz_guess
    sv_start_val  <- if (fix_sv)  sv_fixed else data.to.model$sv_guess

    # ---------- 7) Build the optimizer's starting vector and bounds ----------
    # Each estimated parameter applies to the trials flagged in its tracking column.
    # Its start is the mean guess over those trials (the minimum for t0, which must
    # stay below every RT it applies to). Falls back to the default start if no
    # usable guess exists, and is kept inside the optimizer bounds.
    # t0 is the earliest possible response, so any trial faster than t0 has zero
    # likelihood. Its upper bound is therefore set just below the fastest RT it
    # applies to; otherwise the optimizer can get stuck at a t0 that writes off
    # the fastest trials.
    start_i <- start
    upper_i <- upper_bounds
    for (p in seq_along(param_names)) {
      par_name <- param_names[p]
      par_type <- list_with_split_pars[[p]][1]
      trials <- contains_conditions_list[[p]]
      if (par_type == "t0" && any(trials)) {
        upper_i[p] <- 0.999 * t0_cap(data.to.model$RT[trials], p_outlier)
      }
      guesses <- data.to.model[[paste0(par_type, "_guess")]][trials]
      guesses <- guesses[is.finite(guesses)]
      if (length(guesses) > 0) {
        start_i[par_name] <- if (par_type == "t0") min(guesses) else mean(guesses)
      }
    }
    start_i <- pmin(pmax(start_i, lower_bounds), upper_i)

    # A parameter that applies to no trials of this participant (e.g. a condition the
    # participant has no data for) cannot be estimated. It is left out of this
    # participant's fit and reported as NA, so it neither stays at its starting value nor
    # makes the Hessian singular for the other parameters.
    estimable <- vapply(contains_conditions_list, any, logical(1))
    fit_names <- param_names[estimable]
    fit_splits <- list_with_split_pars[estimable]
    fit_conditions <- contains_conditions_list[estimable]
    lower_i <- lower_bounds[estimable]
    upper_i <- upper_i[estimable]
    start_i <- start_i[estimable]

    ### --- Diagnostic: Print starting values per condition --- ###

    if (print_starting_vals == TRUE) {
      cat("\n======================================\n")
      cat(paste0("Participant ", i, " -- Summary:\n"))
      cat("======================================\n")
      print(cond_summaries)
      cat("\nStarting values passed to the optimizer:\n")
      print(round(start_i, 3))
      cat("======================================\n\n")
    }

    # Trials that share all parameter values (same condition cell and the same tracking
    # pattern), so the likelihood can be computed once per group
    trial_groups <- unname(split(seq_len(nrow(data.to.model)), interaction(
      data.to.model$cond, do.call(paste0, lapply(fit_conditions, as.integer)),
      drop = TRUE)))

    # fit the model from a given starting vector
    run_fit <- function(start_vals) {
      nlminb(start_vals, objective = ll_fun,
             RT = data.to.model[["RT"]], response = data.to.model[["response"]],
             a_start_val = a_start_val, t0_start_val = t0_start_val, z_start_val = z_start_val, v_start_val = v_start_val, st0_start_val = st0_start_val, sz_start_val = sz_start_val, sv_start_val = sv_start_val,
             precision = precision, stop_on_error = stop_on_error,
             param_names = fit_names, contains_conditions_list = fit_conditions, list_with_split_pars = fit_splits,
             groups = trial_groups, p_outlier = p_outlier, w_outlier = w_outlier,
             control = list(iter.max = 10000, eval.max = 10000),
             lower = lower_i,
             upper = upper_i)
    }

    # Alternative starting vector number k (k >= 1), spread over the plausible parameter
    # space with a Halton sequence (one prime base per parameter type, so all conditions of
    # a parameter move together and the EZ pattern across conditions is kept). a and v are
    # scaled by 0.5-2 times their EZ start; the others cover their usual range. Deterministic,
    # so results are reproducible and the user's random seed is untouched.
    par_types <- vapply(fit_splits, `[`, character(1), 1)
    alternative_start <- function(k) {
      s <- start_i
      primes <- c(a = 2, t0 = 3, z = 5, v = 7, st0 = 11, sz = 13, sv = 17)
      for (p in seq_along(fit_names)) {
        type <- par_types[p]
        u <- halton(k, primes[[type]])
        s[p] <- switch(type,
          a   = start_i[p] * 2^(2 * u - 1),
          v   = start_i[p] * 2^(2 * u - 1),
          t0  = upper_i[p] * (0.2 + 0.75 * u),
          z   = 0.3 + 0.4 * u,
          st0 = 0.3 * u,
          sv  = 2.5 * u,
          sz  = {
            # sz must stay below 2 * min(z, a - z) on the absolute scale to be valid
            a_mean <- if (any(par_types == "a")) mean(s[par_types == "a"]) else a_fixed
            z_mean <- if (any(par_types == "z")) mean(s[par_types == "z"]) else z_fixed
            0.9 * u * 2 * min(z_mean, 1 - z_mean) * a_mean
          })
      }
      pmin(pmax(s, lower_i), upper_i)
    }

    # Starting vector that moves sv and sz estimates stuck at 0 to clearly positive values.
    # sv and sz spread drift and starting point symmetrically, so near 0 they change the
    # likelihood only quadratically: its slope at exactly 0 is zero, and the optimizer can
    # report convergence there even when a positive value fits better.
    nudged_start <- function(pars) {
      s <- pars
      s[par_types == "sv" & s <= 1e-3] <- 1
      if (any(par_types == "sz" & s <= 1e-3)) {
        a_mean <- if (any(par_types == "a")) mean(s[par_types == "a"]) else a_fixed
        z_mean <- if (any(par_types == "z")) mean(s[par_types == "z"]) else z_fixed
        s[par_types == "sz" & s <= 1e-3] <- 0.5 * 2 * min(z_mean, 1 - z_mean) * a_mean
      }
      pmin(pmax(s, lower_i), upper_i)
    }
    stuck_at_zero <- function(f) any(par_types %in% c("sv", "sz") & f$par <= 1e-3)

    fit <- run_fit(start_i)
    n_runs <- 1
    replicated <- NA

    # Retry (up to max_retries further optimizer runs) if the fit did not converge or an sv
    # or sz estimate is at 0 (below 0.001; how close the optimizer gets to exactly 0 differs
    # between platforms):
    # 1) if it did not converge, restart the optimizer from where it stopped, which often
    #    completes a fit that stopped early;
    # 2) fit from further starting points (first with sv/sz stuck at 0 moved to positive
    #    values, then alternative starting values) until the best converged solution has
    #    been reached from two independent starting points, which indicates it is the
    #    global optimum rather than a local one.
    # The fit with the highest likelihood is kept.
    if ((fit$convergence != 0 || stuck_at_zero(fit)) && max_retries > 0) {
      if (fit$convergence != 0) {
        restarted <- run_fit(fit$par)
        n_runs <- n_runs + 1
        if (restarted$objective <= fit$objective) fit <- restarted
      }

      start_optima <- fit$objective           # best solution from each independent start
      start_converged <- fit$convergence == 0
      replicated <- FALSE
      tried_nudge <- FALSE
      k <- 0
      while (n_runs <= max_retries && !replicated) {
        if (!tried_nudge && stuck_at_zero(fit)) {
          next_start <- nudged_start(fit$par)
          tried_nudge <- TRUE
        } else {
          k <- k + 1
          next_start <- alternative_start(k)
        }
        candidate <- run_fit(next_start)
        n_runs <- n_runs + 1
        start_optima <- c(start_optima, candidate$objective)
        start_converged <- c(start_converged, candidate$convergence == 0)
        if (candidate$objective < fit$objective) fit <- candidate
        replicated <- fit$convergence == 0 &&
          sum(start_converged & abs(start_optima - fit$objective) < 0.01) >= 2
      }
    }

    # ---------- Standard errors ----------
    # SEs are the square roots of the diagonal of the inverse Hessian of the negative
    # log-likelihood at the estimates (the observed Fisher information), computed numerically
    # with Richardson extrapolation (numDeriv). The Hessian is evaluated with at least
    # precision 4, because the less precise numerical integration at precision 3 noticeably
    # distorts the second derivatives when sv or sz are estimated.
    # Wald SEs require the estimate to lie inside its bounds. Parameters at (within 0.001 of)
    # a bound, e.g. sv = 0 or t0 at the fastest RT, get no SE, and the SEs of the other
    # parameters are computed with them held at their estimates.
    se <- setNames(rep(NA_real_, length(fit_names)), fit_names)
    vcov_i <- NULL
    se_status <- "not computed"
    if (compute_se) {
      at_bound <- (fit$par - lower_i < 1e-3) | (upper_i - fit$par < 1e-3)
      free <- which(!at_bound)
      bound_note <- if (any(at_bound)) {
        paste0("no SE at bound: ", paste(fit_names[at_bound], collapse = ", "))
      } else NULL
      se_status <- if (is.null(bound_note)) "ok" else bound_note
      if (length(free) > 0) {
        nll_free <- function(x) {
          p <- fit$par
          p[free] <- x
          ll_fun(p, RT = data.to.model[["RT"]], response = data.to.model[["response"]],
                 a_start_val = a_start_val, t0_start_val = t0_start_val, z_start_val = z_start_val, v_start_val = v_start_val, st0_start_val = st0_start_val, sz_start_val = sz_start_val, sv_start_val = sv_start_val,
                 precision = max(precision, 4), stop_on_error = stop_on_error,
                 param_names = fit_names, contains_conditions_list = fit_conditions, list_with_split_pars = fit_splits,
                 groups = trial_groups, p_outlier = p_outlier, w_outlier = w_outlier)
        }
        # numDeriv's largest step is d * |x| (d = 0.1 by default). A step past a bound, e.g.
        # t0 beyond the fastest RT, gives trials a likelihood of 0 and a meaningless Hessian,
        # so d is reduced to keep every step at most halfway to the nearest bound.
        x_free <- fit$par[free]
        gap <- pmin(x_free - lower_i[free], upper_i[free] - x_free)
        rel_gap <- gap / pmax(abs(x_free), 1e-4)
        step_d <- min(0.1, 0.5 * min(rel_gap))
        H <- tryCatch(numDeriv::hessian(nll_free, x_free, method.args = list(d = step_d)),
                      error = function(e) NULL)
        if (is.null(H) || any(!is.finite(H))) {
          se_status <- paste(c("Hessian could not be computed", bound_note), collapse = "; ")
        } else {
          H <- (H + t(H)) / 2
          if (min(eigen(H, symmetric = TRUE, only.values = TRUE)$values) <= 0) {
            # not a proper maximum (or a parameter is not identified): no valid SEs
            se_status <- paste(c("Hessian not positive definite", bound_note), collapse = "; ")
          } else {
            vcov_i <- solve(H)
            dimnames(vcov_i) <- list(fit_names[free], fit_names[free])
            se[free] <- sqrt(diag(vcov_i))
          }
        }
      }
    }
    seEsts[loopCount, "ID"] <- i
    seEsts[loopCount, fit_names] <- se
    vcov_list[as.character(i)] <- list(vcov_i)

    fit_info[loopCount, ] <- list(i, n_runs, fit$convergence, replicated, fit$message, se_status)

    # save to empirical_ests
    parms <- fit$par
    paramEsts[loopCount,"ID"] <- i
    for (Est in param_names) {
      paramEsts[loopCount, Est] <- parms[Est]
    }
    for (fixed in names(fixed_par_vals)) {
      paramEsts[loopCount, fixed] <- fixed_par_vals[fixed]
    }
    paramEsts[loopCount,"converge"] <- fit$convergence
    paramEsts[loopCount,"LLE"] <- -1*fit$objective  #LLE
    paramEsts[loopCount,"AIC"] <- length(fit$par)*2-2*(-1*fit$objective)
    paramEsts[loopCount,"BIC"] <- -2*(-1*fit$objective) + length(fit$par)*log(nrow(data.to.model))

    #print(loopCount)
    loopCount = loopCount+1
    setTxtProgressBar(pb, i)
  }

  close(pb)

  # SECTION I: ORGANIZE OUTPUT TABLE AND CREATE SUMMARY TABLE

  # Report participants that needed retries or still did not converge
  n_retried <- sum(fit_info$optimizer_runs > 1)
  n_failed <- sum(fit_info$converge != 0)
  if (n_retried > 0) {
    message(n_retried, " participant(s) needed more than one optimizer run; ", n_failed,
            " did not converge after retries. See results$fit_info for details.")
  }
  names(fit_info)[names(fit_info) == "ID"] <- sub_id
  seEsts$ID <- fit_info[[sub_id]]
  names(seEsts)[names(seEsts) == "ID"] <- sub_id

  # Reorder the columns and rename subject ID if it was specified as custom
  desired_order <- c("ID", col_names, "converge", "LLE", "AIC", "BIC")
  paramEsts <- paramEsts[, desired_order]
  names(paramEsts)[names(paramEsts) == "ID"] <- sub_id

  # Mean and SD across participants of each parameter (fit statistics such as the
  # log-likelihood are reported per participant, not averaged)
  means <- round(colMeans(paramEsts[, col_names, drop = FALSE], na.rm = TRUE), digits = 2)
  sd <- round(apply(paramEsts[, col_names, drop = FALSE], 2, sd, na.rm = TRUE), digits = 2)
  # summary table
  summary <- data.frame(
    Parameter = col_names,
    Mean = means,
    SD = sd,
    Fixed = !(col_names %in% param_names)
  )
  table <- summary %>%
    gt() %>%
    tab_spanner(label = "Model Summary", columns = c("Mean", "SD")) %>%
    fmt_number(columns = c("Mean", "SD"), decimals = 2)
  print(table)

  # Store the processed data (in seconds) with the original condition columns, used by
  # predict() and simulate()
  processed_data <- df.rtdists[, 1:4]
  names(processed_data)[names(processed_data) == "ID"] <- S
  names(processed_data)[names(processed_data) == "response"] <- R
  names(processed_data)[names(processed_data) == "RT"] <- React
  names(processed_data)[names(processed_data) == "cond"] <- ".cell"
  processed_data <- cbind(processed_data, data[conditions])

  # One row per parameter column: its type, whether it is fixed, and the condition values it
  # applies to (NA for condition columns it does not vary across). Used by estimates().
  parameter_info <- data.frame(name = col_names, parameter = unname(param_type[col_names]),
                               fixed = !(col_names %in% param_names), stringsAsFactors = FALSE)
  for (column in conditions) {
    values <- data[[column]][rep(NA_integer_, length(col_names))]
    for (k in seq_along(col_names)) {
      cell <- param_cells[[col_names[k]]]
      if (!is.null(cell[[column]])) values[k] <- cell[[column]]
    }
    parameter_info[[column]] <- values
  }

  # One row per condition cell: its label (.cell) and its condition values
  cells_table <- unique(cbind(data.frame(.cell = cond), data[conditions]))
  rownames(cells_table) <- NULL

  results <- list(
    individual_estimates = paramEsts,
    summary = summary,  # Your summary table
    settings = list(
      varying_parameters = list(
        a = vary_a, t0 = vary_t0, z = vary_z, v = vary_v,
        st0 = vary_st0, sz = vary_sz, sv = vary_sv
      ),
      fixed_parameters = list(
        a = fix_a, t0 = fix_t0, z = fix_z, v = fix_v,
        st0 = fix_st0, sz = fix_sz, sv = fix_sv
      ),
      n_participants = length(RPlist),
      n_conditions = length(unique(df.rtdists$cond)),
      cell_parameters = cell_parameters,
      cells = cells_table,
      parameters = parameter_info,
      model = model_spec,
      condition_columns = conditions,
      columns = list(sub_id = S, response = R, rt = React),
      precision = precision,
      p_outlier = p_outlier,
      w_outlier = w_outlier
    ),
    data = processed_data,
    convergence = table(paramEsts$converge),
    fit_info = fit_info,
    standard_errors = seEsts,
    vcov = vcov_list
  )

  class(results) <- "driftd_results"
  return(results)

}
