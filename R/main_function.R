
#' @title Drift-Diffusion Model Fitting
#' @description Fits a drift-diffusion model (DDM) to response time (RT) and response data for individual participants using maximum likelihood estimation via \code{nlminb}.
#'
#' @param data A data frame containing the experimental data.
#' @param sub_id Character string specifying the column name for subject/participant ID. Default is "sub_id".
#' @param response Character string specifying the column name for responses (coded as 0/1 or 'lower'/'upper'). Default is "response".
#' @param rt Character string specifying the column name for reaction times. Default is "rt".
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
#' @param sz_fixed Numeric: Fixed value for \code{sz} if \code{fix_sz=TRUE}. Default is 0.
#' @param sv_fixed Numeric: Fixed value for \code{sv} if \code{fix_sv=TRUE}. Default is 0.
#' @param precision Numeric: Precision parameter for \code{rtdists::ddiffusion}. Default is 3.
#' @param stop_on_error Logical: Passed to \code{rtdists::ddiffusion}. Default is FALSE.
#' @param print_starting_vals Logical: If TRUE, prints the derived starting parameter estimates for each participant. Default is FALSE.
#'
#' @details
#' The function estimates DDM parameters (\code{a, t0, z, v, st0, sz, sv}) for each participant individually, allowing parameters to vary across experimental conditions specified in the \code{vary_} arguments. It uses the maximum likelihood method with the \code{rtdists::ddiffusion} function. Note that \code{z} is estimated as a proportion (0 to 1) and internally converted to \code{z * a} for \code{ddiffusion}. It includes pre-checks for data quality and suggests initial parameter guesses based on participant-condition RT summaries.
#'
#' @return A list of class \code{driftd_results} containing:
#' \item{individual_estimates}{A data frame with estimated parameters, convergence status, Log-Likelihood (LLE), AIC, and BIC for each participant.}
#' \item{summary}{A data frame with mean and standard deviation of all estimated parameters across participants.}
#' \item{settings}{A list summarizing the model structure (varying and fixed parameters).}
#' \item{data}{The processed data frame used for fitting.}
#' \item{convergence}{A table summarizing the convergence codes from \code{nlminb}.}
#' @export
#'
#' @importFrom rtdists ddiffusion
#' @importFrom dplyr %>% group_by summarise n filter select mutate case_when left_join bind_cols across everything all_of
#' @importFrom rlang sym
#' @importFrom tibble tibble
#' @importFrom gt gt tab_spanner fmt_number
#' @importFrom stats nlminb sd quantile relevel
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

                   print_starting_vals = FALSE

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

  # Validate parameter bounds
  if (a_fixed <= 0) stop("a_fixed must be positive")
  if (t0_fixed < 0) stop("t0_fixed cannot be negative")
  if (z_fixed < 0 || z_fixed > 1) stop("z_fixed must be between 0 and 1")
  if (st0_fixed < 0) stop("st0_fixed cannot be negative")
  if (sz_fixed < 0 || sz_fixed > 1) stop("sz_fixed must be between 0 and 1")
  if (sv_fixed < 0) stop("sv_fixed cannot be negative")

  run_pre_checks <- function(data, sub_id, vary_a, vary_t0, vary_z, vary_v,
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

        for (i in 1:nrow(by_condition)) {
          message("  Condition '", by_condition$temp_cond[i], "': ", by_condition$participants[i])
        }

        message("This may cause estimation problems.")
        warnings_issued <- TRUE
      }

      # Remove temporary column
      data$temp_cond <- NULL
    }

    # Participant checks (total trials per participant)
    participant_counts <- table(data[[sub_id]])
    small_participants <- participant_counts[participant_counts < 20]
    if (length(small_participants) > 0) {
      message("WARNING: Some participants have very few total trials (<20): ",
              paste(names(small_participants), " (n=", small_participants, ")", sep = "", collapse = ", "),
              ". This may cause estimation problems.")
      warnings_issued <- TRUE
    }

    return(warnings_issued)
  }

  # Run the pre-checks
  cat("Running pre-fit checks...\n")
  warnings_issued <- run_pre_checks(data, sub_id, vary_a, vary_t0, vary_z, vary_v,
                                    vary_st0, vary_sz, vary_sv, fix_st0, st0_fixed,
                                    fix_sz, sz_fixed)

  # Ask for confirmation
  if (warnings_issued && interactive()) {
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
  df.rtdists$response <- factor(df.rtdists$response) # define 1 = upper boundary, 0 = lower boundary
  if (length(levels(df.rtdists$response)) > 1) {     # Check if the factor has more than 1 level before releveling
    df.rtdists$response <- relevel(df.rtdists$response, ref="lower")
  }

  # RP List. Will be used to loop over RPs later on
  RPlist <- unique(df.rtdists$ID)
  #RPlist <- RPlist[RPlist != 77]

  # SECTION B: INITIALIZING FUNCTION TO GENERATE NAMES

  #order vary lists alphabetically
  vary_a <- sort(vary_a)
  vary_t0 <- sort(vary_t0)
  vary_z <- sort(vary_z)
  vary_v <- sort(vary_v)
  vary_st0 <- sort(vary_st0)
  vary_sz <- sort(vary_sz)
  vary_sv <- sort(vary_sv)

  # Function to generate names. This function uses the levels list to generate cell names, based on the condition names that are given to it.
  generate_names <- function(levels_list, index = 1, current = NULL) {
    # Base case: if all lists have been processed, return the current combination
    if (index > length(levels_list)) {
      return(paste(current, collapse = "_"))
    }

    # Initialize an empty vector to store combinations for the current list
    list_names <- c()

    # Iterate over each element in the current list
    for (item in levels_list[[index]]) {
      # Recursive call to generate combinations for the next list
      list_names <- c(list_names, generate_names(levels_list, index + 1, c(current, item)))
    }

    return(list_names)
  }

  # SECTION C: GENERATING STARTING PARAMETER LIST THAT WILL BE FED TO NLMINB

  # for each parameter, see which conditions it varies across and generate names for each condition under which parameter will be estimated
  names_a <- generate_names(levels_list[vary_a])
  names_t0 <- generate_names(levels_list[vary_t0])
  names_z <- generate_names(levels_list[vary_z])
  names_v <- generate_names(levels_list[vary_v])
  names_st0 <- generate_names(levels_list[vary_st0])
  names_sz <- generate_names(levels_list[vary_sz])
  names_sv <- generate_names(levels_list[vary_sv])

  # if a parameter doesn't vary across any conditions, simply give it its name (since it will only be estimated once)
  names_a <- ifelse(names_a == "", "a", paste0("a_", names_a))
  names_t0 <- ifelse(names_t0 == "", "t0", paste0("t0_", names_t0))
  names_z <- ifelse(names_z == "", "z", paste0("z_", names_z))
  names_v <- ifelse(names_v == "", "v", paste0("v_", names_v))
  names_st0 <- ifelse(names_st0 == "", "st0", paste0("st0_", names_st0))
  names_sz <- ifelse(names_sz == "", "sz", paste0("sz_", names_sz))
  names_sv <- ifelse(names_sv == "", "sv", paste0("sv_", names_sv))

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

  # take the cond column and split it into the individual strings that were separated by "_". Will be used to create tracking columns
  cond_words <- strsplit(cond, "_")

  split_pars <- function(string) unlist(strsplit(string, "_"))
  list_with_split_pars <- lapply(param_names, split_pars)

  tracking_cols <- lapply(seq_along(list_with_split_pars), function(i) {
    param <- list_with_split_pars[[i]]
    if (length(param) > 1) {
      contains_condition <- sapply(cond_words, function(x) all(param[-1] %in% x))
    } else {
      contains_condition <- rep(TRUE, length(cond_words))
    }
    name_i <- paste0("contains_condition_", param_names[i])
    tibble::as_tibble(setNames(list(contains_condition), name_i))
  })

  df.rtdists <- bind_cols(df.rtdists, tracking_cols)

  # SECTION F: CREATING THE WRAPPER FUNCTION FOR RTDISTS

  # Define wrapper calling rtdists. Wrapper will be given to optim to minimize
  # SECTION F: CREATING THE WRAPPER FUNCTION FOR RTDISTS
  ll_fun <- function(pars, RT, response,
                     a_start_val, t0_start_val, z_start_val, v_start_val, st0_start_val, sz_start_val, sv_start_val,
                     precision, stop_on_error,
                     param_names, contains_conditions_list, list_with_split_pars) {

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

    # Compute likelihoods
    dens <- ddiffusion(
      rt = RT, response = response,
      a = a, v = v, t0 = t0, z = z * a,
      st0 = st0, sz = sz, sv = sv,
      precision = precision, stop_on_error = stop_on_error
    )

    # Handle invalid densities
    if (any(is.na(dens) | dens <= 0)) {
      # Penalize but allow optimizer to move
      dens[dens <= 0 | is.na(dens)] <- 1e-10
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

    current_id <- RPlist[i]

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

    # ---------- 1) RT units: convert to seconds if likely ms ----------
    RT_raw <- data.to.model$RT
    if (mean(RT_raw, na.rm = TRUE) > 10) {
      RT_sec <- RT_raw / 1000
    } else {
      RT_sec <- RT_raw
    }
    data.to.model$RT_sec <- RT_sec

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
        sdRT   = sd(.data$RT_sec, na.rm = TRUE),
        minRT  = min(.data$RT_sec, na.rm = TRUE),
        meanRT_upper = mean(.data$RT_sec[.data$response == "upper"], na.rm = TRUE),
        meanRT_lower = mean(.data$RT_sec[.data$response == "lower"], na.rm = TRUE),
        prop_upper = mean(.data$response == "upper", na.rm = TRUE),
        .groups = "drop"
      )

    # ---------- 4) Compute start guesses (no correctness required) ----------
    # Explanation of v_start:
    #  - prop_upper - 0.5 gives direction and (rough) magnitude of bias
    #  - RT difference (lower - upper) / sdRT gives an RT-based cue (positive -> faster upper)
    #  - combine both with scaling constants, and clamp to [-3,3]
    cond_summaries <- cond_summaries %>%
      dplyr::mutate(
        # fallback handling for sdRT or meanRT_* NA
        # ----- choose sign from proportion (fallback: RT diff sign) -----
        sign_prop = dplyr::case_when(
          is.na(.data$prop_upper) ~ 0,
          .data$prop_upper > 0.5  ~  1,
          .data$prop_upper < 0.5  ~ -1,
          TRUE                    ~  0  # near 0.5
        ),

        sign_rt = sign(.data$meanRT_lower - .data$meanRT_upper), # positive if upper is faster

        # ----- magnitude components -----
        # component from proportion: maps [0.5..1] to [min_prop_mag..max_prop_mag]
        prop_mag = (abs(.data$prop_upper - 0.5)) * (3.0),   # scale_prop

        # component from RT: shorter RTs -> larger magnitude.
        # use inverse meanRT (normalized) as a cue; avoid div by zero
        rt_mag = ifelse(is.na(.data$meanRT) | .data$meanRT <= 0, 0,
                        (1 / .data$meanRT) * 0.4),             # scale_rt

        # combine magnitude (weighted sum)
        v_raw_mag = .data$prop_mag + .data$rt_mag,

        # ----- final sign and clamp -----
        chosen_sign = ifelse(.data$sign_prop == 0, .data$sign_rt, .data$sign_prop),
        v_signed = .data$chosen_sign * .data$v_raw_mag,
        v_guess = pmin(pmax(.data$v_signed, -5), 5),  # clamp to [-5, 5]
        v_guess = .data$v_guess * 0.6,

        # a_start: boundary ~ 2 * sdRT (clamped to reasonable range)
        a_guess = pmin(pmax(.data$sdRT * 2.5, 0.3), 3),
        a_guess = .data$a_guess * 1.2,

        # t0_start: floor near fast RTs (10th percentile is robust)
        t0_guess = pmin(pmax(quantile(.data$RT_sec, probs = 0.1, na.rm = TRUE), 0.05), 1),

        # default small variability starts
        z_guess = 0.5,
        st0_guess = 0.1,
        sz_guess = 0.1,
        sv_guess = 0.1
      ) %>%
      dplyr::select(.data$cond, .data$n_trials, .data$meanRT, .data$sdRT, .data$prop_upper,
                    .data$v_guess, .data$a_guess, .data$t0_guess, .data$z_guess,
                    .data$st0_guess, .data$sz_guess, .data$sv_guess)

    # ---------- 5) If any conditions are very sparse, give a fallback start ----------
    # e.g., if n_trials < 10 or prop_upper is NA, set sensible defaults
    cond_summaries <- cond_summaries %>%
      dplyr::mutate(
        v_guess = ifelse(is.na(.data$v_guess) | .data$n_trials < 8,
                         0.5 * sign(.data$prop_upper - 0.5 + 1e-6),
                         .data$v_guess),
        a_guess = ifelse(is.na(.data$a_guess) | .data$n_trials < 8, 1.0, .data$a_guess),
        t0_guess = ifelse(is.na(.data$t0_guess) | .data$n_trials < 8, 0.25, .data$t0_guess)
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

    ### --- Diagnostic: Print starting values per condition --- ###

    if (print_starting_vals == TRUE) {
      cat("\n======================================\n")
      cat(paste0("Participant ", i, " -- Summary:\n"))
      cat("======================================\n")
      print(cond_summaries)
      cat("======================================\n\n")
    }

    # fit the model
    fit <- nlminb(start, objective = ll_fun,
                  RT = data.to.model[["RT"]], response = data.to.model[["response"]],
                  a_start_val = a_start_val, t0_start_val = t0_start_val, z_start_val = z_start_val, v_start_val = v_start_val, st0_start_val = st0_start_val, sz_start_val = sz_start_val, sv_start_val = sv_start_val,
                  precision = precision, stop_on_error = stop_on_error,
                  param_names = param_names, contains_conditions_list = contains_conditions_list, list_with_split_pars = list_with_split_pars,
                  control = list(iter.max = 10000, eval.max = 10000),
                  lower = lower_bounds,
                  upper = upper_bounds)
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

  # Reorder the columns and rename subject ID if it was specified as custom
  desired_order <- c("ID", col_names, "converge", "LLE", "AIC", "BIC")
  paramEsts <- paramEsts[, desired_order]
  names(paramEsts)[names(paramEsts) == "ID"] <- sub_id

  # calculate means and sd for parameters of model
  means <- colMeans(paramEsts[, -1], na.rm = TRUE)
  means <- round(means, digits=2)
  sd <- apply(paramEsts[, -1], 2, sd, na.rm = TRUE)
  sd <- round(sd, digits=2)
  # summary table
  summary <- data.frame(
    Parameter = names(means),
    Mean = means,
    SD = sd
  )
  table <- summary %>%
    gt() %>%
    tab_spanner(label = "Model Summary", columns = c("Mean", "SD")) %>%
    fmt_number(columns = c("Mean", "SD"), decimals = 2)
  print(table)

  # Store the processed data for posterior predictive checks
  processed_data <- df.rtdists[, 1:4]
  names(processed_data)[names(processed_data) == "ID"] <- S
  names(processed_data)[names(processed_data) == "response"] <- R
  names(processed_data)[names(processed_data) == "RT"] <- React

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
      n_conditions = length(unique(df.rtdists$cond))
    ),
    data = processed_data,
    convergence = table(paramEsts$converge)
  )

  class(results) <- "driftd_results"
  return(results)

}
