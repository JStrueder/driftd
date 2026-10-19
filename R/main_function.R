
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
  # ---- 1. Check the input ------------------------------------------------------------------

  if (!is.data.frame(data)) stop("'data' must be a data frame")
  missing_cols <- setdiff(c(sub_id, response, rt), names(data))
  if (length(missing_cols) > 0) {
    stop("The following specified columns were not found in data: ",
         paste(missing_cols, collapse = ", "))
  }
  if (any(is.na(data[[sub_id]]))) stop("NA values found in subject ID column.")
  if (any(is.na(data[[response]]))) stop("NA values found in response column.")
  if (any(is.na(data[[rt]]))) stop("NA values found in RT column.")
  if (any(data[[rt]] <= 0)) stop("All RT values must be positive")

  # rtdists works in seconds. RTs averaging over 10 are assumed to be milliseconds.
  if (mean(data[[rt]]) > 10) {
    message("RTs appear to be in milliseconds (mean RT = ", round(mean(data[[rt]]), 1),
            "). Converting to seconds; all estimates and simulated RTs are in seconds.")
    data[[rt]] <- data[[rt]] / 1000
  }

  # A parameter can either be fixed to one value or vary across conditions, not both
  vary_list <- list(a = vary_a, t0 = vary_t0, z = vary_z, v = vary_v,
                    st0 = vary_st0, sz = vary_sz, sv = vary_sv)
  fixed <- c(a = fix_a, t0 = fix_t0, z = fix_z, v = fix_v,
             st0 = fix_st0, sz = fix_sz, sv = fix_sv)
  conflicting <- names(fixed)[fixed & lengths(vary_list) > 0]
  if (length(conflicting) > 0) {
    stop("These parameters are both fixed and set to vary across conditions: ",
         paste(conflicting, collapse = ", "), ". Set fix_", conflicting[1],
         " = FALSE to estimate it per condition, or remove vary_", conflicting[1], ".")
  }

  # The condition columns must exist and have no missing values
  vary_cols <- unique(unlist(vary_list))
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

  cat("Running pre-fit checks...\n")
  warnings_issued <- run_pre_checks(data, sub_id, response, rt, vary_list,
                                    fix_st0, st0_fixed, fix_sz, sz_fixed)
  if (warnings_issued && ask && interactive()) {
    cat("\n")
    msg <- readline("Do you want to continue with model fitting? (y/n): ")
    if (tolower(msg) != "y") {
      message("Model fitting cancelled by user.")
      return(invisible(NULL))
    }
    cat("Continuing with model fitting...\n\n")
  }

  # ---- 2. Trial data in rtdists' format ----------------------------------------------------

  # The condition columns any parameter varies across, and their levels
  conditions <- unlist(vary_list)
  conditions <- unique(conditions[!is.na(conditions) & conditions != ""])
  levels_list <- lapply(conditions, function(x) sort(unique(data[[x]])))
  names(levels_list) <- conditions

  # The condition cell of each trial, labelled by its levels joined with "_"
  if (length(conditions) < 1) {
    cond <- rep("NA", nrow(data))
  } else {
    cond <- do.call(paste, c(data[conditions], sep = "_"))
    if (length(unique(cond)) != nrow(unique(data[conditions]))) {
      stop("Some condition levels contain '_' in a way that gives two different conditions ",
           "the same label. Please recode these levels without '_'.")
    }
  }

  resp <- data[[response]]
  if (all(resp %in% c(0, 1))) {
    resp <- ifelse(resp == 1, "upper", "lower")
  } else if (!all(resp %in% c("upper", "lower"))) {
    stop("Response variable must be coded as 0/1 or 'lower'/'upper'")
  }
  trials <- data.frame(ID = data[[sub_id]], cond = cond, response = resp, RT = data[[rt]])
  # rtdists reads a factor response by its level number (1 = lower, 2 = upper), so both
  # levels must always exist, even when every response is at the same boundary
  trials$response <- factor(trials$response, levels = c("lower", "upper"))
  ids <- unique(trials$ID)

  # ---- 3. The parameters and the trials each applies to ------------------------------------

  # Each parameter is estimated once per level combination of the condition columns it
  # varies across ("v_easy", "v_hard", ...), or once if it varies across none ("a").
  # param_cells keeps each parameter's level values, so trials are matched by value.
  vary_list <- lapply(vary_list, sort)
  param_cells <- list()
  param_type <- c()
  for (par in names(vary_list)) {
    if (length(vary_list[[par]]) == 0) {
      param_cells[par] <- list(list())
      param_type[par] <- par
    } else {
      # all level combinations, one row each; the first column varies slowest
      grid <- expand.grid(rev(levels_list[vary_list[[par]]]), KEEP.OUT.ATTRS = FALSE,
                          stringsAsFactors = FALSE)[vary_list[[par]]]
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

  defaults <- c(a = a_fixed, t0 = t0_fixed, z = z_fixed, v = v_fixed,
                st0 = st0_fixed, sz = sz_fixed, sv = sv_fixed)
  col_names <- names(param_type)                         # every parameter, fixed or not
  param_names <- col_names[!fixed[param_type]]           # the estimated parameters
  par_types <- unname(param_type[param_names])
  # default starting values, used where a participant's data give no starting guess
  start <- setNames(defaults[par_types], param_names)
  fixed_par_vals <- defaults[fixed]

  # For each estimated parameter, the trials it applies to: those whose value in every
  # condition column the parameter varies across equals the parameter's level
  applies <- lapply(param_names, function(par_name) {
    cell <- param_cells[[par_name]]
    in_cell <- rep(TRUE, nrow(data))
    for (column in names(cell)) {
      in_cell <- in_cell & data[[column]] == cell[[column]]
    }
    in_cell
  })
  names(applies) <- param_names

  # Which estimated parameter applies in each condition cell (NA if the parameter is
  # fixed). predict() and simulate() use this for each cell.
  cells <- unique(cond)
  first_trial <- match(cells, cond)
  cell_parameters <- data.frame(.cell = cells)
  for (par in names(vary_list)) {
    candidates <- param_names[par_types == par]
    cell_parameters[[par]] <- vapply(first_trial, function(trial) {
      hit <- vapply(candidates, function(name) applies[[name]][trial], logical(1))
      if (any(hit)) candidates[hit][1] else NA_character_
    }, character(1))
  }

  # Optimizer bounds. t0's upper bound is set per participant (see fit_participant()).
  lower_bounds <- c(a = 0.01, t0 = 0, z = 0, v = -Inf, st0 = 0, sz = 0, sv = 0)[par_types]
  upper_bounds <- c(a = Inf, t0 = Inf, z = 1, v = Inf, st0 = Inf, sz = 1, sv = Inf)[par_types]

  # ---- 4. Fit each participant -------------------------------------------------------------

  n_par <- length(param_names)
  paramEsts <- data.frame(matrix(NA, ncol = n_par + length(fixed_par_vals) + 5))
  names(paramEsts) <- c("ID", param_names, names(fixed_par_vals), "converge", "LLE", "AIC", "BIC")
  fit_info <- data.frame(ID = ids, optimizer_runs = NA_integer_,
                         converge = NA_integer_, best_replicated = NA, message = NA_character_,
                         se_status = NA_character_)
  seEsts <- data.frame(matrix(NA_real_, nrow = length(ids), ncol = n_par + 1))
  names(seEsts) <- c("ID", param_names)
  vcov_list <- list()

  settings_fit <- list(fixed = fixed, defaults = defaults, param_names = param_names,
                       par_types = par_types, start = start, lower = unname(lower_bounds),
                       upper = unname(upper_bounds), precision = precision,
                       stop_on_error = stop_on_error, p_outlier = p_outlier, w_outlier = w_outlier,
                       max_retries = max_retries, compute_se = compute_se,
                       print_starting_vals = print_starting_vals)

  cat("Fitting diffusion model for", length(ids), "participants...\n")
  pb <- txtProgressBar(min = 0, max = length(ids), style = 3)

  for (k in seq_along(ids)) {
    # (as a for loop over the IDs would give them: factor IDs as their labels)
    i <- if (is.factor(ids)) as.character(ids[k]) else ids[k]
    rows <- which(trials$ID == i)
    res <- fit_participant(trials[rows, , drop = FALSE], i,
                           lapply(applies, `[`, rows), settings_fit)
    fit <- res$fit

    seEsts[k, "ID"] <- i
    seEsts[k, res$fit_names] <- res$se
    vcov_list[as.character(i)] <- list(res$vcov)
    fit_info[k, ] <- list(i, res$n_runs, fit$convergence, res$replicated, fit$message, res$se_status)

    paramEsts[k, "ID"] <- i
    for (est in param_names) paramEsts[k, est] <- fit$par[est]
    for (fixed_par in names(fixed_par_vals)) paramEsts[k, fixed_par] <- fixed_par_vals[fixed_par]
    paramEsts[k, "converge"] <- fit$convergence
    paramEsts[k, "LLE"] <- -fit$objective
    paramEsts[k, "AIC"] <- 2 * length(fit$par) + 2 * fit$objective
    paramEsts[k, "BIC"] <- 2 * fit$objective + length(fit$par) * log(length(rows))

    setTxtProgressBar(pb, k)
  }
  close(pb)

  # ---- 5. Collect the results --------------------------------------------------------------

  n_retried <- sum(fit_info$optimizer_runs > 1)
  n_failed <- sum(fit_info$converge != 0)
  if (n_retried > 0) {
    message(n_retried, " participant(s) needed more than one optimizer run; ", n_failed,
            " did not converge after retries. See results$fit_info for details.")
  }
  names(fit_info)[names(fit_info) == "ID"] <- sub_id
  seEsts$ID <- fit_info[[sub_id]]
  names(seEsts)[names(seEsts) == "ID"] <- sub_id

  paramEsts <- paramEsts[, c("ID", col_names, "converge", "LLE", "AIC", "BIC")]
  names(paramEsts)[names(paramEsts) == "ID"] <- sub_id

  # Mean and SD across participants of each parameter (fit statistics such as the
  # log-likelihood are reported per participant, not averaged)
  summary <- data.frame(
    Parameter = col_names,
    Mean = round(colMeans(paramEsts[, col_names, drop = FALSE], na.rm = TRUE), digits = 2),
    SD = round(apply(paramEsts[, col_names, drop = FALSE], 2, sd, na.rm = TRUE), digits = 2),
    Fixed = !(col_names %in% param_names)
  )
  table <- summary %>%
    gt() %>%
    tab_spanner(label = "Model Summary", columns = c("Mean", "SD")) %>%
    fmt_number(columns = c("Mean", "SD"), decimals = 2)
  print(table)

  # The processed data (in seconds) with the original condition columns, used by
  # predict() and simulate()
  processed_data <- trials
  names(processed_data) <- c(sub_id, ".cell", response, rt)
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
    summary = summary,
    settings = list(
      varying_parameters = vary_list,
      fixed_parameters = as.list(fixed),
      n_participants = length(ids),
      n_conditions = length(unique(cond)),
      cell_parameters = cell_parameters,
      cells = cells_table,
      parameters = parameter_info,
      model = model_spec,
      condition_columns = conditions,
      columns = list(sub_id = sub_id, response = response, rt = rt),
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
  results
}




# ---- Internal helpers ---------------------------------------------------------------------

# Data-quality checks before fitting. Prints a warning for each problem found and returns
# TRUE if any were found.
#' @noRd
run_pre_checks <- function(data, sub_id, response, rt, vary_list, fix_st0, st0_fixed,
                           fix_sz, sz_fixed) {
  warnings_issued <- FALSE

  if (fix_st0 && st0_fixed > 0.5) {
    message("WARNING: st0_fixed seems large (", st0_fixed, "). This may cause estimation issues.")
    warnings_issued <- TRUE
  }
  if (fix_sz && sz_fixed > 0.5) {
    message("WARNING: sz_fixed seems large (", sz_fixed, "). This may cause estimation issues.")
    warnings_issued <- TRUE
  }

  # Trials per participant and condition cell
  check_conditions <- unique(unlist(vary_list))
  check_conditions <- check_conditions[!is.na(check_conditions) & check_conditions != ""]
  if (length(check_conditions) > 0) {
    data$temp_cond <- do.call(paste, c(data[check_conditions], sep = "_"))
    low_trial_combos <- data %>%
      group_by(!!sym(sub_id), .data$temp_cond) %>%
      summarise(trial_count = n(), .groups = "drop") %>%
      filter(.data$trial_count < 25)

    if (nrow(low_trial_combos) > 0) {
      message("WARNING: Some participant-condition combinations have very few trials (<25):")
      by_condition <- low_trial_combos %>%
        group_by(.data$temp_cond) %>%
        summarise(
          participants = paste(!!sym(sub_id), " (n=", .data$trial_count, ")", sep = "", collapse = ", "),
          .groups = "drop"
        )
      for (i in seq_len(nrow(by_condition))) {
        message("  Condition '", by_condition$temp_cond[i], "': ", by_condition$participants[i])
      }
      message("This may cause estimation problems.")
      warnings_issued <- TRUE
    }
  }

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

  warnings_issued
}

# Negative log-likelihood of one participant's trials. `base` holds, per parameter type,
# its value for every trial (a single value for fixed parameters); each estimated
# parameter in `pars` then replaces the values of the trials it applies to (`applies`).
# Trials in the same group share all parameter values, so ddiffusion is called once per
# group with scalar parameters. (Passing trial-level parameter vectors gives the same
# densities, but rtdists then has to find the unique parameter combinations on every
# call, which dominates the fitting time.)
#' @noRd
ddm_nll <- function(pars, RT, response, base, types, applies, groups, precision,
                    stop_on_error, p_outlier = 0, w_outlier = 0.1) {
  for (p in seq_along(types)) {
    base[[types[p]]][applies[[p]]] <- pars[[p]]
  }
  dens <- numeric(length(RT))
  for (g in groups) {
    k <- g[1]
    at <- function(type) {
      x <- base[[type]]
      if (length(x) == 1) x else x[k]
    }
    dens[g] <- ddiffusion(
      rt = RT[g], response = response[g],
      a = at("a"), v = at("v"), t0 = at("t0"), z = at("z") * at("a"),
      st0 = at("st0"), sz = at("sz"), sv = at("sv"),
      precision = precision, stop_on_error = stop_on_error
    )
  }

  # Outlier mixture as in HDDM: each trial comes from the diffusion model with probability
  # 1 - p_outlier, and otherwise from a uniform outlier distribution with density w_outlier
  if (p_outlier > 0) {
    dens[is.na(dens) | dens < 0] <- 0
    dens <- (1 - p_outlier) * dens + p_outlier * w_outlier
  }

  # A density of 0 usually means the true density is too small to represent (e.g. an RT
  # just above t0), so it is replaced by the smallest positive double. A larger floor (this
  # used to be 1e-10) caps the cost of a poorly fitted trial, and the optimizer can then
  # improve the fit to the remaining trials by effectively discarding one.
  if (any(is.na(dens) | dens <= 0)) {
    dens[dens <= 0 | is.na(dens)] <- .Machine$double.xmin
  }
  -sum(log(dens))
}

# Fit the model to one participant's trials (`d`, with columns cond, response and RT).
# `applies` gives, per estimated parameter, the trials of `d` it applies to; `s` holds the
# model settings built by driftd(). Returns the nlminb fit, the number of optimizer runs,
# whether the best solution was replicated, and the standard errors.
#' @noRd
fit_participant <- function(d, id, applies, s) {
  types <- s$par_types

  # ---- Starting values ----
  # EZ-diffusion (Wagenmakers et al., 2007) gives closed-form estimates of v, a and t0 per
  # condition cell. t0 is kept below the fastest RT (the 10% quantile with an outlier
  # component) so every trial has a non-zero likelihood at the start. z starts unbiased,
  # and the variability parameters at small values.
  cond_summaries <- d %>%
    dplyr::group_by(.data$cond) %>%
    dplyr::summarise(
      n_trials = n(),
      meanRT = mean(.data$RT, na.rm = TRUE),
      varRT  = stats::var(.data$RT, na.rm = TRUE),
      minRT  = t0_cap(.data$RT, s$p_outlier),
      q10RT  = quantile(.data$RT, probs = 0.1, na.rm = TRUE, names = FALSE),
      prop_upper = mean(.data$response == "upper", na.rm = TRUE),
      .groups = "drop"
    )
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
    ) %>%
    # fallbacks where EZ is undefined or implausible; a non-positive EZ t0 (e.g. strongly
    # skewed RTs) falls back to the fast-RT floor
    dplyr::mutate(
      v_guess = ifelse(is.na(.data$v_guess), 0.5 * sign(.data$prop_upper - 0.5 + 1e-6), .data$v_guess),
      a_guess = ifelse(is.na(.data$a_guess), 1.0, .data$a_guess),
      t0_guess = ifelse(is.na(.data$t0_guess) | .data$t0_guess <= 0,
                        pmin(.data$q10RT, 0.9 * .data$minRT), .data$t0_guess)
    )
  guesses <- dplyr::left_join(d, cond_summaries, by = "cond")

  # Per parameter type, the value of every trial before the estimated parameters are
  # filled in: the fixed value, or the starting guesses
  base <- lapply(names(s$fixed), function(type) {
    if (s$fixed[[type]]) s$defaults[[type]] else guesses[[paste0(type, "_guess")]]
  })
  names(base) <- names(s$fixed)

  # Each estimated parameter starts at the mean guess over the trials it applies to (the
  # minimum for t0, which must stay below every RT it applies to), or at its default if
  # there is no usable guess, and is kept inside the bounds. t0 is the earliest possible
  # response, so any trial faster than t0 has zero likelihood; its upper bound is therefore
  # just below the fastest RT it applies to (or the 10% quantile with an outlier component),
  # otherwise the optimizer can get stuck at a t0 that writes off the fastest trials.
  start <- s$start
  upper <- s$upper
  for (p in seq_along(s$param_names)) {
    trials <- applies[[p]]
    if (types[p] == "t0" && any(trials)) {
      upper[p] <- 0.999 * t0_cap(d$RT[trials], s$p_outlier)
    }
    g <- guesses[[paste0(types[p], "_guess")]][trials]
    g <- g[is.finite(g)]
    if (length(g) > 0) start[p] <- if (types[p] == "t0") min(g) else mean(g)
  }
  start <- pmin(pmax(start, s$lower), upper)

  # A parameter that applies to no trials of this participant (e.g. a condition the
  # participant has no data for) cannot be estimated. It is left out of this participant's
  # fit and reported as NA, so it neither stays at its starting value nor makes the Hessian
  # singular for the other parameters.
  estimable <- vapply(applies, any, logical(1))
  fit_names <- s$param_names[estimable]
  types <- types[estimable]
  applies <- applies[estimable]
  lower <- s$lower[estimable]
  upper <- upper[estimable]
  start <- start[estimable]

  if (s$print_starting_vals) {
    cat("\n======================================\n")
    cat(paste0("Participant ", id, " -- Summary:\n"))
    cat("======================================\n")
    print(cond_summaries)
    cat("\nStarting values passed to the optimizer:\n")
    print(round(start, 3))
    cat("======================================\n\n")
  }

  # Trials that share all parameter values (same condition cell and the same pattern of
  # applying parameters), so the likelihood is computed once per group
  groups <- unname(split(seq_len(nrow(d)), interaction(
    d$cond, do.call(paste0, lapply(applies, as.integer)), drop = TRUE)))

  nll <- function(pars, precision = s$precision) {
    ddm_nll(pars, RT = d$RT, response = d$response, base = base, types = types,
            applies = applies, groups = groups, precision = precision,
            stop_on_error = s$stop_on_error, p_outlier = s$p_outlier, w_outlier = s$w_outlier)
  }
  run_fit <- function(start_vals) {
    nlminb(start_vals, objective = nll,
           control = list(iter.max = 10000, eval.max = 10000), lower = lower, upper = upper)
  }

  # ---- Fit, with retries ----
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
  stuck_at_zero <- function(f) any(types %in% c("sv", "sz") & f$par <= 1e-3)
  fit <- run_fit(start)
  n_runs <- 1
  replicated <- NA
  if ((fit$convergence != 0 || stuck_at_zero(fit)) && s$max_retries > 0) {
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
    while (n_runs <= s$max_retries && !replicated) {
      if (!tried_nudge && stuck_at_zero(fit)) {
        next_start <- nudged_start(fit$par, types, lower, upper, s$defaults)
        tried_nudge <- TRUE
      } else {
        k <- k + 1
        next_start <- alternative_start(k, start, types, lower, upper, s$defaults)
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

  # ---- Standard errors ----
  # SEs are the square roots of the diagonal of the inverse Hessian of the negative
  # log-likelihood at the estimates (the observed Fisher information), computed numerically
  # with Richardson extrapolation (numDeriv). The Hessian is evaluated with at least
  # precision 4, because the less precise numerical integration at precision 3 noticeably
  # distorts the second derivatives when sv or sz are estimated.
  # Wald SEs require the estimate to lie inside its bounds. Parameters at (within 0.001 of)
  # a bound, e.g. sv = 0 or t0 at the fastest RT, get no SE, and the SEs of the other
  # parameters are computed with them held at their estimates.
  se <- setNames(rep(NA_real_, length(fit_names)), fit_names)
  vcov <- NULL
  se_status <- "not computed"
  if (s$compute_se) {
    at_bound <- (fit$par - lower < 1e-3) | (upper - fit$par < 1e-3)
    free <- which(!at_bound)
    bound_note <- if (any(at_bound)) {
      paste0("no SE at bound: ", paste(fit_names[at_bound], collapse = ", "))
    } else NULL
    se_status <- if (is.null(bound_note)) "ok" else bound_note
    if (length(free) > 0) {
      nll_free <- function(x) {
        p <- fit$par
        p[free] <- x
        nll(p, precision = max(s$precision, 4))
      }
      # numDeriv's largest step is d * |x| (d = 0.1 by default). A step past a bound, e.g.
      # t0 beyond the fastest RT, gives trials a likelihood of 0 and a meaningless Hessian,
      # so d is reduced to keep every step at most halfway to the nearest bound.
      x_free <- fit$par[free]
      gap <- pmin(x_free - lower[free], upper[free] - x_free)
      step_d <- min(0.1, 0.5 * min(gap / pmax(abs(x_free), 1e-4)))
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
          vcov <- solve(H)
          dimnames(vcov) <- list(fit_names[free], fit_names[free])
          se[free] <- sqrt(diag(vcov))
        }
      }
    }
  }

  list(fit = fit, fit_names = fit_names, n_runs = n_runs, replicated = replicated,
       se = se, vcov = vcov, se_status = se_status)
}

# Alternative starting vector number k (k >= 1), spread over the plausible parameter space
# with a Halton sequence (one prime base per parameter type, so all conditions of a
# parameter move together and the EZ pattern across conditions is kept). a and v are
# scaled by 0.5-2 times their EZ start; the others cover their usual range. Deterministic,
# so results are reproducible and the user's random seed is untouched.
#' @noRd
alternative_start <- function(k, start, types, lower, upper, defaults) {
  s <- start
  primes <- c(a = 2, t0 = 3, z = 5, v = 7, st0 = 11, sz = 13, sv = 17)
  for (p in seq_along(types)) {
    u <- halton(k, primes[[types[p]]])
    s[p] <- switch(types[p],
      a   = start[p] * 2^(2 * u - 1),
      v   = start[p] * 2^(2 * u - 1),
      t0  = upper[p] * (0.2 + 0.75 * u),
      z   = 0.3 + 0.4 * u,
      st0 = 0.3 * u,
      sv  = 2.5 * u,
      # sz must stay below 2 * min(z, a - z) on the absolute scale to be valid
      sz  = 0.9 * u * sz_limit(s, types, defaults))
  }
  pmin(pmax(s, lower), upper)
}

# Starting vector that moves sv and sz estimates stuck at 0 to clearly positive values.
# sv and sz spread drift and starting point symmetrically, so near 0 they change the
# likelihood only quadratically: its slope at exactly 0 is zero, and the optimizer can
# report convergence there even when a positive value fits better.
#' @noRd
nudged_start <- function(pars, types, lower, upper, defaults) {
  s <- pars
  s[types == "sv" & s <= 1e-3] <- 1
  if (any(types == "sz" & s <= 1e-3)) {
    s[types == "sz" & s <= 1e-3] <- 0.5 * sz_limit(s, types, defaults)
  }
  pmin(pmax(s, lower), upper)
}

# The largest valid sz, 2 * min(z, 1 - z) * a, at the mean a and z of a parameter vector
# (or their fixed values)
#' @noRd
sz_limit <- function(pars, types, defaults) {
  a_mean <- if (any(types == "a")) mean(pars[types == "a"]) else defaults[["a"]]
  z_mean <- if (any(types == "z")) mean(pars[types == "z"]) else defaults[["z"]]
  2 * min(z_mean, 1 - z_mean) * a_mean
}
