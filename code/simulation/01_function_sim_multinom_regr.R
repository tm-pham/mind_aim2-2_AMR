# ============================================================================ #
# Project: MInD Aim 2.2
# FUNCTION
# Title: Simulation study for multinomial logistic regression
# Author: Thi Mui Pham, tmpham@hsph.harvard.edu
# ---------------------------------------------------------------------------- #
# This script contains functions to
# 1. Generate synthetic data for multinomial logistic regression with random effects
# 2. Fit the hierarchical multinomial logistic regression using the mclogit package
# ============================================================================ #

# ------------------------------------------------------------------------------
# Function to generate synthetic data
# ------------------------------------------------------------------------------
# Written for simulating S. aureus-like data with three antibiotic classes
# The following AMR profiles are considered: 
# R-R-R, R-R-S, R-S-R, S-R-R, S-S-R, S-S-S
# Structure of the function
# 1. Determine the total number of observations
# 2. Generate time variables (year, day, month) at random
# 3. Generate antibiotic use distributions based on parameters provided 
# 4. Draw the random intercepts for facilities according to N(0, sd)
# 5. Generate the linear predictor for the multinomial outcome
# 6. Simulate the outcome based on the probabilities
generate.data <- function(params, # List of parameters 
                          seed = NULL,
                          abx_params = list(Class_1 = cbind(mean = 4.712055, sd = 0.353727), 
                                            Class_2 = cbind(mean = 4.712055, sd = 0.353727), 
                                            Class_3 = cbind(mean = 4.712055, sd = 0.353727)), # List of parameters for antibiotic use
                          prefix = "simData", 
                          # Covariates for facility-level antibiotic use 
                          ...
) {
  if(is.null(seed)){
    seed <- as.numeric(Sys.time())
  }
  params$seed <- seed
  set.seed(seed)
  # Save values in ... 
  args_list <- list(...)
  
  # Make values in params available in this environment
  list2env(params, envir = environment())
  
  # ----------------------------------------------------------------------------
  # Check whether parameters used in this functions were provided
  necessary_params <- c("sim_id", "time_period", "n_facilities", "n_obs_per_facility", 
                        "beta0", "beta", "gamma", "sigma", "rho", 
                        "random_effect_sd", "antibiograms", "reference")
  if(!all(necessary_params%in% names(params))){
    # Check which parameter is missing
    cat("These parameters are missing", setdiff(necessary_params, names(params)))
    stop("Please provide all necessary parameters.")
  }
  
  # Create results folder 
  output_folder <- here::here("results", "simulation", data_id)
  dir.create(output_folder, showWarnings = FALSE)
  
  sim_id_folder <- here::here("results", "simulation", data_id, sim_id)
  dir.create(sim_id_folder, showWarnings = FALSE)
  setwd(sim_id_folder)
  # ----------------------------------------------------------------------------
  # 1. Total number of observations 
  # = Number of facilites x Number of observations per facility
  n_isolates <- n_facilities * n_obs_per_facility
  
  # ----------------------------------------------------------------------------
  # 2. Generate sequence of dates within the time period years
  dates <- seq.Date(from = as.Date(paste0(time_period[1], "-01-01")),
                    to   = as.Date(paste0(time_period[2], "-12-31")),
                    by   = "day")
  
  # For each isolate, randomly assign a facility and a date
  sim_data <- data.frame(
    facility = sample(1:n_facilities, n_isolates, replace = TRUE),
    date = sample(dates, n_isolates, replace = TRUE)
  )
  
  # Extract day, month, year
  sim_data$day   <- as.numeric(format(sim_data$date, "%d"))
  sim_data$month <- as.numeric(format(sim_data$date, "%m"))
  sim_data$year  <- as.numeric(format(sim_data$date, "%Y"))
  
  # ----------------------------------------------------------------------------
  # Set census region, facility rurality and complexitylevel to 0 if not provided by user
  if(!"census_region" %in% names(params)){
    census_region <- 0
  }
  if(!"facility_rurality" %in% names(params)){
    facility_rurality <- 0
  }
  if(!"complexitylevel" %in% names(params)){
    complexitylevel <- 0
  }
  
  # ----------------------------------------------------------------------------
  # 3. Generate lognormal distribution for antibiotic use
  cat("Facility-level antibiotic use is generated according to a lognormal distribution based on provided or default parameters.\n")
  if(length(abx_params) == 0){
    stop("Please provide distributions for antibiotic use.")
  }
  abx_names <- names(abx_params)

  # For each unique facility-date combination, simulate antibiotic use for each class.
  facility_date <- sim_data %>% distinct(facility, date)
  
  # For each antibiotic class to simulate exposure from log-normal distribution
  # This assumes that for each antibiotic class, the same values will be mapped 
  # to the same facility-date combination. In particular, 
  for(class in abx_names){
    facility_date[[class]] <- rlnorm(nrow(facility_date),
                                     meanlog = abx_params[[class]][, "meanlog"],
                                     sdlog   = abx_params[[class]][, "sdlog"])
  }
  
  # Merge the facility-date exposures back into the isolate-level data
  sim_data <- left_join(sim_data, facility_date, by = c("facility", "date"))
  
  # ----------------------------------------------------------------------------
  # 4. Random intercepts for facilities
  # Normal distribution with mean = 0 and sd = random_effect_sd
  # For each facility and outcome category, simulate a random intercept.
  facility_random_effect <- lapply(antibiograms, function(cat) {
    rnorm(n_facilities, mean = 0, sd = random_effect_sd)
  })
  names(facility_random_effect) <- antibiograms
  
  # ----------------------------------------------------------------------------
  # Convert R-R-R to c(1, 1, 1), R-R-S to c(1, 1, 0), etc.
  # Split each string by the delimiter "-"
  abx_list <- strsplit(antibiograms, split = "-")
  # Convert each element: "R" becomes 1, "S" becomes 0
  ab_outcomes <- lapply(abx_list, function(x) as.numeric(x == "R"))
  names(ab_outcomes) <- antibiograms

  
  # ----------------------------------------------------------------------------
  # 5. Linear predictor for multinomial outcome
  lp <- list()
  constant <- beta[1]*sim_data$year + beta[2]*sim_data$day + beta[3]*sim_data$month + sum(sigma * c(census_region, facility_rurality, complexitylevel))
  for(i in 1:length(antibiograms)){
    ab <- antibiograms[i]
    # The last term multiplies the regression coefficient gamma with the 
    # antibiotic use of the respective antibiotic class for that facility-date
    # combination. This is done for each antibiotic class.
    lp[[ab]] <- beta0[i] + constant + rowSums(sweep(sim_data[, abx_names], 2, gamma * ab_outcomes[[ab]], "*"))
    # Save the computed linear predictor in a new column
    sim_data[[paste0("lp_", ab)]] <- lp[[ab]]
    # Compute exponentiated linear predictors for each non-reference outcome
    sim_data[[paste0("exp_lp_", ab)]] <- exp(sim_data[[paste0("lp_", ab)]])
  }
  
  # Compute the denominator for the multinomial logit model: 
  # 1 (for the reference) + sum of exp(lp) for non-reference outcomes
  exp_lp_matrix <- sapply(antibiograms, function(ab) sim_data[[paste0("exp_lp_", ab)]])
  denom <- 1 + rowSums(exp_lp_matrix)
  
  # Compute probabilities for non-reference outcomes and the reference
  probs_non_ref <- exp_lp_matrix / denom
  prob_ref <- 1 / denom
  
  # Build a matrix of probabilities with column names corresponding to outcomes
  # Order: non-reference outcomes (in the same order as in non_ref) then the reference
  probs_matrix <- cbind(probs_non_ref, prob_ref)
  colnames(probs_matrix) <- c(antibiograms, reference)
  
  # ----------------------------------------------------------------------------
  # 6. Simulate antibiogram membership (outcome) using the probabilities 
  # generated in the previous step. 
  # Apply the sampling function over each row of the probability matrix
  sim_data$outcome <- apply(probs_matrix, 1, function(p) sample(names(p), size = 1, prob = p))
  outcome <- sim_data$outcome <- factor(sim_data$outcome, levels = c(reference, antibiograms))
  cat("Outcome table:\n")
  print(table(outcome))

  # ----------------------------------------------------------------------------
  # Write output 
  # ----------------------------------------------------------------------------
  # Save parameters in text file in the folder 
  # Save time in parameter
  params$time_stamp <- Sys.time()
  list_to_df <- do.call(rbind, lapply(names(params), function(name) {
    data.frame(Variable = name, Value = I(list(params[[name]])), stringsAsFactors = FALSE)
  }))
  list_to_df[list_to_df$Variable=="time_stamp", "Value"][[1]] <- as.character(format(as.POSIXct(list_to_df[list_to_df$Variable=="time_stamp", "Value"][[1]], origin = "1970-01-01"), "%Y-%m-%d %H:%M:%S %Z"))
  
  write.table(list_to_df, file = paste0(prefix, sim_id, "_params.txt"), sep = " = ", 
              row.names = F, col.names = F, quote = F)
  
  # ----------------------------------------------------------------------------
  # Save antibiogram outcome summary 
  write.table(table(outcome), file = paste0(prefix, sim_id, "_outcome_table.txt"), 
              row.names = F, quote = F)
  
  # ----------------------------------------------------------------------------
  # Save simulation data
  simReg <- list(data = sim_data, params = params, abx_params = abx_params)
  save(simReg, args_list, file = paste0(prefix, sim_id, ".RData"))
  
  return(sim_data)
}

# ------------------------------------------------------------------------------
# Function to fit the hierarchical multinomial logistic regression using the 
# mclogit packge (frequentist approach)
# ------------------------------------------------------------------------------
mblogit.fit.model <- function(data, 
                              formula = outcome ~ 1+ year + day + month + fluoroquinolones + lincosamides + macrolides + beta_lactams, 
                              random_effects = ~1|facility) {
  # Load the necessary package
  if (!requireNamespace("mclogit", quietly = TRUE)) {
    install.packages("mclogit")
  }
  
  model <- mclogit::mblogit(formula, 
                            random = random_effects,
                            data = data)
  
  return(model)
}
