# ============================================================================ #
# Project: MInD Aim 2.2
# FUNCTION
# Title: Simulation study for multinomial logistic regression
# Author: Thi Mui Pham, tmpham@hsph.harvard.edu
# ---------------------------------------------------------------------------- #
# Written for simulating S. aureus data with three antibiotic classes
# The following antibiograms are considered: 
# R-R-R
# R-R-S
# R-S-R
# S-R-R
# S-S-R
# S-S-S
# ============================================================================ #
# Function to generate synthetic data
generate.data <- function(params, # List of parameters 
                          seed = NULL,
                          abx_distr = list(), # List of distributions for antibiotic use
                          # Covariates for facility-level antibiotic use 
                          ...
) {
  if(is.null(seed)){
    seed <- as.numeric(Sys.time())
  }
  set.seed(seed)
  # Save values in ... 
  args_list <- list(...)
  
  # Make values in params available in this environment
  list2env(params, envir = environment())
  facilities <- rep(1:n_facilities, each = n_obs_per_facility)
  n <- length(facilities)
  
  ##############################################################################
  # Generate calendar time variables
  year <- sample(time_period, n, replace = TRUE)
  day <- sample(1:31, n, replace = TRUE)
  month <- sample(1:12, n, replace = TRUE)
  
  ##############################################################################
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
  
  ##############################################################################
  # Random intercepts for facilities
  random_intercepts <- rnorm(n_facilities, mean = 0, sd = random_effect_sd)
  facility_effect <- random_intercepts[facilities]
  
  # Convert R-R-R to c(1, 1, 1), R-R-S to c(1, 1, 0), etc.
  ab_outcomes <- lapply(strsplit(antibiograms, "-"), function(x) {
    x <- x[[1]]
    c(1*(x == "R"), 1*(x == "S"), 1*(x == "S"))
  })
  names(ab_outcomes) <- antibiograms
  
  ##############################################################################
  # Linear predictor for multinomial outcome
  lp <- list()
  constant <- beta[2]*year + beta[3]*day + beta[4]*month + facility_effect + sum(sigma * c(census_region, facility_rurality, complexitylevel))
  
  temp <- NULL
  for(i in 1:length(antibiograms)){
    lp[[i]] <- 0
    temp[[i]] <- beta0[i]
    for(j in 1:length(abx_distr)){
      temp[[i]] <- temp[[i]] + gamma[i] * abx_distr[[j]] * ab_outcomes[[i]][j]
    }
    lp[[i]] <- lp[[i]] + temp[[i]] + constant
  }
  
  denominator <- 1 + Reduce("+", lapply(lp, function(x) exp(x)))
  P <- cbind(1/denominator, do.call(cbind, lapply(lapply(lp, exp), function(x) x/denominator)))
  
  # Simulate antibiogram membership using the probabilities 
  outcome <- apply(P, MARGIN = 1, function(x) sample(x = c(reference, antibiograms), size = 1, prob = x))
  cat("Outcome table:\n")
  print(table(outcome))
  
  sim_data <- data.frame(sim_id = sim_id, facility = facilities, year = year, day = day, month = month,
                        outcome = outcome, do.call("cbind", abx_distr))
  sim_data$outcome <- factor(sim_data$outcome, levels = c(reference, antibiograms))
  
  params$time_stamp <- Sys.time()
  
  # Create folder for the simulation dataset
  folder_path <- paste0("data_output/simReg/", org_short, "/", sim_id)
  # Check if the folder exists, if not, create it
  if (!dir.exists(folder_path)) {
    dir.create(folder_path, recursive = TRUE)
    message("Folder created successfully!")
  } else {
    message(paste0("Folder ", folder_path, " already exists."))
  }
  
  # Save parameters in text file in the folder 
  list_to_df <- do.call(rbind, lapply(names(params), function(name) {
    data.frame(Variable = name, Value = I(list(params[[name]])), stringsAsFactors = FALSE)
  }))
  list_to_df[list_to_df$Variable=="time_stamp", "Value"][[1]] <- as.character(format(as.POSIXct(list_to_df[list_to_df$Variable=="time_stamp", "Value"][[1]], origin = "1970-01-01"), "%Y-%m-%d %H:%M:%S %Z"))
  
  write.table(list_to_df, file = paste0(folder_path, "/mind_aim2-2_simReg_", sim_id, "_params.txt"), sep = " = ", 
              row.names = F, col.names = F, quote = F)
  
  # Save some summary statistics 
  write.table(table(outcome), file = paste0(folder_path, "/mind_aim2-2_simReg_", sim_id, "_outcome_table.txt"), 
              row.names = F, quote = F)
  
  # Save simulation data
  simReg <- list(data = sim_data, params = params)
  save(simReg, args_list, file = paste0(folder_path, "/mind_aim2-2_simReg_", sim_id, ".RData"))
  
  return(sim_data)
}

# ---------------------------------------------------------------------------- #
# Function to fit the hierarchical multinomial logistic regression using the 
# mclogit packge (frequentist approach)
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

# ---------------------------------------------------------------------------- #
# Function to fit the hierarchical multinomial logistic regression using the 
# brms package (Bayesian approach)
brms.fit.model <- function(data, 
                           formula, 
                           family = "categorical", 
                           priors, 
                           n_iter = 1000, 
                           n_cores = 3, 
                           n_chains = 3, 
                           seed = 123, 
                           maxTreedepth = 13, 
                           adaptDelta = 0.99, 
                           folder_path_prefix){
  
  sim_id <- unique(data$sim_id)
  
  brms_fit <- brms::brm(formula = formula,
                        data = data, 
                        family = "categorical", 
                        iter = n_iter, 
                        warmup = floor(n_iter/2), 
                        cores = n_cores, 
                        chains = n_chains, 
                        prior = priors, 
                        seed = seed, 
                        control = list(max_treedepth = maxTreedepth, 
                                       adapt_delta = adaptDelta))
  
  # Create folder for the simulation dataset
  folder_path <- paste0(folder_path_prefix, sim_id)
  # Check if the folder exists, if not, create it
  if (!dir.exists(folder_path)) {
    dir.create(folder_path, recursive = TRUE)
    message("Folder created successfully!")
  } else {
    message(paste0("Folder ", folder_path, " already exists."))
  }
  
  params <- list(formula = formula, family = family, n_iter = n_iter, n_cores = n_cores, n_chains = n_chains)
  
  # Save parameters in text file in the folder 
  list_to_df <- do.call(rbind, lapply(names(params), function(name) {
    data.frame(Variable = name, Value = I(list(params[[name]])), stringsAsFactors = FALSE)
  }))
  
  write.table(list_to_df, file = paste0(folder_path, "/mind_aim2-2_simReg_", sim_id, "_brms_params.txt"), sep = " = ", 
              row.names = F, col.names = F, quote = F)
  
  return(brms_fit)
}

# Function to perform simulation study
simulation.study <- function(n_simulations, n_facilities, n_obs_per_facility, beta, gamma, random_effect_sd) {
  results <- replicate(n_simulations, {
    data <- generate_data(n_facilities, n_obs_per_facility, beta, gamma, random_effect_sd)
    model <- fit_model(data)
    coef(summary(model))
  }, simplify = FALSE)
  
  return(results)
}

