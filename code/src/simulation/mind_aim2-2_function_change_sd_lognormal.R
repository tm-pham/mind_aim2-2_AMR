# ============================================================================ #
# Project: MInD Aim 2.2
# RUNFILE
# Title: Change standard deviation of log-normal distribution
# Author: Thi Mui Pham, tmpham@hsph.harvard.edu
# ---------------------------------------------------------------------------- #
# This function changes the standard deviation of a log-normal distribution
# using the following steps: 
# 1. Estimates the parameters of the underlying normal distribution (mu, sigma)
# 2. Reduce the standard deviation by parameter given by user
# 3. Calculate the new sigma_log corresponding to the new log-normal sd
# 4. Adjust mu_log to maintain the same mean in the log-normal distribution
# 5. Generate the new log-normal distribution
# 6. Check the new distribution's mean and standard deviation
# ============================================================================ #
change.sd.lognormal <- function(mean_lognormal, 
                                sd_lognormal, 
                                k = 0.5 # Parameter indicating whether sd should be increased (>1) or decreased (<1)
                                ){
  # Step 1: Estimate the parameters of the underlying normal distribution
  mu_log <- log((mean_lognormal^2) / sqrt(sd_lognormal^2 + mean_lognormal^2))
  sigma_log <- sqrt(log(1 + (sd_lognormal^2 / mean_lognormal^2)))
  
  # Print the estimated parameters of the normal distribution
  cat("Estimated mu_log:", mu_log, "\n")
  cat("Estimated sigma_log:", sigma_log, "\n")
  
  # Step 2: Choose a smaller standard deviation for the log-normal distribution (e.g., half the current sd)
  sd_lognormal_new <- k*sd_lognormal
  
  # Step 3: Calculate the new sigma_log corresponding to the new log-normal sd
  sigma_log_new <- sqrt(log(1 + (sd_lognormal_new^2 / mean_lognormal^2)))
  
  # Step 4: Adjust mu_log to maintain the same mean in the log-normal distribution
  mu_log_new <- log(mean_lognormal) - (sigma_log_new^2) / 2
  
  # Step 5: Generate the new log-normal distribution
  new_data <- rlnorm(1000, meanlog = mu_log_new, sdlog = sigma_log_new)
  
  # Step 6: Check the new distribution's mean and standard deviation
  new_mean <- mean(new_data)
  new_sd <- sd(new_data)
  
  cat("Original log-normal mean:", mean_lognormal, "\n")
  cat("New log-normal mean:", new_mean, "\n")
  cat("Original log-normal sd:", sd_lognormal, "\n")
  cat("New log-normal sd:", new_sd, "\n")
  
  return(list(mean = new_mean, sd = new_sd))
}


