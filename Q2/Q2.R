# Q2_0

set.seed(123)

# Function to price Asian Call Option using Monte Carlo
price_asian_call_MC <- function(S0, K, T, r, q, sigma, m, N) {
  start_time <- proc.time()
  delta_t <- T / m
  sqrt_delta_t <- sqrt(delta_t)
  Z <- matrix(rnorm(N * m), nrow = N, ncol = m) # Generate N x m standard normal random variables
  
  drift <- (r - q - 0.5 * sigma^2) * delta_t
  diffusion <- sigma * sqrt_delta_t * Z
  W <- t(apply(diffusion, 1, cumsum)) # Calculate cumulative sum for Brownian motion
  
  log_S <- log(S0) + outer(rep(1, N), (r - q - 0.5 * sigma^2) * (delta_t * (1:m))) + W
  S <- exp(log_S)
  S_bar <- rowMeans(S)  # Calculate arithmetic average S_bar_T
  
  payoffs <- pmax(S_bar - K, 0)
  discounted_payoffs <- exp(-r * T) * payoffs
  
  option_price <- mean(discounted_payoffs)

  std_error <- sd(discounted_payoffs) / sqrt(N)
  
  # Calculate 95% confidence interval
  CI_lower <- option_price - 1.96 * std_error
  CI_upper <- option_price + 1.96 * std_error
  
  end_time <- proc.time()
  comp_time <- (end_time - start_time)[["elapsed"]]

  return(list(
    N = N,
    Option_Price = round(option_price, 2),
    Standard_Error = round(std_error, 5),
    Confidence_Interval = c(round(CI_lower, 2), round(CI_upper, 2)),
    Computation_Time_sec = round(comp_time, 4)
  ))
}

S0 <- 100      # Initial asset price
K <- 100       # Strike price
T <- 1         # Time to maturity (in years)
r <- 0.10      # Risk-free rate
q <- 0         # Dividend yield
sigma <- 0.20  # Volatility
m <- 50        # Number of monitoring points

sample_sizes <- c(1000, 4000, 16000, 64000, 256000)

# Initialize a data frame to store results
results <- data.frame(
  Sample_Size = numeric(),
  Option_Price = numeric(),
  Standard_Error = numeric(),
  CI_Lower = numeric(),
  CI_Upper = numeric(),
  Computation_Time_sec = numeric(),
  stringsAsFactors = FALSE
)

for (N in sample_sizes) {
  cat("Processing N =", N, "...\n")
  res <- price_asian_call_MC(S0, K, T, r, q, sigma, m, N)

  results <- rbind(results, data.frame(
    Sample_Size = res$N,
    Option_Price = res$Option_Price,
    Standard_Error = res$Standard_Error,
    CI_Lower = res$Confidence_Interval[1],
    CI_Upper = res$Confidence_Interval[2],
    Computation_Time_sec = res$Computation_Time_sec
  ))
}

print(results)

# result
#  Sample_Size Option_Price Standard_Error CI_Lower CI_Upper Computation_Time_sec
# 1        1000         7.02        0.26869     6.49     7.55                 0.00
# 2        4000         7.39        0.14233     7.11     7.67                 0.03
# 3       16000         7.10        0.06805     6.96     7.23                 0.14
# 4       64000         7.16        0.03433     7.09     7.23                 0.66
# 5      128000         7.16        0.02433     7.11     7.20                 1.55



#Q2_2
set.seed(123)

# function to calculate the analytical price of Geometric Asian Call Option
price_geometric_asian_call <- function(S0, K, T, r, q, sigma, m){
  delta_t <- T/m
  sigma_sq <- sigma^2 
  sigma_z_sq <- (sigma_sq) * (m + 1) * (2 * m + 1) / (6 * m^2) # the effective volatility of the log of the geometric average
  drift <-(r - q - 0.5*sigma_sq) * (m + 1) / (2 * m) + 0.5 * sigma_z_sq
  sigma_z <- sqrt(sigma_z_sq)
  d1 <- (log(S0/K) + (drift + 0.5 * sigma_z_sq) * T) / (sigma_z * sqrt(T))
  d2 <- (log(S0/K) + (drift - 0.5 * sigma_z_sq) * T) / (sigma_z * sqrt(T))
  geo_asian_price <- exp(-r*T) * (S0* exp(drift * T) * pnorm(d1) - K * pnorm(d2))
  return(geo_asian_price)
}

# function to price asian call option using monte carlo with controla variate
price_asian_call_MC_control_variate <- function(S0, K, T, r, q, sigma, m, N) {
  start_time <- proc.time()

  delta_t <- T / m
  sqrt_delta_t <- sqrt(delta_t)
  
  # Generate N x m standard normal random variables
  Z <- matrix(rnorm(N * m), nrow = N, ncol = m)
  
  drift <- (r - q - 0.5 * sigma^2) * delta_t
  diffusion <- sigma * sqrt_delta_t * Z
  
  # Calculate cumulative sum for Brownian motion, rows-wisely
  W <- t(apply(diffusion, 1, cumsum))

  log_S <- log(S0) + outer(rep(1, N), (r - q - 0.5 * sigma^2) * (delta_t * (1:m))) + W 
  # 2nd term: a matrix of size N*m, where each row contains the drift term values for one path over all m time steps
  S <- exp(log_S)
  
  # Calculate arithmetic average S_bar_T
  S_bar_arith <- rowMeans(S)
  
  # Calculate geometric average S_bar_G
  S_bar_geo <- exp(rowMeans(log(S)))
  
  # Calculate payoffs for arithmetic Asian call
  payoffs_arith <- pmax(S_bar_arith - K, 0)
  
  # Calculate payoffs for geometric Asian call
  payoffs_geo <- pmax(S_bar_geo - K, 0)
  
  # Discount payoffs to present value
  discounted_payoffs_arith <- exp(-r * T) * payoffs_arith
  discounted_payoffs_geo <- exp(-r * T) * payoffs_geo
  
  # Calculate analytical price of Geometric Asian Call
  geo_asian_price <- price_geometric_asian_call(S0, K, T, r, q, sigma, m)
  
  # Calculate covariance and variance for Control Variate
  cov_xy <- cov(discounted_payoffs_arith, discounted_payoffs_geo)
  var_y <- var(discounted_payoffs_geo)
  
  # Optimal coefficient
  theta <- cov_xy / var_y
  
  # Calculate Control Variate estimator
  control_variate_estimator <- discounted_payoffs_arith - theta * (discounted_payoffs_geo - geo_asian_price)
  
  # Calculate option price using Control Variate
  option_price_cv <- mean(control_variate_estimator)

  std_error_cv <- sd(control_variate_estimator) / sqrt(N)
  
  CI_lower_cv <- option_price_cv - 1.96 * std_error_cv
  CI_upper_cv <- option_price_cv + 1.96 * std_error_cv
  
  end_time <- proc.time()
  comp_time <- (end_time - start_time)[["elapsed"]]
  
  return(list(
    N = N,
    Option_Price_CV = round(option_price_cv, 2),
    Standard_Error_CV = round(std_error_cv, 5),
    Confidence_Interval_CV = c(round(CI_lower_cv, 2), round(CI_upper_cv, 2)),
    Computation_Time_sec = round(comp_time, 4)
  ))
}


# Function to perform standard Monte Carlo pricing
price_asian_call_MC_standard <- function(S0, K, T, r, q, sigma, m, N) {
  start_time <- proc.time()

  delta_t <- T / m
  sqrt_delta_t <- sqrt(delta_t)
  
  # Generate N x m standard normal random variables
  Z <- matrix(rnorm(N * m), nrow = N, ncol = m)
  
  drift <- (r - q - 0.5 * sigma^2) * delta_t
  diffusion <- sigma * sqrt_delta_t * Z
  
  # Calculate cumulative sum for Brownian motion
  W <- t(apply(diffusion, 1, cumsum))
  
  # Calculate log S_ti
  log_S <- log(S0) + outer(rep(1, N), (r - q - 0.5 * sigma^2) * (delta_t * (1:m))) + W
  
  S <- exp(log_S)
  
  # Calculate arithmetic average S_bar_T
  S_bar <- rowMeans(S)
  
  # Calculate payoffs
  payoffs <- pmax(S_bar - K, 0)
  
  discounted_payoffs <- exp(-r * T) * payoffs
  
  option_price <- mean(discounted_payoffs)
  
  std_error <- sd(discounted_payoffs) / sqrt(N)
  
  CI_lower <- option_price - 1.96 * std_error
  CI_upper <- option_price + 1.96 * std_error
  
  end_time <- proc.time()
  comp_time <- (end_time - start_time)[["elapsed"]]
  
  return(list(
    N = N,
    Option_Price = round(option_price, 2),
    Standard_Error = round(std_error, 5),
    Confidence_Interval = c(round(CI_lower, 2), round(CI_upper, 2)),
    Computation_Time_sec = round(comp_time, 4)
  ))
}


S0 <- 100      # Initial asset price
K <- 100       # Strike price
T <- 1         # Time to maturity (in years)
r <- 0.10      # Risk-free rate
q <- 0         # Dividend yield
sigma <- 0.20  # Volatility
m <- 50        # Number of monitoring points

sample_sizes <- c(1000, 4000, 16000, 64000, 256000)

results_standard <- data.frame(
  Sample_Size = numeric(),
  Option_Price = numeric(),
  Standard_Error = numeric(),
  CI_Lower = numeric(),
  CI_Upper = numeric(),
  Computation_Time_sec = numeric(),
  stringsAsFactors = FALSE
)

results_cv <- data.frame(
  Sample_Size = numeric(),
  Option_Price_CV = numeric(),
  Standard_Error_CV = numeric(),
  CI_Lower_CV = numeric(),
  CI_Upper_CV = numeric(),
  Computation_Time_sec = numeric(),
  stringsAsFactors = FALSE
)


for (N in sample_sizes) {
  cat("Processing N =", N, "...\n")
  
  # Standard Monte Carlo
  res_std <- price_asian_call_MC_standard(S0, K, T, r, q, sigma, m, N)
  results_standard <- rbind(results_standard, data.frame(
    Sample_Size = res_std$N,
    Option_Price = res_std$Option_Price,
    Standard_Error = res_std$Standard_Error,
    CI_Lower = res_std$Confidence_Interval[1],
    CI_Upper = res_std$Confidence_Interval[2],
    Computation_Time_sec = res_std$Computation_Time_sec
  ))
  
  # Monte Carlo with Control Variate
  res_cv <- price_asian_call_MC_control_variate(S0, K, T, r, q, sigma, m, N)
  results_cv <- rbind(results_cv, data.frame(
    Sample_Size = res_cv$N,
    Option_Price_CV = res_cv$Option_Price_CV,
    Standard_Error_CV = res_cv$Standard_Error_CV,
    CI_Lower_CV = res_cv$Confidence_Interval_CV[1],
    CI_Upper_CV = res_cv$Confidence_Interval_CV[2],
    Computation_Time_sec = res_cv$Computation_Time_sec
  ))
}

# Combine and display results
final_results <- merge(results_standard, results_cv, by = "Sample_Size")
print(final_results)

# normal monte carlo method
#  Sample_Size Option_Price Standard_Error CI_Lower CI_Upper Computation_Time_sec.x 
# 1        1000         7.55        0.28892     6.98     8.11                   0.01            
# 2        4000         7.12        0.13784     6.85     7.39                   0.06            
# 3       16000         7.04        0.06748     6.91     7.18                   0.35            
# 4       64000         7.15        0.03440     7.08     7.22                   1.12            
# 5      128000         7.17        0.02433     7.13     7.22                   1.55         
# with control variate
#  Sample_Size  Option_Price_CV  Standard_Error_CV CI_Lower_CV CI_Upper_CV Computation_Time_sec.y
# 1        1000         7.17        0.00877        7.15        7.19                   0.01
# 2        4000         7.17        0.00414        7.16        7.18                   0.08
# 3       16000         7.17        0.00200        7.16        7.17                   0.34
# 4       64000         7.17        0.00101        7.16        7.17                   0.88
# 5      128000         7.17        0.00072        7.16        7.17                   1.69
# Save the combined results into CSV files
# write.csv(final_results, "E:\\IE 522 stats\\statsFinalQ3\\statsFinalQ2\\option_pricing_results.csv", row.names = FALSE)
write.csv(results_standard, "E:\\IE 522 stats\\statsFinalQ3\\statsFinalQ2\\standard_monte_carlo_results.csv", row.names = FALSE)
write.csv(results_cv, "E:\\IE 522 stats\\statsFinalQ3\\statsFinalQ2\\control_variate_results.csv", row.names = FALSE)

# Prepare data for plotting
sample_sizes <- final_results$Sample_Size

# Standard Errors
standard_error_std <- final_results$Standard_Error
standard_error_cv <- final_results$Standard_Error_CV

# Computation Times
comp_time_std <- final_results$Computation_Time_sec.x
comp_time_cv <- final_results$Computation_Time_sec.y

# Plot Standard Error vs. Sample Size
par(mfrow = c(1, 2))  # Set up plotting area for two plots side by side

# Plot Standard Error
plot(sample_sizes, standard_error_std, type = "b", col = "blue", lwd = 2, log = "x",
     xlab = "Sample Size (log scale)", ylab = "Standard Error",
     main = "Standard Error vs. Sample Size",
     ylim = range(c(standard_error_std, standard_error_cv)))

lines(sample_sizes, standard_error_cv, type = "b", col = "red", lwd = 2)

legend("topright",
       legend = c("Standard Monte Carlo", "Control Variate Monte Carlo"),
       col = c("blue", "red"),
       lty = 1, lwd = 2)

# Plot Computation Time vs. Sample Size
plot(sample_sizes, comp_time_std, type = "b", col = "blue", lwd = 2, log = "x",
     xlab = "Sample Size (log scale)", ylab = "Computation Time (seconds)",
     main = "Computation Time vs. Sample Size",
     ylim = range(c(comp_time_std, comp_time_cv)))

lines(sample_sizes, comp_time_cv, type = "b", col = "red", lwd = 2)

legend("topleft",
       legend = c("Standard Monte Carlo", "Control Variate Monte Carlo"),
       col = c("blue", "red"),
       lty = 1, lwd = 2)

# Reset plotting area
par(mfrow = c(1, 1))


















# 
# # Further study: stratified MC with control variates
# price_asian_call_MC_control_variate_stratified<-function(S0, K, T, r, q, sigma, m, N, strata = 10){
#   start_time <- proc.time()
#   if (N %% strata != 0){
#     stop('Sample size N must be divisible by the number of strata') # ensure N is divisible by strata
#   }
#   N_stratum <- N/strata # number of samples in each strata
#   delta_t <- T/m
#   sqrt_delta_t <-sqrt(delta_t)
#   
#   # initialize vectors to store payoffs
#   discounted_payoffs_arith<- numeric(N)
#   discounted_payoffs_geo<- numeric(N)
#   
#   # Analytical price of Geometric Asian Call
#   geo_asian_price <- price_geometric_asian_call(S0, K, T, r, q, sigma, m)
#   
#   # Stratified Sampling over the first time step
#   # For simplicity, we stratify based on the first random variable Z1
#   # Divide the standard normal distribution into 'strata' equal probability intervals
#   quantiles <- qnorm(seq(0, 1, length.out = strata + 1))
#   
#   for (i in 1:strata){
#     # the interval for the current stratum
#     lower <- quantiles[i]
#     upper <- quantiles[i+1]
#     
#     # Sample uniformly within the interval and then apply inverse CDF to get stratified Z1
#     U <- runif(N_stratum, 0, 1)
#     Z1 <- qnorm(pnorm(lower) + U*(pnorm(upper) - pnorm(lower)))
#     
#     # generate the reamining m-1 standard normal variables for each path
#     Z_rest <- matrix(rnorm(N_stratum * (m-1)), nrow = N_stratum, ncol = m-1)
#     Z <- cbind(Z1, Z_rest)
#     
#     drift <- (r - q - 0.5 * sigma^2) * delta_t
#     diffusion <- sigma * sqrt_delta_t * Z
#     
#     W <- t(apply(diffusion, 1, cumsum))
#     log_S <- log(S0) + outer(rep(1, N_stratum), (r -q - 0.5 * sigma^2)*(delta_t*(1:m))) + W
#     
#     S <- exp(log_S)
#     
#     # Calculate arithmetic and geometric averages
#     S_bar_arith <- rowMeans(S)
#     S_bar_geo <- exp(rowMeans(log(S)))
#     
#     # Calculate payoffs
#     payoffs_arith <- pmax(S_bar_arith - K, 0)
#     payoffs_geo <- pmax(S_bar_geo - K, 0)
#     
#     discounted_payoffs_arith[((i - 1) * N_stratum + 1):(i * N_stratum)] <- exp(-r * T) * payoffs_arith
#     discounted_payoffs_geo[((i - 1) * N_stratum + 1):(i * N_stratum)] <- exp(-r * T) * payoffs_geo
#   }
#   
#   cov_xy <- cov(discounted_payoffs_arith, discounted_payoffs_geo) 
#   var_y <- var(discounted_payoffs_geo)
#   
#   theta <- cov_xy/ var_y
#   
#   # Calculate Control Variate estimator
#   control_variate_estimator <- discounted_payoffs_arith - theta * (discounted_payoffs_geo - geo_asian_price)
#   
#   # Calculate option price using Control Variate
#   option_price_cv_strat <- mean(control_variate_estimator)
#   
#   std_error_cv_strat <- sd(control_variate_estimator) / sqrt(N)
#   
#   CI_lower_cv_strat <- option_price_cv_strat - 1.96 * std_error_cv_strat
#   CI_upper_cv_strat <- option_price_cv_strat + 1.96 * std_error_cv_strat
#   
#   end_time <- proc.time()
#   comp_time <- (end_time - start_time)[["elapsed"]]
#   
#   return(list(
#     Method = "Control Variate + Stratified Sampling",
#     N = N,
#     Option_Price = round(option_price_cv_strat, 2),
#     Standard_Error = round(std_error_cv_strat, 5),
#     Confidence_Interval = c(round(CI_lower_cv_strat, 2), round(CI_upper_cv_strat, 2)),
#     Computation_Time_sec = round(comp_time, 4)
#   ))
# }
# 
# 
# # Function to price Asian Call Option using Monte Carlo with Control Variate and Brownian Bridge
# price_asian_call_MC_control_variate_brownian_bridge <- function(S0, K, T, r, q, sigma, m, N) {
#   start_time <- proc.time()
#   
#   delta_t <- T / m
#   sqrt_delta_t <- sqrt(delta_t)
#   
#   # Generate N x m standard normal random variables
#   Z <- matrix(rnorm(N * m), nrow = N, ncol = m)
#   
#   # Apply Brownian Bridge
#   # Adjust the sampling to better estimate the Brownian paths: reorder the Z to simulate the Brownian bridge
#   # For simplicity, we will implement a simple Brownian bridge technique where we first simulate the endpoints and then fill in the midpoints.
#   
#   # Initialize Brownian paths matrix
#   W_T <- rnorm(N) * sqrt(T)
#   W <- matrix(0, nrow = N, ncol = m)
#   
#   # Set the first and last points of the Brownian paths
#   W[,1] <- 0         # W_0 = 0
#   W[,m] <- W_T       # W_T set explicitly
#   
#   # Simulate intermediate points using the Brownian Bridge
#   for (j in 2:(m-1)) {
#     t_j <- (j-1) * delta_t       # Time at the j-th step
#     W[,j] <- (t_j / T) * W_T + rnorm(N) * sqrt((T - t_j) * t_j / T)
#   }
#   
#   # Calculate log S_ti using Brownian motion paths
#   log_S <- log(S0) + outer(rep(1, N), (r - q - 0.5 * sigma^2) * (delta_t * (1:m))) + sigma * sqrt_delta_t * W
#   
#   # Calculate S_ti
#   S <- exp(log_S)
#   
#   # Calculate arithmetic average S_bar_T
#   S_bar_arith <- rowMeans(S)
#   
#   # Calculate geometric average S_bar_G
#   S_bar_geo <- exp(rowMeans(log(S)))
#   
#   # Calculate payoffs for arithmetic Asian call
#   payoffs_arith <- pmax(S_bar_arith - K, 0)
#   
#   # Calculate payoffs for geometric Asian call
#   payoffs_geo <- pmax(S_bar_geo - K, 0)
#   
#   # Discount payoffs to present value
#   discounted_payoffs_arith <- exp(-r * T) * payoffs_arith
#   discounted_payoffs_geo <- exp(-r * T) * payoffs_geo
#   
#   # Calculate analytical price of Geometric Asian Call
#   geo_asian_price <- price_geometric_asian_call(S0, K, T, r, q, sigma, m)
#   
#   # Compute covariance and variance
#   cov_xy <- cov(discounted_payoffs_arith, discounted_payoffs_geo)  # Covariance of X and Y
#   var_y <- var(discounted_payoffs_geo)                             # Variance of Y
#   
#   # Optimal coefficient
#   theta <- cov_xy / var_y
#   
#   # Calculate Control Variate estimator
#   control_variate_estimator <- discounted_payoffs_arith - theta * (discounted_payoffs_geo - geo_asian_price)
#   
#   # Calculate option price using Control Variate
#   option_price_cv_bb <- mean(control_variate_estimator)
#   
#   # Calculate standard error
#   std_error_cv_bb <- sd(control_variate_estimator) / sqrt(N)
#   
#   # Calculate 95% confidence interval
#   CI_lower_cv_bb <- option_price_cv_bb - 1.96 * std_error_cv_bb
#   CI_upper_cv_bb <- option_price_cv_bb + 1.96 * std_error_cv_bb
#   
#   # End time
#   end_time <- proc.time()
#   comp_time <- (end_time - start_time)[["elapsed"]]
# 
#   return(list(
#     Method = "Control Variate + Brownian Bridge",
#     N = N,
#     Option_Price = round(option_price_cv_bb, 2),
#     Standard_Error = round(std_error_cv_bb, 5),
#     Confidence_Interval = c(round(CI_lower_cv_bb, 2), round(CI_upper_cv_bb, 2)),
#     Computation_Time_sec = round(comp_time, 4)
#   ))
# }
# 
# # Main Parameters
# S0 <- 100      # Initial asset price
# K <- 100       # Strike price
# T <- 1         # Time to maturity (in years)
# r <- 0.10      # Risk-free rate
# q <- 0         # Dividend yield
# sigma <- 0.20  # Volatility
# m <- 50        # Number of monitoring points
# 
# # Sample sizes to investigate convergence
# sample_sizes <- c(1000, 4000, 16000, 64000, 128000)
# 
# # Initialize data frames to store results
# results_standard <- data.frame(
#   Method = character(),
#   Sample_Size = numeric(),
#   Option_Price = numeric(),
#   Standard_Error = numeric(),
#   CI_Lower = numeric(),
#   CI_Upper = numeric(),
#   Computation_Time_sec = numeric(),
#   stringsAsFactors = FALSE
# )
# 
# results_cv <- data.frame(
#   Method = character(),
#   Sample_Size = numeric(),
#   Option_Price = numeric(),
#   Standard_Error = numeric(),
#   CI_Lower = numeric(),
#   CI_Upper = numeric(),
#   Computation_Time_sec = numeric(),
#   stringsAsFactors = FALSE
# )
# 
# results_cv_strat <- data.frame(
#   Method = character(),
#   Sample_Size = numeric(),
#   Option_Price = numeric(),
#   Standard_Error = numeric(),
#   CI_Lower = numeric(),
#   CI_Upper = numeric(),
#   Computation_Time_sec = numeric(),
#   stringsAsFactors = FALSE
# )
# 
# results_cv_bb <- data.frame(
#   Method = character(),
#   Sample_Size = numeric(),
#   Option_Price = numeric(),
#   Standard_Error = numeric(),
#   CI_Lower = numeric(),
#   CI_Upper = numeric(),
#   Computation_Time_sec = numeric(),
#   stringsAsFactors = FALSE
# )
# 
# # Loop over different sample sizes
# for (N in sample_sizes) {
#   cat("Processing N =", N, "...\n")
#   
#   # Monte Carlo with Control Variate and Stratified Sampling
#   # Note: Stratified Sampling requires N to be divisible by strata (default strata=10)
#   if (N %% 10 == 0) {
#     res_cv_strat <- price_asian_call_MC_control_variate_stratified(S0, K, T, r, q, sigma, m, N, strata = 10)
#     results_cv_strat <- rbind(results_cv_strat, data.frame(
#       Method = res_cv_strat$Method,
#       Sample_Size = res_cv_strat$N,
#       Option_Price = res_cv_strat$Option_Price,
#       Standard_Error = res_cv_strat$Standard_Error,
#       CI_Lower = res_cv_strat$Confidence_Interval[1],
#       CI_Upper = res_cv_strat$Confidence_Interval[2],
#       Computation_Time_sec = res_cv_strat$Computation_Time_sec
#     ))
#   }
#   
#   # Monte Carlo with Control Variate and Brownian Bridge
#   res_cv_bb <- price_asian_call_MC_control_variate_brownian_bridge(S0, K, T, r, q, sigma, m, N)
#   results_cv_bb <- rbind(results_cv_bb, data.frame(
#     Method = res_cv_bb$Method,
#     Sample_Size = res_cv_bb$N,
#     Option_Price = res_cv_bb$Option_Price,
#     Standard_Error = res_cv_bb$Standard_Error,
#     CI_Lower = res_cv_bb$Confidence_Interval[1],
#     CI_Upper = res_cv_bb$Confidence_Interval[2],
#     Computation_Time_sec = res_cv_bb$Computation_Time_sec
#   ))
# }
# 
# final_results <- merge(results_cv_strat, results_cv_bb, by = "Sample_Size")
# print(final_results)
# 
# # Order the results by Sample_Size and Method
# # final_results <- final_results[order(final_results$Sample_Size, final_results$Method), ]
# 
# # Reset row names
# rownames(final_results) <- NULL
# 
# # Print the final results
# print(final_results)
# 
# #     Sample_Size                              Method.x Option_Price.x Standard_Error.x CI_Lower.x CI_Upper.x Computation_Time_sec.x
# # 1        1000 Control Variate + Stratified Sampling           7.17          0.00936       7.15       7.19                   0.01
# # 2        4000 Control Variate + Stratified Sampling           7.16          0.00411       7.16       7.17                   0.05
# # 3       16000 Control Variate + Stratified Sampling           7.17          0.00201       7.16       7.17                   0.17
# # 4       64000 Control Variate + Stratified Sampling           7.17          0.00102       7.16       7.17                   0.64
# # 5      128000 Control Variate + Stratified Sampling           7.16          0.00071       7.16       7.17                   1.95
# #                                              Method.y Option_Price.y Standard_Error.y CI_Lower.y CI_Upper.y Computation_Time_sec.y
# # 1        1000 Control Variate + Brownian Bridge           7.17          0.00823       7.15       7.18                   0.01
# # 2        4000 Control Variate + Brownian Bridge           7.16          0.00395       7.15       7.17                   0.03
# # 3       16000 Control Variate + Brownian Bridge           7.16          0.00202       7.16       7.17                   0.14
# # 4       64000 Control Variate + Brownian Bridge           7.17          0.00102       7.16       7.17                   0.41
# # 5      128000 Control Variate + Brownian Bridge           7.17          0.00072       7.16       7.17                   1.24
# 
# # Q3
# 
# 
# 
# install.packages("knitr")
# install.packages("kableExtra")
# 

















