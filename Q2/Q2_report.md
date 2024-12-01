---
title: "Q2_main"
date: "2024-10-30"
output: pdf_document
---


The table below compares the results of pricing Asian call options using the **Standard Monte Carlo Method** and the **Control Variate Method** for various sample sizes.


\section*{2.Monte Carlo on Asian Call Option}

\textbf{2.1 Introduction}


An **Asian Option** is a type of exotic option where the payoff depends on the average price of the underlying asset over a certain period, rather than its price at a specific point in time. This average can be computed as either an arithmetic average or a geometric average. Asian options are particularly useful in markets where the underlying asset's price is highly volatile, as they reduce the impact of price manipulation or extreme price movements near expiration.

For an Asian call or put option, the payoff depends on whether the average price of the underlying asset is greater than (for calls) or less than (for puts) the strike price \( K \). The payoff formulas are listed as follows.

The payoff for a continuous arithmetic average Asian call or put option is:

\[
\Phi(S) = \max\left(\frac{1}{T} \int_0^T S(t) dt - K, 0\right)
\quad \text{or} \quad
\Phi(S) = \max\left(K - \frac{1}{T} \int_0^T S(t) dt, 0\right).
\]

The payoff for a continuous geometric average Asian call or put option is:

\[
\Phi(S) = \max\left(e^{\frac{1}{T} \int_0^T \log S(t) dt} - K, 0\right)
\quad \text{or} \quad
\Phi(S) = \max\left(K - e^{\frac{1}{T} \int_0^T \log S(t) dt}, 0\right).
\]

For discrete monitoring, the arithmetic average Asian call or put option has the following payoff:

\[
\Phi(S) = \max\left(\frac{1}{m+1} \sum_{i=0}^m S\left(\frac{iT}{m}\right) - K, 0\right)
\quad \text{or} \quad
\Phi(S) = \max\left(K - \frac{1}{m+1} \sum_{i=0}^m S\left(\frac{iT}{m}\right), 0\right).
\]

Similarly, for discrete monitoring, the geometric average Asian call or put option has the payoff:

\[
\Phi(S) = \max\left(e^{\frac{1}{m+1} \sum_{i=0}^m \log S\left(\frac{iT}{m}\right)} - K, 0\right)
\quad \text{or} \quad
\Phi(S) = \max\left(K - e^{\frac{1}{m+1} \sum_{i=0}^m \log S\left(\frac{iT}{m}\right)}, 0\right).
\]

Asian options are widely used in financial markets for energy and commodity trading. Their structure reduces the risk of price manipulation near expiration and smooths out price volatility through averaging. 

Despite their practical utility, pricing Asian options analytically is challenging because the arithmetic average does not lead to a closed-form solution for the option price. However, the geometric average Asian option can be priced analytically under the Black-Scholes framework, providing a crucial benchmark for numerical methods.

In the following part, we will focus on:

1. Pricing an arithmetic average Asian call option using the **Monte Carlo simulation** method.
2. Enhancing the efficiency of the Monte Carlo method by implementing **variance reduction techniques**, including the **Control Variate Method**.
3. Exploring further variance reduction methods, such as combining control variates with **stratified sampling** and the **Brownian bridge**.

By examining these approaches, we will highlights the trade-offs between computational cost and accuracy in Asian option pricing.

\textbf{2.2 Asian Option Pricing-- Monte Carlo Method}

Monte Carlo (MC) simulation is a widely used numerical method for pricing financial derivatives, including Asian options. The MC method involves generating a large number of simulated paths for the underlying asset price and computing the average payoff over these paths to estimate the option price. Its flexibility makes it suitable for pricing Asian options, as it can handle the arithmetic averaging of the underlying prices, which does not have a closed-form solution.


The key steps in pricing an Asian call option using the Monte Carlo method are as follows:

1. **Simulate Asset Price Paths**: 
   The underlying asset price \( S_t \) is modeled as a geometric Brownian motion (GBM) under the risk-neutral measure:
   \[
   dS_t = (r - q) S_t \, dt + \sigma S_t \, dW_t,
   \]
   where:
   - \( r \): Risk-free interest rate,
   - \( q \): Dividend yield,
   - \( \sigma \): Volatility of the asset,
   - \( W_t \): Standard Brownian motion.

   Discretizing this SDE using the Euler-Maruyama scheme over \( m \) monitoring points gives:
   \[
   S_{t_{i+1}} = S_{t_i} \exp\left(\left(r - q - \frac{\sigma^2}{2}\right) \Delta t + \sigma \sqrt{\Delta t} Z_i\right),
   \]
   where \( Z_i \) are independent standard normal random variables, and \( \Delta t = T / m \) is the time step.

2. **Calculate Average Prices**: 
   For each simulated path, compute the arithmetic average of the asset prices:
   \[
   \bar{S}_{\text{arith}} = \frac{1}{m} \sum_{i=1}^m S_{t_i}.
   \]

3. **Determine Payoff**: 
   Calculate the payoff of the option for each path as:
   \[
   \text{Payoff}_i = \max(\bar{S}_{\text{arith}} - K, 0),
   \]
   where \( K \) is the strike price.

4. **Discount Payoffs**: 
   Discount the payoff to present value using the risk-free rate:
   \[
   \text{Discounted Payoff}_i = e^{-rT} \cdot \text{Payoff}_i.
   \]

5. **Estimate Option Price**: 
   The Monte Carlo estimate of the Asian option price is the average of the discounted payoffs:
   \[
   \text{Option Price} = \frac{1}{N} \sum_{i=1}^N \text{Discounted Payoff}_i,
   \]
   where \( N \) is the number of simulated paths.

6. **Measure Uncertainty**: 
   The standard error of the Monte Carlo estimate is given by:
   \[
   \text{Standard Error} = \frac{\text{Standard Deviation of Payoffs}}{\sqrt{N}}.
   \]

Below is an R function to compute the Asian call option price using the Monte Carlo method:

```r
price_asian_call_MC <- function(S0, K, T, r, q, sigma, m, N) {
  delta_t <- T / m
  sqrt_delta_t <- sqrt(delta_t)
  
  # Generate N x m standard normal random variables
  Z <- matrix(rnorm(N * m), nrow = N, ncol = m)
  
  drift <- (r - q - 0.5 * sigma^2) * delta_t
  diffusion <- sigma * sqrt_delta_t * Z
  
  # Simulate log prices and then transform to prices
  W <- t(apply(diffusion, 1, cumsum))
  log_S <- log(S0) + outer(rep(1, N), drift * (1:m)) + W
  S <- exp(log_S)
  
  # Calculate arithmetic average and payoffs
  S_bar <- rowMeans(S)
  payoffs <- pmax(S_bar - K, 0)
  discounted_payoffs <- exp(-r * T) * payoffs
  
  # Estimate option price
  option_price <- mean(discounted_payoffs)
  std_error <- sd(discounted_payoffs) / sqrt(N)
  
  return(list(
    Option_Price = option_price,
    Standard_Error = std_error
  ))
}
```

The accuracy of Monte Carlo simulation is governed by the **Central Limit Theorem**. The error in the estimated option price decreases as the number of simulated paths \( N \) increases. Specifically, the standard error of the Monte Carlo estimate is proportional to \( \mathcal{O}(1/\sqrt{N}) \), where \( N \) is the number of simulations.

This relationship implies a relatively slow convergence rate compared to other numerical methods, such as finite difference methods or quadrature techniques. Achieving higher accuracy requires a significant increase in the number of simulated paths, leading to higher computational cost.

To halve the error, the number of simulations must be increased by a factor of four. While the Monte Carlo method is flexible and robust, this slow convergence rate necessitates the use of **variance reduction techniques** to improve efficiency.

\textbf{2.3 Variance Reduction Techniques-- Control Variate Method}

We want to estimate the expected value \( \theta = \mathbb{E}[Y] \), where \( Y = g(X) \) is a function of some random variable \( X \). If we can find another random variable \( Z \), for which the expected value \( \mathbb{E}[Z] \) is known, we can construct alternative estimators for \( \theta \). These include:

1. The standard Monte Carlo estimator:
   \[
   \hat{\theta} = Y
   \]

2. The control variate estimator:
   \[
   \hat{\theta}_c = Y + c \cdot (Z - \mathbb{E}[Z]),
   \]
   where \( c \) is a constant.

It can be shown that the control variate estimator \( \hat{\theta}_c \) is unbiased, as:
\[
\mathbb{E}[\hat{\theta}_c] = \mathbb{E}[Y] + c \cdot (\mathbb{E}[Z] - \mathbb{E}[Z]) = \mathbb{E}[Y] = \theta.
\]

To minimize the variance of \( \hat{\theta}_c \), we compute its variance:
\[
\text{Var}(\hat{\theta}_c) = \text{Var}(Y) + c^2 \cdot \text{Var}(Z) + 2c \cdot \text{Cov}(Y, Z).
\]

Using calculus, the optimal \( c \) that minimizes the variance is:
\[
c_{\text{opt}} = -\frac{\text{Cov}(Y, Z)}{\text{Var}(Z)}.
\]

Substituting \( c_{\text{opt}} \) into the variance formula gives:
\[
\text{Var}(\hat{\theta}_{c_{\text{opt}}}) = \text{Var}(Y) - \frac{\text{Cov}(Y, Z)^2}{\text{Var}(Z)}.
\]

This demonstrates that the control variate method reduces variance by leveraging the correlation between \( Y \) and \( Z \), particularly when they are highly correlated.

Using geometric Asian Call as a control variant is efficient, since the correlation between the arithmetic average and the geometric average is very high. Also, the counterpart geometric average Asian option price has a closed-form solution

\[
\begin{aligned}
\sigma_z^2 &= \sigma^2 \cdot \frac{(m + 1)(2m + 1)}{6 m^2} \\
\mu &= \left( r - q - \frac{1}{2} \sigma^2 \right) \cdot \frac{m + 1}{2 m} + \frac{1}{2} \sigma_z^2 \\
d_1 &= \frac{\ln\left( \frac{S_0}{K} \right) + \left( \mu + \frac{1}{2} \sigma_z^2 \right) T}{\sigma_z \sqrt{T}} \\
d_2 &= \frac{\ln\left( \frac{S_0}{K} \right) + \left( \mu - \frac{1}{2} \sigma_z^2 \right) T}{\sigma_z \sqrt{T}} \\
\text{Geometric Asian Call Price} &= e^{-r T} \left[ S_0 \cdot e^{\mu T} \cdot N(d_1) - K \cdot N(d_2) \right]
\end{aligned}
\]


Below is the R implementation of the control variate method:

```r
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
```


\textbf{2.4 Results and Discussion}

The following tables summarize the results obtained from pricing the Asian call option using the **standard Monte Carlo method** and the **control variate method**. 


From both method, we can see that as the sample size increase 4 times, the standard errors decrease half, confirmed that the standard error of the Monte Carlo estimate is proportional to \( \mathcal{O}(1/\sqrt{N}) \), where \( N \) is the number of simulations.

Considering standard error, the standard Monte Carlo method produces significantly higher standard errors compared to the control variate method for the same sample sizes. This difference in precision is evident across all sample sizes. For instance, with 1,000 samples, the standard error for the standard Monte Carlo method is 0.28892, whereas the control variate method achieves a much smaller standard error of 0.00877. Even with a larger sample size of 128,000, the standard Monte Carlo method’s standard error reduces to 0.02433, which is still considerably higher than the control variate method’s 0.00072. This demonstrates the remarkable variance reduction achieved by the control variate method.

While both methods exhibit similar computation times, the control variate method incurs a slight increase in runtime due to the additional calculations required for covariance, variance, and the optimal coefficient \( \theta \). For instance, with 128,000 samples, the standard Monte Carlo method requires 1.55 seconds, whereas the control variate method takes 1.69 seconds. Despite this small increase, the control variate method's significant reduction in standard error far outweighs the marginally higher computation time, making it an efficient approach.

By leveraging the high correlation between the arithmetic and geometric Asian call options, it achieves remarkable variance reduction with minimal additional computational cost.

The following tables summarize the results obtained from pricing the Asian call option using the **standard Monte Carlo method** and the **control variate method**. 


```{r setup, include=FALSE}
# Load necessary libraries
library(knitr)
library(kableExtra)

results_standard_MC <- data.frame(
  Sample_Size = c(1000, 4000, 16000, 64000, 256000),
  Option_Price = c(7.02, 7.39, 7.09, 7.15, 7.18),
  Standard_Error = c(0.26869, 0.13957, 0.06741, 0.03433, 0.01717),
  CI_Lower = c(6.49, 7.12, 6.96, 7.08, 7.14),
  CI_Upper = c(7.55, 7.66, 7.22, 7.22, 7.21),
  Computation_Time_sec = c(0.00, 0.03, 0.14, 0.65, 2.80)
)

# Generate the table with kable
kable(
  results_standard_MC, 
  caption = "Results for standard Monte Carlo Method",
  booktabs = TRUE
) %>%
  kable_styling(latex_options = c("striped", "hold_position", "scale_down"))



results_cv_MC <- data.frame(
  Sample_Size = c(1000, 4000, 16000, 64000, 256000),
  Option_Price_CV = c(7.18, 7.16, 7.17, 7.17, 7.17),
  Standard_Error_CV = c(0.00844, 0.00400, 0.00203, 0.00101, 0.00051),
  CI_Lower_CV = c(7.16, 7.16, 7.16, 7.16, 7.16),
  CI_Upper_CV = c(7.19, 7.17, 7.17, 7.17, 7.17),
  Computation_Time_sec = c(0.02, 0.06, 0.21, 0.72, 3.22)
)

# Generate the table with kable
kable(
  results_cv_MC, 
  caption = "Results for Monte Carlo Method with Control Variate",
  booktabs = TRUE
) %>%
  kable_styling(latex_options = c("striped", "hold_position", "scale_down"))
```

```{r setup, include=FALSE}
# Load necessary libraries
library(ggplot2)
# Data for both methods
sample_size <- c(1000, 4000, 16000, 64000, 256000)

# Standard Monte Carlo Results
option_prices_mc <- c(7.02, 7.39, 7.09, 7.15, 7.18)
standard_error_mc <- c(0.26869, 0.13957, 0.06741, 0.03433, 0.01717)
computation_time_mc <- c(0.00, 0.03, 0.14, 0.65, 2.80)

# Control Variate Results
option_prices_cv <- c(7.18, 7.16, 7.17, 7.17, 7.17)
standard_error_cv <- c(0.00844, 0.00400, 0.00203, 0.00101, 0.00051)
computation_time_cv <- c(0.02, 0.06, 0.21, 0.72, 3.22)

# Combine into data frames for plotting
data_option_prices <- data.frame(
  Sample_Size = rep(sample_size, 2),
  Option_Price = c(option_prices_mc, option_prices_cv),
  Method = rep(c("Standard Monte Carlo", "Control Variate"), each = 5)
)

data_standard_error <- data.frame(
  Sample_Size = rep(sample_size, 2),
  Standard_Error = c(standard_error_mc, standard_error_cv),
  Method = rep(c("Standard Monte Carlo", "Control Variate"), each = 5)
)

data_computation_time <- data.frame(
  Sample_Size = rep(sample_size, 2),
  Computation_Time = c(computation_time_mc, computation_time_cv),
  Method = rep(c("Standard Monte Carlo", "Control Variate"), each = 5)
)


ggplot(data_option_prices, aes(x = Sample_Size, y = Option_Price, color = Method)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  scale_x_log10() +
  labs(
    title = "Option Prices vs. Sample Size",
    x = "Sample Size (Log Scale)",
    y = "Option Price"
  ) +
  theme_minimal()


ggplot(data_standard_error, aes(x = Sample_Size, y = Standard_Error, color = Method)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  scale_x_log10() +
  labs(
    title = "Standard Errors vs. Sample Size",
    x = "Sample Size (Log Scale)",
    y = "Standard Error"
  ) +
  theme_minimal()

ggplot(data_computation_time, aes(x = Sample_Size, y = Computation_Time, color = Method)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  scale_x_log10() +
  labs(
    title = "Computation Time vs. Sample Size",
    x = "Sample Size (Log Scale)",
    y = "Computation Time (s)"
  ) +
  theme_minimal()

```

\textbf{2.5 Further Discussion on Variance Reduction Techniques}

In addition to the control variate method, variance reduction can be further enhanced by combining it with other techniques such as stratified sampling and the Brownian bridge construction. These methods introduce additional refinements to the Monte Carlo simulation process, addressing specific sources of randomness and variability in the simulation.

Stratified sampling divides the sample space into non-overlapping strata and ensures that an equal number of samples are drawn from each stratum. This method reduces variance by minimizing randomness in the sampling process.

Mathematically, the variance of the stratified sampling estimator can be expressed as:
\[
\text{Var}(\hat{\theta}_{\text{stratified}}) = \sum_{k=1}^L \frac{w_k^2}{n_k} \text{Var}(\theta_k),
\]
where:
- \( L \) is the number of strata,
- \( w_k \) is the weight of the \( k \)-th stratum,
- \( n_k \) is the number of samples in the \( k \)-th stratum,
- \( \text{Var}(\theta_k) \) is the variance within the \( k \)-th stratum.

By ensuring that \( n_k \) is proportional to \( w_k \), stratified sampling can achieve a lower variance compared to direct sampling. A common application of stratified sampling in Monte Carlo simulations is in the context of Brownian motion, where stratification is often applied to the first step of the simulation, \( Z_1 \), which introduces the most variability in the path. This step simplifies the implementation and ensures that the variability in the simulation is minimized from the outset.However, stratified sampling has its limitations. It is computationally more expensive than standard Monte Carlo sampling because it requires the sample space to be divided and the simulation to be carefully managed within each stratum. 

The Brownian bridge construction is another technique that can be used in conjunction with the control variate method to reduce variance. The basic idea of the Brownian bridge is to refine the simulation of Brownian paths by starting with the known values of the process at the beginning and end points and then interpolating the intermediate points. This approach introduces a conditional dependency among the simulated points, ensuring that the generated paths are more consistent with the overall process. The Brownian bridge is particularly useful in scenarios where the terminal value of the asset is known to dominate the option payoff.

One advantage of the Brownian bridge is that it is computationally efficient for high-dimensional problems, as it reduces the number of random variables that need to be generated independently. Additionally, it complements the control variate method by improving the accuracy of the simulated paths, which further enhances the correlation between the control variate and the target variable.

However, the implementation of the Brownian bridge introduces additional complexity. It requires careful interpolation of the intermediate points and additional adjustments to the simulation framework. Furthermore, the variance reduction achieved by the Brownian bridge depends on the specific characteristics of the option being priced and the structure of the underlying asset paths. In some cases, the benefits may be marginal compared to the added implementation effort.







