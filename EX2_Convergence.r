# Exercise2_Convergence_Study.R
# Path-dependent simulation and convergence analysis

library(ggplot2)
library(dplyr)

set.seed(42)

# Function to simulate GBM path using Euler scheme
simulate_gbm_path <- function(S0, r, sigma, T, m, M) {
  dt <- T / m
  paths <- matrix(0, nrow = M, ncol = m + 1)
  paths[, 1] <- S0
  
  for (k in 1:m) {
    Z <- rnorm(M)
    paths[, k + 1] <- paths[, k] * exp((r - 0.5 * sigma^2) * dt + sigma * sqrt(dt) * Z)
  }
  
  return(paths)
}

# Function to price European call using Euler scheme
price_euler_call <- function(S0, K, r, sigma, T, m, M) {
  paths <- simulate_gbm_path(S0, r, sigma, T, m, M)
  ST <- paths[, m + 1]
  payoffs <- pmax(ST - K, 0)
  discount_factor <- exp(-r * T)
  price <- discount_factor * mean(payoffs)
  se <- discount_factor * sd(payoffs) / sqrt(M)
  return(list(price = price, se = se))
}

# Parameters
S0 <- 100
K <- 100
r <- 0.05
sigma <- 0.2
T <- 1

# Black-Scholes price for comparison
d1 <- (log(S0/K) + (r + 0.5 * sigma^2) * T) / (sigma * sqrt(T))
d2 <- d1 - sigma * sqrt(T)
BS_price <- S0 * pnorm(d1) - K * exp(-r * T) * pnorm(d2)

# Part 1: Fix m = 100, vary M
M_values <- c(100, 1000, 10000, 100000)
m_fixed <- 100

part1_results <- data.frame()

for (M in M_values) {
  result <- price_euler_call(S0, K, r, sigma, T, m_fixed, M)
  part1_results <- rbind(part1_results, data.frame(
    M = M,
    m = m_fixed,
    Price = result$price,
    SE = result$se,
    Error = result$price - BS_price
  ))
}

# Part 2: Fix M = 100000, vary m
M_fixed <- 100000
m_values <- c(10, 50, 100, 500, 1000)

part2_results <- data.frame()

for (m in m_values) {
  result <- price_euler_call(S0, K, r, sigma, T, m, M_fixed)
  part2_results <- rbind(part2_results, data.frame(
    M = M_fixed,
    m = m,
    Price = result$price,
    SE = result$se,
    Error = result$price - BS_price
  ))
}

# Part 3: Complete convergence table
part3_results <- data.frame()

for (M in M_values) {
  for (m in m_values) {
    result <- price_euler_call(S0, K, r, sigma, T, m, M)
    part3_results <- rbind(part3_results, data.frame(
      M = M,
      m = m,
      Price = result$price,
      SE = result$se,
      Error = result$price - BS_price,
      Abs_Error = abs(result$price - BS_price)
    ))
  }
}

# Part 4: Optimal allocation for budget C = 10^7
# Euler scheme has weak order gamma = 1
# Optimal: m* ~ C^(1/3), M* = C/m*
C_budget <- 1e7
m_optimal <- round(C_budget^(1/3))
M_optimal <- round(C_budget / m_optimal)

cat("\n======= Exercise 2 Results =======\n")
cat("\nBlack-Scholes Price:", round(BS_price, 6))

cat("\n\nPart 1: Fixed m = 100, Varying M")
print(part1_results, digits = 6)

cat("\n\nPart 2: Fixed M = 100000, Varying m")
print(part2_results, digits = 6)

cat("\n\nPart 3: Complete Convergence Table")
print(part3_results, digits = 6)

cat("\n\nPart 4: Optimal Allocation")
cat(sprintf("\nBudget C = %.0f", C_budget))
cat(sprintf("\nOptimal time steps (m*): %d", m_optimal))
cat(sprintf("\nOptimal samples (M*): %d", M_optimal))
cat(sprintf("\nTotal operations: %.0f", m_optimal * M_optimal))

# Plots
# Plot 1: SE vs M (log-log)
p1 <- ggplot(part1_results, aes(x = M, y = SE)) +
  geom_point(size = 3) +
  geom_line() +
  geom_smooth(method = "lm", se = FALSE, linetype = "dashed") +
  scale_x_log10() +
  scale_y_log10() +
  labs(title = "Standard Error vs Number of Simulations",
       subtitle = "Fixed m = 100",
       x = "M (log scale)",
       y = "Standard Error (log scale)") +
  theme_minimal()
print(p1)

# Plot 2: Bias vs 1/m
part2_results$inv_m <- 1 / part2_results$m
p2 <- ggplot(part2_results, aes(x = inv_m, y = Price)) +
  geom_point(size = 3) +
  geom_line() +
  geom_hline(yintercept = BS_price, linetype = "dashed", color = "red") +
  labs(title = "Discretization Bias Analysis",
       subtitle = "Fixed M = 100000",
       x = "1/m (time step size)",
       y = "Option Price") +
  theme_minimal()
print(p2)

# Plot 3: Heatmap of absolute error
p3 <- ggplot(part3_results, aes(x = factor(M), y = factor(m), fill = Abs_Error)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "red") +
  labs(title = "Absolute Error Heatmap",
       x = "Number of Simulations (M)",
       y = "Number of Time Steps (m)",
       fill = "Absolute Error") +
  theme_minimal()
print(p3)

# Save results
saveRDS(list(
  part1 = part1_results,
  part2 = part2_results,
  part3 = part3_results,
  optimal = data.frame(m = m_optimal, M = M_optimal)
), "exercise2_results.rds")