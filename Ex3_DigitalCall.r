# Exercise3_Digital_Call_Option.R
# Monte Carlo Pricing of Digital Call Option

library(ggplot2)
library(dplyr)

set.seed(42)

# Digital call price functions
digital_call_BS <- function(S0, K, r, sigma, T) {
  d2 <- (log(S0/K) + (r - 0.5 * sigma^2) * T) / (sigma * sqrt(T))
  return(exp(-r * T) * pnorm(d2))
}

digital_call_MC <- function(S0, K, r, sigma, T, M) {
  Z <- rnorm(M)
  ST <- S0 * exp((r - 0.5 * sigma^2) * T + sigma * sqrt(T) * Z)
  payoffs <- as.numeric(ST > K)
  discount_factor <- exp(-r * T)
  price <- discount_factor * mean(payoffs)
  se <- discount_factor * sd(payoffs) / sqrt(M)
  return(list(price = price, se = se))
}

vanilla_call_MC <- function(S0, K, r, sigma, T, M) {
  Z <- rnorm(M)
  ST <- S0 * exp((r - 0.5 * sigma^2) * T + sigma * sqrt(T) * Z)
  payoffs <- pmax(ST - K, 0)
  discount_factor <- exp(-r * T)
  price <- discount_factor * mean(payoffs)
  se <- discount_factor * sd(payoffs) / sqrt(M)
  return(list(price = price, se = se))
}

# Vanilla call Black-Scholes
vanilla_call_BS <- function(S0, K, r, sigma, T) {
  d1 <- (log(S0/K) + (r + 0.5 * sigma^2) * T) / (sigma * sqrt(T))
  d2 <- d1 - sigma * sqrt(T)
  return(S0 * pnorm(d1) - K * exp(-r * T) * pnorm(d2))
}

# Parameters (same as from Section 5.2)
S0 <- 5
K <- 5
r <- 0.06
sigma <- 0.3
T <- 1

# Analytical prices
BS_digital <- digital_call_BS(S0, K, r, sigma, T)
BS_vanilla <- vanilla_call_BS(S0, K, r, sigma, T)

# Sample sizes
M_values <- c(1000, 10000, 100000, 1000000)

# Results storage
digital_results <- data.frame()
vanilla_results <- data.frame()

# Run simulations
for (M in M_values) {
  # Digital call
  dig <- digital_call_MC(S0, K, r, sigma, T, M)
  digital_results <- rbind(digital_results, data.frame(
    M = M,
    Price = dig$price,
    SE = dig$se,
    Error = dig$price - BS_digital,
    Abs_Error = abs(dig$price - BS_digital)
  ))
  
  # Vanilla call
  van <- vanilla_call_MC(S0, K, r, sigma, T, M)
  vanilla_results <- rbind(vanilla_results, data.frame(
    M = M,
    Price = van$price,
    SE = van$se,
    Error = van$price - BS_vanilla,
    Abs_Error = abs(van$price - BS_vanilla)
  ))
}

# Print results
cat("\n======= Exercise 3 Results =======\n")
cat("\nParameters:")
cat(sprintf("\nS0 = %.2f, K = %.2f, r = %.2f, sigma = %.2f, T = %.2f", 
            S0, K, r, sigma, T))

cat("\n\nAnalytical Prices:")
cat(sprintf("\nDigital Call: %.6f", BS_digital))
cat(sprintf("\nVanilla Call: %.6f", BS_vanilla))

cat("\n\nDigital Call Monte Carlo Results:")
print(digital_results, digits = 6)

cat("\n\nVanilla Call Monte Carlo Results:")
print(vanilla_results, digits = 6)

cat("\n\nComparison:")
digital_results$Option <- "Digital"
vanilla_results$Option <- "Vanilla"
comparison <- rbind(digital_results[, c("M", "SE", "Option")],
                    vanilla_results[, c("M", "SE", "Option")])

# Plot comparison
p1 <- ggplot(comparison, aes(x = M, y = SE, color = Option)) +
  geom_point(size = 3) +
  geom_line() +
  scale_x_log10() +
  scale_y_log10() +
  labs(title = "Standard Error Comparison: Digital vs Vanilla Call",
       subtitle = "Digital option has higher variance due to discontinuous payoff",
       x = "Number of Simulations (M)",
       y = "Standard Error (log scale)") +
  theme_minimal() +
  theme(legend.position = "bottom")
print(p1)

# Plot convergence with confidence bands
digital_results$lower <- digital_results$Price - 1.96 * digital_results$SE
digital_results$upper <- digital_results$Price + 1.96 * digital_results$SE
vanilla_results$lower <- vanilla_results$Price - 1.96 * vanilla_results$SE
vanilla_results$upper <- vanilla_results$Price + 1.96 * vanilla_results$SE

# Digital call convergence plot
p2 <- ggplot(digital_results, aes(x = M, y = Price)) +
  geom_point(size = 3, color = "blue") +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2, alpha = 0.5) +
  geom_hline(yintercept = BS_digital, linetype = "dashed", color = "red", size = 1) +
  scale_x_log10() +
  labs(title = "Digital Call: Monte Carlo Convergence",
       x = "Number of Simulations (M)",
       y = "Option Price") +
  annotate("text", x = max(M_values), y = BS_digital + 0.01, 
           label = "BS Price", hjust = 1) +
  theme_minimal()
print(p2)

# Combined convergence plot
digital_results$Type <- "Digital"
vanilla_results$Type <- "Vanilla"
combined <- rbind(
  digital_results[, c("M", "Price", "SE", "Type")],
  vanilla_results[, c("M", "Price", "SE", "Type")]
)

p3 <- ggplot(combined, aes(x = M, y = Price, color = Type)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = Price - 1.96 * SE, ymax = Price + 1.96 * SE), 
                width = 0.2, alpha = 0.5) +
  geom_hline(yintercept = BS_digital, linetype = "dashed", color = "blue", size = 1) +
  geom_hline(yintercept = BS_vanilla, linetype = "dashed", color = "red", size = 1) +
  scale_x_log10() +
  labs(title = "Convergence Comparison",
       subtitle = "Digital vs Vanilla Call Options",
       x = "Number of Simulations (M)",
       y = "Option Price") +
  theme_minimal() +
  theme(legend.position = "bottom")
print(p3)

# Save results
saveRDS(list(
  digital = digital_results,
  vanilla = vanilla_results,
  BS = data.frame(Digital = BS_digital, Vanilla = BS_vanilla)
), "exercise3_results.rds")

# Additional analysis: Variance comparison
cat("\n\nVariance Comparison:")
cat(sprintf("\nFor M = 1000000:"))
cat(sprintf("\nDigital Option SE: %.6f", digital_results$SE[4]))
cat(sprintf("\nVanilla Option SE: %.6f", vanilla_results$SE[4]))
cat(sprintf("\nSE Ratio (Digital/Vanilla): %.2f", 
            digital_results$SE[4] / vanilla_results$SE[4]))
cat("\n\nThe digital call has larger standard error because:")
cat("\n1. Payoff is discontinuous (indicator function)")
cat("\n2. Bernoulli variable has variance p(1-p) which is maximized at p=0.5")
cat("\n3. The vanilla call payoff is continuous and smoother")