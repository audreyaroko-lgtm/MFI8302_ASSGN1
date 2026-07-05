# Exercise1_European_Options_Put_Call_Parity.R
# Monte Carlo Pricing of European Call and Put Options

# Load required libraries
library(ggplot2)
library(dplyr)

# Set seed for reproducibility
set.seed(42)

# Function to price European call and put using Monte Carlo
price_european_options <- function(S0, K, r, sigma, T, M) {
  # Generate standard normal random variables
  Z <- rnorm(M)
  
  # Exact GBM simulation - terminal stock price
  ST <- S0 * exp((r - 0.5 * sigma^2) * T + sigma * sqrt(T) * Z)
  
  # Calculate payoffs
  call_payoffs <- pmax(ST - K, 0)
  put_payoffs <- pmax(K - ST, 0)
  
  # Discounted prices
  discount_factor <- exp(-r * T)
  call_price <- discount_factor * mean(call_payoffs)
  put_price <- discount_factor * mean(put_payoffs)
  
  # Standard errors
  call_se <- discount_factor * sd(call_payoffs) / sqrt(M)
  put_se <- discount_factor * sd(put_payoffs) / sqrt(M)
  
  return(list(
    call_price = call_price,
    call_se = call_se,
    put_price = put_price,
    put_se = put_se,
    call_payoffs = call_payoffs,
    put_payoffs = put_payoffs
  ))
}

# Parameters
S0 <- 100
K <- 100
r <- 0.05
sigma <- 0.2
T <- 1

# Theoretical values (Black-Scholes)
d1 <- (log(S0/K) + (r + 0.5 * sigma^2) * T) / (sigma * sqrt(T))
d2 <- d1 - sigma * sqrt(T)
BS_call <- S0 * pnorm(d1) - K * exp(-r * T) * pnorm(d2)
BS_put <- K * exp(-r * T) * pnorm(-d2) - S0 * pnorm(-d1)

# Put-call parity theoretical
put_call_parity_theoretical <- S0 - K * exp(-r * T)

# Sample sizes
M_values <- c(1000, 10000, 100000, 1000000)

# Results storage
results <- data.frame()

# Run simulations for each M
for (M in M_values) {
  result <- price_european_options(S0, K, r, sigma, T, M)
  
  # Calculate put-call parity difference
  pc_diff <- result$call_price - result$put_price
  pc_diff_theoretical <- put_call_parity_theoretical
  pc_error <- pc_diff - pc_diff_theoretical
  
  # Standard error of the difference
  # Var(C - P) = Var(C) + Var(P) - 2*Cov(C,P)
  # Approximate using individual SEs (assuming independence)
  combined_se <- sqrt(result$call_se^2 + result$put_se^2)
  
  # Check if within 2 SE
  within_2se <- abs(pc_error) < 2 * combined_se
  
  results <- rbind(results, data.frame(
    M = M,
    Call_Price = result$call_price,
    Call_SE = result$call_se,
    Put_Price = result$put_price,
    Put_SE = result$put_se,
    PC_Diff = pc_diff,
    PC_Theoretical = pc_diff_theoretical,
    PC_Error = pc_error,
    Combined_SE = combined_se,
    Within_2SE = within_2se
  ))
}

# Print results
cat("\n======= Exercise 1 Results =======\n")
cat("\nParameters:")
cat(sprintf("\nS0 = %.2f, K = %.2f, r = %.2f, sigma = %.2f, T = %.2f", 
            S0, K, r, sigma, T))
cat("\n\nBlack-Scholes Prices:")
cat(sprintf("\nCall: %.4f", BS_call))
cat(sprintf("\nPut: %.4f", BS_put))
cat(sprintf("\nPut-Call Parity (S0 - K*e^(-rT)): %.4f\n", put_call_parity_theoretical))

cat("\n\nMonte Carlo Results:")
print(results, digits = 6)

# Convergence plot
convergence_data <- data.frame(
  M = M_values,
  Call_Price = results$Call_Price,
  Call_SE = results$Call_SE,
  Put_Price = results$Put_Price,
  Put_SE = results$Put_SE
)

# Plot call and put convergence
p1 <- ggplot(convergence_data, aes(x = M)) +
  geom_point(aes(y = Call_Price, color = "Call")) +
  geom_errorbar(aes(ymin = Call_Price - 1.96 * Call_SE, 
                    ymax = Call_Price + 1.96 * Call_SE), 
                color = "blue", width = 0.2) +
  geom_point(aes(y = Put_Price, color = "Put")) +
  geom_errorbar(aes(ymin = Put_Price - 1.96 * Put_SE, 
                    ymax = Put_Price + 1.96 * Put_SE), 
                color = "red", width = 0.2) +
  geom_hline(yintercept = BS_call, linetype = "dashed", color = "blue") +
  geom_hline(yintercept = BS_put, linetype = "dashed", color = "red") +
  scale_x_log10() +
  labs(title = "Monte Carlo Convergence for European Options",
       x = "Number of Simulations (M)",
       y = "Option Price",
       color = "Option Type") +
  theme_minimal()

print(p1)

# Plot put-call parity check
p2 <- ggplot(results, aes(x = M)) +
  geom_point(aes(y = PC_Error)) +
  geom_hline(yintercept = 0, linetype = "solid") +
  geom_hline(yintercept = c(-2, 2) * results$Combined_SE, 
             linetype = "dashed", color = "red") +
  scale_x_log10() +
  labs(title = "Put-Call Parity Error",
       x = "Number of Simulations (M)",
       y = "Error (MC - Theoretical)") +
  theme_minimal()

print(p2)

cat("\n\nConvergence Rate Analysis:")
cat("\nError decreases as O(1/sqrt(M))")
cat("\nThe standard errors decrease by approximately sqrt(M) when M increases by factor 10\n")

# Save results
saveRDS(results, "exercise1_results.rds")
