# Monte Carlo Methods in Finance - Solutions

## Overview

This repository contains R implementations for the Monte Carlo methods exercises from the "Part III: Introduction to Monte Carlo Methods in Finance" course material.

## Exercises

### Exercise 1: European Options and Put-Call Parity

**File:** `Exercise1_European_Options_Put_Call_Parity.R`

**Description:**
- Implements plain Monte Carlo estimator using exact GBM simulation
- Prices both European call and put options
- Verifies put-call parity across different sample sizes
- Analyzes convergence rate

**Parameters:**
- S0 = 100, K = 100, r = 0.05, σ = 0.2, T = 1

**Outputs:**
- Monte Carlo price estimates with standard errors
- Put-call parity verification
- Convergence plots

### Exercise 2: Convergence Study

**File:** `Exercise2_Convergence_Study.R`

**Description:**
- Implements path-dependent simulation using Euler scheme
- Studies convergence with respect to M (number of simulations)
- Studies discretization bias with respect to m (time steps)
- Creates comprehensive convergence table
- Determines optimal (M*, m*) allocation for given budget

**Parameters:**
- S0 = 100, K = 100, r = 0.05, σ = 0.2, T = 1

**Outputs:**
- Standard error vs M analysis
- Bias vs 1/m analysis
- Convergence heatmap
- Optimal allocation recommendation

### Exercise 3: Digital Call Option

**File:** `Exercise3_Digital_Call_Option.R`

**Description:**
- Computes analytical Black-Scholes price for digital call
- Implements Monte Carlo pricing for digital call
- Compares convergence with vanilla call
- Plots convergence curves with confidence bands

**Parameters:**
- S0 = 5, K = 5, r = 0.06, σ = 0.3, T = 1

**Outputs:**
- Digital call Monte Carlo estimates
- Comparison with vanilla call convergence
- Confidence band plots

## Prerequisites

### R Packages Required

```r
install.packages(c("ggplot2", "dplyr"))
