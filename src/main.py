"""Main module."""
import numpy as np
import saw_monte_carlo as smc

# Use the deterministic algorithm to count SAWs of length L
L = 10
print("Number of SAWs of length", L, "=", smc.count_saws(L))

# Estimate using naive Monte Carlo
est_cL = smc.estimate_cL(L, N=200000)
print("Naive MC estimate:", est_cL)

# PERM example
smc.init_statistics(L)
trials = 200000
sum_est = 0
for _ in range(trials):
    visited = {(0,0)}
    sum_est += smc.perm_grow(0, 0, 0, visited, weight=1.0, L_max=L)
print("PERM estimate:", sum_est / trials)

# Pivot MCMC
mu_est = smc.run_pivot_get_mu_estimate(n=100, pivot_attempts=200000, burn_in=20000)
print("Pivot mu estimate for n=100:", mu_est)