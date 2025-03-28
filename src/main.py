"""Main module."""
import numpy as np
import saw_monte_carlo as smc
from saw_monte_carlo.utils import get_deviation

# Use the deterministic algorithm to count SAWs of length L
L = 10
trials = 500000
# print("Number of SAWs of length", L, "=", smc.count_saws(L))

# Estimate using naive Monte Carlo
cL_est_MC = smc.estimate_cL(L, N=trials)
print(f"Naive MC estimate for c_{L} ≈ {cL_est_MC}")
print(f"Naive MC estimate for mu ≈ {cL_est_MC**(1/L)}")
print("Naive MC estimate deviation:", get_deviation(cL_est_MC, L))

# Rosenbluth example
cL_est_rosenbluth = smc.rosenbluth_estimate_cL(L, trials=trials)
print(f"Vanilla Rosenbluth estimate for c_{L} ≈ {cL_est_rosenbluth}")
print(f"Vanilla Rosenbluth estimate for mu ≈ {cL_est_rosenbluth**(1/L)}")
print("Rosenbluth estimate deviation:", get_deviation(cL_est_rosenbluth, L))

# PERM example
cL_est_perm = smc.perm_estimate_cL(L, trials, c_minus=0.2, c_plus=3.0)
print(f"PERM estimate c_{L} ≈ {cL_est_perm}")
print(f"PERM estimate mu ≈ {cL_est_perm**(1/L)}")
print("PERM estimate deviation:", get_deviation(cL_est_perm, L))

# Pivot MCMC
mu_est = smc.run_pivot_get_mu_estimate(n=100, pivot_attempts=200000, burn_in=20000)
print("Pivot mu estimate for n=100:", mu_est)