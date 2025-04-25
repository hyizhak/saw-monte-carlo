===============
saw-monte-carlo
===============

.. .. image:: https://img.shields.io/pypi/v/saw_monte_carlo.svg
..    :target: https://pypi.python.org/pypi/saw_monte_carlo

.. .. image:: https://img.shields.io/travis/hyizhak/saw_monte_carlo.svg
..    :target: https://travis-ci.com/hyizhak/saw_monte_carlo

.. .. image:: https://readthedocs.org/projects/saw-monte-carlo/badge/?version=latest
..    :target: https://saw-monte-carlo.readthedocs.io/en/latest/?version=latest
..    :alt: Documentation Status

A Python library for estimating properties of **Self-Avoiding Walks (SAWs)** on a 2D square lattice. This includes:

- Estimating the number of SAWs of a given length (``c_L``).
- Estimating the connectivity constant (:math:`\mu`) via **Pivot Monte Carlo**.
- Enumerations, naive Monte Carlo, and advanced methods like **PERM** (Pruned-Enriched Rosenbluth Method).

**Free software**: MIT license  

.. **Documentation**: https://saw-monte-carlo.readthedocs.io

Features
--------

- **Exact enumeration** (recursive counting) for smaller lengths.
- **Naive Monte Carlo** to estimate :math:`c_L`.
- **Rosenbluth** method for estimating :math:`c_L`.
- **PERM** (Pruned-Enriched Rosenbluth Method) for more sophisticated sampling.
- **Pivot MCMC** for estimating the connectivity constant :math:`\mu`.
- Written in pure Python with minimal dependencies (``numpy``).

Installation
------------

You can install ``saw_monte_carlo`` using pip:

.. code-block:: bash

   pip install -e .

Usage Example
-------------

Below is a brief example demonstrating how to use the package’s key functions:

.. code-block:: python

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
   L = 50
   cL_est_rosenbluth = smc.rosenbluth_estimate_cL(L, trials=trials)
   print(f"Vanilla Rosenbluth estimate for c_{L} ≈ {cL_est_rosenbluth}")
   print(f"Vanilla Rosenbluth estimate for mu ≈ {cL_est_rosenbluth**(1/L)}")
   print("Rosenbluth estimate deviation:", get_deviation(cL_est_rosenbluth, L))

   # PERM example
   L = 200
   cL_est_perm = smc.perm_estimate_cL(L, trials, c_minus=0.2, c_plus=3.0)
   print(f"PERM estimate c_{L} ≈ {cL_est_perm}")
   print(f"PERM estimate mu ≈ {cL_est_perm**(1/L)}")
   print("PERM estimate deviation:", get_deviation(cL_est_perm, L))

   # Pivot MCMC example
   # 1. Choose a set of chain–lengths to run:
   ns = [100, 200, 400, 800, 1600]

   # 2. Run pivot+atmosphere on each and collect mu_n estimates:
   mu_estimates = []
   for n in ns:
      burn_in = int(10 * n * np.log(n))          # O(n log n) burn‐in
      pivot_attempts = n * 10000
      mu_n = smc.run_pivot_get_mu_estimate(
         n=n,
         pivot_attempts=pivot_attempts,
         burn_in=burn_in,
         seed=42
      )
      print(f"n={n:4d}, mu_n = {mu_n:.6f}")
      mu_estimates.append(mu_n)

   # 3. Fit μ̂ₙ = μ + A / n via linear regression in (1/n, μ̂ₙ):
   x = 1.0 / np.array(ns)          # independent variable
   y = np.array(mu_estimates)      # observed μ̂ₙ
   A, mu_inf = np.polyfit(x, y, 1) # slope, intercept

   print()
   print(f"Extrapolated mu ≈ {mu_inf:.6f}")

Roadmap
-------

- Parallel implementations.
   - I do have tried this for pivot method but I realized that as every walk needs different pivot attempts at different points,
     it is not possible to parallelize this to achieve any efficiency gain.
   - There are some commented code for numpy vectorization of the pivot method, but I did not observe any significant speedup.
- More advanced MCMC analyses (e.g., radius of gyration, end-to-end distance).

Credits
-------

This package was created with Cookiecutter_ and the `audreyr/cookiecutter-pypackage`_ project template.

.. _Cookiecutter: https://github.com/audreyr/cookiecutter
.. _`audreyr/cookiecutter-pypackage`: https://github.com/audreyr/cookiecutter-pypackage