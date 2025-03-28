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

   import saw_monte_carlo as smc

   # 1) Exact counting of SAWs (small L)
   L = 10
   num_saws = smc.count_saws(L)
   print(f"Exact number of SAWs of length {L}: {num_saws}")

   # 2) Naive Monte Carlo estimate of c_L
   N = 200000
   estimate_cL = smc.estimate_cL(L, N=N)
   print(f"Naive Monte Carlo estimate of c_{L}: {estimate_cL}")

   # 3) PERM
   smc.init_statistics(L)
   trials = 200000
   sum_estimates = 0.0
   for _ in range(trials):
       visited = {(0, 0)}
       # Grow up to length L using PERM
       sum_estimates += smc.perm_grow(0, 0, 0, visited, weight=1.0, L_max=L)
   perm_estimate = sum_estimates / trials
   print(f"PERM estimate of c_{L}: {perm_estimate}")

   # 4) Pivot MCMC for mu estimate
   mu_est = smc.run_pivot_get_mu_estimate(n=100, pivot_attempts=200000, burn_in=20000)
   print(f"Pivot-based mu estimate: {mu_est}")

Roadmap
-------

- Parallel implementations.
- More advanced MCMC analyses (e.g., radius of gyration, end-to-end distance).

Credits
-------

This package was created with Cookiecutter_ and the `audreyr/cookiecutter-pypackage`_ project template.

.. _Cookiecutter: https://github.com/audreyr/cookiecutter
.. _`audreyr/cookiecutter-pypackage`: https://github.com/audreyr/cookiecutter-pypackage