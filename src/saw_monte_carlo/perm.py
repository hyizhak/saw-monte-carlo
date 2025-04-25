"""
PERM (Pruned-Enriched Rosenbluth Method) for SAWs.
"""

import numpy as np
import sys

def perm_recursive(n, L, weight, pos, visited, rng, stats, c_plus, c_minus, results):
    """
    Recursively grow a self-avoiding walk (SAW) using the PERM algorithm.
    
    Parameters
    ----------
    n : int
        Current length of the walk (number of steps taken so far).
    L : int
        Target length of the walk.
    weight : float
        Current weight of the walk.
    pos : tuple of int
        Current lattice position (x, y).
    visited : set of tuple
        Set of lattice sites already visited by the walk.
    rng : np.random.Generator
        Random number generator for reproducibility.
    stats : dict
        Dictionary holding statistics for each step (keys "weight_sum" and "tours").
        These are used to compute local average weights for pruning/enrichment.
    c_plus : float
        Enrichment threshold multiplier.
    c_minus : float
        Pruning threshold multiplier.
    results : list
        List to collect the weights of completed walks (of length L).
    """
    # Update statistics for the current length.
    stats["tours"][n] += 1
    stats["weight_sum"][n] += weight

    # If the target length is reached, record the completed walk's weight.
    if n == L:
        results.append(weight)
        return

    # Find all allowed moves (neighbors) from the current position.
    x, y = pos
    neighbors = []
    for dx, dy in [(1, 0), (-1, 0), (0, 1), (0, -1)]:
        new_pos = (x + dx, y + dy)
        if new_pos not in visited:
            neighbors.append(new_pos)
    
    # If no moves are available, the walk is trapped.
    if not neighbors:
        return

    # Multiply the weight by the number of available moves (Rosenbluth factor).
    a = len(neighbors)
    new_weight = weight * a

    # Determine the average weight for the next step (if any walks have reached that length).
    if stats["tours"][n + 1] > 0:
        avg = stats["weight_sum"][n + 1] / stats["tours"][n + 1]
    else:
        avg = new_weight  # Default if no data exists yet.

    # Set thresholds for enrichment and pruning.
    W_plus = c_plus * avg
    W_minus = c_minus * avg

    # Enrichment: if the weight is too high, create multiple copies.
    if new_weight > W_plus:
        # Number of clones (at least one)
        m = max(int(new_weight / W_plus), 1)
        # Divide the weight equally among the clones.
        new_weight /= m
        for _ in range(m):
            # Choose one neighbor uniformly at random.
            chosen = neighbors[rng.integers(len(neighbors))]
            new_visited = visited.copy()
            new_visited.add(chosen)
            perm_recursive(n + 1, L, new_weight, chosen, new_visited, rng, stats, c_plus, c_minus, results)
    # Pruning: if the weight is too low, kill the walk with probability 1/2.
    elif new_weight < W_minus:
        if rng.random() < 0.5:
            return
        else:
            new_weight *= 2
            chosen = neighbors[rng.integers(len(neighbors))]
            new_visited = visited.copy()
            new_visited.add(chosen)
            perm_recursive(n + 1, L, new_weight, chosen, new_visited, rng, stats, c_plus, c_minus, results)
    # Otherwise, continue with normal growth.
    else:
        chosen = neighbors[rng.integers(len(neighbors))]
        new_visited = visited.copy()
        new_visited.add(chosen)
        perm_recursive(n + 1, L, new_weight, chosen, new_visited, rng, stats, c_plus, c_minus, results)


def perm_estimate_cL(L, iterations=1000, c_minus=0.2, c_plus=3.0, seed=42):
    """
    Estimate the number of SAWs (c_L) for a given length L using the PERM algorithm.
    
    Parameters
    ----------
    L : int
        The target walk length (number of steps).
    iterations : int
        Number of independent PERM runs (each run can produce one or more complete walks).
    c_minus : float
        Pruning threshold multiplier.
    c_plus : float
        Enrichment threshold multiplier.
    seed : int or None
        Seed for reproducibility (optional).
    
    Returns
    -------
    float
        Estimate of c_L (the total weight of complete walks averaged over the iterations).
    """
    sys.setrecursionlimit(max(10_000, L*2))

    rng = np.random.default_rng(seed)
    results = []
    # Initialize statistics for each step (from 0 to L).
    stats = {
        "weight_sum": [0.0] * (L + 1),
        "tours": [0] * (L + 1)
    }
    
    # Run the PERM growth process several times.
    for _ in range(iterations):
        perm_recursive(
            n=0,
            L=L,
            weight=1.0,
            pos=(0, 0),
            visited={(0, 0)},
            rng=rng,
            stats=stats,
            c_plus=c_plus,
            c_minus=c_minus,
            results=results
        )
    
    # The estimator for c_L is the average weight of completed walks per iteration.
    if len(results) == 0:
        return 0.0
    return sum(results) / iterations


if __name__ == "__main__":
    # Example usage
    L = 200
    iterations = 200000  
    cL_est = perm_estimate_cL(L, iterations=iterations, c_minus=0.2, c_plus=3.0, seed=42)
    print(f"PERM estimate for c_{L} ≈ {cL_est}")
    print(f"PERM estimate for mu ≈ {cL_est**(1/L)}")