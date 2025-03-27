"""
PERM (Pruned-Enriched Rosenbluth Method) for SAWs.
"""

import numpy as np

# Global accumulators for partial-walk statistics
global_w = []
global_s = []

def init_statistics(L_max):
    """
    Initialize global statistics arrays for lengths 0..L_max.
    """
    global global_w, global_s
    global_w = [0.0] * (L_max + 1)
    global_s = [0] * (L_max + 1)


def perm_grow(n, x, y, visited, weight, L_max):
    """
    Recursively grow a self-avoiding walk using PERM.

    Parameters
    ----------
    n : int
        Current length of the walk.
    x, y : int
        Current position on the lattice.
    visited : set
        Set of visited coordinates.
    weight : float
        Current weight for this partial walk.
    L_max : int
        Target maximum length of the walk.

    Returns
    -------
    float
        Contribution to the final c_{L_max} estimate from this branch.
    """
    global global_w, global_s

    # Record partial walk stats
    global_s[n] += 1
    global_w[n] += weight

    if n == L_max:
        # Reached the full length
        return weight  # This weight contributes to c_{L_max}

    # Determine allowed moves
    moves = [(1, 0), (-1, 0), (0, 1), (0, -1)]
    neighbors = []
    for dx, dy in moves:
        nx, ny = x + dx, y + dy
        if (nx, ny) not in visited:
            neighbors.append((nx, ny))
    if not neighbors:
        # Dead end
        return 0.0

    a = len(neighbors)  # number of allowed moves
    # Current average weight at length n (or fallback to current weight)
    if global_s[n] > 0:
        target_weight = global_w[n] / global_s[n]
    else:
        target_weight = weight

    # Compute ratio
    R = weight / target_weight
    copies = int(R)
    if np.random.rand() < (R - copies):
        copies += 1
    # Pruning if copies=0
    if copies == 0:
        return 0.0

    # Set new weight for continuing
    new_weight = target_weight
    total_weight_contrib = 0.0

    for _ in range(copies):
        # Choose a random neighbor
        nx, ny = neighbors[np.random.randint(len(neighbors))]
        visited.add((nx, ny))
        branch_contrib = perm_grow(n + 1, nx, ny, visited, new_weight * a, L_max)
        total_weight_contrib += branch_contrib
        visited.remove((nx, ny))

    return total_weight_contrib


if __name__ == "__main__":
    # Example usage
    L_max = 20
    init_statistics(L_max)
    trials = 100000  # number of independent PERM growth attempts
    sum_estimates = 0.0
    for _ in range(trials):
        visited = {(0, 0)}
        sum_estimates += perm_grow(0, 0, 0, visited, weight=1.0, L_max=L_max)

    c_L_est = sum_estimates / trials
    print(f"PERM estimate c_{L_max} ≈ {c_L_est}")