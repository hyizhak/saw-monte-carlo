"""
Pivot MCMC method for SAWs.
"""

import numpy as np
from saw_monte_carlo.utils import SYM_FUNCTIONS

def generate_initial_saw(n):
    """
    Generate a trivial n-step SAW (straight line along x-axis).
    """
    return [(i, 0) for i in range(n + 1)]


def walk_to_set(walk):
    """
    Return a set of coordinates in the given walk for fast occupancy checks.
    """
    return set(walk)


def attempt_pivot_optimized(walk, occupied_set, rng=None):
    """
    Perform a pivot move with partial subwalk updates.

    Parameters
    ----------
    walk : list of (int, int)
        The current walk coordinates.
    occupied_set : set of (int, int)
        Occupied positions for quick membership checks.
    rng : np.random.Generator, optional
        Random number generator for reproducibility.

    Returns
    -------
    (new_walk, new_occupied_set, accepted)
        new_walk : list of (int, int)
            The updated walk if accepted, or the old walk if rejected.
        new_occupied_set : set of (int, int)
            Updated occupied set if accepted, else old set.
        accepted : bool
            Whether the pivot was accepted.
    """
    if rng is None:
        rng = np.random.default_rng(42)

    n = len(walk) - 1
    p = rng.integers(n + 1)
    sym = SYM_FUNCTIONS[rng.integers(len(SYM_FUNCTIONS))]

    if p == 0 or p == n:
        # Pivoting around an endpoint typically doesn't change anything
        return walk, occupied_set, False

    pivot_x, pivot_y = walk[p]
    left_size = p
    right_size = n - p

    # Prepare new copies
    new_walk = list(walk)
    new_occ = set(occupied_set)

    def transform(px, py, x, y):
        rx, ry = x - px, y - py
        rx_t, ry_t = sym(rx, ry)
        return (rx_t + px, ry_t + py)

    if left_size < right_size:
        # Transform the subwalk 0..(p-1)
        for i in range(p):
            old_coord = new_walk[i]
            new_occ.remove(old_coord)
        for i in range(p):
            ox, oy = new_walk[i]
            new_walk[i] = transform(pivot_x, pivot_y, ox, oy)
        for i in range(p):
            if new_walk[i] in new_occ:
                return walk, occupied_set, False
            new_occ.add(new_walk[i])
    else:
        # Transform the subwalk (p+1)..n
        for i in range(p + 1, n + 1):
            old_coord = new_walk[i]
            new_occ.remove(old_coord)
        for i in range(p + 1, n + 1):
            ox, oy = new_walk[i]
            new_walk[i] = transform(pivot_x, pivot_y, ox, oy)
        for i in range(p + 1, n + 1):
            if new_walk[i] in new_occ:
                return walk, occupied_set, False
            new_occ.add(new_walk[i])

    return new_walk, new_occ, True


def count_free_forward_moves(walk):
    """
    For a 2D walk, check how many of the 3 'forward directions' are free
    (excluding the immediate back-step direction).
    """
    n = len(walk) - 1
    x_n, y_n = walk[-1]
    x_nm1, y_nm1 = walk[-2]

    dx = x_n - x_nm1
    dy = y_n - y_nm1
    all_moves = [(1, 0), (-1, 0), (0, 1), (0, -1)]
    back_move = (-dx, -dy)  # move that goes back onto the old site
    forward_moves = [mv for mv in all_moves if mv != back_move]

    occupied = set(walk)
    free_count = 0
    for fx, fy in forward_moves:
        new_site = (x_n + fx, y_n + fy)
        if new_site not in occupied:
            free_count += 1
    return free_count


def run_pivot_get_mu_estimate(n=100, pivot_attempts=20000, burn_in=2000, seed=42):
    """
    Runs pivot MCMC on an n-step SAW and estimates mu by measuring 
    the average number of free forward moves at the chain end.

    Parameters
    ----------
    n : int
        Number of steps in the SAW.
    pivot_attempts : int
        Total number of pivot attempts to run.
    burn_in : int
        Number of initial pivot attempts to discard for equilibration.
    seed : int, optional
        Seed for the random number generator.

    Returns
    -------
    float
        Estimate of mu based on the average forward move count.
    """
    rng = np.random.default_rng(seed)
    walk = generate_initial_saw(n)
    occupied_set = walk_to_set(walk)

    accepted = 0
    sum_free_moves = 0
    samples = 0

    for step in range(pivot_attempts):
        new_walk, new_occ, ok = attempt_pivot_optimized(walk, occupied_set, rng=rng)
        if ok:
            walk, occupied_set = new_walk, new_occ
            accepted += 1

        if step >= burn_in:
            # measure extension after burn_in
            free_count = count_free_forward_moves(walk)
            sum_free_moves += free_count
            samples += 1

    avg_free = sum_free_moves / samples if samples > 0 else 0.0
    mu_est = avg_free

    # print(f"Pivot MCMC for n={n}")
    # print(f" - Acceptance rate: {accepted / pivot_attempts:.3f}")
    # print(f" - mu estimate ~ {mu_est:.6f}")
    return mu_est


if __name__ == "__main__":
    # Example usage
    for test_n in [50, 100, 200, 400, 800, 1600]:
        mu_est = run_pivot_get_mu_estimate(n=test_n, pivot_attempts=100000, burn_in=10000)
        print(f'Pivot mu estimate for n={test_n}: {mu_est}')