"""
Basic SAW methods and Monte Carlo estimators.
"""

import numpy as np

def count_saws(L, pos=(0, 0), visited=None):
    """
    Recursively count the number of self-avoiding walks (SAWs) of length L
    on a 2D square lattice, starting from the origin.

    Parameters
    ----------
    L : int
        Remaining number of steps to take in the walk.
    pos : tuple of int
        Current (x, y) position on the lattice.
    visited : set of (int, int)
        Set of positions already visited.

    Returns
    -------
    int
        Number of SAWs of length L starting from 'pos'.
    """
    if visited is None:
        visited = {(0, 0)}

    if L == 0:
        return 1

    count = 0
    # Possible moves on the 2D square lattice
    moves = [(1, 0), (-1, 0), (0, 1), (0, -1)]
    for dx, dy in moves:
        new_pos = (pos[0] + dx, pos[1] + dy)
        if new_pos not in visited:
            visited.add(new_pos)
            count += count_saws(L - 1, new_pos, visited)
            visited.remove(new_pos)
    return count


def generate_rw(L, rng=None):
    """
    Generate a random walk (not necessarily self-avoiding) of length L 
    on a 2D square lattice.

    Parameters
    ----------
    L : int
        The desired length of the random walk.

    Returns
    -------
    list of (int, int)
        A list of coordinates representing the walk.
    """
    if rng is None:
        rng = np.random.default_rng(42)

    walk = [(0, 0)]
    moves = [(1, 0), (-1, 0), (0, 1), (0, -1)]
    for _ in range(L):
        x, y = walk[-1]
        dx, dy = moves[rng.integers(len(moves))]
        walk.append((x + dx, y + dy))
    return walk


def estimate_cL(L, N=100000, seed=None):
    """
    Estimate the average number of self-avoiding walks (SAWs) of length L
    on a 2D square lattice using naive Monte Carlo.

    Parameters
    ----------
    L : int
        The length of the walks to consider.
    N : int
        The number of Monte Carlo trials.

    Returns
    -------
    float
        The estimated average number of SAWs of length L.
    """
    rng = np.random.default_rng(seed)
    
    count_saw = 0
    for _ in range(N):
        walk = generate_rw(L, rng=rng)
        # Check if all positions in walk are unique
        if len(walk) == len(set(walk)):
            count_saw += 1
    return count_saw / N * (4 ** L)
    

if __name__ == "__main__":
    L = 10
    print(f"Number of SAWs of length {L}: {count_saws(L)}")
    N = 100000
    est = estimate_cL(L, N, seed=42)
    print(f"Estimated number of SAWs of length {L} with {N} trials: {est}")
    print(f"Estimated mu for L={L}: {est**(1 / L)}")