'''
Rosenbluth method for SAWs
'''
import numpy as np

def rosenbluth_single_walk(L, rng=None):
    """
    Perform a single Rosenbluth-guided self-avoiding walk (SAW) of length L on a 2D square lattice.
    Returns the Rosenbluth weight for this walk. If the walk gets stuck before length L, the weight is 0.
    
    Parameters
    ----------
    L : int
        The target length of the walk (number of steps).
    rng : np.random.Generator or None
        Optional random generator for reproducibility. If None, use np.random default.
    
    Returns
    -------
    float
        The Rosenbluth weight for the generated walk (0 if trapped before reaching L steps).
    """
    if rng is None:
        rng = np.random.default_rng(42)

    # Initialize
    x, y = 0, 0
    visited = {(x, y)}
    weight = 1.0

    for step in range(L):
        # Find all possible neighbor moves (up, down, left, right)
        neighbors = []
        for dx, dy in [(1,0), (-1,0), (0,1), (0,-1)]:
            nx, ny = x + dx, y + dy
            if (nx, ny) not in visited:
                neighbors.append((nx, ny))
        
        # If no available moves, walk is trapped
        if not neighbors:
            return 0.0

        # The Rosenbluth factor: multiply weight by number of allowed moves
        a = len(neighbors)
        weight *= a

        # Choose one neighbor uniformly at random
        chosen = neighbors[rng.integers(0, a)]
        x, y = chosen
        visited.add((x, y))

    return weight


def rosenbluth_estimate_cL(L, trials=100000, seed=None):
    """
    Estimate the number of SAWs (c_L) for length L using a vanilla Rosenbluth method.
    
    Parameters
    ----------
    L : int
        The target walk length.
    trials : int
        Number of independent walk attempts.
    seed : int or None
        Seed for reproducibility (optional).
    
    Returns
    -------
    float
        Estimate of c_L.
    """
    rng = np.random.default_rng(seed)
    total_weight = 0.0

    for _ in range(trials):
        w = rosenbluth_single_walk(L, rng)
        total_weight += w
    
    # Average of the weights is the unbiased estimator for c_L
    return total_weight / trials


if __name__ == "__main__":
    # Example usage
    L = 20
    trials = 200000
    cL_est = rosenbluth_estimate_cL(L, trials=trials, seed=42)
    print(f"Vanilla Rosenbluth estimate for c_{L} ≈ {cL_est}")
    print(f"Vanilla Rosenbluth estimate for mu ≈ {cL_est**(1/L)}")