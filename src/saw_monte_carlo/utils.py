"""
Utilities for self-avoiding walk (SAW) simulations.
"""

def rotate90(x, y):
    """Rotate (x, y) by 90 degrees counterclockwise."""
    return -y, x

def rotate180(x, y):
    """Rotate (x, y) by 180 degrees."""
    return -x, -y

def rotate270(x, y):
    """Rotate (x, y) by 270 degrees counterclockwise."""
    return y, -x

def reflect_x(x, y):
    """Reflect (x, y) across the X-axis."""
    return (x, -y)

def reflect_y(x, y):
    """Reflect (x, y) across the Y-axis."""
    return (-x, y)

def reflect_diag(x, y):
    """Reflect (x, y) across the line y = x."""
    return (y, x)

def reflect_antidiag(x, y):
    """Reflect (x, y) across the line y = -x."""
    return (-y, -x)

SYM_FUNCTIONS = [
    lambda x, y: (x, y),
    rotate90,
    rotate180,
    rotate270,
    reflect_x,
    reflect_y,
    reflect_diag,
    reflect_antidiag
]