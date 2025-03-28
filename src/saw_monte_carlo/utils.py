"""
Utilities for self-avoiding walk (SAW) simulations.
"""

EXACT_VALUES = {
    0: 1,
    1: 4,
    2: 12,
    3: 36,
    4: 100,
    5: 284,
    6: 780,
    7: 2172,
    8: 5916,
    9: 16268,
    10: 44100,
    11: 120292,
    12: 324932,
    13: 881500,
    14: 2374444,
    15: 6416596,
    16: 17245332,
    17: 46466676,
    18: 124658732,
    19: 335116620,
    20: 897697164,
    21: 2408806028,
    22: 6444560484,
    23: 17266613812,
    24: 46146397316,
    25: 123481354908,
    26: 329712786220,
    27: 881317491628,
    28: 2351378582244,
    29: 6279396229332,
    30: 16741957935348,
    31: 44673816630956,
    32: 119034997913020,
    33: 317406598267076,
    34: 845279074648708,
    35: 2252534077759844,
    36: 5995740499124412,
    37: 15968852281708724,
    38: 42486750758210044,
    39: 113101676587853932,
    40: 300798249248474268,
    41: 800381032599158340,
    42: 2127870238872271828,
    43: 5659667057165209612,
    44: 15041631638016155884,
    45: 39992704986620915140,
    46: 106255762193816523332,
    47: 282417882500511560972,
    48: 750139547395987948108,
    49: 1993185460468062845836,
    50: 5292794668724837206644,
    51: 14059415980606050644844,
    52: 37325046962536847970116,
    53: 99121668912462180162908,
    54: 263090298246050489804708,
    55: 698501700277581954674604,
    56: 1853589151789474253830500,
    57: 4920146075313000860596140,
    58: 13053884641516572778155044,
    59: 34642792634590824499672196,
    60: 91895836025056214634047716,
    61: 243828023293849420839513468,
    62: 646684752476890688940276172,
    63: 1715538780705298093042635884,
    64: 4549252727304405545665901684,
    65: 12066271136346725726547810652,
    66: 31992427160420423715150496804,
    67: 84841788997462209800131419244,
    68: 224916973773967421352838735684,
    69: 596373847126147985434982575724,
    70: 1580784678250571882017480243636,
    71: 4190893020903935054619120005916
}

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

def get_deviation(estimation, L_max):
    """
    Compute the deviation of the estimation from the exact value for c_L.
    
    Parameters
    ----------
    estimation : float
        The estimated value of c_L.
    L_max : int
        The length of the walks.
    
    Returns
    -------
    float
        The deviation of the estimation from the exact value.
    """
    reference = EXACT_VALUES.get(L_max, None)
    if reference is None:
        raise ValueError(f"No exact value available for c_{L_max}")
    return abs(estimation - reference) / reference