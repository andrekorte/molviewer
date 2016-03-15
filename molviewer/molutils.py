# -*- coding: utf-8 -*-


# This function is used to find the symmetry of a molecule in
# function find_symmetry.
# For other comparisons use np.isclose.
def isclose(i, j, tolerance=0.3):
    '''
    Returns true if |i - j| is within tolerance

    '''
    return abs(i - j) <= tolerance
