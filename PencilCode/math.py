#!/usr/bin/python3
# Last Modification: $Id$
#=======================================================================
# math.py
#
# Facilities for mathematical operations.
#
# Chao-Chin Yang, 2016-06-29
#=======================================================================
def _deriv1(a, axis=-1):
    """Takes the first derivative over the specified axis, assuming the
    grid spacing is unity.

    Positional Argument
        a
            Array to be operated on.

    Keyword Arguments
        axis
            Axis of a to be opearted on.
    """
    b = a.view() if axis == -1 else a.swapaxes(-1, axis)
    db = 0.75 * (b[...,4:-2] - b[...,2:-4]) - 0.15 * (b[...,5:-1] - b[...,1:-5]) + (1/60) * (b[...,6:] - b[...,:-6])
    return db if axis == -1 else db.swapaxes(axis,-1)
