# simple_centered.py
#
"""
Contains simple centered differential operators.
"""


def simple_centered(x, y):
    """
    Calculate dy by center differencing using array slices
    """

    import numpy as np

    try:
        # We know it will be this size.
        dy = np.zeros(y.shape, np.float)
        dy[1:-1] = (y[2:] - y[0:-2]) / (x[2:] - x[0:-2])
        dy[0] = (y[1] - y[0]) / (x[1] - x[0])
        dy[-1] = (y[-1] - y[-2]) / (x[-1] - x[-2])
    except:
        print("!! ERROR: Something is wrong with your input:")
        print("!!        Shape x: " + str(x.shape))
        print("!!        Shape y: " + str(y.shape))

    return dy
