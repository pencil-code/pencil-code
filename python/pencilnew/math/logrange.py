def logrange(a, b, include_ab=False, finer_factor=1):
    """ Gives back a log-space from a to b.
    	e.g.: logrange(4,40) = [4,5,6,7,8,9,10,20,30,40]
    """
    import math
    import numpy as np
    ten_range = np.arange(1,10,1./finer_factor)
    low_exp = int(math.log10(a))
    high_exp = int(math.log10(b))
    exp_range = range(low_exp,high_exp+1)

    ## check consistency
    if a == b: return a

    if include_ab:
        logrange = [a]
    else:
        logrange = []
    for i in exp_range:
        logrange = logrange + [10.0**i * j for j in ten_range]

    if include_ab: logrange = logrange + [b]

    ## cleanup to requested range
    seen = set()
    return np.array([c for c in logrange if c >= a and c <= b and c not in seen and not seen.add(c)])
