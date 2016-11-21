def interpolate_nans(array):
    """Interpolates nans by nearest neighbor method to fill gaps in arrays.
    Beware this method does not invoke bondaries correctly!

    Args:
        - array:        array with nans inside

    Returns:
        array without nans inside.
    """

    import numpy as np
    # from pencilnew.math import interpolate_nans
    from scipy.interpolate import interp1d
    import pencilnew

    array = np.array(array)
    shape = array.shape

    if len(shape) == 1:                             # do interpolation
        xrange = np.arange(len(array))

        interpos = []
        knownvalues = []
        for x, y in zip(xrange, array):
            if np.isnan(y):
                interpos.append(x)
            else:
                knownvalues.append(np.array([x,y]))

        knownvalues = np.array(knownvalues).T
        # pencilnew.io.debug_breakpoint()
        f = interp1d(knownvalues[0], knownvalues[1], kind='cubic')

        for x in interpos:
            if x < knownvalues[0][0]: array[x] = knownvalues[1][0]
            elif x > knownvalues[0][-1]: array[x] = knownvalues[1][-1]
            else: array[x] = f(x)

    else:
        for i, subarray in enumerate(array):        # go a level deeper
            array[i] = interpolate_nans(subarray)

    return array
