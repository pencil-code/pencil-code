def fill_gaps_in_grid(array, key='nan', method='cubic'):
    """Interpolates nans by nearest neighbor method to fill gaps in arrays.
    Beware this method does not invoke bondaries correctly!

    Args:
        - array:        array with nans inside
        - key:          indicates array entries that will be interpolated
        - method:       as specified in scipy.interpolate.griddata,
                        i.e. use 'linear' for 3D data

    Returns:
        array without nans inside.

    Example:
        [[0,  1,  0],                 [[0, 1, 0,],
         [1, nan, 1],   -interpol->    [1, 2, 1],
         [0,  1,  0]]                  [0, 1, 0]]

    """


    import numpy as np

    ######## type cast array into np.array
    array = np.array(array)

    grid_x, grid_y, grid_z = np.mgrid[0:]


    ######## get indexes where array entries match key
    if key == 'nan' or np.isnan(key):
        idx_bad = np.isnan(array)
        coords_bad = np.argwhere(idx_bad)

        idx_good = ~np.isnan(array)
        coords_good = np.argwhere(idx_good)

    else:
        idx_bad = np.argwhere(array == key)
        idx_good = array != key


    ######## replace each of these entries from indx will interpolated values
    # A does not work :(
    # from scipy.interpolate import griddata
    # array[idx_bad] = griddata(points=coords_good, values=array[idx_good], xi=coords_bad, method=method)

    # B
    from scipy.interpolate import interpn
    newArray = interpn(points=coords_good, values=array[idx_good], xi=coords_bad)

    if key in newArray: print('! WARNING: Still '+str(key)+' found inside interpolated array!')

    return array
