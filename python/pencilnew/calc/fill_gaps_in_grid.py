def fill_gaps_in_grid(array, key='nan', steps=['interpolate', 'extrapolate'], order=1, DEBUG=False):
    """Interpolates nans by nearest neighbor method to fill gaps in arrays.
    Beware this method does not invoke bondaries correctly!
    This method does not work with array beeing a vector field! Hence, use componentwise
    and stitch together manually!

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

    def get_columns_with_cell_idx(matrix, cell_idx):
        slicing = [[e if i!=j else slice(None) for j,e in enumerate(cell_idx)] for i in range(np.size(cell_idx))]
        return [matrix[c] for c in slicing]

    ######## type cast array into np.array
    array = np.array(array)

    ######## Two iteration steps: first interpolate, then extrapolate
    for step in steps:
        ######## get indexes where array entries match key
        if key == 'nan' or np.isnan(key):
            idx_bad = np.isnan(array)
            coords_bad = np.argwhere(idx_bad)

            idx_good = ~np.isnan(array)
            coords_good = np.argwhere(idx_good)

        else:
            idx_bad = np.argwhere(array == key)
            idx_good = array != key

        ####### Step 1: Interpolate
        if step == 'interpolate':
            from scipy.interpolate import LinearNDInterpolator

            f = LinearNDInterpolator(coords_good, array[idx_good])
            for x_bad in coords_bad:
                y_bad = f(*x_bad)
                if DEBUG: print('Interpolating '+str(x_bad)+' entry to '+str(y_bad))
                array.__setitem__(tuple(x_bad), y_bad)

        ####### Step 2: Extrapolate
        elif step == 'extrapolate':
            from scipy.interpolate import InterpolatedUnivariateSpline

            for x_bad in coords_bad:
                ## get columns in all axis having a certain nan position inside
                columns = get_columns_with_cell_idx(array, x_bad)

                ## extrapolate at x_bad for each column and store in val
                vals = []
                for pos, column in zip(x_bad, columns):
                    xi=[]; yi=[]

                    ## clean up all additional locations with key or nan inside of column
                    for x,y in enumerate(column):
                        if not ( (y == key) or (np.isnan(y)) ):
                            xi.append(x); yi.append(y)

                    s = InterpolatedUnivariateSpline(xi, yi, k=order)
                    vals.append(s(pos))

                y_bad = np.mean(vals)
                if DEBUG: print('Extrapolating '+str(x_bad)+' entry to '+str(y_bad))
                array.__setitem__(tuple(x_bad), y_bad)      ## replace in array with mean extrapolated value

        else:
            print('! ERROR: Could not identify operation step "'+step+'"!')

    if key in array: print('! WARNING: Still '+str(key)+' found inside interpolated array! Something went wrong...')

    return array
