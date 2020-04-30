def part_to_grid(xp, yp, zp=False, quantity=False, Nbins=[1024,1024,1024], sim=False, extent=False, fill_gaps=False):
    """Bins quantity based on position data xp and yp to 1024^2 bins like a histrogram.
    This method is not using TSC.

    Args:
        - xp, yp:       array of x and y positions
        - zp:           specify if 3D run, set False to have zp == 0
        - quantity:     array of same shape as xp and yp, but with quantity to bin, set it False to count number of occurrences/histrogram2d
        - Nbins:        number of histrogram bins for each direction. if 2d only the first two entries in Nbins are used
        - sim:          to extract extent from by loading ghostzone-free grid
        - extent:       [[xmin, xmax],[ymin, ymax]] or set false and instead give a sim
                        set extent manually e.g. if you want to include ghost zones
        - fill_gaps     interpolate empty grid cells

    Returns: arr, xgrid, ygrid
        - arr:          2d array with binned values
        - x-/ygrid:     linspace of used x/y grid
        - zgrid:        if zp != False

    Example:
        vpx = part_to_grid_2d(pvar.xp, pvar.yp, pvar.vpx), notice that this will execute pc.get_sim() internally to get the extent
    """

    import numpy as np
    from .. import get_sim
    from ..calc import fill_gaps_in_grid

    if not xp.shape == yp.shape and yp.shape == xp.shape and yp.shape == zp.shape:
        print('! ERROR: Shape of xp, yp, zp and quantity needs to be equal!')

    if extent == False and sim == False:
        sim = get_sim()

    if quantity == False:
        quantity = xp/xp

    if extent == False:
        grid = sim.grid
        if type(zp) == type(False) and zp == False:
            extent = [[grid.x[0]-grid.dx/2, grid.x[-1]+grid.dx/2],
                      [grid.y[0]-grid.dy/2, grid.y[-1]+grid.dy/2]]
        else:
            extent = [[grid.x[0]-grid.dx/2, grid.x[-1]+grid.dx/2],
                      [grid.y[0]-grid.dy/2, grid.y[-1]+grid.dy/2],
                      [grid.z[0]-grid.dz/2, grid.z[-1]+grid.dz/2]]


    if type(zp) == type(False) and zp == False:
        arr = np.zeros((2, Nbins[0], Nbins[1]))
    else:
        arr = np.zeros((3, Nbins[0], Nbins[1], Nbins[2]))

    arr[:] = np.NAN
    xgrid = (np.linspace(extent[0][0], extent[0][1], num=Nbins[0]+1)[:-1]+np.linspace(extent[0][0], extent[0][1], num=Nbins[0]+1)[1:])/2
    ygrid = (np.linspace(extent[1][0], extent[1][1], num=Nbins[1]+1)[:-1]+np.linspace(extent[1][0], extent[1][1], num=Nbins[1]+1)[1:])/2
    if type(zp) == type(False) and zp == False:
        for x, y, q in zip(xp, yp, quantity):
            idx = np.argmin(np.abs(x-xgrid))
            idy = np.argmin(np.abs(y-ygrid))
            if np.isnan(arr[0, idx, idy]): arr[0, idx, idy] = 0; arr[1, idx, idy] = 0
            arr[0, idx, idy] += q
            arr[1, idx, idy] += 1

    else:
        zgrid = (np.linspace(extent[2][0], extent[2][1], num=Nbins[2]+1)[:-1]+np.linspace(extent[2][0], extent[2][1], num=Nbins[2]+1)[1:])/2
        for x, y, z, q in zip(xp, yp, zp, quantity):
            idx = np.argmin(np.abs(x-xgrid))
            idy = np.argmin(np.abs(y-ygrid))
            idz = np.argmin(np.abs(z-zgrid))
            if np.isnan(arr[0, idx, idy, idz]): arr[0, idx, idy, idz] = 0; arr[1, idx, idy, idz] = 0
            arr[0, idx, idy, idz] += q
            arr[1, idx, idy, idz] += 1


    arr = arr[0]/arr[1]

    if fill_gaps == True: arr = fill_gaps_in_grid(arr, key=np.NAN)

    if type(zp) == type(False) and zp == False:
        return arr, xgrid, ygrid
    else:
        return arr, xgrid, ygrid, zgrid
