def part_to_grid_2d(xp, yp, quantity=False, Nbins=[1024,1024], sim=False, extent=False):
    """Bins quantity based on position data xp and yp to 1024^2 bins like a histrogram.
    This method is not using TSC.

    Args:
        - xp, yp:       array of x and y positions
        - quantity:     array of same shape as xp and yp, but with quanittiy to bin, set ti False to count number of occurrences/histrogram2d
        - Nbins:        number of histrogram bins for each direction
        - sim:          to extract extent from by loading ghostzone-free grid
        - extent:       [[xmin, xmax],[ymin, ymax]] or set false and instead give a sim
                        set extent manually e.g. if you want to include ghost zones

    Returns: arr, xgrid, ygrid
        - arr:          2d array with binned values
        - x-/ygrid:     linspace of used x/y grid

    Example:
        vpx = part_to_grid_2d(pvar.xp, pvar.yp, pvar.vpx), notice that this will execute pencilnew.get_sim() internally to get the extent
    """

    import numpy as np
    from pencilnew import get_sim

    if not xp.shape == yp.shape and yp.shape == quantity.shape:
        print('! ERROR: Shape of xp, yp and quantity needs to be equal!')

    if extent == False and sim == False:
        sim = get_sim()

    if extent == False:
        grid = sim.grid
        extent = [[grid.x[0]-grid.dx/2, grid.x[-1]+grid.dx/2],
                  [grid.y[0]-grid.dy/2, grid.y[-1]+grid.dy/2]]


    arr = np.zeros((2, Nbins[0], Nbins[1]))
    xgrid = np.linspace(extent[0][0], extent[0][1], num=Nbins[0])
    ygrid = np.linspace(extent[1][0], extent[1][1], num=Nbins[1])

    for x, y, q in zip(xp, yp, quantity):
        idx = np.argmin(np.abs(x-xgrid))
        idy = np.argmin(np.abs(y-ygrid))
        arr[0, idx, idy] += q
        arr[1, idx, idy] += 1

    arr[0] = arr[0]/arr[1]

    return arr[0], xgrid, ygrid
