.. _IDL_to_python:

IDL to Python Translation Guide for Pencil Code
=================================================

Overview
--------

This comprehensive guide explains how to translate analysis scripts from IDL to Python using the Pencil Code analysis libraries. Both IDL (in ``pencil-code/idl/``) and Python (in ``pencil-code/python/pencil``) libraries provide access to simulation data, but with different syntax and capabilities.

**Key Concepts**:

- The Pencil Code Python library is the modern, recommended approach for new analysis scripts
- Both languages can read the same simulation output formats (HDF5, NetCDF, binary)
- Python offers superior data analysis, visualization, and integration with scientific computing ecosystem
- IDL remains useful for legacy code and specific visualization needs
- Simulation object caching in Python provides performance advantages over direct function calls

**When to use each language**:

- **Python**: New scripts, data analysis, numerical computing, interactive Jupyter notebooks
- **IDL**: Legacy code maintenance, existing IDL-based workflows, specific IDL library functions


Architecture Comparison
----------------------------

+------------------------+----------------------------------+------------------------------------+
| Feature                | IDL                              | Python                             |
+========================+==================================+====================================+
| Library location       | ``pencil-code/idl/``             | ``pencil-code/python/pencil``      |
+------------------------+----------------------------------+------------------------------------+
| Main namespace         | ``pc_*`` procedures              | ``pencil`` module                  |
+------------------------+----------------------------------+------------------------------------+
| Data structures        | IDL objects (``obj``)            | Python objects, dictionaries       |
+------------------------+----------------------------------+------------------------------------+
| Array operations       | Built-in IDL functions           | NumPy library                      |
+------------------------+----------------------------------+------------------------------------+
| Plotting               | IDL graphics, custom routines    | Matplotlib, Jupyter integration    |
+------------------------+----------------------------------+------------------------------------+
| Caching                | Manual (file writes)             | Automatic via ``simulation`` object|
+------------------------+----------------------------------+------------------------------------+
| Parameter access       | Direct attributes                | Dictionary with ``.get()``         |
+------------------------+----------------------------------+------------------------------------+

Reading Simulation Data
--------------------------

Option 1: Direct Functions (Original Approach)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Both languages offer direct function calls to read data:

**IDL**::

    ; IDL - Direct function approach
    pc_read_param, obj=par, datadir='data/'
    pc_read_dim, obj=dim, datadir='data/'
    pc_read_ts, obj=ts, datadir='data/'
    pc_read_var, obj=var, datadir='data/'
    pc_read_xyaver, obj=aver, datadir='data/'
    pc_read_slices, obj=slices, datadir='data/'

**Python**::

    # Python - Direct function approach
    import pencil as pc
    
    par = pc.read.param(datadir='data/')
    dim = pc.read.dim(datadir='data/')
    ts = pc.read.ts(datadir='data/')
    var = pc.read.var(datadir='data/')
    aver = pc.read.aver(datadir='data/')
    slices = pc.read.slices(datadir='data/')

**Advantages**:
  - Simple, direct function calls
  - No object initialization needed
  - Works with any data directory path

**Disadvantages**:
  - Must specify datadir for each call
  - Reads from disk every time (slower for repeated access)
  - Data is not cached between calls
  - More verbose code


Option 2: Simulation Object (Python Recommended)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Python's modern recommended approach:

**Python**::

    import pencil as pc
    
    # Initialize simulation object (caches data on first access)
    sim = pc.sim.simulation('data/')
    
    # Access data (cached automatically)
    par = sim.param           # Parameters dictionary
    dim = sim.dim             # Dimensions object
    ts = sim.ts               # Time series data
    var = sim.get_var(step=1) # Variable data at specific step
    aver = sim.get_aver()     # XY-averaged data
    slices = sim.get_slices() # Slice data

**Advantages**:
  - Data is cached (faster repeated access)
  - Cleaner, more intuitive code
  - Automatic path resolution
  - Object-oriented design
  - Single initialization point

**Disadvantages**:
  - Slight memory overhead for caching
  - Need to know simulation directory path

**Recommendation**: Use simulation object (Option 2) for all new Python scripts.


Core Data Access Functions
-----------------------------

Parameters
^^^^^^^^^^^^^^^

**IDL**::

    pc_read_param, obj=par, datadir='data/'
    cp = par.cp
    gamma = par.gamma
    rho0 = par.rho0

**Python (Direct)**::

    par = pc.read.param(datadir='data/')
    cp = par.get('cp', 1.0)
    gamma = par.get('gamma', 1.667)
    rho0 = par.get('rho0', 1.0)

**Python (Object)**::

    sim = pc.sim.simulation('data/')
    par = sim.param
    cp = par.get('cp', 1.0)
    gamma = par.get('gamma', 1.667)
    rho0 = par.get('rho0', 1.0)

**Key Differences**:
  - Python parameters are stored in a dictionary
  - Must use ``.get()`` method with optional defaults
  - Parameter names are strings (case-sensitive)
  - No direct attribute access like IDL's ``par.cp``


Dimensions and Grid
^^^^^^^^^^^^^^^^^^^^^^

**IDL**::

    pc_read_dim, obj=dim, datadir='data/'
    nx = dim.nx
    ny = dim.ny
    nz = dim.nz
    x = dim.x      ; 1D coordinate array
    y = dim.y
    z = dim.z

**Python (Direct)**::

    dim = pc.read.dim(datadir='data/')
    nx = dim.nx
    ny = dim.ny
    nz = dim.nz
    x = dim.x
    y = dim.y
    z = dim.z

**Python (Object)**::

    sim = pc.sim.simulation('data/')
    dim = sim.dim
    nx = getattr(dim, 'nx', None)  # Safe access with default
    ny = getattr(dim, 'ny', None)
    nz = getattr(dim, 'nz', None)
    x = dim.x
    y = dim.y
    z = dim.z


Time Series
^^^^^^^^^^^^^^

**IDL**::

    pc_read_ts, obj=ts, datadir='data/'
    time = ts.t
    urms = ts.urms
    umax = ts.umax

**Python**::

    # Direct
    ts = pc.read.ts(datadir='data/')
    
    # Object
    sim = pc.sim.simulation('data/')
    ts = sim.ts
    
    # Access data
    time = ts.t
    urms = ts.urms
    umax = ts.umax


Variable Data (Full 3D)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**IDL**::

    pc_read_var, obj=var, datadir='data/', ivar=1
    uu = var.uu       ; Shape: (nx, ny, nz, 3) or similar
    aa = var.aa
    rho = var.rho

**Python**::

    # Direct
    var = pc.read.var(datadir='data/', varfile='var.dat')
    
    # Object
    sim = pc.sim.simulation('data/')
    var = sim.get_var(varint=1)
    
    # Access data (returns dictionaries)
    uu = var.uu       ; Velocity field
    aa = var.aa       ; Magnetic potential
    rho = var.rho     ; Density


Averaged Data (XY-averaged, YZ-averaged, etc.)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**IDL**::

    pc_read_xyaver, obj=aver, datadir='data/'
    z = aver.z
    frad = aver.fradz
    fkin = aver.fkinz
    time = aver.t

**Python**::

    # Direct
    aver = pc.read.aver(datadir='data/', plane='xy')
    
    # Object
    sim = pc.sim.simulation('data/')
    aver = sim.get_aver(plane='xy')
    
    # Access data
    z = aver.z
    frad = aver.fradz
    fkin = aver.fkinz
    time = aver.t


Slice Data (2D cross-sections)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**IDL**::

    pc_read_slices, obj=slices, datadir='data/'
    uu_xy = slices.uu_xy    ; Velocity at xy plane
    tt_z = slices.tt_z      ; Temperature at z plane
    bb_xz = slices.bb_xz    ; Magnetic field at xz plane

**Python**::

    # Direct
    slices = pc.read.slices(datadir='data/')
    
    # Object
    sim = pc.sim.simulation('data/')
    slices = sim.get_slices()
    
    # Access data
    uu_xy = slices.uu_xy
    tt_z = slices.tt_z
    bb_xz = slices.bb_xz


Grid Information
^^^^^^^^^^^^^^^^^^^^^

**IDL**::

    pc_read_grid, obj=grid, datadir='data/'
    x = grid.x
    y = grid.y
    z = grid.z

**Python**::

    # Direct
    grid = pc.read.grid(datadir='data/')
    x = grid.x
    y = grid.y
    z = grid.z
    
    # Object
    sim = pc.sim.simulation('data/')
    grid = sim.get_grid()


IDL to Python Array Operations
---------------------------------

Basic Operations
^^^^^^^^^^^^^^^^^^^^^

+----------------------------------+--------------------------------------+
| IDL                              | Python                               |
+==================================+======================================+
| ``where(arr > 0)``               | ``np.where(arr > 0)[0]``             |
+----------------------------------+--------------------------------------+
| Boolean indexing: ``arr[ind]``   | ``arr[mask]`` (direct)               |
+----------------------------------+--------------------------------------+
| ``arr[*]`` (all elements)        | ``arr`` or ``arr[:]``                |
+----------------------------------+--------------------------------------+
| ``arr[0:10]`` (slice)            | ``arr[0:10]`` (same)                 |
+----------------------------------+--------------------------------------+
| Size: ``size(arr)``              | ``len(arr)`` or ``arr.size``         |
+----------------------------------+--------------------------------------+
| Shape: ``size(arr, /dim)``       | ``arr.shape``                        |
+----------------------------------+--------------------------------------+
| Transpose: ``transpose(arr)``    | ``arr.T`` or ``np.transpose(arr)``   |
+----------------------------------+--------------------------------------+
| Reshape: ``reform(arr, nx, ny)`` | ``arr.reshape(nx, ny)``              |
+----------------------------------+--------------------------------------+

Mathematical Functions
^^^^^^^^^^^^^^^^^^^^^^^^

+----------------------------------+--------------------------------------+
| IDL                              | Python                               |
+==================================+======================================+
| ``alog(x)``                      | ``np.log(x)``                        |
+----------------------------------+--------------------------------------+
| ``alog10(x)``                    | ``np.log10(x)``                      |
+----------------------------------+--------------------------------------+
| ``exp(x)``                       | ``np.exp(x)``                        |
+----------------------------------+--------------------------------------+
| ``sqrt(x)``                      | ``np.sqrt(x)``                       |
+----------------------------------+--------------------------------------+
| ``abs(x)``                       | ``np.abs(x)``                        |
+----------------------------------+--------------------------------------+
| ``sin(x)``                       | ``np.sin(x)``                        |
+----------------------------------+--------------------------------------+
| ``cos(x)``                       | ``np.cos(x)``                        |
+----------------------------------+--------------------------------------+
| ``tan(x)``                       | ``np.tan(x)``                        |
+----------------------------------+--------------------------------------+
| ``min(arr)``                     | ``np.min(arr)``                      |
+----------------------------------+--------------------------------------+
| ``max(arr)``                     | ``np.max(arr)``                      |
+----------------------------------+--------------------------------------+
| ``total(arr)``                   | ``np.sum(arr)``                      |
+----------------------------------+--------------------------------------+
| ``mean(arr)``                    | ``np.mean(arr)``                     |
+----------------------------------+--------------------------------------+
| ``stddev(arr)``                  | ``np.std(arr)``                      |
+----------------------------------+--------------------------------------+
| ``median(arr)``                  | ``np.median(arr)``                   |
+----------------------------------+--------------------------------------+
| ``fix(x)`` (truncate)            | ``np.trunc(x)``                      |
+----------------------------------+--------------------------------------+
| ``floor(x)``                     | ``np.floor(x)``                      |
+----------------------------------+--------------------------------------+
| ``ceil(x)``                      | ``np.ceil(x)``                       |
+----------------------------------+--------------------------------------+
| ``round(x)``                     | ``np.round(x)``                      |
+----------------------------------+--------------------------------------+

Array Construction
^^^^^^^^^^^^^^^^^^^^^^

+----------------------------------+--------------------------------------+
| IDL                              | Python                               |
+==================================+======================================+
| ``fltarr(10)``                   | ``np.zeros(10)``                     |
+----------------------------------+--------------------------------------+
| ``intarr(10)``                   | ``np.zeros(10, dtype=int)``          |
+----------------------------------+--------------------------------------+
| ``fltarr(10,20)``                | ``np.zeros((10, 20))``               |
+----------------------------------+--------------------------------------+
| ``findgen(10)`` (0 to 9)         | ``np.arange(10)``                    |
+----------------------------------+--------------------------------------+
| ``indgen(10)`` (integers 0-9)    | ``np.arange(10)``                    |
+----------------------------------+--------------------------------------+
| ``linspace(0, 1, 100)``          | ``np.linspace(0, 1, 100)``           |
+----------------------------------+--------------------------------------+
| ``replicate(val, 10)``           | ``np.full(10, val)``                 |
+----------------------------------+--------------------------------------+
| ``[arr1, arr2]`` (concatenate)   | ``np.concatenate([arr1, arr2])``     |
+----------------------------------+--------------------------------------+

Reduction Operations
^^^^^^^^^^^^^^^^^^^^^^^^^

+----------------------------------+--------------------------------------+
| IDL                              | Python                               |
+==================================+======================================+
| ``total(arr, /cumulative)``      | ``np.cumsum(arr)``                   |
+----------------------------------+--------------------------------------+
| ``min(arr, max=m)``              | ``m = np.max(arr)``                  |
+----------------------------------+--------------------------------------+
| Sum along dimension: ``total()`` | ``np.sum(arr, axis=2)``              |
+----------------------------------+--------------------------------------+

Integration
-----------

**IDL** (custom function or library):

.. code:: 

    ; IDL - Simple integration
    result = integr(data, x=z)  ; Custom function


**Python**

.. code:: python

    import numpy as np
    
    # Trapezoidal rule
    result = np.trapz(data, z)
    
    # Simpson's rule (more accurate)
    from scipy.integrate import simps
    result = simps(data, z)

**Example: Integrate cooling flux vertically**::

.. code:: python

    # Python
    import numpy as np
    import pencil as pc
    
    sim = pc.sim.simulation('data/')
    aver = sim.get_aver()
    
    # Integrate fcoolz over all heights
    fcool_total = np.trapz(aver.fcoolz, aver.z)
    
    # Integrate only in unstable layer
    mask = aver.z <= 0.5
    fcool_unstable = np.trapz(aver.fcoolz[mask], aver.z[mask])


Boolean Masking and Conditional Selection
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**IDL**

.. code::

    ; Find where z >= z2
    ind = where(z >= z2, count)
    if count gt 0 then begin
        subset = data[ind]
    endif

**Python**

.. code:: python

    # Create boolean mask
    mask = z >= z2
    subset = data[mask]
    
    # Or inline
    subset = data[z >= z2]
    
    # Check if any elements match
    if np.any(mask):
        subset = data[mask]
    
    # Count matching elements
    count = np.sum(mask)


Working with Parameters in Depth
-----------------------------------

Common Parameter Access Patterns
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Single parameter with default**

.. code:: python

    # Python
    rho0 = par.get('rho0', 1.0)

**Multiple parameters safely**

.. code:: python

    # Python
    params = {
        'cp': par.get('cp', 1.0),
        'gamma': par.get('gamma', 1.667),
        'gravity': par.get('gravz', -1.0),
        'viscosity': par.get('nu', 1.0e-3)
    }

**Check if parameter exists**

.. code:: python


    if 'reference_state' in par:
        ref_state = par['reference_state']
    else:
        ref_state = 'polytropic'

**Iterate over all parameters**

.. code:: python

    # Print all parameters
    for key in sorted(par.keys()):
        value = par[key]
        print(f"{key}: {value}")

**Get all parameters matching a pattern**

.. code:: python

    # All flux-related parameters
    flux_pars = {k: v for k, v in par.items() if 'flux' in k}
    
    # All entropy-related parameters
    entropy_pars = {k: v for k, v in par.items() if 'entropy' in k.lower() or 's' == k}


Advanced Data Access Techniques
---------------------------------

Reading Multiple Variable Snapshots
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**IDL**

.. code:: 

    ; Read multiple snapshots
    for ivar = 1, nvar do begin
        pc_read_var, obj=var, ivar=ivar, datadir='data/'
        ; Process var...
    endfor

**Python**

.. code:: python

    import pencil as pc
    
    # Get list of available variable files
    sim = pc.sim.simulation('data/')
    
    # Read multiple snapshots
    for step in range(1, 10):
        var = sim.get_var(varint=step)
        # Process var...
        uu = var.uu
        rho = var.rho


Time Series Analysis
^^^^^^^^^^^^^^^^^^^^^^^

**IDL**

.. code::

    pc_read_ts, obj=ts, datadir='data/'
    
    ; Find maximum velocity
    imax = (where(ts.urms eq max(ts.urms)))[0]
    t_max = ts.t[imax]

**Python**

.. code:: python

    import pencil as pc
    import numpy as np
    
    sim = pc.sim.simulation('data/')
    ts = sim.ts
    
    # Find maximum velocity
    imax = np.argmax(ts.urms)
    t_max = ts.t[imax]
    urms_max = ts.urms[imax]


Vertical Profile Analysis
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**IDL**

.. code::

    ; Compute profile quantities
    for i=0, nz-1 do begin
        profile[i] = mean(field[*, *, i])
    endfor

**Python**

.. code:: python

    import numpy as np
    
    # Python - Much simpler
    profile = np.mean(field, axis=(0, 1))  # Average over x,y
    
    # Or compute at specific heights
    z_indices = [0, nz//4, nz//2, 3*nz//4, nz-1]
    profiles = field[:, :, z_indices]


Plotting and Visualization
---------------------------------

Basic Plot Function Pattern
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**IDL**

.. code:: 

    ; IDL procedure
    pro plot_temperature, temp, z, title=title
        plot, z, temp, xtitle='z', ytitle='T', title=title
        if keyword_set(savefig) then begin
            device, /close
        endif
    end

**Python**

.. code:: python

    import matplotlib.pyplot as plt
    import os
    
    def plot_temperature(z, temp, title="", savefig=False):
        """Plot temperature profile.
        
        Parameters
        ----------
        z : array_like
            Height coordinates
        temp : array_like
            Temperature values
        title : str
            Plot title
        savefig : bool
            If True, save figure
        """
        fig, ax = plt.subplots(figsize=(10, 6))
        ax.plot(z, temp, 'b-', linewidth=2, label='Temperature')
        ax.set_xlabel('Height (z)')
        ax.set_ylabel('Temperature (T)')
        ax.set_title(title)
        ax.grid(True, alpha=0.3)
        ax.legend()
        
        if savefig:
            figdir = 'fig'
            if not os.path.isdir(figdir):
                os.makedirs(figdir)
            figpath = os.path.join(figdir, 'temperature.png')
            fig.savefig(figpath, dpi=150, bbox_inches='tight')
            print(f"Figure saved to {figpath}")
        
        return fig, ax


Multi-Panel Plotting
^^^^^^^^^^^^^^^^^^^^^

**Python**

.. code:: python

    import matplotlib.pyplot as plt
    import numpy as np
    
    def plot_three_profiles(z, temp, rho, entropy):
        """Create 3-panel plot of thermodynamic profiles."""
        
        fig, axes = plt.subplots(1, 3, figsize=(15, 5))
        
        # Temperature
        axes[0].plot(z, temp, 'r-', linewidth=2)
        axes[0].set_xlabel('Height (z)')
        axes[0].set_ylabel('Temperature (T)')
        axes[0].set_title('Temperature Profile')
        axes[0].grid(True, alpha=0.3)
        
        # Density
        axes[1].plot(z, rho, 'g-', linewidth=2)
        axes[1].set_xlabel('Height (z)')
        axes[1].set_ylabel('Density (ρ)')
        axes[1].set_title('Density Profile')
        axes[1].grid(True, alpha=0.3)
        
        # Entropy
        axes[2].plot(z, entropy, 'b-', linewidth=2)
        axes[2].set_xlabel('Height (z)')
        axes[2].set_ylabel('Entropy (s)')
        axes[2].set_title('Entropy Profile')
        axes[2].grid(True, alpha=0.3)
        
        plt.tight_layout()
        return fig, axes


Time Evolution Plots
^^^^^^^^^^^^^^^^^^^^^

**Python**

.. code:: python

    import matplotlib.pyplot as plt
    import pencil as pc
    
    sim = pc.sim.simulation('data/')
    ts = sim.ts
    
    # Plot energy evolution
    fig, axes = plt.subplots(2, 2, figsize=(12, 8))
    
    axes[0, 0].plot(ts.t, ts.urms)
    axes[0, 0].set_ylabel('RMS Velocity')
    axes[0, 0].grid(True)
    
    axes[0, 1].plot(ts.t, ts.umax)
    axes[0, 1].set_ylabel('Max Velocity')
    axes[0, 1].grid(True)
    
    axes[1, 0].plot(ts.t, ts.rhom)
    axes[1, 0].set_ylabel('Mean Density')
    axes[1, 0].grid(True)
    
    axes[1, 1].plot(ts.t, ts.ssm)
    axes[1, 1].set_ylabel('Mean Entropy')
    axes[1, 1].grid(True)
    
    axes[1, 0].set_xlabel('Time')
    axes[1, 1].set_xlabel('Time')
    
    plt.suptitle('Evolution of Thermodynamic Quantities')
    plt.tight_layout()
    plt.show()


Interactive Visualization
^^^^^^^^^^^^^^^^^^^^^^^^^^

**Python with Pencil Code tools**

.. code:: python

    import pencil as pc
    import matplotlib.pyplot as plt
    
    sim = pc.sim.simulation('data/')
    
    # Interactive slice viewer (modern approach)
    pc.visu.animate_interactive(
        sim.get_slices().uu_xy,
        sim.ts.t,
        title='Velocity Field (XY plane)',
        cmap='RdBu_r'
    )

**Advantages**:
  - Real-time slider control
  - Automatic scaling
  - Built-in colorbar
  - Better than writing custom animation code


Contour and Heatmap Plots
^^^^^^^^^^^^^^^^^^^^^^^^^^

**Python**

.. code:: python

    import matplotlib.pyplot as plt
    import numpy as np
    
    def plot_field_2d(field, x, z, title="", cmap='viridis'):
        """Plot 2D field as contour/heatmap."""
        
        fig, ax = plt.subplots(figsize=(10, 6))
        
        # Heatmap
        im = ax.contourf(x, z, field.T, levels=20, cmap=cmap)
        cbar = fig.colorbar(im, ax=ax)
        cbar.set_label('Value')
        
        # Contours
        cs = ax.contour(x, z, field.T, levels=10, colors='black', alpha=0.3, linewidths=0.5)
        ax.clabel(cs, inline=True, fontsize=8)
        
        ax.set_xlabel('x')
        ax.set_ylabel('z')
        ax.set_title(title)
        
        return fig, ax


Data I/O Patterns
------------------

Reading Text Data
^^^^^^^^^^^^^^^^^^

**IDL**

.. code:: 

    ; IDL
    data = read_ascii('data.txt')

**Python**

.. code:: python

    import numpy as np
    
    # Simple 1D data
    data = np.loadtxt('data.txt')
    
    # With headers/comments
    data = np.loadtxt('data.txt', comments='#')
    
    # Multiple columns into separate arrays
    x, y, z = np.loadtxt('data.txt', unpack=True)


Saving Results
^^^^^^^^^^^^^^

**Python**

.. code:: python

    import numpy as np
    
    # Save 1D array
    np.savetxt('result.txt', array_1d, fmt='%.6e')
    
    # Save 2D array
    np.savetxt('result.txt', array_2d, fmt='%.6e')
    
    # Save with header
    np.savetxt('result.txt', array, fmt='%.6e', 
               header='x  y  z')
    
    # CSV format
    np.savetxt('result.csv', array, fmt='%.6e', delimiter=',')


Working with HDF5 Data
^^^^^^^^^^^^^^^^^^^^^^^

**Python**

.. code:: python

    import h5py
    import numpy as np
    
    # Open HDF5 file
    with h5py.File('simulation.h5', 'r') as f:
        # List groups
        print(list(f.keys()))
        
        # Read dataset
        velocity = f['velocity'][:]
        
        # Read attributes
        metadata = dict(f.attrs)


Error Handling and Debugging
-----------------------------

Common Issues and Solutions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Problem: Missing Parameter**

.. code:: python

    KeyError: 'unknown_parameter'

**Solution**

.. code:: python

    # Use .get() with default
    value = par.get('unknown_parameter', 0.0)
    
    # Or check existence
    if 'unknown_parameter' in par:
        value = par['unknown_parameter']

**Problem: Dimension Mismatch**

.. code:: python

    ValueError: operands could not be broadcast together

**Solution**

.. code:: python

    # Check shapes
    print(f"Shape of var1: {var1.shape}")
    print(f"Shape of var2: {var2.shape}")
    
    # Explicit slicing
    var1_2d = var1[:, :, 0]
    var2_expanded = var2[:, :, np.newaxis]


**Problem: File Not Found**

.. code:: python

    FileNotFoundError: [...]/xyaverages.dat

**Solution**

.. code:: python

    import os
    
    sdir = 'data/'
    data_dir = os.path.join(sdir, 'data')
    
    if os.path.isdir(data_dir):
        print(f"Data directory found: {data_dir}")
    else:
        print(f"Available items in {sdir}:")
        for item in os.listdir(sdir):
            print(f"  {item}")


**Problem: Type Errors**

.. code:: python

    TypeError: unsupported operand type(s)

**Solution**

.. code:: python

    # Ensure correct data types
    arr_float = np.array(arr, dtype=float)
    result = arr_float + 1.0  # Now works


Inspecting Data Objects
^^^^^^^^^^^^^^^^^^^^^^^^

**List all parameters**

.. code:: python

    import pencil as pc
    
    sim = pc.sim.simulation('data/')
    par = sim.param
    
    print(f"Total parameters: {len(par)}")
    for key in sorted(par.keys()):
        print(f"  {key}: {par[key]}")

**Inspect dimension object**

.. code:: python

    dim = sim.dim
    
    print("Dimension attributes:")
    for attr in dir(dim):
        if not attr.startswith('_'):
            try:
                val = getattr(dim, attr)
                if not callable(val):
                    if isinstance(val, np.ndarray):
                        print(f"  {attr}: array shape {val.shape}")
                    else:
                        print(f"  {attr}: {val}")
            except:
                pass

**Inspect variable data**

.. code:: python

    var = sim.get_var(varint=1)
    
    print("Available variables:")
    for key in dir(var):
        if not key.startswith('_'):
            try:
                val = getattr(var, key)
                if isinstance(val, np.ndarray):
                    print(f"  {key}: shape {val.shape}, dtype {val.dtype}")
            except:
                pass


Best Practices
---------------

1. **Use Python Simulation Object in New Scripts**
   
   Always use ``sim = pc.sim.simulation(...)`` for new Python scripts.
   
   .. code:: python
   
       import pencil as pc
       sim = pc.sim.simulation('data/')
       par = sim.param

2. **Safe Parameter Access**
   
   Always use ``.get()`` with default values.
   
   .. code:: python
   
       value = par.get('parameter', default_value)

3. **Function-Based Organization**
   
   Encapsulate plotting and computation into reusable functions.
   
   .. code:: python
   
       def compute_profile(field, z_axis=0):
           """Compute vertical profile by averaging."""
           return np.mean(field, axis=tuple(i for i in range(3) if i != z_axis))
       
       def plot_profile(z, profile, title=""):
           """Plot vertical profile."""
           fig, ax = plt.subplots()
           ax.plot(z, profile)
           return fig, ax

4. **Directory Management**
   
   Always check and create output directories.
   
   .. code:: python
   
       import os
       
       figdir = 'figures'
       if not os.path.isdir(figdir):
           os.makedirs(figdir)
           print(f"Created directory {figdir}")

5. **Error Handling**
   
   Use try-except for I/O operations.
   
   .. code:: python
   
       import os
       
       try:
           data = np.loadtxt('data.txt')
       except FileNotFoundError:
           print("File not found. Creating dummy data...")
           data = np.zeros((100,))

6. **Documentation**
   
   Always include docstrings and comments.
   
   .. code:: python
   
       def analyze_convection(sim):
           """
           Analyze convection statistics from simulation.
           
           Parameters
           ----------
           sim : pencil.sim.simulation
               Simulation object with loaded data
           
           Returns
           -------
           stats : dict
               Dictionary with statistics
           """
           par = sim.param
           ts = sim.ts
           
           stats = {
               'mean_urms': np.mean(ts.urms),
               'max_urms': np.max(ts.urms),
               'Ra': compute_rayleigh_number(par)
           }
           return stats

7. **Testing and Validation**
   
   Write test code at module level.
   
   .. code:: python
   
       if __name__ == '__main__':
           # Test with example data
           sim = pc.sim.simulation('examples/conv-slab')
           result = main(sim)
           print(f"Analysis complete: {result}")


Complete Example: From IDL to Python
-------------------------------------

**Original IDL Script**

.. code::

    pro analyze_convection, datadir=datadir
        ; Read data
        pc_read_param, obj=par, datadir=datadir
        pc_read_ts, obj=ts, datadir=datadir
        pc_read_xyaver, obj=aver, datadir=datadir
        
        ; Compute diagnostics
        urms_mean = mean(ts.urms)
        urms_max = max(ts.urms)
        
        ; Plot
        plot, ts.t, ts.urms, xtitle='Time', ytitle='RMS Velocity'
        
        ; Save result
        printf, 1, 'Mean RMS velocity: ', urms_mean
    end

**Equivalent Python Script**

.. code:: python

    import numpy as np
    import matplotlib.pyplot as plt
    import pencil as pc
    import os
    
    def analyze_convection(datadir='data/'):
        """
        Analyze convection simulation.
        
        Parameters
        ----------
        datadir : str
            Path to simulation data directory
        """
        # Initialize
        sim = pc.sim.simulation(datadir)
        par = sim.param
        ts = sim.ts
        aver = sim.get_aver()
        
        # Compute diagnostics
        urms_mean = np.mean(ts.urms)
        urms_max = np.max(ts.urms)
        
        # Plot
        fig, ax = plt.subplots(figsize=(10, 6))
        ax.plot(ts.t, ts.urms, 'b-', linewidth=2)
        ax.set_xlabel('Time')
        ax.set_ylabel('RMS Velocity')
        ax.set_title('Convection Evolution')
        ax.grid(True, alpha=0.3)
        
        # Save figure
        figdir = 'figures'
        if not os.path.isdir(figdir):
            os.makedirs(figdir)
        figpath = os.path.join(figdir, 'convection.png')
        fig.savefig(figpath, dpi=150, bbox_inches='tight')
        
        # Print results
        print(f"Mean RMS velocity: {urms_mean:.6e}")
        print(f"Max RMS velocity:  {urms_max:.6e}")
        print(f"Figure saved to: {figpath}")
        
        return urms_mean, urms_max
    
    if __name__ == '__main__':
        analyze_convection('data/')


Summary of Key Differences
---------------------------

+------------------------+-------------------------------------+-------------------------------------+
| Aspect                 | IDL                                 | Python                              |
+========================+=====================================+=====================================+
| Parameters             | Object attributes (par.cp)          | Dictionary (.get('cp', default))    |
+------------------------+-------------------------------------+-------------------------------------+
| Arrays                 | Built-in operations                 | NumPy arrays                        |
+------------------------+-------------------------------------+-------------------------------------+
| Plotting               | IDL graphics / custom routines      | Matplotlib (more powerful)          |
+------------------------+-------------------------------------+-------------------------------------+
| Integration            | Custom functions or libraries       | SciPy (trapz, simps, etc.)          |
+------------------------+-------------------------------------+-------------------------------------+
| Data caching           | Manual (file writes)                | Automatic (simulation object)       |
+------------------------+-------------------------------------+-------------------------------------+
| Syntax                 | Procedural, IDL-specific            | Pythonic, scientific stack          |
+------------------------+-------------------------------------+-------------------------------------+
| Learning curve         | IDL knowledge required              | Python knowledge transferable       |
+------------------------+-------------------------------------+-------------------------------------+

---

**Last Updated**: January 2026

**Format**: ReStructuredText

**Status**: Complete with comprehensive examples



Acknowledgments
----------------

This tutorial has been written with the aid of GitHub Copilot.


