IDL to Python Translation Guide
================================

Overview
--------

This guide explains how to translate analysis scripts from IDL to Python using the Pencil Code Python library. Two completed examples are provided: ``thermo.py`` and ``pc_flux.py``, translated from ``thermo.pro`` and ``pc_fluxz.pro``.

**Key Concepts**:

- Use ``pc.sim.simulation()`` to initialize simulation object (caches parameters)
- Access parameters via ``sim.param`` (dictionary) instead of ``pc.read.param()``
- Access dimensions via ``sim.dim`` (object) instead of ``pc.read.dim()``
- Use NumPy for array operations instead of IDL's built-in functions
- Use Matplotlib for visualization (often superior to IDL's plotting)
- Organize code into functions for reusability and testing


Reading Simulation Data
======================

Option 1: Direct Functions (Original Approach)
----------------------------------------------

This is the traditional method and still works::

    import pencil as pc
    
    # Read parameters
    par = pc.read.param()  # Returns: dict
    
    # Read dimensions  
    dim = pc.read.dim()    # Returns: object
    
    # Read data
    ts = pc.read.ts()
    var = pc.read.var()
    aver = pc.read.aver()

**Advantages**:
  - Simple, direct function calls
  - No object initialization needed

**Disadvantages**:
  - Must specify datadir for each call (or rely on environment variable)
  - Reads from disk every time (slower for repeated access)
  - Data is not cached between calls


Option 2: Simulation Object (Recommended)
------------------------------------------

This is the modern approach and preferred method::

    import pencil as pc
    
    # Initialize simulation object (caches data)
    sim = pc.sim.simulation('illa/test/conv-slab')
    
    # Access parameters (cached dictionary)
    par = sim.param  # Returns: dict
    
    # Access dimensions (cached object)
    dim = sim.dim    # Returns: object
    
    # Access time series (cached)
    ts = sim.ts
    
    # Access variable data (specify step)
    var = sim.get_var(varint=1)
    
    # Access xyaverage data
    aver = sim.get_aver()

**Advantages**:
  - Data is cached (faster repeated access)
  - Cleaner code (no datadir arguments)
  - Automatic path resolution
  - Object-oriented design
  - Single initialization point

**Disadvantages**:
  - Slight memory overhead for caching
  - Need to know simulation directory path


Comparison Table: Parameter Access
===================================

+-------------------------+-------------------+------------------------------------+
| Task                    | Direct Function   | Simulation Object                  |
+=========================+===================+====================================+
| Read single param       | ``pc.read.param`` | ``sim.param.get('cp', 1.0)``      |
+-------------------------+-------------------+------------------------------------+
| Fast repeated access    | Slower (disk I/O) | Faster (cached)                    |
+-------------------------+-------------------+------------------------------------+
| Initialize             | No setup needed   | ``sim = pc.sim.simulation(...)``   |
+-------------------------+-------------------+------------------------------------+
| Specify directory      | ``datadir='...'`` | Automatic in simulation object     |
+-------------------------+-------------------+------------------------------------+
| Get with default       | Manual if logic   | ``par.get('key', default)``        |
+-------------------------+-------------------+------------------------------------+

**Recommendation**: Use simulation object (Option 2) for all new scripts.


Working with Parameters
=======================

Reading Parameters
------------------

**Old IDL approach**::

    ; IDL
    pc_read_param, obj=par
    cp = par.cp
    gamma = par.gamma
    reference_state = par.reference_state

**New Python approach**::

    # Python
    sim = pc.sim.simulation('illa/test/conv-slab')
    par = sim.param
    
    cp = par.get('cp', 1.0)           # With default
    gamma = par.get('gamma', 1.667)   # With default
    reference_state = par.get('reference_state', 'polytropic')

**Key Differences**:
  - Dictionary uses ``.get()`` method with optional default
  - Must handle missing parameters explicitly
  - Parameter names are strings (dictionary keys)


Common Parameter Access Patterns
---------------------------------

**Single parameter with default**::

    rho0 = par.get('rho0', 1.0)

**Multiple parameters**::

    cp = par.get('cp', 1.0)
    gamma = par.get('gamma', 1.667)
    gravity_profile = par.get('gravity_profile', 'const')

**Check if parameter exists**::

    if 'reference_state' in par:
        reference_state = par['reference_state']
    else:
        reference_state = 'polytropic'

**Iterate over all parameters**::

    for key, value in par.items():
        print(f"{key}: {value}")

**Get all keys matching pattern**::

    flux_pars = {k: v for k, v in par.items() if 'flux' in k}


Grid Dimensions
===============

Option 1: Direct Function
--------------------------

::

    # Python (direct)
    dim = pc.read.dim()  # Returns object
    
    nx = dim.nx  # Direct attribute access
    ny = dim.ny
    nz = dim.nz
    
    x = dim.x    # 1D arrays
    y = dim.y
    z = dim.z


Option 2: Simulation Object
----------------------------

::

    # Python (recommended)
    sim = pc.sim.simulation('illa/test/conv-slab')
    dim = sim.dim  # Returns cached object
    
    # Safe attribute access with defaults
    nx = getattr(dim, 'nx', None)
    ny = getattr(dim, 'ny', None)
    nz = getattr(dim, 'nz', None)
    
    # Arrays (these exist)
    z = dim.z
    zGhost = dim.zGhost  # Include ghost cells


Comparison Table: Dimension Access
===================================

+---------------------+------------------+--------------------------------+
| Task                | Direct Function  | Simulation Object              |
+=====================+==================+================================+
| Grid size (nx)      | ``dim.nx``       | ``getattr(dim, 'nx', None)``   |
+---------------------+------------------+--------------------------------+
| Coordinate array    | ``dim.z``        | ``dim.z`` (same)               |
+---------------------+------------------+--------------------------------+
| With ghost cells    | ``dim.zGhost``   | ``dim.zGhost``                 |
+---------------------+------------------+--------------------------------+
| Safe access         | Direct (errors)  | ``getattr()`` (safe)           |
+---------------------+------------------+--------------------------------+

**Recommendation**: Use simulation object with ``getattr()`` for safe access.


Key Pencil Functions
====================

Time Series Data
----------------

**IDL**::

    ; IDL
    pc_read_ts, obj=ts
    time = ts.t
    urms = ts.urms

**Python**::

    # Python
    ts = pc.read.ts(datadir='illa/test/conv-slab')
    time = ts.t
    urms = ts.urms
    
    # Or with simulation object:
    sim = pc.sim.simulation('illa/test/conv-slab')
    ts = sim.ts
    time = ts.t
    urms = ts.urms


Variable Data
-------------

**IDL**::

    ; IDL
    pc_read_var, obj=var, datadir='...'
    print, var.uu[*,*,*, 0]  ; uu component at first point

**Python**::

    # Python
    var = pc.read.var(datadir='illa/test/conv-slab')
    uu = var.uu  # Shape: (nx, ny, nz)
    
    # Or with simulation object:
    sim = pc.sim.simulation('illa/test/conv-slab')
    var = sim.get_var(varint=1)  # varint = index
    uu = var.uu


XY-Averaged Data
----------------

**IDL**::

    ; IDL
    pc_read_xyaver, obj=aver, datadir='...'
    z = aver.z
    frad = aver.frad  ; Radiative flux

**Python**::

    # Python
    aver = pc.read.aver(datadir='illa/test/conv-slab')
    z = aver.z
    frad = aver.frad
    
    # Or with simulation object:
    sim = pc.sim.simulation('illa/test/conv-slab')
    aver = sim.get_aver()
    z = aver.z
    frad = aver.frad


Slice Data
----------

**IDL**::

    ; IDL
    pc_read_slices, obj=slices, datadir='...'
    uu_xy = slices.uu_xy

**Python**::

    # Python
    slices = pc.read.slices(datadir='illa/test/conv-slab')
    uu_xy = slices.uu_xy
    
    # Or:
    sim = pc.sim.simulation('illa/test/conv-slab')
    slices = sim.get_slices()
    uu_xy = slices.uu_xy


Grid Data
---------

**IDL**::

    ; IDL
    pc_read_grid, obj=grid, datadir='...'

**Python**::

    # Python
    grid = pc.read.grid(datadir='illa/test/conv-slab')
    
    # Or:
    sim = pc.sim.simulation('illa/test/conv-slab')
    x = grid.x
    y = grid.y
    z = grid.z


IDL to Python Function Patterns
================================

Array Operations
----------------

+-------------------------------+--------------------------------+
| IDL                           | Python                         |
+===============================+================================+
| ``where(arr > 0)``            | ``np.where(arr > 0)``          |
+-------------------------------+--------------------------------+
| Boolean indexing: ``arr[ind]``| ``arr[mask]`` (direct indexing)|
+-------------------------------+--------------------------------+
| ``arr[*]`` (all elements)     | ``arr`` or ``arr[:]``          |
+-------------------------------+--------------------------------+
| Size: ``size(arr)``           | ``len(arr)`` or ``arr.size``   |
+-------------------------------+--------------------------------+
| Transpose: ``transpose(arr)`` | ``arr.T`` or ``np.transpose()``|
+-------------------------------+--------------------------------+

Mathematical Operations
-----------------------

+-------------------------------+--------------------------------+
| IDL                           | Python                         |
+===============================+================================+
| ``alog(x)``                   | ``np.log(x)``                  |
+-------------------------------+--------------------------------+
| ``alog10(x)``                 | ``np.log10(x)``                |
+-------------------------------+--------------------------------+
| ``exp(x)``                    | ``np.exp(x)``                  |
+-------------------------------+--------------------------------+
| ``sqrt(x)``                   | ``np.sqrt(x)``                 |
+-------------------------------+--------------------------------+
| ``abs(x)``                    | ``np.abs(x)``                  |
+-------------------------------+--------------------------------+
| ``min(arr)``                  | ``np.min(arr)``                |
+-------------------------------+--------------------------------+
| ``max(arr)``                  | ``np.max(arr)``                |
+-------------------------------+--------------------------------+
| ``total(arr)``                | ``np.sum(arr)``                |
+-------------------------------+--------------------------------+
| ``mean(arr)``                 | ``np.mean(arr)``               |
+-------------------------------+--------------------------------+

Integration
-----------

**IDL**: ``integr(data, x=z)`` (custom function)

**Python**: ``np.trapz(data, z)`` (trapezoidal integration)

**Example**::

    # Integrate cooling flux vertically
    fcool_z = np.trapz(fcoolz, z)  # Returns scalar (total)


Boolean Masking
---------------

**IDL**::

    ; Find where z >= z2
    ind = where(z >= z2)
    subset = data[ind]

**Python**::

    # Create boolean mask
    mask = z >= z2
    subset = data[mask]
    
    # Or one-liner
    subset = data[z >= z2]


Example: thermo.py Pattern
===========================

**Original IDL**: ``thermo.pro`` (IDL procedure)

**Python translation**::

    import numpy as np
    import matplotlib.pyplot as plt
    import pencil as pc
    import os
    
    # Initialize
    sdir = 'illa/test/conv-slab'
    sim = pc.sim.simulation(sdir)
    par = sim.param
    dim = sim.dim
    
    # Read parameters
    cp = par.get('cp', 1.0)
    gamma = par.get('gamma', 1.667)
    g0 = par.get('g0', 1.0)
    
    # Create grid
    z = dim.z
    nz = dim.nz
    
    # Compute profiles
    temperature = compute_temperature(z, par)
    lnrho = compute_density(z, par)
    ss = compute_entropy(z, par)
    
    # Plot
    fig, ax = plot_profiles(temperature, lnrho, ss, savefig=True)
    plt.show()

**Key features**:
  - Use simulation object for data initialization
  - Access parameters with ``.get()`` and defaults
  - Organize computation into functions
  - Plotting function with optional savefig parameter


Example: pc_flux.py Pattern
============================

**Original IDL**: ``pc_fluxz.pro`` (IDL procedure)

**Python translation**::

    import numpy as np
    import matplotlib.pyplot as plt
    import pencil as pc
    
    # Initialize
    sdir = 'illa/test/conv-slab'
    sim = pc.sim.simulation(sdir)
    dim = sim.dim
    
    # Read averages
    aver = pc.read.aver(datadir=sdir)
    z = aver.z
    t = aver.t
    
    # Extract flux components
    fkinz = aver.fkinz   # Kinetic flux
    fradz = aver.fradz   # Radiative flux
    fconvz = aver.fconvz # Convective flux
    
    # Compute cooling flux (integrate)
    fcoolz = compute_cooling_flux(aver)
    
    # Plot
    fig1, ax1 = plot_fluxes_final_time(z, fluxes, t, savefig=True)
    fig2, ax2 = plot_fluxes_evolution(z, fluxes, t, savefig=True)


Data I/O Patterns
=================

Reading NetCDF Slice Data
-------------------------

**Python**::

    import numpy as np
    import pencil as pc
    
    # Read slices
    slices = pc.read.slices(datadir='illa/test/conv-slab', quiet=True)
    
    # Access specific slice
    uu_xy = slices.uu_xy  # Velocity at xy plane
    tt_z = slices.tt_z    # Temperature at z plane
    
    # Convert to numpy array
    uu_xy_array = np.array(uu_xy)


Saving Data to Text
--------------------

**Python**::

    import numpy as np
    
    # Save 1D array
    np.savetxt('data.txt', array_1d, fmt='%.6e')
    
    # Save 2D array
    np.savetxt('data.txt', array_2d, fmt='%.6e')
    
    # Save with header
    np.savetxt('data.txt', array, fmt='%.6e', 
               header='Column 1 | Column 2')


Plotting Patterns
=================

Basic Plot Function
--------------------

**Pattern**::

    def plot_something(x, y, title="", savefig=False):
        """Create and optionally save figure."""
        fig, ax = plt.subplots(figsize=(10, 6))
        ax.plot(x, y, 'b-', linewidth=2)
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_title(title)
        ax.grid(True, alpha=0.3)
        
        if savefig:
            sdir = 'illa/test/conv-slab'
            figdir = os.path.join(sdir, 'fig')
            if not os.path.isdir(figdir):
                os.makedirs(figdir)
                print(f"Created directory 'fig' inside {sdir}")
            
            figpath = os.path.join(figdir, 'something.png')
            fig.savefig(figpath, format='png', dpi=150, bbox_inches='tight')
            print(f"Figure saved as {figpath}")
        
        return fig, ax


Multi-Panel Plot Function
-------------------------

**Pattern**::

    def plot_multipanel(data_list, titles, savefig=False):
        """Create grid of subplots."""
        nplots = len(data_list)
        nrows, ncols = 2, 2
        
        fig, axes = plt.subplots(nrows, ncols, figsize=(12, 10))
        axes = axes.flatten()  # Convert to 1D array
        
        for idx, (data, title) in enumerate(zip(data_list, titles)):
            ax = axes[idx]
            ax.plot(data['x'], data['y'])
            ax.set_title(title)
        
        plt.tight_layout()
        
        if savefig:
            figdir = os.path.join(sdir, 'fig')
            if not os.path.isdir(figdir):
                os.makedirs(figdir)
            figpath = os.path.join(figdir, 'multipanel.png')
            fig.savefig(figpath, dpi=150)
        
        return fig, axes


Interactive Visualization
--------------------------

**Using Pencil Code built-in visualization**::

    import pencil as pc
    
    sim = pc.sim.simulation('illa/test/conv-slab')
    
    # Interactive slice viewer
    pc.visu.animate_interactive(
        sim,
        varfile='uu_xy',  # Velocity field
        plane='xy',
        step_ratio=10  # Plot every 10th frame
    )

**Advantages**:
  - Real-time slider control
  - Automatic scaling
  - Built-in colorbar
  - No need to write custom animation code


Debugging and Error Handling
=============================

Common Error: Missing Parameter
--------------------------------

**Problem**::

    KeyError: 'unknown_parameter'

**Solution 1: Use .get() with default**::

    value = par.get('unknown_parameter', 0.0)

**Solution 2: Check if key exists**::

    if 'unknown_parameter' in par:
        value = par['unknown_parameter']

**Solution 3: Print all available keys**::

    print("Available parameters:")
    for key in sorted(par.keys()):
        print(f"  {key}")


Common Error: Dimension Mismatch
--------------------------------

**Problem**::

    ValueError: operands could not be broadcast together

**Solution: Check array shapes**::

    print(f"Shape of var1: {var1.shape}")
    print(f"Shape of var2: {var2.shape}")
    
    # Use explicit indexing
    var1_2d = var1[:, :, 0]  # Take first z slice
    result = var1_2d + var2_2d


Common Error: File Not Found
----------------------------

**Problem**::

    FileNotFoundError: [...]/xyaverages.dat

**Solution: Check data directory**::

    import os
    sdir = 'illa/test/conv-slab'
    data_dir = os.path.join(sdir, 'data')
    
    if os.path.isdir(data_dir):
        print(f"Data directory found: {data_dir}")
    else:
        print(f"Data directory NOT found: {data_dir}")
        print(f"Available files in {sdir}:")
        for item in os.listdir(sdir):
            print(f"  {item}")


Verifying Parameter Values
---------------------------

**Print simulation parameters**::

    sim = pc.sim.simulation('illa/test/conv-slab')
    par = sim.param
    
    important_pars = ['cp', 'gamma', 'g0', 'rho0', 'Nz']
    for key in important_pars:
        value = par.get(key, 'NOT FOUND')
        print(f"{key:15s} = {value}")


Inspecting Dimension Object
----------------------------

**Print all attributes**::

    import pencil as pc
    
    dim = pc.read.dim(datadir='illa/test/conv-slab')
    
    print("Dimension attributes:")
    for attr in dir(dim):
        if not attr.startswith('_'):
            try:
                val = getattr(dim, attr)
                if not callable(val):
                    print(f"  {attr}: {type(val).__name__}")
            except:
                pass


Running Analysis Scripts
========================

Basic Execution
---------------

**Run from project directory**::

    cd /home/darkside/pencil-code
    python illa/test/conv-slab/python/thermo.py

**Run with output redirection**::

    python illa/test/conv-slab/python/thermo.py > thermo.log 2>&1

**Run in background**::

    nohup python illa/test/conv-slab/python/thermo.py &


From Within Python
------------------

**Interactive testing**::

    import sys
    sys.path.insert(0, '/home/darkside/pencil-code')
    
    # Now import the script as module
    import illa.test.conv_slab.python.thermo as thermo
    
    # Call main function
    thermo.main()


Best Practices
==============

1. **Use Simulation Object**
   
   Always use ``sim = pc.sim.simulation(...)`` for data loading in new scripts.

2. **Safe Parameter Access**
   
   Always use ``.get()`` with default values::
   
       value = par.get('parameter', default_value)

3. **Function-Based Organization**
   
   Encapsulate plotting and computation into reusable functions.

4. **Directory Management**
   
   Always check and create output directories::
   
       if not os.path.isdir(figdir):
           os.makedirs(figdir)
           print(f"Created directory...")

5. **Error Handling**
   
   Use try-except for file I/O operations::
   
       try:
           data = np.loadtxt('file.txt')
       except FileNotFoundError:
           print("File not found. Creating default data...")
           data = np.zeros((10,))

6. **Documentation**
   
   Always include docstrings::
   
       def my_function(x, y):
           """
           Brief description.
           
           Parameters
           ----------
           x : array_like
               Description of x
           y : float
               Description of y
           
           Returns
           -------
           result : ndarray
               Description of result
           """
           pass

7. **Testing**
   
   Write simple test code at module level::
   
       if __name__ == '__main__':
           main()

---

**Last Updated**: January 2026  
**Format**: ReStructuredText  
**Status**: Complete with working examples
