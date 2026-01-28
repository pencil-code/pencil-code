Summary: IDL to Python Translation for conv-slab Analysis
===========================================================

Completed Work
==============

1. ✅ **thermo.py** - Thermodynamic Analysis
   
   **Status**: Complete and updated
   
   **Features**:
   
   - Reads simulation parameters using ``sim.param`` (dictionary)
   - Reads grid dimensions using ``sim.dim`` 
   - Computes initial stratification profiles (temperature, density, entropy)
   - Includes plotting function for visualization
   - Generates figure with configurable directory structure
   
   **Key Functions**:
   
   - ``initialize_stratification()`` - Computes polytropic stratification
   - ``plot_profiles()`` - Creates 3-panel comparison plot with optional save
   
   **Usage**::
   
      python thermo.py
   
   **Output**:
   
   - Figures saved to ``fig/stratification_profiles.png``


2. ✅ **pc_flux.py** - Energy Flux Visualization
   
   **Status**: New, fully translated from IDL
   
   **Features**:
   
   - Reads xy-averaged data using ``pc.read.aver()``
   - Extracts flux components: kinetic, radiative, convective, turbulent, cooling
   - Computes cooling flux by vertical integration
   - Multiple visualization functions with configurable saving
   
   **Key Functions**:
   
   - ``plot_fluxes_final_time()`` - All fluxes at final timestep (with savefig option)
   - ``plot_fluxes_time_evolution()`` - Grid showing flux evolution (with savefig option)
   - ``plot_total_flux_time_series()`` - Total flux time series (with savefig option)
   
   **Output Figures**:
   
   - ``fig/flux_final_time.png`` - Flux profiles at t_final
   - ``fig/flux_evolution.png`` - 6-panel time evolution
   - ``fig/total_flux_timeseries.png`` - Time series at different heights
   
   **Usage**::
   
      python pc_flux.py


3. ✅ **TRANSLATION_GUIDE.md** - Complete Documentation
   
   **Status**: Comprehensive and updated
   
   **Contents**:
   
   - Two methods for reading parameters: ``pc.read.param()`` vs ``sim.param``
   - Two methods for grid dimensions: ``pc.read.dim()`` vs ``sim.dim``
   - Detailed comparison tables
   - All key Pencil functions with examples
   - IDL ↔ Python pattern mappings
   - Debugging tips and error solutions
   - Running instructions


Key Design Decisions
====================

Parameter Access: ``sim.param`` (Dictionary)
---------------------------------------------

**Recommended approach**::

   sim = pc.sim.simulation('illa/test/conv-slab')
   par = sim.param  # Returns: dict
   cp = par.get('cp', 1.0)  # Use .get() with default

**Why**: 
   - Faster (cached in simulation object)
   - Cleaner syntax with ``.get()``
   - No need for ``getattr()``

Grid Dimensions: ``sim.dim`` (Object)
--------------------------------------

**Recommended approach**::

   dim = sim.dim  # Returns: object, cached
   nx = getattr(dim, 'nx', None)  # Still an object
   z = dim.z  # Access coordinate array

Plotting: Functions Over Scripts
---------------------------------

All plotting logic is encapsulated in functions:

- Reusable and testable
- Can be called from other scripts
- Easy to modify and extend
- Professional code organization
- Integrated figure directory management with ``savefig`` option

**Example**::

   def plot_fluxes_final_time(z, fluxes, t, savefig=False):
       """docstring"""
       fig, ax = plt.subplots(...)
       # plotting code
       
       if savefig:
           figdir = os.path.join(sdir, 'fig')
           if not os.path.isdir(figdir):
               os.makedirs(figdir)
           # save figure
       
       return fig, ax


Translation Patterns Used
=========================

+-----------------------------------+------------------------------------+-------------------------------------------+
| IDL                               | Python                             | Example                               |
+===================================+====================================+===================================+
| ``where(z >= z2)``                | Boolean mask                       | ``mask = z >= z2; data[mask]``        |
+-----------------------------------+------------------------------------+-------------------------------------------+
| ``alog(x)``                       | ``np.log(x)``                      | ``np.log(rho0)``                      |
+-----------------------------------+------------------------------------+-------------------------------------------+
| ``exp(x)``                        | ``np.exp(x)``                      | ``np.exp(gamma*ss/cp)``               |
+-----------------------------------+------------------------------------+-------------------------------------------+
| ``integr(data, x=z)``             | ``np.trapz(data, z)``              | Integration over z                    |
+-----------------------------------+------------------------------------+-------------------------------------------+
| ``make_array()``                  | ``np.zeros_like()``                | Array initialization                  |
+-----------------------------------+------------------------------------+-------------------------------------------+
| ``pc_read_param``                 | ``pc.read.param()`` or ``sim.param``| Parameter loading                     |
+-----------------------------------+------------------------------------+-------------------------------------------+
| ``pc_read_xyaver``                | ``pc.read.aver()``                 | Averages loading                      |
+-----------------------------------+------------------------------------+-------------------------------------------+


File Locations
==============

::

   pencil-code/
   ├── illa/test/conv-slab/
   │   ├── python/
   │   │   ├── thermo.py                    (✅ Updated)
   │   │   ├── pc_flux.py                   (✅ Updated)
   │   │   ├── TRANSLATION_GUIDE.md         (✅ Updated)
   │   │   ├── ptvsurms.py                  (existing)
   │   │   ├── pstability.py                (existing)
   │   │   └── ... other scripts
   │   ├── idl/
   │   │   ├── thermo.pro                   (reference)
   │   │   └── pc_fluxz.pro                 (reference)
   │   └── tutorial-conv-slab.rst           (tutorial document)


Testing the Scripts
===================

Test 1: Check Thermo Script
---------------------------

::

   cd /home/darkside/pencil-code
   python illa/test/conv-slab/python/thermo.py
   
   # Expected output:
   # - Console: parameter values and ranges
   # - File: fig/stratification_profiles.png


Test 2: Check Flux Script
--------------------------

::

   python illa/test/conv-slab/python/pc_flux.py
   
   # Expected output:
   # - Console: flux statistics
   # - Files: fig/flux_final_time.png, fig/flux_evolution.png, fig/total_flux_timeseries.png


Test 3: Verify Translation Guide
---------------------------------

::

   # View in text editor or markdown viewer
   cat illa/test/conv-slab/python/TRANSLATION_GUIDE.md


Next Steps (Optional)
=====================

1. Create additional analysis scripts for other IDL routines
2. Add Jupyter notebook versions for interactive analysis
3. Create data export routines (e.g., to HDF5 or CSV)
4. Implement spectral analysis for frequency content
5. Add animation generation for time-series data


Advantages of Python Translation
=================================

+-------------------+-------------------+--------------------------------------+
| Aspect            | IDL               | Python                             |
+===================+===================+======================================+
| Speed             | Compiled; faster  | Interpreted; often fast enough     |
+-------------------+-------------------+--------------------------------------+
| Learning          | Proprietary       | Standard, widely used              |
+-------------------+-------------------+--------------------------------------+
| Integration       | Limited           | Excellent (numpy, scipy, etc.)     |
+-------------------+-------------------+--------------------------------------+
| Plotting          | Limited           | Excellent (matplotlib, plotly)     |
+-------------------+-------------------+--------------------------------------+
| Open-source       | No                | Yes                                |
+-------------------+-------------------+--------------------------------------+
| Community         | Small             | Huge                               |
+-------------------+-------------------+--------------------------------------+
| Maintenance       | Vendor-dependent  | Community-driven                   |
+-------------------+-------------------+--------------------------------------+
| Cost              | $                 | Free                               |
+-------------------+-------------------+--------------------------------------+


References
==========

- Pencil Code Python library: ``pencil-code/python/pencil/``
- Tutorial document: ``illa/test/conv-slab/tutorial-conv-slab.rst``
- IDL reference scripts: ``illa/test/conv-slab/idl/``
- NumPy documentation: https://numpy.org/doc/
- Matplotlib guide: https://matplotlib.org/

---

**Author**: AI Assistant  
**Date**: January 2026  
**Status**: Complete and tested
