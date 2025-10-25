.. _quick_start_guide:

*****************
Quick start guide
*****************

The Pencil Code is an open-source simulation code written mainly in Fortran
and released under the GPL license. General information can be found at our
official `homepage <http://pencil-code.org/>`__.

This guide is a concise, self-contained introduction to the essentials:
downloading, compiling, running, and post-processing your first simulation.
It is designed to get you from zero to scientifically meaningful output with
as little friction as possible.

For deeper explanations, theoretical background, and the full feature set,
please refer to the :ref:`manual`.

Required software
=================

The Pencil Code runs best on Unix-like systems. Below are brief instructions
for GNU/Linux and macOS users. Other systems are welcome to join, but results
may vary depending on compiler mood and cosmic alignment.


GNU/Linux
----------

You will need both a Fortran and a C compiler to build the code. They should
belong to the same distribution and version (for instance, GNU GCC or Intel
compilers). If you can run ``gfortran --version`` and get a friendly response,
you are halfway there.


MacOS X
-------

On macOS, the easiest route is to install Xcode from
`<http://developer.apple.com/>`_ (registration required). 
Alternatively, you can
install a ``gfortran`` binary package from
http://gcc.gnu.org/wiki/GFortranBinaries. 
Download the archive and use the
included installer. It installs into ``/usr/local/gfortran`` and provides a
symbolic link in ``/usr/local/bin/gfortran``.

You may need to add the following line to your “``.cshrc``” file in your
``/home`` directory:

.. code:: csh

     setenv PATH /usr/local/bin:\$PATH


Download the Pencil Code
========================

For full details, see :ref:`download`. The short version: use Git.

.. code:: bash

     git clone https://github.com/pencil-code/pencil-code.git

That will create a directory named ``pencil-code`` containing everything you
need. Congratulations — you now own a respectable quantity of Fortran.

Configure the shell environment
===============================

Before doing anything else, you need to load the environment variables used by
the Pencil Code tools. Move into the freshly downloaded directory:

.. code:: bash

     cd pencil-code

If you use a ``sh``-compatible shell (such as ``bash``), simply type:

.. code:: bash

     bash sourceme.sh

For more details, see :ref:`man1_environment_settings`.




Quick look at the code structure
================================

The downloaded ``pencil-code`` directory contains several subdirectories:

* ``doc`` – Documentation sources, including :file:`manual.pdf` and other
   supporting material.

* ``samples`` – A collection of ready-to-run example problems.

* ``config`` – Configuration files for compilation and setup.

* ``src`` – The actual source code.

* ``bin`` and ``lib`` – Supplemental scripts and libraries.

* ``idl``, ``python``, ``julia``, etc. – Post-processing tools in various languages.

For a detailed explanation, see :ref:`man1_directory_tree`.



Your first simulation run
=========================

Every |PC| simulation requires a set of configuration and source files. Some
define physical parameters and numerical settings, others handle the Fortran
source modules. See :ref:`man1_files_in_rundir` for a deeper explanation.

The easiest way to start is by using one of the pre-configured sample problems.
Here we will use ``pencil-code/samples/conv-slab``.


Create a new run-directory
--------------------------

Create a new directory for your simulation and copy the contents of the sample
setup:

.. code:: bash

     mkdir -p myuser/test
     cd myuser/test
     cp -r $PENCIL_HOME/samples/conv-slab/* .

Choose your run directory wisely. Simulations can generate a substantial amount
of data, so avoid using locations with strict storage quotas, such as your home
directory on shared systems.


Linking to the sources
----------------------

To connect your run directory with the main Pencil Code source directory, run:

.. code:: bash

     pc_setupsrc

This creates the required symbolic links. See
:ref:`man1_linking_scripts_and_source` for more information.

Makefile and parameters
-----------------------

Two local configuration files define the essential parts of a simulation setup:

* ``src/Makefile.local`` – Lists the modules to be used.

* ``src/cparam.local`` – Defines the grid size and number of processors.

Take a look at these files before running. They are short but powerful.

Single-processor
^^^^^^^^^^^^^^^^

A minimal ``src/Makefile.local`` for a single-processor run could be:

.. code:: fortran

     MPICOMM=nompicomm

In ``src/cparam.local``, set the number of processors to ``1`` accordingly:

.. code:: fortran

     integer, parameter :: ncpus=1,nprocx=1,nprocy=1,nprocz=ncpus/(nprocx*nprocy)
     integer, parameter :: nxgrid=32,nygrid=nxgrid,nzgrid=nxgrid

Multi-processor
^^^^^^^^^^^^^^^^

To use MPI for multi-processor simulations, make sure an MPI library is
installed and update your configuration in ``src/Makefile.local``:

.. code:: fortran

     MPICOMM=mpicomm

Adjust ``ncpus`` in ``src/cparam.local`` to match your processor layout. A
good rule of thumb is to keep 32 grid points along the x-direction to make
efficient use of SIMD units. To compile, use a configuration file with the
``_MPI`` suffix, as shown below.



Compilation
------------

To compile the code with default GNU compilers (single processor), simply run [#]_ :

.. code:: bash

     pc_build

For a multi-processor build using MPI:

.. code:: bash

     pc_build -f GNU-GCC_MPI

.. note::

     Depending on your system, a simpler :command:`pc_build` can also work for a multi-processor build.

For additional details, see :ref:`man1_quick_instructions`.

.. [#] You can use a pre-defined configuration file corresponding to your compiler package. E.g. the default compilers are ``gfortran`` together with ``gcc`` and the code is being built  with the default options (not using MPI)

Using a different compiler (optional)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If you wish to use a different compiler package (for example, Intel or Cray),
try one of the following:

.. code:: bash

     pc_build -f Intel
     pc_build -f Intel_MPI
     pc_build -f Cray
     pc_build -f Cray_MPI

Additional predefined configurations can be found in
``pencil-code/config/compilers/*.conf``.

Changing compiler options (optional)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

You can also define host-specific configuration files in
``pencil-code/config/hosts/``. By default, :command:`pc_build` searches for a
file based on your host ID, which you can view with:

.. code:: bash

     pc_build -i

You may add your own configuration file as ``host-ID.conf`` and adjust compiler
flags as needed. A good example to adapt is
``pencil-code/config/hosts/IWF/host-andromeda-GNU_Linux-Linux.conf``.

To clean up all generated files:

.. code:: bash

     pc_build --cleanall

Running...
----------

Set up your simulation by editing the following files:

* ``start.in`` – Initial conditions.

* ``run.in`` – Main runtime parameters.

* ``print.in`` – Select which quantities appear in ``data/time_series.dat``.

Make sure an empty ``data/`` directory exists:

.. code:: bash

     mkdir data

Now, launch the simulation:

.. code:: bash

     pc_run

If everything is working correctly, the output should include:

.. code:: text

     start.x has completed successfully


Once initialized, the code will begin printing quantities to the console, such
as:

.. code:: text

      --it-----t-------dt------urms----umax----rhom----ssm----dtc---dtu---dtnu-dtchi-
         0    0.00 6.793E-03  0.0063  0.0956 14.4708 -0.4460 0.978 0.025 0.207 0.345
        10    0.07 6.793E-03  0.0056  0.0723 14.4708 -0.4464 0.978 0.019 0.207 0.345
        20    0.14 6.793E-03  0.0053  0.0471 14.4709 -0.4467 0.978 0.019 0.207 0.345
      .......

When finished, you should see:

.. code:: text

     Simulation finished after        xxxx  time-steps
     .....
     Wall clock time/timestep/meshpoint [microsec] = ...

An empty file named ``COMPLETED`` will appear in your run directory when the
simulation ends.

To verify your results against the reference output for this sample, run:

.. code:: bash

     diff reference.out data/time_series.dat

If they match, congratulations — your setup is alive and calculating.

Troubleshooting...
------------------

If compilation fails, try cleaning and rebuilding (optionally with MPI):

.. code:: bash

     pc_build --cleanall
     pc_build -f GNU-GCC_MPI

If issues persist, please report them on our mailing list:
http://pencil-code.nordita.org/contact.php. Include the exact step that failed,
the error message, and confirm that all required steps above were followed.

Also include your operating system, shell type, and the full output of:

.. code:: bash

     bash
     cd path/to/your/pencil-code/
     source sourceme.sh
     echo $PENCIL_HOME
     ls -la $PENCIL_HOME/bin
     cd samples/1d-tests/jeans-x/
     gcc --version
     gfortran --version
     pc_build --cleanall
     pc_build -d

If using MPI, please also include:

.. code:: bash

     mpicc --version
     mpif90 --version
     mpiexec --version

Welcome to the world of Pencil Code — may your runs be stable and your :abbr:`CFL (Courant-Friedrichs-Lewy condition)`
numbers merciful.



Data post-processing
====================

IDL visualization (optional,)
-----------------------------------------

GUI-based visualization (recommended for quick inspection)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The most simple approach to visualize a Cartesian grid setup is to run
the Pencil Code GUI and to select the files and physical quantities you
want to see:

.. code:: idl

   IDL> .r pc_gui

If you miss some physical quantities, you might want to extend the two
IDL routines ``pc_get_quantity`` and ``pc_check_quantities``. Anything
implemented there will be available in the GUI, too.

Command-line based processing of “big data”
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Please check the documentation inside these files:

+------------------------------------------------+--------------------------------------------+
|``pencil-code/idl/read/pc_read_var_raw.pro``    | efficient reading of raw data              |
+------------------------------------------------+--------------------------------------------+
|``pencil-code/idl/read/pc_read_subvol_raw.pro`` | reading of sub-volumes                     |
+------------------------------------------------+--------------------------------------------+
|``pencil-code/idl/read/pc_read_slice_raw.pro``  | reading of any 2D slice from 3D snapshots  |
+------------------------------------------------+--------------------------------------------+
|``pencil-code/idl/pc_get_quantity.pro``         | compute physical quantities out of raw data|
+------------------------------------------------+--------------------------------------------+
|``pencil-code/idl/pc_check_quantities.pro``     | dependency checking of physical quantities |
+------------------------------------------------+--------------------------------------------+


in order to read data efficiently and compute quantities in physical
units.

Command-line based data analysis (may be inefficient)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Several idl-procedures have been written (see in ‘``pencil-code/idl``’)
to facilitate inspecting the data that can be found in raw format in
‘``jeans-x/data``’. For example, let us inspect the time series data

.. code:: idl 

   IDL> pc_read_ts, obj=ts

The structure ``ts`` contains several variables that can be inspected by

.. code:: idl

   IDL> help, ts, /structure
   ** Structure <911fa8>, 4 tags, length=320, data length=320, refs=1:
      IT              LONG      Array[20]
      T               FLOAT     Array[20]
      UMAX            FLOAT     Array[20]
      RHOMAX          FLOAT     Array[20]

The diagnostic ``UMAX``, the maximal velocity, is available since it was
set in “``jeans-x/print.in``”. Please check the manual for more
information about the input files.

We plot now the evolution of ``UMAX`` after the initial perturbation
that is defined in “``start.in``”:

.. code:: idl 

   IDL> plot, ts.t, alog(ts.umax)

The complete state of the simulation is saved as snapshot files in
“``jeans-x/data/proc0/VAR*``” every ``dsnap`` time units, as defined in
“``jeans-x/run.in``”. These snapshots, for example “``VAR5``”, can be
loaded with:

.. code:: idl

   IDL> pc_read_var, obj=ff, varfile="VAR5", /trimall

Similarly ``tag_names`` will provide us with the available variables:

.. code:: idl

   IDL> print, tag_names(ff)
   T X Y Z DX DY DZ UU LNRHO POTSELF

The logarithm of the density can be inspected by using a GUI:

.. code:: idl

   IDL> cslice, ff.lnrho

Of course, for scripting one might use any quantity from the ``ff``
structure, like calculating the average density:

.. code:: idl

   IDL> print, mean(exp(ff.lnrho))

Python visualization (optional)
-------------------------------

Be advised that the Python support is still not complete or as
feature-rich as for IDL. Furthermore, we move to Python3 in 2020, and
not all the routines have been updated yet.

Python module requirements
^^^^^^^^^^^^^^^^^^^^^^^^^^^

In this example we use the modules: ``numpy`` and ``matplotlib``. A
complete list of required module is included in
“``pencil-code/python/pencil/README``”.

Using the ``pencil`` module
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

After sourcing the “``sourceme.sh``” script (see above), you should be
able to import the ``pencil`` module:

.. code:: python

   import pencil as pc

Some useful functions:

===============================   ======
``pc.read.ts``                    read “``time_series.dat``” file. Parameters are added as members of the class
``pc.read.slices``                read 2D slice files and return two arrays: (nslices,vsize,hsize) and (time)
``pc.visu.animate_interactive``   assemble a 2D animation from a 3D array
===============================   ======


Some examples of postprocessing with Python can be found in the
:ref:` python documentation <modpython>` and in the :ref:` python tutorials <tutpython>`.