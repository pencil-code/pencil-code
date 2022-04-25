.. tutmathematica:

****************************
Pencil Mathematica Tutorials
****************************

Here you can find some tutorials on using Mathematica for post-processing.


Loading the package
===================

We need to modify ``init.m`` so that the path to the package is automatically added to the ``$Path`` variable in Mathematica.
First, type

.. code::

  FileNameJoin[{$UserBaseDirectory, "Kernel", "init.m"}]

in Mathematica to locate this ``init.m`` file.
Then, add line

.. code::

  AppendTo[$Path, "your/pencil/home/mathematica"]

in this file and save it. In general the path will be ``$PENCIL_HOME/mathematica/``, but of course you may put it somewhere else.
Mathematica will not search in subdirectories, so make sure the package in right in the folder.

After updating ``init.m``, restart the Mathematica kernel (``Evaluation`` -> ``Quit Kernel``).
To use the package, call ``Needs["pc`"]`` and then ``pcInitialize[]`` in a notebook or a script.

To use the package on subkernels, call ``pcParallelize[n]``.
This will launch ``n`` subkernels and load the package on each of them.
Then you can do things like ``ParallelTable[readTS[...],...]``.
Only loading the package on the master kernel is not enough.
See the discussions `here <https://mathematica.stackexchange.com/questions/11595/package-found-with-needs-but-not-with-parallelneeds>`_, and the 'Possible issues' section `here <https://reference.wolfram.com/language/ref/ParallelNeeds.html>`_.


Each time you have updated the data, remember to do ``pcInitialize[]`` and ``pcParallelize[n]`` again.
These two functions remove some persistent variables defined.


Pencil Code Commands in General
===============================

For a list of all Pencil Code commands, load the package and type ``pcFunction[]``.
To access the help of any command just type '?' followed by the command, e.g. ``?readTS``.
You can also check the full definition of the command by typing '??' followed by the command.


Reading and Plotting Time Series
================================

To read the time series, type

.. code::

  data = readTS[sim,var1,var2,...]

where ``var1``, ``var2`` etc. are entries in the time series.
The return of the right side is a ``List`` object, with its elements corresponding to ``var1``, ``var2`` etc.
You can then access, for example, the time series of ``var2`` through ``data[[2]]`` (indexing in Mathematica starts from ``1``).

Alternatively, you may also put a ``List`` object on the left side so that ``var1``, ``var2`` etc. will be assigned to each of its elements.
For example,

.. code ::

  {t,urms} = readTS[sim,"t","urms"]

Make sure that ``Length`` of the left side is equal to the number of ``var``; otherwise Mathematica will complain.

To plot the data, you can say

.. code ::

  fig = ListPlot[Transpose[{t,urms}],Joined->True]

or, in a one-line command,

.. code ::

  (** same as ListPlot[Transpose[readTS[sim,"t","urms"]]] **)
  fig = readTS[sim,"t","urms"]//Transpose//ListPlot

A few options for some internal plotting functions have been reset by the package.
For details check ``??pcLabelStyle`` and ``??pcPlotStyle``.

To export the figure in ``.eps`` format,

.. code ::

  Export["directory/to/export/figure.eps",fig]


Reading VAR files
================================

VAR files can be read using

.. code ::

  data = readVARN[sim,iVAR]

Here ``iVAR`` is the index of the VAR file and starts from 0.

By default, ghost zones will not be trimmed.
You can do it using the option ``"ltrim"->True``, or ``"ltrim"->More``;
the latter will trim ``2*nghost`` cells on each boundary.

To compute "magic" variables, you can use

.. code ::

  data = readVARN[sim,iVAR,{"oo","bb","jj"}]

Here "oo", "bb", "jj" refer to vorticity, magnetic, and current fields, respectively.

The return of ``readVARN`` is an ``Association`` object (i.e., ``Head[data]=Association``).
You can obtain all of its keys by ``Keys[data]``. Here is an example of its return:

.. code ::

  data = readVARN[sim,iVAR,{"oo","bb","jj"}];
  Keys[data]
  (* {"t", "dx", "dy", "dz", "deltay", "lx", "ly", "lz", "x", "y", "z",
     "uu1", "uu2", "uu3", "lnrho", "ooo1", "ooo2", "ooo3", "bbb1", "bbb2",
     "bbb3", "jjj1", "jjj2", "jjj3"} *)

Magic variables are named using triple characters, to avoid shadowing the auxilliary ones
written by the code (which will be "oo1" etc.).

The ``x`` coordinates of the mesh points is then ``data["x"]``, which will have length
``(16+6)^3`` if the resolutoin is ``16^3`` and ``nghost=3``.
One can form a three-dimensional map of ``uu1`` using

.. code ::

  uu1 = Transpose[ data/@{"x","y","z","uu1"} ];
  (* {{x1,y1,z1,f1},{x2,y2,z2,f2},...} *)

Sometimes the following method is also useful:

.. code ::

  Clear[uu1]
  grid = Transpose[ data/@{"x","y","z"} ];
  uu1 = Association[ Thread[ grid->data["uu1"] ] ];

Then ``uu1`` becomes a "function" and its value at ``{x1,y1,z1}`` is simply ``uu1[{x1,y1,z1}]``.

Visualizing slices from VAR files
================================

A quick way to make a density plot from ``data`` is

.. code ::

  showSlice[data, "uu1", {"z", 8}]

Here ``{"z",8}`` instructs to plot the 8th slice in the ``z`` direction.

For vector fields one can also use

.. code ::

  showSliceVector[data, "uu", {"z", 8}]

Notice the second argument is just ``"uu"`` with no index.
The function then makes a density plot of the out-of-plane component of (here ``"uu3"``),
and a superposed vector plot of the in-plane components (here ``"uu1"`` and ``"uu2"``).

Reading video files
================================

To read video or slice files, one uses

.. code ::

  {slices,times,position}=readSlice[sim,"uu1","xy2"]

The returned ``slices`` variable is a ``List`` of all slices at different times, and can
be visualized by, say, ``DensityPlot[ slices[[1]] ]``.
``position`` tells you the spatial coordinate of the slices.

Here is an example to make a video:

.. code ::

  Clear[makeFrame]
  makeFrame[ slice_,time_ ] := DensityPlot[ slice, PlotLabel->"t="<>ToString@time]
  frames = MapThread[ makeFrame, {slices,times} ];
  (* to view the video in the notebook; can be slow if too many frames*)
  ListAnimate[ frame, AnimationRunning->False ]
  (* output to a movie file *)
  Export[ "your/output/directory/video.mov", frames, FrameRate->24 ]

One can also visualize variables in a 3D box.
For more information see the comments of ``makeBox`` and ``makeBoxes``.

Running on supercomputers
================================

First, make sure Mathematica is available on the machine.
You can check this by saying ``which wolfram`` in the terminal.
If it is not installed, contact your administrator to see if it can be loaded.

Once you have loaded the Mathematica module, try ``wolfram`` in the terminal.
It should bring you to the text-based interface of Mathematica.
You can then follow the steps in the previous sections to set up the package.

There is a sample script in the directory ``$PENCIL_HOME/mathematica/sample_script.wls``.
Modify its first line according to where your ``wolfram`` is.
Remember to include the ``-script`` option.

To run a script, use ``wolframscript your_script.wls``.











