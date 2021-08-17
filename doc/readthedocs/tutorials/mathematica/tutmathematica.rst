.. tutmathematica:

****************************
Pencil Mathematica Tutorials
****************************

Here you can find some tutorials on using Mathematica for post-processing.

.. warning::

        This page is still under development.


Loading the package
===================

There are two ways of telling Mathematica where the package is:

1. Modifying ``init.m`` so that the path to the package is automatically added to the ``$Path`` variable in Mathematica.
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
To use the package, call ``Needs["pc`"]`` in a notebook or a script.

2. Alternatively, if you don't want to modify ``init.m``, you may also call

.. code::

 Needs["pc`","path/to/this/package"]

each time.


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
















