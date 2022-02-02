.. tutpython:

***********************
Pencil Python Tutorials
***********************

Here you can find some tutorials on how to modify/contribute to the Python Code 
using the Coding style :ref:`pythonstyle` and how to use the code for post-processing :ref:`pythongeneral`.



.. _pythonstyle: 

Python Coding Style
===================

Good coding style greatly improves the readability of the code. Similar
to the guidelines for the Fortran routines, it is strongly recommended
to follow some basic style rules for the python routines. These are some
recommendations extracted from `PEP 008 <https://www.python.org/dev/peps/pep-0008/>`_ and 
`Google Python Style Guide
<https://google-styleguide.googlecode.com/svn/trunk/pyguide.html>`_.


General Style Guide for Python
------------------------------

Indentation and Spaces
~~~~~~~~~~~~~~~~~~~~~~

-  Use 4 spaces per indentation level.
-  Use hanging indent for function calls over multiple lines:

   .. code:: python


        # Aligned with opening delimiter.
        foo = long_function_name(var_one, var_two,
                                 var_three, var_four)


-  Wildcard imports ( from import \* ) should be avoided, as they make
   it unclear which names are present in the namespace, confusing both
   readers and many automated tools.
-  More than one space around an assignment (or other) operator to align
   it with another should be avoided. **No**:

   .. code:: python

      x             = 1
      y             = 2
      long_variable = 3

   **Yes**:

   .. code:: python

      x = 1
      y = 2
      long_variable = 3

-  Always surround these binary operators with a single space on either
   side: assignment ( = ), augmented assignment ( += , -= etc.),
   comparisons ( == , < , > , != , <> , <= , >= , in , not in , is , is
   not ), Booleans ( and , or , not ).
-  If operators with different priorities are used, consider adding
   whitespace around the operators with the lowest priority(ies).
   
   **Yes**:

   .. code:: python

      i = i + 1
      submitted += 1
      x = x*2 - 1

   **No**:

   .. code:: python

      
      i=i+1
      submitted +=1
      x = x * 2 - 1
      
-  Don’t use spaces around the = sign when used to indicate a keyword
   argument or a default parameter value. 
   
   **Yes**:

   .. code:: python

      def complex(real, imag=0.0):
            return magic(r=real, i=imag)
      

   **No**:

   .. code:: python

      def complex(real, imag = 0.0):
            return magic(r = real, i = imag)
     
Comments
~~~~~~~~

-  Comments should be complete sentences.
-  Block comments generally apply to some (or all) code that follows
   them, and are indented to the same level as that code. Each line of a
   block comment starts with a # and a single space (unless it is
   indented text inside the comment). Paragraphs inside a block comment
   are separated by a line containing a single # .

Docstrings
~~~~~~~~~~

Always use docstrings for classes and functions which can be accessed by
the user. 

We are now working with read the docs and sphinx to create automatic documentation for the code, hence we have updated the style guide for creating docstrings.

We are using Numpy docstring style, and require the following fields in the docstring:

- General description of the Class/function
- Signature: how the function can be called
- Parameters: list of parameters of the class/function
- Returns: type of variable the function returns
- Examples: at least one example of usage
- Notes (ptional): any further comments to the function


.. code:: python

   def complex(real=0.0, imag=0.0):
        """
        Form a complex number.

        Signature
        ---------
        complex(real, imag)

        Parameters
        ----------
         *real*: float
             the real part (default 0.0)
         *imag*: float
             the imaginary part (default 0.0)

        Returns
        -------
        complex number with real and imaginary part

        Examples 
        --------
        Define two complex numbers:
        >>> a = complex(3, 5)
        >>> b = complex(4, 7)
        >>> print(a)
        (3+5j)
        >>> a + b
        (7+12j)
        """
  
Naming Convention
~~~~~~~~~~~~~~~~~

module_name, package_name, ClassName, method_name, ExceptionName,
function_name, GLOBAL_CONSTANT_NAME, global_var_name, instance_var_name,
function_parameter_name, local_var_name

Exceptions for >our< code: datadir, varfile, varfiles, …

pylint
~~~~~~

Run pylint over your code. pylint is a tool for finding bugs and style
problems in Python source code. It finds problems that are typically
caught by a compiler for less dynamic languages like C and C++.

black
~~~~~~

Run black over your code for automatic formatting.
This makes sure that all the above criteria (apart from the doc string)
are fullfilled.

Default Function Arguments
~~~~~~~~~~~~~~~~~~~~~~~~~~

Do not use mutable objects as default values in the function or method
definition. 

**Yes**:

.. code:: python

   def foo(a, b=None):
       if b is None:
           b = []

**No**: 

.. code:: python

    def foo(a, b=[]):


Private Methods
~~~~~~~~~~~~~~~

Python does not know any private methods or class member. In order to
somewhat hide such methods use two underscores in the function
definition: ``def __magicAttributes(self, param):``.

Others
~~~~~~

-  Use ``''.startswith()`` and ``''.endswith()`` instead of string
   slicing to check for prefixes or suffixes. startswith() and
   endswith() are cleaner and less error prone. For example: **Yes**:
   ``if foo.startswith('bar'):`` **No**: ``if foo[:3] == 'bar':``
-  For sequences, (strings, lists, tuples), use the fact that empty
   sequences are false. 

   **Yes**:

   .. code:: python
     
      if not seq:
      if seq:
      

   **No**:

   .. code:: python
      
      if len(seq)
      if not len(seq)
      

-  Don’t compare boolean values to True or False using == . 

   **Yes**:

   .. code:: python

      if greeting:

   **No**:

   .. code:: python

      if greeting == True:

-  Check if a variable has a particular type by using ``isinstance``,
   e.g.: ``isinstance(my_variable, list)``.


Pencil Code Specific Style
--------------------------

Classes/Objects
~~~~~~~~~~~~~~~

Use classes as much as possible. When you write a function try to embed
it into a class as **init** function which should return the desired
result. This has the advantage of adding methods to the returned object
which can modify the data. Read-methods always give back objects
containing the whole information (container philosophy). Therefore we
use classes if possible.

Data Directory
~~~~~~~~~~~~~~

The default data directory is always ‘./data’ and not ‘data’.

File Headers
~~~~~~~~~~~~

Start each file with the file ID and  a short
description of the routines.
(The authors' list is no longer required since it can be easily accesed through git history.)

.. code:: python

   
   # varfile.py
   #
   # Read VAR files. Based on the read_var.pro IDL script.
   #
   # NB: the f array returned is C-ordered: f[nvar,nz,ny,nx]
   #     NOT Fortran as in Pencil (& IDL):  f[nx,ny,nz,nvar]
   
  

Import Libraries
~~~~~~~~~~~~~~~~

-  Import numpy as *np* instead of *N*.
-  Import pylab as *plt* instead of *P*.

If you need to access libraries in some routines in your module, import
them in the routine, rather than the head of the module. That way they
are not visible by the user.

**Yes**:

.. code:: python

   # my_module.py

   class MyClass(object):
       """
       Some documentation.
       """

       def __init__(self):
           import numpy as np

           self.pi = np.pi

**No**:

.. code:: python

        # my_module.py
        import numpy as np

        class MyClass(object):
        """
        Some documentation.
        """

        def __init__(self):
                self.pi = np.pi</pre>


Further Reading
---------------

`<https://www.python.org/dev/peps/pep-0008/#tabs-or-spaces>`_

`<https://google-styleguide.googlecode.com/svn/trunk/pyguide.html>`_



.. _pythongeneral: 

Pencil Code Commands in General
===============================

For a list of all Pencil Code commands start IPython and type ``pc. <TAB>`` (as with auto completion).
To access the help of any command just type the command followed by a '?' (no spaces), e.g.:

.. code:: 

        pc.math.dot?
        Type:       function
        String Form:<function dot at 0x7f9d96cb0cf8>
        File:       ~/pencil-code/python/pencil/math/vector_multiplication.py
        Definition: pc.math.dot(a, b)
        Docstring:
        take dot product of two pencil-code vectors a & b with shape

        a.shape = (3, mz, my, mx)
        
You can also use ``help(pc.math.dot)`` for a more complete documentation of the command.

There are various reading routines for the Pencil Code data. All of them return an object with the data. To store the data into a user defined variable type e.g.

.. code:: python

        ts = pc.read.ts()

Most commands take some arguments. For most of them there is a default value, e.g.

.. code:: python

        pc.read.ts(file_name='time_series.dat', datadir='data')

You can change the values by simply typing e.g.


.. code:: python

        pc.read.ts(datadir='other_run/data')


Reading and Plotting Time Series
================================

Reading the time series file is very easy. Simply type

.. code:: python

        ts = pc.read.ts()

and python stores the data in the variable ``ts``. 
The physical quantities are members of the object ``ts`` and can be accessed accordingly, e.g. ``ts.t, ts.emag``. 
To check which other variables are stored simply do the tab auto completion ``ts. <TAB>``.

 Plot the data with the matplotlib commands:

.. code:: python

        plt.plot(ts.t, ts.emag)


The standard plots are not perfect and need a little polishing. See further down about making pretty plots.
You can save the plot into a file using the GUI or with

.. code:: python

        plt.savefig('plot.eps')

Reading and Plotting VAR files and slice files
==============================================

Read var files:

.. code:: python

        var = pc.read.var()

Read slice files:

.. code:: python

        slices = pc.read.slices(field='bb1', extension='xy')

This returns an object ``slices`` with members ``t`` and ``xy``. 
The last contains the additional member ``xy``.


If you want to plot e.g. the x-component of the magnetic field at the central plane simply type:

.. code:: python
        
        plt.imshow(var.bb[0, 128, :, :].T, origin='lower', extent=[-4, 4, -4, 4], interpolation='nearest', cmap='hot')

For a complete list of arguments of ``plt.imshow`` refer to its documentation.

For a more interactive function plot use:

.. code:: python

        pc.visu.animate_interactive(slices.xy.bb, slices.t)

.. warning::

        arrays from the reading routines are ordered ``f[nvar, mz, my, mx]``, i.e. reversed to IDL. 
        This affects reading var files and slice files.

Create a custom VAR0 or var.dat
===============================

With the functionality of writing snapshots directly into ``VAR*`` or ``var.dat`` the user can now generate an initial condition directly from a numpy array or modify the last snapshot and continue running. The function to be used is in ``python/pencil/io/snapshot.py`` and is called ``write_snapshot``. Here we outline how to generate an initial condition. For modifying the ``var.dat`` only the last steps are necessary.

First we need an empty run. For this let us use ``samples/kin-dynamo``


.. code:: python

        cd pencil-code/samples/kin-dynamo
        pc_setupsrc

In principle we can use any initial condition, as we are going to over write it. But it is cleaner to use

.. code::

        INITIAL_CONDITION = noinitial_condition

in ``src/Makefile.local``. Compile and start:

.. code:: bash

        make
        pc_start

This generates a ``VAR0`` and ``var.dat`` in every proc directory.

Our snapshot writing routine needs to know the cpu structure. Furthermore, we need to know the indices of the primary variables. The first can be obtained from ``src/cparam.local``, while the latter can be read from the newly generated ``data/index.pro``. The numpy arrays that are written need to have the shape [nvar, nz, ny, nz] with the correct order of variables and no ghost zones. Optionally, the number of ghost zones, which is usually 3, can be specified.

Putting it all together our python routine would look something like this:

.. code:: python

        import numpy as np
        import pencil as pc

        # Read the data to obtain the shape of the arrays, rather than the actual data.
        var = pc.read.var(trimall=True)

        # Modify the data.
        var.aa += np.random.random(var.aa.shape)

        # Write the new VAR0 and var.dat files.
        pc.io.write_snapshot(var.aa, file_name='VAR0', nprocx=1, nprocy=1, nprocz=1)
        pc.io.write_snapshot(var.aa, file_name='var.dat', nprocx=1, nprocy=1, nprocz=1)


Examples
========

Standard plots with any plotting library are not the prettiest ones. The same is true for matplotlib. Here are a few pretty examples of plots where the default style is changed. You can add your commands into a script e.g. ``plot_results.py`` and execute it in IPython with ``execfile('plot_results.py')``.

Simple plot:

.. code:: python

        import pencil as pc
        import numpy as np
        import pylab as plt

        # Read the time_series.dat.
        ts = pc.read.ts()

        # Prepare the plot.
        # Set the size and margins.
        width = 8
        height = 6
        plt.rc('text', usetex=True)
        plt.rc('font', family='arial')
        plt.rc("figure.subplot", left=0.2)
        plt.rc("figure.subplot", right=0.95)
        plt.rc("figure.subplot", bottom=0.15)
        plt.rc("figure.subplot", top=0.90)
        figure = plt.figure(figsize=(width, height))
        axes = plt.subplot(111)

        # Make the actual plot.
        plt.semilogy(ts.t, ts.brms, linestyle='-', linewidth=2, color='black', label=r'$\langle\bar{B}\rangle$')
        plt.semilogy(ts.t, ts.jrms, linestyle='--', linewidth=2, color='blue', label=r'$\langle\bar{J}\rangle$')
        plt.semilogy(ts.t, ts.jmax, linestyle=':', linewidth=2, color='red', label=r'$J_{\rm max}$')

        plt.xlabel(r'$t$', fontsize=25)
        plt.ylabel(r'$\langle\bar{B}\rangle, \langle\bar{J}\rangle, J_{\rm max}$', fontsize=25)
        plt.title('various quantities', fontsize=25, family='serif')

        # Prepare the legend.
        plt.legend(loc=1, shadow=False, fancybox=False, numpoints=1)
        leg = plt.gca().get_legend()
        # Change the font size of the legend.
        ltext = leg.get_texts() # all the text.Text instance in the legend
        for k in range(len(ltext)):
                legLine = ltext[k]
                legLine.set_fontsize(25)
        frame = leg.get_frame()
        frame.set_facecolor('1.0')
        leg.draw_frame(False)

        # Make plot pretty.
        plt.xticks(fontsize=20, family='serif')
        plt.yticks(fontsize=20, family='serif')
        axes.tick_params(axis='both', which='major', length=8)
        axes.tick_params(axis='both', which='minor', length=4)

        # Create an offset between the xylabels and the axes.
        for label in axes.xaxis.get_ticklabels():
                label.set_position((0, -0.03))
        for label in axes.yaxis.get_ticklabels():
                label.set_position((-0.03, 0))


Simple 2d plot:

.. code:: python

        import pencil as pc
        import numpy as np
        import pylab as plt

        # Read the slices.
        slices = pc.read.slices()

        # Read the grid size.
        grid = pc.read.grid()
        x0 = grid.x[3]
        x1 = grid.x[-4]
        y0 = grid.y[3]
        y1 = grid.y[-4]

        # Prepare the plot.
        # Set the size and margins.
        width = 8
        height = 6
        plt.rc('text', usetex=True)
        plt.rc('font', family='arial')
        plt.rc("figure.subplot", left=0.15)
        plt.rc("figure.subplot", right=0.95)
        plt.rc("figure.subplot", bottom=0.15)
        plt.rc("figure.subplot", top=0.95)
        figure = plt.figure(figsize=(width, height))
        axes = plt.subplot(111)

        # Make the actual plot.
        plt.imshow(slices.xy.bb1[0, :, :].T, origin='lower', interpolation='nearest', cmap='hot', extent=[x0, x1, y0, y1])
        plt.xlabel(r'$x$', fontsize=25)
        plt.ylabel(r'$y$', fontsize=25)

        # Set the colorbar.
        cb = plt.colorbar()
        cb.set_label(r'$B_{x}(x,y,z=0)$', fontsize=25)
        cbytick_obj = plt.getp(cb.ax.axes, 'yticklabels')
        plt.setp(cbytick_obj, fontsize=15, family='serif')

        # Make plot pretty.
        plt.xticks(fontsize=20, family='serif')
        plt.yticks(fontsize=20, family='serif')
        axes.tick_params(axis='both', which='major', length=8)
        axes.tick_params(axis='both', which='minor', length=4)

        # Create an offset between the xylabels and the axes.
        for label in axes.xaxis.get_ticklabels():
                label.set_position((0, -0.03))
        for label in axes.yaxis.get_ticklabels():
                label.set_position((-0.03, 0))


IDL to Python guide
===================

A large array of idl scripts have been developed over the years, and many of them served their purpose at the time, but there are many others
of general purpose. Below is a small selection of examples of idl call sequences along with their python counterparts.

Here are the links to a few potentially useful sites:

1. `IDL to Python bridge <https://www.l3harrisgeospatial.com/docs/IDLToPython.html>`_

2. `IDL commands in numerical Python <http://mathesaurus.sourceforge.net/idl-python-xref.pdf>`_

===============================   ======
IDL                               Python
===============================   ======
pc_read_var,obj=var,/trimall      var = pc.read.var(var_file = 'var.dat', trimall = True, sim = SIM)    
help,var                          help(var)       
pc_read_param,obj=param           pc.read.param()
===============================   ======
