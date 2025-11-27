.. _devpcpython:

*****************************
Python Package Development
*****************************

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
