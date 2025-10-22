.. _fortranautodochowto:
 
**********************************************
How to generate autodocumentation for Fortran
**********************************************


Sphinx (and by extension, ReadTheDocs) speaks fluent :code:`Python`. Python
modules can be documented out-of-the-box, with minimal effort.

Fortran, on the other hand, is not a native Sphinx citizen. To bring our
|PC| Fortran modules into the documentation realm, we rely on community-developed Sphinx extensions, with a few tweaks of our own.


Extensions
===============

The key Sphinx extensions for Fortran are:

* `sphinxfortran.fortran_autodoc
  <https://sphinx-fortran.readthedocs.io/en/latest/lib.autodoc.html#module-sphinxfortran.fortran_autodoc>`_
* `sphinxfortran.fortran_domain
  <https://sphinx-fortran.readthedocs.io/en/latest/lib.domain.html>`_

Both are documented in the `Sphinx Fortran <https://sphinx-fortran.readthedocs.io/en/latest/index.html>`__ ReadTheDocs webpage.

After trial and error, we realized we needed to adapt these libraries for
|PC|. Our modified versions live in the :file:`readthedocs/_ext` directory:

* :file:`_ext/fortran_autodoc.py`
* :file:`_ext/fortran_domain.py`


Configuration
==============


Generating Fortran autodocs requires more than just adding an extension
to :file:`conf.py`. Here’s how we configure Sphinx for |PC|.


Configuration of ``conf.py``
------------------------------

1. **Add the extensions**

   Include ``fortran_domain`` and ``fortran_autodoc`` in the extensions
   list. The relevant section of :file:`conf.py` now looks like this:


    .. code:: python

        # Add any Sphinx extension module names here, as strings. They can be
        # extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
        # ones.
        extensions = [
            "autoapi.extension",
            "sphinx.ext.autosectionlabel",
            "sphinx.ext.autodoc",
            "sphinx.ext.autosummary",
            "sphinx.ext.doctest",
            "sphinx.ext.todo",
            "sphinx.ext.mathjax",
            "sphinx.ext.ifconfig",
            "sphinx.ext.viewcode",
            "sphinx.ext.githubpages",
            "sphinx.ext.napoleon",
            "sphinxcontrib.images",
            "fortran_domain",
            "fortran_autodoc",
            "sphinx.ext.intersphinx",
        ]

2. **Add extra Fortran configuration**

    .. code:: python

        ## Fortran configuration
        fortran_src = create_fortran_modules_rst("../../src")
        fortran_ext = ["f90"]



Extra files
-------------


The function :command:`create_fortran_modules_rst` is implemented in
:file:`fortran_rst_generator.py`, located in the :file:`readthedocs`
directory.

This script:

* Processes each ``.f90`` file to extract module-level comments.

* Generates nice tables for the :ref:`fortran_modules`.

* Defines a ``FILES_THAT_DONT_WORK`` list to exclude problematic files.

* Creates individual :file:`.rst` files containing the ``.. f:autosrcfile::`` Sphinx directive,
  which during the Sphinx build create one html page with automatically-generated documentation per ``.f90`` file.

* It handles the creation of the table-of-contents entries in the sidebar.



How to document your Fortran code
===================================


Why is this important?
-----------------------

Once configured, autodocumentation worked for the majority of modules.
Out of 541 files, only 24 caused issues — not bad at all!

Here is the preliminary list:



.. code:: python

    FILES_THAT_DONT_WORK = [
        "src/chemistry.f90", # 'ascii' codec can't decode byte 0xcb in position 7543: ordinal not in range(128).
        "src/diagnostics.f90", # ø character in file
        "src/diagnostics_outlog.f90", # ø character in file
        "src/forcing.f90", # ± character in file
        "src/fourier_fftpack.f90", # § character in file
        "src/initcond.f90", # 'ascii' codec can't decode byte 0xe2 in position 6496: ordinal not in range(128).
        "src/nosolid_cells.f90", # ø character in file
        "src/particles_chemistry.f90", # 'ascii' codec can't decode byte 0xe2 in position 613: ordinal not in range(128).
        "src/particles_dust.f90", # 'ascii' codec can't decode byte 0xe2 in position 2129: ordinal not in range(128).
        "src/solid_cells.f90", # ø character in file
        "src/solid_cells_ogrid.f90", # ø, é characters in file
        "src/solid_cells_ogrid_mpicomm.f90", # é character in file
        "src/timestep_rkf_lowsto.f90", # í character in file
        "initial_condition/1D_loop_init.f90", # 'ascii' codec can't decode byte 0xc4 in position 6280: ordinal not in range(128).
        "initial_condition/alfven_wave.f90", # 'ascii' codec can't decode byte 0xc3 in position 39: ordinal not in range(128).
        "initial_condition/coronae_init.f90", # 'ascii' codec can't decode byte 0xc4 in position 3033: ordinal not in range(128).
        "special/streamfunction_fullmultigrid.f90", # 'ascii' codec can't decode byte 0xcf in position 5737: ordinal not in range   (128).
        "special/streamfunction_multigrid.f90", # 'ascii' codec can't decode byte 0xe2 in position 4021: ordinal not in range(128).
        "test_methods/testfield_xy.f90", # 'ascii' codec can't decode byte 0xc3 in position 361: ordinal not in range(128).
        "src/io_dist.f90", # Found non-(space,digit) char in the first column. Are you sure that this code is in fix form?  line='kloop:do kk=kka,kke ')
        "src/hydro.f90", # UNKNOWN ERROR
        "src/polynomialroots.f90", # CRITICAL: Unexpected section title or transition.
        "src/slices.f90", # (exception: '=')
        "src/sub.f90", # (exception: expected string or bytes-like object, got 'NoneType')
    ]


Most issues are caused by non-ASCII characters, which are easily avoidable.
However, this illustrates why having clear guidelines and some preprocessing
is useful before attempting automatic Fortran documentation.


Removing ascii problems
^^^^^^^^^^^^^^^^^^^^^^^^

Parsing the file to get the line and specific problem:

.. code:: bash

    ~/pencil-code/src$ grep --color='auto' -nP '[^\x00-\x7F]' chemistry.f90 
    186:  real :: conc_sat_spec_cgs=1e-8 !units of mol/cmˆ3

Change line.

Or use the provided Python script: :file:`find_nonascii.py` to find the files, line, and columns of the non-ascii characters.


After solving ascii issues, the list of non-working files was reduced considerably:

.. code:: python

    FILES_THAT_DONT_WORK = [
        "src/io_dist.f90", # Found non-(space,digit) char in the first column. Are you sure that this code is in fix form?  line='kloop:do kk=kka,kke ')
        "src/hydro.f90", # UNKNOWN ERROR
        "src/polynomialroots.f90", # CRITICAL: Unexpected section title or transition.
        "src/slices.f90", # (exception: '=')
        "src/sub.f90", # (exception: expected string or bytes-like object, got 'NoneType')  
    ]

