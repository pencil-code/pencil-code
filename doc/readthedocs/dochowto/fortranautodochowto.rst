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


Besides the coding style rules described in :ref:`pencilstyleguide`, it is
important to follow a few additional conventions to ensure that the automatic
Fortran documentation builds correctly and consistently. Think of these as
rules for keeping the autodoc engine from screaming into the void.


Basic rules
------------

#. **Use only ASCII characters.**

   Non-ASCII symbols (like `ø`, `é`, or `±`) may look innocent, but they tend
   to break Sphinx’s parsing and the Fortran autodoc extensions. If you need a
   special symbol, spell it out in plain text (for example, write ``sigma``
   instead of `σ`). You can find and fix these characters using the helper
   script :file:`find_nonascii_dir.py`.

#. **Avoid reStructuredText headers in comments.**

   Avoid long runs of ``====`` or ``----`` in your Fortran comments; Sphinx may
   mistake them for section headers and crash dramatically. A few dashes or
   equals signs are fine, but if you can fence them off with extra comment
   markers or shorten them, even better. Check :ref:`rst-headings`.

#. **Do not indent doc comments more than the code.**

   The autodoc parser associates a comment with the following entity (module,
   function, variable) only if the indentation matches.

See `Optimize the process <https://sphinx-fortran.readthedocs.io/en/latest/user.autodoc.html#optimize-the-process>` for further reference. 



Following these conventions will keep your documentation pipeline smooth and
prevent mysterious parsing errors that look like temporal anomalies.


Why is this important?
-----------------------

Once properly configured, autodocumentation worked for the majority of modules.
Out of **541** Fortran files, only **24** caused issues — not bad at all!

Here’s the preliminary list of problematic files:



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


Most issues were caused by **non-ASCII characters**, which are easily avoidable.
Still, this illustrates why having **clear documentation rules** and some
**preprocessing scripts** is invaluable before attempting automatic
Fortran documentation.

Removing ascii problems
^^^^^^^^^^^^^^^^^^^^^^^^

To locate non-ASCII characters within a file, you can use:

.. code:: bash

    ~/pencil-code/src$ grep --color='auto' -nP '[^\x00-\x7F]' chemistry.f90
    186:  real :: conc_sat_spec_cgs=1e-8 !units of mol/cmˆ3

Then simply edit the offending line to replace or remove the problematic symbol.

Alternatively, use the provided Python script
:file:`find_nonascii.py`, which lists all files containing non-ASCII characters
along with their exact line and column positions.

After removing the non-ASCII issues, the list of problematic files
was reduced considerably:

.. code:: python

    FILES_THAT_DONT_WORK = [
        "src/io_dist.f90",
        "src/hydro.f90",
        "src/polynomialroots.f90",
        "src/slices.f90",
        "src/sub.f90",
    ]

Common Remaining Issues
^^^^^^^^^^^^^^^^^^^^^^^^


``src/io_dist.f90``
""""""""""""""""""""""

*Error:*

.. code::

    | Handler <function fortran_parse at 0x76498181acb0> for event 'builder-inited' threw an exception (exception: readfortrancode: Found non-(space,digit) char in the first column.
    pencil-sphinx  |        Are you sure that this code is in fix form?
    pencil-sphinx  |        line='kloop:do kk=kka,kke ')
    pencil-sphinx  | parsing fortran sources.../app/pencil-code/src/io_dist.f90

*Fix:*

The problem was actually in the included file ``src/io_dist.h``.

Add a **leading space** before label-starting lines, for example:

.. code::

    "kloop:do kk=kka,kke"  →  " kloop:do kk=kka,kke"
    "jloop: do jj=jja,jje" →  " jloop: do jj=jja,jje"
    "iloop: do ii=iia,iie" →  " iloop: do ii=iia,iie"

This prevents the autodoc parser from interpreting the file as
fixed-form Fortran (where column 1 must be blank or numeric).
Since ``io_dist.h`` is a free-form include, this small change
does not alter semantics but avoids parser confusion.


``src/polynomialroots.f90``
""""""""""""""""""""""""""""


*Error:*



.. code:: text


    10: CRITICAL: Missing matching underline for section title overline.
    108: CRITICAL: Unexpected section title or transition.
    111: WARNING: Block quote ends without a blank line; unexpected unindent. [docutils]
    122: CRITICAL: Unexpected section title or transition.
    125: WARNING: Block quote ends without a blank line; unexpected unindent. [docutils]
    137: CRITICAL: Unexpected section title or transition.
    152: CRITICAL: Unexpected section title or transition.

*Fix:*

Replace all comment-only separator lines (e.g. long ``-----`` or ``====``)
with shorter or descriptive comments, for example:

Sphinx interprets long punctuation lines in the generated RST
as section overlines or transitions. If the lengths mismatch or
the structure is ambiguous, parsing fails.  
Shorter, plain comments prevent RST misinterpretation.

``src/slices.f90``
""""""""""""""""""""""""

*Error:*


.. code:: text
    
    src/sub.f90
    pencil-sphinx  | =type(1:len(trim(type))-1) None
    pencil-sphinx  | ::real3 <re.Match object; span=(0, 7), match='::real3'>
    pencil-sphinx  | make: *** [Makefile:41: html] Error 2



*Fix:*  

Not yet identified. Possibly caused by a malformed type declaration or
nonstandard Fortran syntax that confuses the parser.




``src/hydro.f90``
""""""""""""""""""""

*Error:*


.. code:: text

    Extension error (fortran_autodoc):
    pencil-sphinx  | Handler <function fortran_parse at 0x718c73c0ec20> for event 'builder-inited' threw an exception (exception: (2, '(lsh*(lsh', '(lsh*(lsh'))
    pencil-sphinx  | parsing fortran sources.../app/pencil-code/src/hydro.f90
    pencil-sphinx  | make: *** [Makefile:41: html] Error 2

*Fix:*  

Still under investigation.



``src/sub.f90``
""""""""""""""""""


*Error:*


.. code:: text

    src/sub.f90
    pencil-sphinx | =type(1:len(trim(type))-1) None
    pencil-sphinx | ::real3 <re.Match object; span=(0, 7), match='::real3'>
    pencil-sphinx | make: *** [Makefile:41: html] Error 2

*Fix:*  
Not yet resolved. 


Even with these few stubborn cases, the autodoc coverage rate is above **99%**.
In other words, most of the universe is in order — only a few
files still resist documentation, like rogue planets evading a gravity well.
