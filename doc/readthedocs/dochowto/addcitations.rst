.. _addcitations:
 
********************************************
How to add scientific citations to the |PC|
********************************************

Scientific references are an essential part of the |PC| documentation. They link each
piece of physics implemented in the code to the original research papers where it was
first described or applied.  

To make these links automatic, we embed short BibTeX-style citations directly inside the
Fortran source files. During the documentation build, these references are collected and
used to generate a table that connects every module to its related scientific literature
(:doc:`scientific references per module <code/tables/papers>`).


How to add Scientic citations to the code
===========================================

To associate a new reference with a specific part of the |PC| code, you need to edit two
files:

* The Fortran source file (:file:`.f90`) corresponding to the module you want to cite,
  located in :file:`pencil-code/src/`.

* The bibliography file that contains all reference entries in BibTeX format:
  :file:`pencil-code/doc/citations/ref.bib`.


Fortran source
---------------

Inside the Fortran module, citations are added as commented lines following a specific
pattern. Each reference block must begin and end with the header markers shown below.

For example, in :file:`special/gravitational_waves_hTXk.f90`, the relevant section looks
like this:

.. code::

   !** AUTOMATIC REFERENCE-LINK.TEX GENERATION ********************
   ! Declare relevant citations from pencil-code/doc/citations/ref.bib for this  module.
   ! The entries are taken from pencil-code/doc/citations/notes.tex
   !
   ! 2020GApFD.114..130R,%RoperPol+ "The timestep constraint in solving the   gravitational wave equations sourced by hydromagnetic ..."
   ! 2020PhRvD.102h3512R,%RoperPol+ "Numerical simulations of gravitational   waves from early-universe ..."
   !
   !***************************************************************

To add your own citations, follow the same format:

.. code::

   !** AUTOMATIC REFERENCE-LINK.TEX GENERATION ********************
   ! Declare relevant citations from pencil-code/doc/citations/ref.bib for this  module.
   ! The entries are taken from pencil-code/doc/citations/notes.tex
   !
   ! bibcode "title"
   ! bibcode "title"
   !
   !***************************************************************


List of references
-------------------

Each citation used in the source code must also exist in the central bibliography file
:file:`pencil-code/doc/citations/ref.bib`, written in standard :command:`bibtex` format.

This file gathers all papers referenced across the |PC| source. When adding a new
citation, make sure the BibTeX entry includes the same *bibcode* you used in your
Fortran comments.

For example, the first reference in the previous code snippet corresponds to the
following entry:



.. code::

   @ARTICLE{2020GApFD.114..130R,
       author = {{Roper Pol}, Alberto and {Brandenburg}, Axel and {Kahniashvili}, Tina and
         {Kosowsky}, Arthur and {Mandal}, Sayan},
        title = "{The timestep constraint in solving the gravitational wave equations sourced by hydromagnetic turbulence}",
      journal = {Geophys. Astrophys. Fluid Dynam.},
     keywords = {Gravitational waves, early universe, aeroacoustics, turbulence},
         year = 2020,
        month = mar,
       volume = {114},
        pages = {130-161},
          doi = {10.1080/03091929.2019.1653460},
       adsurl = {https://ui.adsabs.harvard.edu/abs/2020GApFD.114..130R},
      adsnote = {Provided by the SAO/NASA Astrophysics Data System}
   }

Be careful to preserve the same *bibcode* identifier (`2020GApFD.114..130R` in this
example) between the Fortran file and the entry in :file:`ref.bib`. Otherwise, the link
between code and citation will not be recognized during the automatic documentation
build.


How does it work?
=================

When Sphinx builds the documentation, a Python helper script automatically parses all
Fortran files contained in :file:`pencil-code/src` (including subdirectories) and
searches for the marker line:

.. code::

   !** AUTOMATIC REFERENCE-LINK.TEX GENERATION ********************

Once the marker is found, the script skips three lines and starts reading the references
listed in the comment block, until it encounters another header line starting with
``!**``. Each line is expected to contain the *bibcode* of a paper declared in
:file:`pencil-code/doc/citations/ref.bib`.

The Python function that performs this operation is called
:func:`process_papers`. Its logic can be summarized as follows:

#. Use the UNIX command :command:`grep -rl "AUTOMATIC REFERENCE-LINK.TEX GENERATION"`
   to find all Fortran files that contain the citation marker.

#. Read each file, extract the bibcodes declared in the reference comment block, and
   store them in a dictionary where the key is the Fortran file name and the value is a
   list of bibcodes.

#. Parse the BibTeX database :file:`pencil-code/doc/citations/ref.bib` using
   :mod:`bibtexparser`.

#. For each bibcode found in the Fortran sources, match it against the BibTeX database.
   If a match exists, format it using the internal function :func:`format_paper`, which
   converts the BibTeX entry into a readable citation containing the author list, year,
   title, DOI, and ADS link.

#. Generate the final reStructuredText file
   :file:`code/tables/papers.rst` with a table that links every Fortran module to its
   corresponding scientific references.

Each entry in the table contains the name of the Fortran file (for example,
:file:`special/gravitational_waves_hTXk.f90`) and its associated references. The result
is then included automatically in the documentation as
:doc:`code/tables/papers`.
