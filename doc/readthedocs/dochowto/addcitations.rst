.. _addcitations:
 
********************************************
How to add scientific citations to the |PC|
********************************************


Scientific references are a fundamental part of the |PC| documentation.  
They establish a clear link between the implemented physical models and the
published research that introduced or validated them.  
Maintaining accurate citations ensures transparency, reproducibility, and proper
scientific credit.

.. important::

   All references used in the |PC| are centralized in a shared BibTeX file:
   :file:`pencil-code/doc/citations/ref.bib`

   This file is automatically used by the documentation system and by the
   Fortran module listings. 


Scientific usage of the |PC|
=============================

The recommended way to add a scientific citation that can be referenced both from
:ref:`_scientific_usage` and the :doc:`Fortran modules </code/sourceFortran/index>`
is to edit the central BibTeX file directly.

Follow these steps:

#. Open the BibTeX file: 
   :file:`pencil-code/doc/citations/ref.bib`


#. Add your paper:
   Copy the BibTeX entry from a reliable source (e.g., NASA ADS) and paste it
   at the end of the file.


   .. important::

      Before adding a new BibTeX entry, always check if the paper already exists
      in :file:`pencil-code/doc/citations/ref.bib`.  
      Duplicate entries with different keys cause confusion and may lead to broken
      cross-references in the documentation.

      Use consistent and descriptive BibTeX keys, ideally following the ADS bibcode
      format (for example: ``2020PhRvD.102h3512R``).  
      This ensures automatic linking between Fortran modules and the corresponding
      scientific references.

   .. important::

      If you initially add an arXiv version of your paper and later want to replace
      it with the published peer-reviewed version, **update the existing entry**
      in :file:`pencil-code/doc/citations/ref.bib` instead of adding a new one.  

      Creating duplicate entries for the same paper (arXiv and published) can break
      cross-references and lead to inconsistent citations throughout the
      documentation.

#. Assign a topical tag: 
   Add the custom :command:`pcsection` field at the end of the entry to classify
   the paper into one or more topics or subtopics listed in :ref:`papers_by_topic`.

   Example with one topic:

   .. code::

      @ARTICLE{2004IJAsB...3..209B,
            author = {{Brandenburg}, A. and {Multam{\"a}ki}, T.},
             title = "{How long can left and right handed life forms coexist?}",
           journal = {Int. J. Astrobiology},
            eprint = {q-bio/0407008},
          keywords = {exobiology, homochirality, origin of life.},
              year = 2004,
             month = jul,
            volume = 3,
             pages = {209-219},
               doi = {10.1017/S1473550404001983},
            adsurl = {http://adsabs.harvard.edu/abs/2004IJAsB...3..209B},
           adsnote = {Provided by the SAO/NASA Astrophysics Data System},
         pcsection = {Turbulent combustion, front propagation, radiation \& ionization}
      }

   Example with multiple topics:

   .. code::

      @INPROCEEDINGS{2004IAUS..223...57B,
            author = {{Brandenburg}, A. and {Sandin}, C. and {K   {\"a}pyl{\"a}}, P.~J.
         	},
             title = "{Helical coronal ejections and their role   in the solar cycle}",
         booktitle = {Multi-Wavelength Investigations of Solar    Activity},
              year = 2004,
            series = {IAU Symp.},
            volume = 223,
            eprint = {astro-ph/0407598},
            editor = {{Stepanov}, A.~V. and {Benevolenskaya}, E.  ~E. and {Kosovichev}, A.~G.
         	},
             pages = {57-64},
               doi = {10.1017/S1743921304005101},
            adsurl = {http://adsabs.harvard.edu/abs/2004IAUS.. 223...57B},
           adsnote = {Provided by the SAO/NASA Astrophysics Data  System},
         pcsection = {{Large-scale dynamos, helical turbulence,   and catastrophic quenching}, {Hydrodynamic and MHD         instabilities}}
      }

#. Commit and push your changes:
   Once your update is pushed to the main repository, the
   :ref:`_scientific_usage` section will be automatically rebuilt and include
   your new paper.

List of tags 
------------

The complete list of available topics and subtopics is defined in the Python
file :file:`fortran_rst_generator.py`.  
These entries determine the categories used by the automatic documentation.



.. code:: python

       all_topics = {
        "Interstellar and intercluster medium as well as early Universe": [
            "Interstellar and intercluster medium",
            "Small-scale dynamos and reconnection",
            "Primordial magnetic fields and decaying turbulence",
            "Relic gravitational waves \\& axions"
        ],
        "Planet formation and inertial particles": [
            "Planet formation",
            "Inertial, tracer particles, \\& passive scalars"
        ],
        "Accretion discs and shear flows": [
            "Accretion discs and shear flows",
            "Shear flows",

        ],
        "Solar physics": [
            "Coronal heating and coronal mass ejections",
            "Large-scale dynamos, helical turbulence, and catastrophic quenching",
            "Helioseismology",
            "Strongly stratified MHD turbulence and NEMPI",
            "Convection in Cartesian domains",
            "Global convection and dynamo simulations"
        ],
        "Miscellanea": [
            "Turbulent transport and test-field method",
            "Hydrodynamic and MHD instabilities",
            "Chiral MHD",
            "Hydrodynamic and MHD turbulence",
            "Turbulent combustion, front propagation, radiation \\& ionization",
            "Code development, GPU etc"
        ]
    }

You can add new topics or subtopics by editing this dictionary as needed.



Adding citations anywhere in the documentation
-----------------------------------------------

All citations in :file:`ref.bib` are available to every document in the
|PC| documentation.  
They can be referenced using the same syntax as in LaTeX, with the
roles provided by the `sphinxcontrib-bibtex <https://sphinxcontrib-bibtex.readthedocs.io/en/latest/usage.html#roles-and-directives>`_ extension.

* Parenthetical citation

   .. code::

      :cite:p:`2004IJAsB...3..209B`

   → produces :cite:p:`2004IJAsB...3..209B`

* Textual citation

   .. code::

      :cite:t:`2004IJAsB...3..209B`

   → produces :cite:t:`2004IJAsB...3..209B`


How the scientific usage pages are updated
---------------------------------------------------------------------

The automatic update of the :ref:`scientific_usage` section is handled by the
Python function :func:`process_papers`, defined in the file
:file:`fortran_rst_generator.py`, available in
:file:`~/pencil-code/doc/readthedocs`.

This function collects all scientific references related to the |PC| and
organizes them into several reStructuredText (``.rst``) files for inclusion in
the documentation. It is called during the documentation build process to ensure
that all citation pages remain up to date.

Specifically, :func:`process_papers` performs the following steps:

#. **Generates publication statistics and plots.**  
   It first calls :func:`plot_papers_per_year`, which analyses all BibTeX entries
   from :file:`doc/citations/ref.bib` and counts the number of |PC|-related papers
   per year.  
   The function distinguishes between:
   
   * papers **co-authored by Brandenburg** (black line),

   * papers **without Brandenburg** (blue line),

   * papers citing the |PC| for **code comparison or reference** (green line), and

   * **unclassified** papers that cannot be assigned to any category (red line).

   The results are visualized as a step plot saved to
   :file:`_images/papers_per_year.png`.  
   This plot is automatically included in the documentation (see
   :numref:`Fpapers`) and helps track the yearly evolution of |PC|-related
   publications.

#. **Creates yearly citation summaries.**  
   Using the same data, :func:`process_papers` writes the file
   :file:`code/papers_by_year.rst`, listing the number of papers per year and
   their corresponding :rst:role:`cite:t` references.  
   The statistics also report the total number of papers and their percentage
   breakdown by category.

#. **Groups papers by research topic.**  
   The function generates :file:`code/papers_by_topic.rst`, where publications are
   organized into predefined research themes such as *Solar physics*, *Planet
   formation*, and *Accretion discs and shear flows*.  
   Each topic includes alphabetically labelled subtopics, listing all relevant
   citations sorted by publication year.

#. **Lists code-comparison and unclassified papers.**  
   Two additional files are created:
   :file:`code/papers_code.rst`, which contains papers citing the |PC| for
   benchmarking or code-comparison purposes, and
   :file:`code/papers_unclassified.rst`, which collects all entries that could not
   be automatically assigned to a topic.

The resulting files are included directly in the documentation using the
:rst:dir:`include` directive. This ensures that the
:ref:`_scientific_usage` section always reflects the most recent and complete set
of publications referencing the |PC|, without the need for manual editing.



Scientific citations to individual Fortran modules
===================================================

In addition to tracking how the |PC| is used across the scientific literature
(:ref:`_scientific_usage`), each individual Fortran module can also be linked to
its corresponding research papers.  
These **module-level citations** provide a direct connection between the physics
implemented in the code and the original scientific work that introduced or
validated it.

To make these links automatic, short BibTeX-style citation keys are embedded
directly inside the Fortran source files. During the documentation build, these
references are extracted and combined with the entries in
:file:`pencil-code/doc/citations/ref.bib` to populate the third column of the
:doc:`Fortran modules </code/sourceFortran/index>` table.  
This table lists every module together with its description and the associated
scientific references, ensuring transparent traceability between code and
literature.

To associate a new reference with a specific part of the |PC| code, you need to
edit two files:

* The Fortran source file (:file:`.f90`) corresponding to the module you want to
  cite, located in :file:`pencil-code/src/`.
  
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
------------------

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
:func:`process_file`. Its logic can be summarized as follows:

#. First process all Fortran files in the :file:`pencil-code/src` and extracts file names 

#. Add the module description for the second column of the table 

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
is then included automatically in the documentation as the third column in 
:doc:`/code/sourceFortran/index`.
