.. _usingrst:

**********************
Using reStructuredText
**********************

All documentation available on this page is written in the reStructuredText
(reST) markup language.


.. contents:: 
    :depth: 3



About reST
==========

reST is a simple markup language for plain text files. It is used to
semantically mark parts of a document (e.g., section headings, bold text, lists,
source code, etc.) for further processing and uniform rendering.

Sphinx is a documentation generator that, in our case, achieves two goals:

  - processes the reST documents and renders them to HTML and PDF;
  - autogenerates the code documentation for Python, IDL and Fortran projects
    (i.e., it extracts and formats lists of classes, functions, etc.,
    each with descriptions based on comments found in the source code).

While the `reST <https://docutils.sourceforge.io/rst.html>`_
(and `Sphinx <https://www.sphinx-doc.org/en/master/contents.html>`_)
official documentation pages are exhaustive, they are perhaps not recommended
for a beginner, as they necessarily contain a lot of information that is
not relevant for our documentation page.
We suggest starting with https://rest-sphinx-memo.readthedocs.io/en/latest/ReST.html,
which is a quick reference for reST and Sphinx that was specifically created
to cover a small subset of features that are likely to be used on a daily basis.

Another good source of reST commands, is the CheatSheet writen by Thomas Cokelaer
https://thomas-cokelaer.info/tutorials/sphinx/rest_syntax.html.

Converting existing documents
=============================

The utility *pandoc* can be used to convert a variety of formats, including
Microsoft Word (*doc*, *docx*), Libre Office (*odt*), LaTeX,
HTML, XML, and if all else fails even PDF, to reST.

The syntax of the command is:

.. code:: bash

   pandoc input.doc -o output.rst

where *input.doc* is your input document, in any format other than *rst*.

Style guideline
===============

Rules
-----

Basic rules:

* reST is sentitive to identation,
* reST requires blanck lines between paragraphs,
* reST is sentitive to spaces and blank lines.

Special characters
------------------

Special characters used to format text:

=========== ================================== ==============================
usage          syntax                           HTML rendering
=========== ================================== ==============================
italic      `*italic*`                         *italic*
bold        `**bold**`                         **bold**
link        ```python <www.python.org>`_``     `python <www.python.org>`_
verbatim    ````*````                               ``*``
code        ````text````                       ``code``
=========== ================================== ==============================

The docuble backquote is used to enter in verbatim mode, and can be sue as the escaping character. 

Restrictions about  the ``*`` and  the `````` syntax:

   * cannot not be nested,
   * cannot start with space: ``* text*`` is wrong,
   * must be separated from surrounding text by a space.

Other semantic markup¶
^^^^^^^^^^^^^^^^^^^^^^

You can use additional reST roles to provide semantic meaning to the text. 

Abbreviations
   ``:abbr:`LSST (Legacy Survey of Space and Time)``` → :abbr:`LSST (Legacy Survey of Space and Time)` (a tool tip exposes the definition)

Filenames and paths
   ``:file:`repos.yaml``` → :file:`repos.yaml`

Shell commands
   ``:command:`git rebase -i main``` → :command:`git rebase -i main`

User interface labels
   ``:guilabel:`New Pull Request``` → :guilabel:`New Pull Request`. This markup can be used for button labels, menus, or even text labels in interactive shell programs.

Keyboard commands
   ``:kbd:`Control-a s``` → :kbd:`Control-a s`. Spell out the keys rather than using Emacs short hand, such as ``C-x``.

To semantically markup Python or C++ code objects, refer to the section on :ref:`rst-code-link`.


Headings
--------

In reST, headings are marked by underlining them with the same character, at least as long as the text:

.. code:: rst

   This is a heading
   =================

In this documentation, the following markers should be used, following the convection used in `Python Developer's Guide for documenting <https://devguide.python.org/documentation/markup/#sections>`_

.. code:: rst

   ####                   
   Part                      
   ####     

   ***********************
   Chapter (title of page)
   *********************** 
  
   Section
   =======

   Subsection
   ----------

   Sub-subsection
   ^^^^^^^^^^^^^^


   Paragraph
   """""""""



You should not use further levels of headings, as it would prevent optimal
rendering of the table of contents in the left-hand sidebar. You can structure
your document further by using the ``.. rubric::`` directive.

Do not use ``###`` and ``***``, as they are already used for higher-level headings
(e.g., on the main landing page).

Directives
----------

ReST syntax provides directives to include formatted text. 


Insert a code block
^^^^^^^^^^^^^^^^^^^

Example: insert an example of ReST code:

.. code:: rst

   .. code:: rst

      *bold*

produces:

.. code:: rst

   *bold*

In this case *rst* is an argument telling that the code is ReST. 

Instead of *code*, other directives are admonitions, or insert image, video, etc.



Admonitions
-----------

The use of admonition directives can greatly enhance the user experience by
presenting tips, warnings, important notes, etc. in a way that stands out from
the rest of the document.

The following admonition directives are available: ``attention``,
``caution``, ``danger``, ``error``, ``hint``, ``important``, ``note``, ``tip``, ``todo``,
``warning``, ``seealso``.

Any of the previous values can be used as follows:

.. code:: rst

   .. note::

      This is a note.

producing the following output:

.. note::

   This is a note.

Keep in mind that overuse of admonitions will detract from the
document flow too much, and consequently worsen the user experience.
**Use them sparingly.**

There is also the possibility of definining your own admonition:

.. code:: rst

   .. admonition:: Do not forget!!

      something important that I don't rememeber anymore

produces:

.. admonition:: Do not forget!!

      something important that I don't rememeber anymore

Internal and external links
---------------------------

External links
^^^^^^^^^^^^^^

* Simple link to a website:

   .. code:: rst

      `<http://www.python.org/>`_

   produces:

   `<http://www.python.org/>`_

* Link with label:

   .. code:: rst

      `Python <http://www.python.org/>`_

   produces:

   `Python <http://www.python.org/>`_

.. note::

   If you have an underscore, escape it with '\\'.


Implicit links to titles
^^^^^^^^^^^^^^^^^^^^^^^^

All titles are considered hyperlinks. You can link to any title by suing the same quotes and final underscore as above:

.. code:: rst

   `Internal and external links`_

produces:

`Internal and external links`_

But only works if the title and the link are in the same ReST file. Otherwise, you need to use the `Explicit links`_ format.

Explicit links
^^^^^^^^^^^^^^

You can create explicit links between reST files by creating a label, like:

.. code:: rst

   .. _usingrst:

(This label was added at the beginning of this document).

And then refer to the label using one of the following methods:

#. Method:

   .. code:: rst

      usingrst_

   produces: usingrst_

#. Method: 

   .. code:: rst

      :ref:`usingrst`
   
   uses the first title's name after the link, so here you will see: :ref:`usingrst`. 

   You can only use this method if the link is found in an external reST file. 


Footnotes
---------


.. code::

   This is a line.\ [#label]_

   .. [#label] This is the footnote content.



List and bullets
----------------

* Bulleted list: 
   
   the code:

   .. code:: rst

      * Bulleted list

      * also in the list

      - continues the list
      - another item

         - sublist
         * or this

   produces:

   * Bulleted list

   * also in the list

   - continues the list

   - another item

      - sublist

      * or this

* Numbered list:

   the code:

   .. code:: rst
      
      1. first item

      2. second item

      #. also works

      #. n-th item

         #. sublist
         #. sublist 

   produces:

   1. first item
   2. second item

   #. also works
   #. n-th item

      #. sublist
      #. sublist 

The only important part to distinguish between different levels on the lists is the identation, not the blanck lines in between.


Tables
------

There are several ways to write tables, but their rendering depends on the CSS/HTML style, not on sphinx itself.

#. Simple table

   .. code:: rst

      +---------+---------+---------+
      | Title1  | Title2  | Title3  |
      +=========+=========+=========+
      | 1       |  2      |  3      |
      +---------+---------+---------+

   producing the output:

   +---------+---------+---------+
   | Title1  | Title2  | Title3  |
   +=========+=========+=========+
   | 1       |  2      |  3      |
   +---------+---------+---------+


#. Adjusting the size of the cells:

   .. code:: rst

      +----------------+---------+-------+
      | Title1         | Title2  | Title3|
      +================+=========+=======+
      | 1              |       2 |  3    |
      +----------------+---------+-------+

   producing the output:

   +----------------+---------+-------+
   | Title1         | Title2  | Title3|
   +================+=========+=======+
   | 1              |       2 |  3    |
   +----------------+---------+-------+

#. A simplified version with multiple cells:

   .. code:: rst

      ====  ====  =======  =======
      Title1     Title2   Title3 
      ----------  -------  -------
      A      B    
      ====  ====  ======  ======
      1a    1b    2       3 
      ====  ====  ======  ======

   producing the output:

   ====  ====  =======  =======
     Title1    Title2   Title3 
   ----------  -------  -------
   A      B    
   ====  ====  =======  =======
   1a    1b    2        3 
   ====  ====  =======  =======

#. The previous formats may give problems with LaTeX, since the column width is difficult to compute automatically. Use the following directive if you are outputing LaTeX documents:

   .. code:: rst

      .. tabularcolumns:: column spec

   Example:

   .. code:: rst

      .. tabularcolumns:: |l|c|p{5cm}|

      +--------------+---+-----------+
      |  simple text | 2 | 3         |
      +--------------+---+-----------+

   which produces:

   .. tabularcolumns:: |l|c|p{5cm}|

   +--------------+---+-----------+
   |  simple text | 2 | 3         |
   +--------------+---+-----------+

#. Another option is the ``list-table`` directive, which creates a table from data in a uniform two-level bullet list. “Uniform” means that each sublist (second-level list) must contain the same number of list items.

   Example:

   .. code:: rst

      .. list-table:: Frozen Delights!
         :widths: 15 10 30
         :header-rows: 1

         * - Treat
           - Quantity
           - Description
         * - Albatross
           -  2.99
           - On a stick!
         * - Crunchy Frog
           - 1.49
           -  If we took the bones out, it wouldn't be crunchy, now would it?
         * - Gannet Ripple
           - 1.99
           - On a stick!

To add a label for reference, use the normal reference before the table: ``.. _labeloftable:`` 


Images and figures
----------------------


Three different directives allow for the addition images in the documentation.
Please, see `this guide <https://docutils.sourceforge.io/docs/ref/rst/directives.html#images>`_ 
for a full description.

#. The simplest one is the ``image`` directive:

   .. code:: rst

      .. image:: pics/myimage.png

   Accepted options for the directive are the  width and alternative text for screen readers:

   .. code:: rst

      .. image:: pics/myimage.png
         :width: 400
         :height: 100px
         :scale: 50 %
         :alt: alternate text
         :align: right
      

#. The ``figure`` directive supports all the options of the ``image`` directive and  allows for adding a caption to the figure:
   
   .. code:: rst

      .. figure:: pics/myimage.png
         :scale: 50 %
         :alt: Flow patterns in the Sun

         This is the caption of the figure (a simple paragraph).

         This is the legend of the figure, which can include a table:

         +-----------------------+-----------------------+
         | Symbol                | Meaning               |
         +=======================+=======================+
         | .. image:: arrow.png  | Magnetic field lines |
         +-----------------------+-----------------------+
         | .. image:: lines.png  | Velocity lines        |
         +-----------------------+-----------------------+
   
   There must be blank lines before the caption paragraph and before the legend. 
   To specify a legend without a caption, use an empty comment (“..”) in place of the caption.
 
   If you want to add a label to the figure, just use the option :command:`:name:`

   .. code:: rst

      .. figure:: pics/myimage.png
         :name: solarimage
         :scale: 50 %
         :alt: Flow patterns in the Sun

         And write the caption to the figure.


#. The ``thumbnail`` directive allows you expand the image by clicking on it:

   .. code:: rst
   
      .. thumbnail:: pics/myimage.png
         :width: 500px


Videos
------

You can add short movies to your documentation by using the ``.. video::``
directive. Any video that works inside an HTML5 *video* tag can be used (i.e.,
mp4, webm, ogg). Follow these steps to add your video:

- Add the ``.. video:: <video_url>`` directive in your rst file,
  where you want the video to be rendered.
- It is not necessary to specify any options (height, width, etc.), but if
  you want to have a look at the documentation of the extension:
  https://github.com/sphinx-contrib/video

This is the recommended way of adding videos, since they should not
be committed to the *ingdoc* git repository, but rather stored on a
separate server.

However, if you absolutely need to store the video with the documentation,
follow these steps instead:

- Copy the video file to the directory ``_static``. This is necessary at the
  moment, since we have not found a way (yet) for Sphinx to deploy the file
  otherwise.
- Add the ``.. video:: <relative_path_to_video>`` directive in your rst file,
  where you want the video to be rendered. The path is relative to your rst file,
  so it will probably look similar to ``../_static/video.mp4``.

Linking to External Docs (Intersphinx)
===============================================


Sometimes, our documentation needs to refer to things that live outside the Pencil Code universe — like Python functions, NumPy arrays, or SciPy routines.  
Instead of manually typing full URLs (and then forgetting to update them when the external docs change), Sphinx offers a smarter way: **Intersphinx**.

This allows our documentation to automatically link to objects in other Sphinx-based projects.

What it Does
------------

Intersphinx creates a bridge between your documentation and another project’s documentation.  
When you write something like

.. code:: rst

   :class:numpy.ndarray

Sphinx looks up ``numpy.ndarray`` in a remote *inventory file* (``objects.inv``) provided by the NumPy documentation,  
and automatically turns it into a working hyperlink.

No manual linking. No URL maintenance. No fuss.

How to Enable it
----------------

Intersphinx is already enable for this documentation. 

The extension was added to ``conf.py`` file:

.. code-block:: python

   extensions = [
       'sphinx.ext.intersphinx',
       # other extensions...
   ]

Then, it is defined where the external documentation lives using the ``intersphinx_mapping`` dictionary:

.. code-block:: python

   intersphinx_mapping = {
       'python': ('https://docs.python.org/3', None),
       'numpy': ('https://numpy.org/doc/stable/', None),
       'scipy': ('https://docs.scipy.org/doc/scipy/', None),
       'matplotlib': ('https://matplotlib.org/stable/', None),
   }

Each key (like ``numpy`` or ``scipy``) becomes the *prefix* for cross-references.


Example usage
-------------

Now, inside your ``.rst`` files, you can simply write:

.. code-block:: rst

   The simulation data are stored as :class:`numpy.ndarray` objects,
   and can be manipulated using :func:`numpy.mean` or :func:`scipy.signal.convolve`.

When the documentation is built, these references automatically link to the correct pages in the NumPy and SciPy manuals, generating:


The simulation data are stored as :class:`numpy.ndarray` objects, and can be manipulated using :func:`numpy.mean` or :func:`scipy.signal.convolve`.


Advanced Tips
--------------

* You can link between **your own documentation projects**, too.  
  Just copy the generated ``objects.inv`` file from one project to another and add it to the mapping.

* Sphinx caches the inventories locally, so you don’t need internet access for every build.

* To link to a specific domain, use the full role syntax, for example:
  ``:py:func:`numpy.mean``` or ``:py:class:`numpy.ndarray```.


Intersphinx is your documentation’s long-distance calling plan:
it connects your project to the rest of the Python (and Pencil Code) ecosystem,
so your references stay alive, up-to-date, and perfectly linked — even across galaxies of documentation.



Writing equations in LaTeX format
=================================

Luckily, writing equations in LaTeX format is supported natively in Sphinx, see the `official Sphinx Math Documentation <https://www.sphinx-doc.org/en/master/usage/restructuredtext/directives.html#math>`_ , using the directive `math`:

.. code:: rst

   .. math::

   \frac{ \sum_{t=0}^{N}f(t,k) }{N}

creates:

.. math::

   \frac{ \sum_{t=0}^{N}f(t,k) }{N}


Or you can use the same directive inline:

.. code:: rst

   the equation :math:`\frac{ \sum_{t=0}^{N}f(t,k) }{N}` gives blabla

which generates:


the equation :math:`\frac{ \sum_{t=0}^{N}f(t,k) }{N}` gives blabla


.. important:: 

   Do not forget the delimiting backticks \`\`.


Some useful options are:

.. code:: rst

   : name: label (text)

An implicit target name that can be referenced using ``ref``.

.. code:: rst

   : label: label (text)

With this option the equation will get a number, and the equation can be referenced by its number using ``:math:numref:`euler```
The default is that equations are not numbered.


Example:

.. code:: rst

   .. math:: e^{i\pi} + 1 = 0
   :label: euler

   Euler's identity, equation :math:numref:`euler`, was elected one of the
   most beautiful mathematical formulas.

Produces

.. math:: e^{i\pi} + 1 = 0
   :label: euler

Euler's identity, equation :math:numref:`euler`, was elected one of the
most beautiful mathematical formulas.


.. Note::

   If you are using latex notation in rst doc-strings, you need to double escape the backslash using a double-backslash fro the math elements, so then write ``\\frac`` and **not** ``\frac``.

Using LaTeX Commands Defined in ``manual.tex``
----------------------------------------------

The Pencil Code community traditionally uses a set of LaTeX shortcuts, most of them originally defined in the file ``manual.tex``.  
These macros make writing math expressions and physical equations much cleaner and more readable.

For consistency across the documentation, these definitions have been included in ``conf.py`` so they are available **everywhere** in the Sphinx build — not just in the manual itself.

This means you can now use the same LaTeX commands directly inside reStructuredText pages.  
For example:


.. code:: bash

   :math:`\pderiv{u}{x}` or :math:`\grad \cdot \uv = 0`

and get: 

:math:`\pderiv{u}{x}` or :math:`\grad \cdot \uv = 0`


Notes
^^^^^

* All math commands defined in ``manual.tex`` are automatically loaded — no need to redefine them.  
* You can use them anywhere inside ``:math:`` or ``.. math::`` blocks.  
* This keeps the documentation consistent with the notation used in the scientific papers and source comments.  


Example
^^^^^^^

.. code:: rst

   The induction equation is written as :math:`\pderiv{\Bv}{t} = \curl (\Uv \times \Bv - \eta \Jv)`.

Which renders as:

:math:`\pderiv{\Bv}{t} = \curl (\Uv \times \Bv - \eta \Jv)`

How It Works Internally
^^^^^^^^^^^^^^^^^^^^^^^


Below are the LaTeX macros registered in ``conf.py`` and grouped by type.  
Each entry shows the macro name (left) and the TeX definition that MathJax/Sphinx will use (right).


**Operators**

.. list-table::
   :header-rows: 1
   :widths: 25 75

   * - **Macro**
     - **TeX definition**
   * - ``\de``
     - ``\mathrm{d}``
   * - ``\De``
     - ``\mathrm{D}``
   * - ``\const``
     - ``\mathrm{const}``
   * - ``\erf``
     - ``\operatorname{erf}``
   * - ``\erfc``
     - ``\operatorname{erfc}``
   * - ``\grad``
     - ``\boldsymbol{\nabla}``
   * - ``\Div``
     - ``\boldsymbol{\nabla}\!\cdot``
   * - ``\curl``
     - ``\boldsymbol{\nabla}\!\times``
   * - ``\Laplace``
     - ``\nabla^2``
   * - ``\rot``
     - ``\boldsymbol{\nabla}\!\times``

**Derivatives**

.. list-table::
   :header-rows: 1
   :widths: 25 75

   * - **Macro**
     - **TeX definition**
   * - ``\pderiv{#1}{#2}``
     - ``\frac{\partial #1}{\partial #2}``
   * - ``\pderivn{#1}{#2}{#3}``
     - ``\frac{{\partial{}}^{#3} #1}{{\partial #2}^{#3}}``

**Vector notation**

.. list-table::
   :header-rows: 1
   :widths: 25 75

   * - **Macro**
     - **TeX definition**
   * - ``\vec{#1}``
     - ``\boldsymbol{#1}``
   * - ``\vcs{#1}``
     - ``\boldsymbol{\scriptstyle{#1}}``

**Common vectors**

.. list-table::
   :header-rows: 1
   :widths: 25 75

   * - **Macro**
     - **TeX definition**
   * - ``\Av``
     - ``\boldsymbol{A}``
   * - ``\Bv``
     - ``\boldsymbol{B}``
   * - ``\Jv``
     - ``\boldsymbol{J}``
   * - ``\Uv``
     - ``\boldsymbol{U}``
   * - ``\Wv``
     - ``\boldsymbol{W}``
   * - ``\Ev``
     - ``\boldsymbol{E}``
   * - ``\Fv``
     - ``\boldsymbol{F}``
   * - ``\fv``
     - ``\boldsymbol{f}``
   * - ``\gv``
     - ``\boldsymbol{g}``
   * - ``\jv``
     - ``\boldsymbol{j}``
   * - ``\kv``
     - ``\boldsymbol{k}``
   * - ``\ov``
     - ``\boldsymbol{\omega}``
   * - ``\uv``
     - ``\boldsymbol{u}``
   * - ``\vv``
     - ``\boldsymbol{v}``
   * - ``\bv``
     - ``\boldsymbol{b}``
   * - ``\xv``
     - ``\boldsymbol{x}``
   * - ``\zerovect``
     - ``\boldsymbol{0}``
   * - ``\omv``
     - ``\boldsymbol{\omega}``
   * - ``\Bhat``
     - ``\hat{B}``
   * - ``\BBhat``
     - ``\hat{\boldsymbol{B}}``

**Physics symbols**

.. list-table::
   :header-rows: 1
   :widths: 25 75

   * - **Macro**
     - **TeX definition**
   * - ``\Ra``
     - ``\mathrm{Ra}``
   * - ``\Reynolds``
     - ``\mathrm{Re}``
   * - ``\Rm``
     - ``\mathrm{Rm}``
   * - ``\vA``
     - ``v_{\mathrm{A}}``
   * - ``\cs``
     - ``c_{\mathrm{s}}``
   * - ``\csnull``
     - ``c_{{\mathrm{s}},0}``
   * - ``\Heat``
     - ``\mathcal{H}``
   * - ``\Cool``
     - ``\mathcal{C}``
   * - ``\Heavi``
     - ``\theta``
   * - ``\Strain``
     - ``\boldsymbol{\mathsf{S}}``

**Brackets & notation**

.. list-table::
   :header-rows: 1
   :widths: 25 75

   * - **Macro**
     - **TeX definition**
   * - ``\bra{#1}``
     - ``\langle #1\rangle``
   * - ``\Eq{#1}``
     - ``Eq.~(\ref{#1})``
   * - ``\Fig{#1}``
     - ``Fig.~\ref{#1}``

**Exponent helpers**

.. list-table::
   :header-rows: 1
   :widths: 25 75

   * - **Macro**
     - **TeX definition**
   * - ``\EE{n}``
     - ``\times 10^{n}``
   * - ``\ttimes{n}``
     - ``10^{n}``
   * - ``\xtimes{a}{b}``
     - ``a \times 10^{b}``

**Inequality symbols**

.. list-table::
   :header-rows: 1
   :widths: 25 75

   * - **Macro**
     - **TeX definition**
   * - ``\la``
     - ``\lesssim``
   * - ``\ga``
     - ``\gtrsim``


Adding Your Own Macros
^^^^^^^^^^^^^^^^^^^^^^^

If you need to add more macros — for instance, to document a new physical quantity or operator — you can do so by editing the ``mathjax3_config`` section of your ``conf.py``:

.. code:: python

   mathjax3_config = {
       "tex": {
           "macros": {
               "newmacro": [r"\mathrm{NewSymbol}", 0],
           }
       }
   }

The number at the end (here ``0``) specifies how many arguments your macro expects.  
After saving, rebuild the documentation — your new macro will be instantly available across all pages.


Pro tip
^^^^^^^^

If you find yourself defining too many macros, you might be trying to start your own dialect of LaTeX.  
That’s fine — just remember to share it with the rest of the fellowship so everyone’s equations still render beautifully.


Using ``_substitutions.rst``
=============================

In reStructuredText, you often find yourself repeating small bits of text or math across multiple pages — things like the project name, URLs, or symbols. That’s where ``_substitutions.rst`` comes to the rescue. It acts as a central dictionary of short text macros that can be reused anywhere in the documentation.

How it works
------------

A substitution maps a short token (wrapped in pipes) to replacement text. Put definitions in ``_substitutions.rst`` like this:

.. code-block:: rst

   .. |PC| replace:: Pencil Code
   .. |ver| replace:: v2025.1
   .. |repo| replace:: `GitHub repository <https://github.com/pencil-code/pencil-code>`__

Then use the substitution in any page:

.. code-block:: rst

   |PC| (|ver|) is available on |repo|.

At build time, Sphinx replaces each occurrence of ``|PC|``, ``|ver|``, etc., with the corresponding definition.


Notes
-----

* The leading underscore in ``_substitutions.rst`` is a convention signalling this is a helper file (not a standalone page).

* Keep substitutions short and stable — they are intended for tiny inline fragments (names, URLs, symbols), not long paragraphs.

Making the substitutions available
-----------------------------------

You can include the file manually in pages that need it:

.. code-block:: rst

   .. include:: /_substitutions.rst

Or make it global by adding the following to ``conf.py``:

.. code-block:: python

   rst_prolog = """
   .. include:: /_substitutions.rst
   """

A bit more robust way of including the substitutions is:


.. code-block:: python

   substitutions_path = os.path.join(os.path.dirname(__file__), '_substitutions.rst')
   with open(substitutions_path, encoding='utf-8') as f:
      substitutions_content = f.read()

   rst_prolog = substitutions_content

When placed in ``rst_prolog``, the substitutions are implicitly available on every page (no need to include manually).

.. important::

   In the current documentation, ``substitutions.rst`` is globally available by default. 


Typical substitutions (example table)
---------------------------------------

.. list-table::
   :header-rows: 1
   :widths: 20 80

   * - **Token**
     - **Replacement**
   * - ``|PC|``
     - Pencil Code
   * - ``|repo|``
     - `https://github.com/pencil-code/pencil-code <https://github.com/pencil-code/pencil-code>`__
   * - ``|grad|``
     - :math:`\boldsymbol{\nabla}`
   * - ``|div|``
     - :math:`\boldsymbol{\nabla}\cdot`
   * - ``|curl|``
     - :math:`\boldsymbol{\nabla}\times`



Adding math substitutions
^^^^^^^^^^^^^^^^^^^^^^^^^

Substitutions may contain inline math roles. Example entry for ``_substitutions.rst``:

.. code-block:: rst

   .. |grad| replace:: :math:`\boldsymbol{\nabla}`
   .. |Re| replace:: :math:`\mathrm{Re}`

These will render as math wherever you use the substitution.


Tips & conventions
-------------------

* Use short, descriptive tokens wrapped in pipes (e.g., ``|proj_name|``) — avoid cryptic names.  
* Group related substitutions together in ``_substitutions.rst`` (URLs, math, badges, short notices).  
* If you add many math macros, prefer adding them to ``conf.py`` as MathJax macros (so they work directly inside ``:math:`` without substitutions).  
* Document any substitution that is not self-explanatory at the top of the file.

Example: _substitutions.rst (minimal)
-------------------------------------

.. code-block:: rst

   .. |PC| replace:: Pencil Code
   .. |repo| replace:: `GitHub repository <https://github.com/pencil-code/pencil-code>`__
   .. |grad| replace:: :math:`\boldsymbol{\nabla}`