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



Images
------


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