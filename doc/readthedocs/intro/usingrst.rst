Using reStructuredText
======================

All documentation available on this page is written in the reStructuredText
(reST) markup language.

About reST
----------

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

Converting existing documents
-----------------------------

The utility *pandoc* can be used to convert a variety of formats, including
Microsoft Word (*doc*, *docx*), Libre Office (*odt*), LaTeX,
HTML, XML, and if all else fails even PDF, to reST.

The syntax of the command is:

.. code:: bash

   pandoc input.doc -o output.rst

where *input.doc* is your input document, in any format other than *rst*.

Style guideline
---------------

Headings
~~~~~~~~

In reST, headings are marked by underlining them with the same character:

.. code:: rst

   This is a heading
   =================

In the The Pencil Code, the following markers should be used, in this order:

.. code:: rst

   Title of your page
   ==================

   Section
   -------

   Subsection
   ~~~~~~~~~~

   Sub-subsection
   ++++++++++++++

You should not use further levels of headings, as it would prevent optimal
rendering of the table of contents in the left-hand sidebar. You can structure
your document further by using the ``.. rubric::`` directive.

Do not use ``###`` and ``***``, as they are already used for higher-level headings
(e.g., on the main landing page).

Admonitions
~~~~~~~~~~~

The use of admonition directives can greatly enhance the user experience by
presenting tips, warnings, important notes, etc. in a way that stands out from
the rest of the document.

The following admonition directives are available for the Pencil Code: *attention*,
*caution*, *danger*, *error*, *hint*, *important*, *note*, *tip*, *todo*,
*warning*.

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


Images
~~~~~~


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
         | .. image:: arrow.png   | Magnetic field lines |
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
~~~~~~

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

