.. _python-decisions:

Decisions and Outcomes
======================
This page collects all decisions made by the Python postprocessing group,
organized by topic. Each decision includes a link to the meeting where it
was taken. These represent the current consensus and coding policies of
the group.



.. note::

   The decisions listed here summarize the outcome of discussions held in
   our development meetings. They serve as the current working guidelines
   for contributors to the Pencil Code Python postprocessing tools.

   Decisions may evolve as the codebase and tools mature.
   Proposed changes should be discussed in the upcoming meetings or through
   the Python discussion mailing list before being incorporated here.

.. rubric:: Last updated
   
2025 Nov 6


.. rubric:: Contents


.. contents::
   :local:
   :depth: 2





Code Formatting
---------------


* **Use of Black**

  All Python code in the project must be formatted using `Black` before being merged.
  This ensures consistency and readability across all contributions.

  - Add this rule to the style guide.
  - Contributors should run Black before submitting changes.
  - Regular reminders will be sent to the Python group mailing list.

  *Related meeting:* :ref:`2020 Nov 23 <python-pastmeetings>` and :ref:`2021 Apr 14 <python-pastmeetings>`

* **Formatting consistency**

  Formatting differences between contributors are resolved by using Black as the default.
  All new and existing files should be progressively reformatted to match this standard.

  *Related meeting:* :ref:`2021 Sept 17 <python-pastmeetings>`


Documentation
-------------

* **Documentation system**

  The documentation for the Python postprocessing tools is generated using
  `Sphinx`, hosted on `ReadTheDocs` and integrated with GitHub.

  *Related meeting:* :ref:`2021 Feb 24 <python-pastmeetings>` and :ref:`2020 Nov 23 <python-pastmeetings>`

* **Docstring format**

  Every function must follow a unified docstring structure containing the following
  sections:

  - General definition
  - Signature
  - Parameters
  - Returns
  - Return Type
  - Examples
  - Notes (optional)

  This structure must be described in the style guide.

  *Related meeting:* :ref:`2021 Sept 17 <python-pastmeetings>`

* **File-level documentation**

  Each module should include a short description of its submodules in its
  respective ``__init__.py`` file, instead of keeping all descriptions in
  ``pencil/__init__.py``.

  *Related meeting:* :ref:`2021 Apr 14 <python-pastmeetings>`


Imports and Package Structure
-----------------------------

* **Import location**

  According to PEP8, imports should appear at the top of each file.
  However, the group acknowledged that for interactive use (e.g. in IPython),
  importing within functions may still be convenient. The convention remains open,
  with both approaches acceptable depending on context.

  *Related meeting:* :ref:`2020 Nov 23 <python-pastmeetings>`

* **Avoiding cyclic imports**

  Cyclic dependencies should be avoided. A dedicated script
  (``python/utils/cyclic-imports``) was developed to detect them.
  Core utility functions such as ``is_sim_dir()`` were moved to independent
  modules (``pencil/util``) to prevent dependency loops.

  *Related meeting:* :ref:`2021 Apr 14 <python-pastmeetings>`

* **Import statements**

  Use absolute imports instead of relative ones, for example:
  ``from pencil.sim import simulations`` instead of
  ``from ..sim import simulations``.

  *Related meeting:* :ref:`2021 Apr 14 <python-pastmeetings>`

* **Package naming**

  The internal package ``backpack`` should be renamed to ``third_party`` to
  better reflect its purpose as a collection of external or optional libraries.

  *Related meeting:* :ref:`2021 Apr 14 <python-pastmeetings>`

* **Interactive vs. library mode**

  To improve import speed and interactivity, a distinction between “interactive”
  and “library” usage modes was proposed. Lighter imports should be possible
  for exploratory work in IPython or notebooks.

  *Related meeting:* :ref:`2021 Apr 14 <python-pastmeetings>`


Testing
-------

* **Testing framework**

  The group agreed to extend the existing test suite located at
  ``pc/python/pencil/tests``. New routines should include appropriate test
  functions whenever possible.

  *Related meeting:* :ref:`2021 Sept 17 <python-pastmeetings>`

* **General approach**

  New test functions will help ensure stability as functionality expands.
  Regular testing was encouraged as part of the development workflow.

  *Related meeting:* :ref:`2021 Apr 14 <python-pastmeetings>`


Tutorials and Examples
----------------------

* **Tutorial formats**

  Tutorials should exist in both **Jupyter notebook** and **plain script**
  formats. Jupyter notebooks can be converted to scripts using ``jupyter nbconvert --to script``
  or ``pandoc``.

  *Related meeting:* :ref:`2020 Nov 23 <python-pastmeetings>`

* **Directory structure**

  Tutorial directories should be renamed and organized by content relevance.
  For instance, include a "Getting started" directory and other thematic
  folders for specific physics or analysis topics.

  *Related meeting:* :ref:`2020 Nov 23 <python-pastmeetings>`


Deprecated or Removed Components
--------------------------------

* **npfile.py**

  The deprecated ``npfile.py`` module was removed as part of cleanup and modernization.
  Its functionality was replaced in other parts of the code.

  *Related meeting:* :ref:`2021 Apr 14 <python-pastmeetings>` and :ref:`2021 Feb 24 <python-pastmeetings>`


Communication and Community
---------------------------

* **Mailing list**

  A dedicated mailing list was established for discussion and coordination.
  All contributors are encouraged to join.

  *Related meeting:* :ref:`2021 Apr 14 <python-pastmeetings>`

* **Newsletter**

  A regular newsletter will be sent to share project updates, highlight recent
  improvements, and encourage participation.

  *Related meeting:* :ref:`2021 Apr 14 <python-pastmeetings>`


General Policies
----------------

* **Adding new routines**

  New routines and extensions to existing ones should follow consistent
  testing and documentation practices. Each addition should ideally include:

  - Corresponding test functions
  - Updated docstrings in the unified format
  - Example usage if applicable

  *Related meeting:* :ref:`2021 Sept 17 <python-pastmeetings>`

* **pv_* visualization routines**

  Contact the authors of the existing ``visu/pv_*`` routines before extending
  or modifying them.

  *Related meeting:* :ref:`2021 Sept 17 <python-pastmeetings>`

* **Automatic reload**

  Due to issues with Python’s reloading mechanisms (e.g. function names
  overlapping with filenames), reloading behavior remains a known limitation
  when working in interactive sessions.

  *Related meeting:* :ref:`2020 Nov 23 <python-pastmeetings>`