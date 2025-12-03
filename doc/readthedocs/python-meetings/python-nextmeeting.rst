
.. _python-nextmeeting:

Next Meeting
=============

**Date:** TBA. January 2026. 
Please fill out your prefered dates `here <https://www.when2meet.com/?33748149-PptLx>`_.


**Proposed agenda:**

- Tutorials

- Check :file:`pencil-code/python/pencil/TO_DO_LIST.txt`. The last update is from Jun 29, 2024, is there any further update?

- Changing behaviour of ``.keys()`` methods (of ``Averages``, ``Timeseries`` etc.). Currently it just prints the keys, but it would be more useful to return a list
- Decide on a convention for extra debug output (controlled by the ``quiet`` keyword in functions like ``pc.read.var`` and ``pc.read.grid``). Currently, whether it defaults to ``False`` or ``True`` varies from function to function, which is confusing and annoying. Ideally we would control such output by a module-wide flag (such that something like ``pc.shut_up = True`` at the beginning of a script would have the effect of setting ``quiet=True`` in all the ``pc`` functions).
- Changing the behaviour of ``pc.io.get_value_from_file``. Currently, whether it has been successful of not is indicated by a boolean return value (as if it is a C program). I think it makes more sense to explicitly raise an error message. I presume most users would want to know that the function is failing, rather than dealing with crytic errors later on.

- Supported Python versions: is it enough to go by the versions currently supported by upstream (<https://devguide.python.org/versions/>), or do we want to support Python versions older than 5 years?

- The item **Interactive vs. library mode** in :file:`pencil-code/doc/readthedocs/python-meetings/python-decisions.rst` is now obsolete due to the implementation of lazy loading.


To propose topics or dates, please edit this file.

.. _discussion-topics:

Discussion Topics
=================

Ongoing or pending topics for future meetings.




Feel free to add more — this section is a living to-do list.
If an item you want to discuss is already in the list, please add a ``+1`` to give the item priority according to the number of people interested.

.. _ideas-and-proposals:

Ideas and Proposals
===================

If you have a concept or feature suggestion, add it here following this format:

.. code:: text

   Title
   ---------

   Brief description of the idea

   **Motivation:** Why this is useful
   **Proposal:** How to implement it
   **Reference:** (optional links or related work)

Example:


Example: Postprocessing for shock detection
---------------------------------------------

**Motivation:** Identify shocks automatically in 3D data

**Proposal:** Implement gradient-based detection in `diagnostics.py`

**Reference:** Smith et al. (2023), J. Comp. Phys. 492, 112334

See also:
----------

* :ref:`past-meetings` — for notes and outcomes of previous meetings
* :ref:`decisions` — for items already agreed upon
