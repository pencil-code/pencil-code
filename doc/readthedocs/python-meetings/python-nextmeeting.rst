
.. _python-nextmeeting:

Next Meeting
=============

**Date:**  January 29th 2026 at 3pm CET.


**Proposed agenda:**

- Tutorials

- Check :file:`pencil-code/python/pencil/TO_DO_LIST.txt`. The last update is from Jun 29, 2024, is there any further update?

- Changing behaviour of ``.keys()`` methods (of ``Averages``, ``Timeseries`` etc.). Currently it just prints the keys, but it would be more useful to return a list

- Decide on a convention for extra debug output (controlled by the ``quiet`` keyword in functions like ``pc.read.var`` and ``pc.read.grid``). Currently, whether it defaults to ``False`` or ``True`` varies from function to function, which is confusing and annoying. Ideally we would control such output by a module-wide flag (such that something like ``pc.shut_up = True`` at the beginning of a script would have the effect of setting ``quiet=True`` in all the ``pc`` functions).

- Changing the behaviour of ``pc.io.get_value_from_file``. Currently, whether it has been successful of not is indicated by a boolean return value (as if it is a C program). I think it makes more sense to explicitly raise an error message. I presume most users would want to know that the function is failing, rather than dealing with crytic errors later on.

- Rename ``pc.sim.__Simulation__`` to ``pc.sim.Simulation``. Leading underscores usually signify that something is not meant to be publicly accessed; why do we not want users to directly use this class?

- Supported Python versions: is it enough to go by the versions currently supported by upstream (<https://devguide.python.org/versions/>), or do we want to support Python versions older than 5 years?

- The item **Interactive vs. library mode** in :file:`pencil-code/doc/readthedocs/python-meetings/python-decisions.rst` is now obsolete due to the implementation of lazy loading.

- Folders in `python/` (like `meetings` and `docs`) are detected as Python modules if `$PENCIL_HOME/python` is added to `$PYTHONPATH`. Options to avoid this:

  * move the Python modules to `python/src` and instead add this to `$PYTHONPATH`
  * ??

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



See also:
----------

* :ref:`past-meetings` — for notes and outcomes of previous meetings
* :ref:`decisions` — for items already agreed upon
