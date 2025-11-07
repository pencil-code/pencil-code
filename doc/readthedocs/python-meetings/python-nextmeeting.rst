
.. _python-nextmeeting:

Next Meeting
=============

**Tentative date:** TBD

Please fill up your prefered dates `here <https://www.when2meet.com/?33402888-Q69GB>`__ before November 14, 2025. 

If the dates and hours do not work for you, please add your constraints here.


**Proposed agenda:**

- Summary of previous meetings

- Meeting format


To propose topics or dates, please edit this file.

.. _discussion-topics:

Discussion Topics
=================

Ongoing or pending topics for future meetings.

- Deprecate usage of non standard python modules like ``eqtools`` and try to use ``numpy`` instead when possible.
- Changing behaviour of ``.keys()`` methods (of ``Averages``, ``Timeseries`` etc.). Currently it just prints the keys, but it would be more useful to return a list
- Decide on a convention for extra debug output (controlled by the ``quiet`` keyword in functions like ``pc.read.var`` and ``pc.read.grid``). Currently, whether it defaults to ``False`` or ``True`` varies from function to function, which is confusing and annoying. Ideally we would control such output by a module-wide flag (such that something like ``pc.shut_up = True`` at the beginning of a script would have the effect of setting ``quiet=True`` in all the ``pc`` functions).
- Changing the behaviour of ``pc.read.aver`` for ``yaver`` and ``zaver``; see <https://groups.google.com/g/pencil-code-python/c/a6eu61yOMuk>


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
