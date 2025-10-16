.. _manualAppendixA:

**********************************************************
Appendix A: Timings
**********************************************************

Timings
=======

In the following table we list the results of timings of the code on
different machines. Shown is (among other quantities) the wall clock time
per mesh point (excluding the ghost zones) and per full 3-stage time step,
a quantity that is printed by the code at the end of a run [#]_


As these results were assembled during the development phase of the code
(that hasn't really finished yet, ...), you may not get the same numbers,
but they should give some orientation of what to expect for your specific
application on your specific hardware.

The code will output the timing (in microseconds per grid point per time-step)
at the end of a run. You can also specify ``walltime`` in ``print.in`` to
have the code continuously output the physical time it took to reach the
time-steps where diagnostics is done. The time-dependent code speed can then
be calculated by differentiating, e.g., in IDL:

.. code:: bash

    IDL> pc_read_ts, obj=ts
    IDL> plot, ts.it, 1/nw*deriv(ts.it, ts.walltime/1.0e-6), psym=2

where ``nw = nx*ny*nz``.

.. [#] Note that when using ``nompicomm.f90``, the timer currently used may overflow on some machines, so you should not blindly trust the timings
   given by the code.


.. include:: appendixA_table.rst