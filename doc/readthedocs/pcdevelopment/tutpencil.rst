.. tutpencil:

*************************
Modifying the Pencil Code
*************************

The Pencil Code is written in Fortran90, and is hosted at github 
`<https://github.com/pencil-code/pencil-code>`_.


.. note::

        Adapted from the `Pencil Code github wiki page <https://github.com/pencil-code/pencil-code/wiki/>`_

Commit rights for the Pencil Code are given out quite liberally, but
they come with responsibility: 

- only commit code that is meaningful and necessary 
- avoid committing code that breaks any of the auto tests 
- discuss major changes with other developers first (e.g. at the
pencil-code-discuss mailing list) 
- follow the existing :ref:`pencilstyleguide`.

When developing the Pencil Code, we often communicate via commit
messages. A typical example is:

.. code::

        Fixed a bug in hydro.do_something().
        @Paul, you wrote the original subroutine, can you check whether my
        fix is OK?

and ideally a few commits down the line, we would have

.. code::

        Improved Peter's fix of the hydro.do_something() subroutine.
        Thanks for finding and analyzing this.

For this mode of communication to work, we **must** be able to rely on
our co-developers to read the commit messages. If they do not, their
contributions have a strongly reduced value for the code, and they
probably should not have commit rights in the first place.

By far the easiest way of reading the commit messages is to subscribe to
``pencil-code-commits@googlegroups.com`` and read at least superficially
through the commit messages as they appear in the mailbox. Individual
developers may prefer other ways of keeping up to date, but you should
be aware that if you do not follow what is going on with the code, it is
likely that parts of the code that you contributed will eventually get
moved around, altered, or even removed.

If you want to become a member of any of the groups without being a
committer, and if your name is not already known to us (or by googling)
we would appreciate a brief email to us explaining your interest. This
is to prevent spammers entering our lists. The pencil-code-core list, on
the other hand, is reserved for project owners only.

