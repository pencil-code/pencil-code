.. tutpencil:

***********************
Pencil Code Tutorials
***********************

The Pencil Code is written in Fortran90, and is hosted at github 
`<https://github.com/pencil-code/pencil-code>`_.



Modifying the Pencil Code
=========================

.. note::

        Adapted from this `github wiki page <https://github.com/pencil-code/pencil-code/wiki/>`_

Commit rights for the Pencil Code are given out quite liberally, but
they come with responsibility: - only commit code that is meaningful and
necessary - avoid committing code that breaks any of the auto tests -
discuss major changes with other developers first (e.g. at the
pencil-code-discuss mailing list) - follow the existing :ref:`pencilstyleguide`.

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

.. _pencilstyleguide:

Coding Style guide
==================

.. note::

        Adapted from this `github wiki page <https://github.com/pencil-code/pencil-code/wiki/CodingStyle>`_

We describe here the Pencil Code coding style and best practice to write
code in form of a checklist that should applied to each of the given
items. Of course, no rules can be hammered in stone and always need some
reality check to ensure these rules help improving the code and not the
opposite.

Is a module…
------------

-  named well in the sense that its name describes its purpose?
-  abstract enough to form a separate module?
-  consistent within the existing modules scheme?
-  interface obviously revealing how the module should be used?
-  interface abstract enough so that the module can be used without
   thinking about how the services are implemented?


Is a subroutine…
----------------

-  name revealing exactly what this subroutine is doing?
-  implementing only one task that is well defined?
-  containing only code parts that would not better be put in a separate
   subroutine?
-  interface obviously revealing how the subroutine should be used?
-  interface abstract enough so that the subroutine can be used without
   thinking about how it is implemented in detail?


Is a data type…
---------------

-  named well so that its name describes its content type?
-  descriptive so that it helps to document its variable declaration?
-  simple so that it minimizes complexity?
-  that needs to be complex operated only through access subroutines
   (set_XY/get_XY)?

Is a variable…
--------------

-  really necessary to be variable (and not a constant)?
-  given all possible attributes (like “save”, “parameter”,
   “intent(in/out)”)?
-  not some kind of a “magic number” or “magic string” that should be
   converted to a named constant?
-  not redundant and distinct to all other variables that are otherwise
   available?
-  used only for the single purpose that it was intended for by its
   name?
-  additionally defined and used, if the clarity of the code is
   significantly improved?

Is a variable name…
-------------------

-  chosen well so that it describes its content (i.p. not its data
   type)?
-  readable and following the style
   “lower_case_words_connected_with_underscores”?
-  of a boolean indicating its flag-like behavior (currently “lflag” for
   a flag)?
-  of a loop counter more informative than just i, j, k, l, m, n? ([i,j]
   should be reserved only for general matrix and vector operations,
   while [l,m,n] are reserved for the m-n-loop, for handling the f-array
   or p-pencil indices).

Is a logical statement…
-----------------------

-  using only simple boolean expressions?
-  better stored it into an additional boolean variable or put into a
   boolean function, if it is to be reused?
-  not using double negations? (oops!)

Is an ``if``-``else``-construct…
--------------------------------

-  consisting of small blocks where both blocks are of similar size?
-  written so that the “normal” case appears first?
-  used to minimize complexity?

Is a loop…
----------

-  performing exactly one well-defined function, as a subroutine would?
-  implemented using the best matching type:
   ``do``/``while``/``repeat``?
-  nested in another loop only if necessary?

Is the code…
------------

-  representing its own logical structure?
-  nominal path or calling sequence clear and easy to follow?
-  organized so that related statements are grouped together?
-  free of relatively independent code blocks that could stay in
   subroutines?
-  hiding implementation details to the greatest extent?
-  written in respect to the problem solution and not in terms of
   programming-language requirements?
-  initializing all variables outside any conditional statements?
   (e.g. code inside ``if``-``elseif``-``elseif``-constructs might never
   be executed!)
-  compiling without any compiler warnings? (yes, they do have a serious
   background even if it is sometimes not obvious!)

Is a commit message…
--------------------

-  not “new settings” or “some corrections” or “minor changes” or
   similar?
-  telling which code block is affected?
-  containing all relevant major changes?
-  informative about the consequences of previous bugs that are now
   fixed?
-  giving hints what to search or where to start reading if anyone is
   interested in details?

FORTRAN formatting
------------------

-  use the Fortran95 language standard (F95)
-  use spaces for indentation, two spaces represent one indentation
   level
-  use spaces for formatting of tabular data or comments
-  no spaces at the end of a line
-  one empty line is enough to split code blocks
-  two empty lines can be used to split distinct parts of code
-  in-code comments should be indented together with the code
-  block-like comments (e.g. function headers) start at the beginning of
   a line
-  use spaces around operators, where applicable
-  line-breaks are required after 130 characters (F95)
-  line-breaks can be used to significantly improve readability of long
   code lines

Typical rule-breaker and its solution
-------------------------------------

-  ``goto`` => implement a loop or an ``if``-``else``-construct
-  ``entry`` => implement an interface or split into distinct
   subroutines
-  ``format`` => put the format string inside each ``write`` statement
-  hard-coded file units => use a named constant
-  hard-coded string length => use pre-defined global constants

Recommended further reading
---------------------------

-  Kernighan, Brian, and Plauger: “The Elements of Programming Style”,
   2nd ed., McGraw-Hill, New York, 1978
-  Kernighan, Brian, and Pike: “The Practice of Programming”, Addison
   Wesley, Reading (Massachusetts), 1999
-  McConnell: “Code Complete”, 2nd ed., Microsoft Press, Redmont
   (Washington), 2004
-  Hunt and Thomas: “The Pragmatic Programmer”, Addison Wesley, Reading
   (Massachusetts), 1999


Social Rules?
=============

.. note::

        Adapted from this `github wiki page <https://github.com/pencil-code/pencil-code/wiki/SocialRules>`_

Hi guys,

I discussed with a co-developer of a code that is developed pretty much
like Pencil, open-source, a team, version control, etc. Talking to him
about code development, I asked if there were cases of flame fights or
heated arguments in the code community. He mentioned a couple of cases,
and pointed me to **books** on open source development where such stuff
is discussed. Not surprisingly, it is quite a common occurrence.

Chapter 6 of the first link, from 102 on (“Difficult People”), is
particularly relevant.

`<http://producingoss.com/>`_ 

`<http://artofcommunityonline.org/>`_ 

Wlad.

Good electronic communication
-----------------------------

For an electronic discussion, there is no such thing as a meta-level of
information transfer. Therefore, every good electronic communicator just
stays with the facts. And *if* an interpretation needs to be done, one
chooses the interpretation that assumes *best* motives of your opponent.
Only then, one has a chance to understand the opponent right. And
without understanding an opponent *fully*, one has no right to answer.
(Philippe)


Code of Conduct
---------------

Although spaces may feel informal at times, we want to remind ourselves
that this is a professional space. As such, the Pencil Code community
adheres to a code of conduct adapted from the Contributor Covenant
(`<https://www.contributor-covenant.org/>`_) code of conduct. All
contributors will be required to confirm they have read our 
`code of conduct <https://github.com/pencil-code/pencil-code/blob/master/license/CODE_OF_CONDUCT.md>`_,
and are expected to adhere to it in all Pencil Code spaces and
associated interactions.
