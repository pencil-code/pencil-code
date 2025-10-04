.. _pencilstyleguide:

************************
Pencil Code Coding Style
************************

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
