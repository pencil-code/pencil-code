.. _manualAppendixB:

**********************************************************
Appendix B: Coding standard
**********************************************************


The numerous elements that make up the |PC| are written
in a consistent style that has evolved since it was first created.
Many people have contributed their knowledge and experience with
in this and the result is what we believe is and extremely
readable and manageable code.

As well as improving the readability of the code, by having some
naming conventions for example aids greatly in understanding what
the code does.

There is a standard for all aspects of the code, be it Fortran source,
shell scripts, Perl scripts, LaTeX source, Makefiles, or otherwise.
Where nothing has been explicitly stated it is recommended that
similar existing examples found in the code are used as a template.


.. _manApBCodingStandards:

Fortran Code
------------


The code should remain fully compatible with the :code:`Fortran90` standard.
This ensures that the code will run on all platforms.
Indeed, an important aspect of |PC| philosophy is to be maximally flexible.
This also means that useful non-standard extensions to the code should be hidden in and
be made accessible through suitable non-default modules.

Fortran is not case-sensitive but in almost all instances we prescribe
some form of capitalization for readability.

In general all Fortran code including keywords, variable names etc. are
written in lowercase.
Some of the coding standard has already been discussed in
Sect.~:ref:`ProgrammingStyle`.
Here we discuss and amplify some remaining matters.

Indenting and whitespace
^^^^^^^^^^^^^^^^^^^^^^^



Whitespace should be removed from the end of lines.

Blank lines are kept to a minimum, and when occurring in subroutines
or functions are replaced by a single :code:`!'` in the first column.


Tab characters are not used anywhere in the code.  Tab characters are
not in fact allowed by the :code:`Fortran` standard and compilers that accept
them do so as an extension.

All lines are kept to be not more than 80 characters long.
Where lines are longer they must be explicitly wrapped using the
Fortran continuation character :code:`\&`.
Longer lines (up to 132 characters) and additional spaces are allowed
in cases where the readability of the code is enhanced, e.g., when one
line is followed by a similar one with minor differences in some places.

Code in syntactic blocks such as :code:`if -- endif`, :code:`do -- enddo`,
:code:`subroutine -- endsubroutine` etc. is always indented by precisely two spaces.
The exception to this is that nested loops where only the innermost loop
contains executable code should be written with the :code:`do-- enddo`
pairs at the same level of indentation, using the following pattern:

.. code:: fortran

    do n=n1,n2
    do m=m1,m2
      [...]
    enddo
    enddo

Alternatively nested loops may be written on a single line:

.. code:: fortran

    do n=n1,n2; do m=m1,m2
      [...]
    enddo; enddo

Comments
^^^^^^^^



Descriptive comments are written on their own lines unless there is a strong
reason to do otherwise. Comments are never indented and the :code:`!'` should appear
in the first column followed by two spaces and then the text of the comment.
Extremely short comments may follow at the end of a line of code, provided
there is space.

Comments also must not exceed the 78 character line length and should be
wrapped onto more lines as needed.

Typically comments should appear with a blank commented line above and below
the wrapped text of the comment.

All subroutine/functions begin with a standard comment block describing what
they do, when and by whom they were created and when and by whom any
non-trivial modifications were made.

Comments should be written in sentences using the usual capitalization and
punctuation of English, similar to how text is formatted in an e-mail or a
journal article.

For example:

.. code:: console

      some fortran code
      some more fortran code
    !
    !  A descriptive comment explaining what the following few lines
    !  of code do.
    !
      the fortran code being described
      the fortran code being described
      ...
    !
    !  A final detail described here.
    !
      the final fortran code
      the final fortran code
      ...

Subroutines and functions are started with a comment block describing
what they do, when and by whom they were created and when and by whom any
non-trivial modifications were made. The layout of this comment block
is a standard, for example:

.. code:: fortran

    !***********************************************************************
        subroutine initialize_density(f,lstarting)
    !
    !  Perform any post-parameter-read initialization i.e. calculate derived
    !  parameters.
    !
    !  For compatibility with other applications, we keep the possibility
    !  of giving diffrho units of dxmin*cs0, but cs0 is not well defined general.
    !
    !  24-nov-02/tony: coded
    !   1-aug-03/axel: normally, diffrho should be given in absolute units
    !

where dates are written in dd-mmm-yy format as shown and names appearing
after the :code:`/'` are either the users cvs login name or, where such exists
amongst the |PC| community, the accepted short form (:math:`\approx 4`
characters) of the authors name.

Module names
^^^^^^^^^^^^

The names of modules are written with initial letter capitalization
of each word and the multiple words written consecutively without
any separator.

Variable names
^^^^^^^^^^^^^^

Variables are given short but meaningful names and written
in all lowercase. Single character names are avoided except
for commonly used loop indices and the two code data structures
of the |PC|: :code:`f` the main state array (see :ref:`f-array`) and :code:`p` the pencil case
structure (see :ref:`pencil-case`).

Quantities commonly represented by a particular single character
in mathematics are typically given names formed by repeating the
character (usually in lowercase), e.g., the velocity :math:`u` becomes :code:`uu`,
specific entropy :math:`s` becomes :code:`ss` etc.

Temperature in variable names is denoted with a capital T so as not to
be confused with time as represented by a lowercase t. Note however that since
Fortran is not case sensitive the variables, for example :code:`TT` and :code:`tt`, are the same
so distinct names must be used. For this reason time is usually represented
by a single :code:`t` contrary to the above guideline.

The natural log of a quantity is represented by adding :code:`ln` to its
name, for example log of temperature would be :code:`lnTT`.

There are some standard prefixes used to help identify the type and nature
of variables as follows:

*   ``i``       -- Denotes integer variables typically used as array indices.

*    ``i_``      -- Denotes pencil case array indices.

*    ``idiag_``  -- Denotes diagnostic indices.

*    ``l``       -- Denotes logical/boolean flags

*    ``cdt``     -- Denotes timestep constraint parameters.

*   ``unit_``   -- Denotes conversion code/physics unit conversion parameters.

Emacs settings
^^^^^^^^^^^^^^

Here are some settings from wd's :file:`~/.emacs` file:

.. code:: emacs-lisp

    ;;; ~/.f90.emacs
    ;;; Set up indentation and similar things for coding the |PC|.
    ;;; Most of this can probably be set through Emacs' Customize interface
    ;;; as well.
    ;;; To automatically load this file, put the lines
    ;;;   (if (file-readable-p "~/.f90.emacs")
    ;;;       (load-file "~/.f90.emacs"))
    ;;; into your ~/.emacs file.

    ;; F90-mode indentation widths
    (setq f90-beginning-ampersand nil) ; no 2nd ampersand at continuation line
    (setq f90-do-indent           2)
    (setq f90-if-indent           2)
    (setq f90-type-indent         2)
    (setq f90-continuation-indent 4)

    ;; Don't use any tabs for indentation (with TAB key).
    ;; This is actually already set for F90-mode.
    (setq-default indent-tabs-mode nil)

    ;; Ensure Emacs uses F90-mode (and not Fortran-mode) for F90 files:
    (setq auto-mode-alist
          (append
           '(
             ("\\.[fF]90$"   . f90-mode)
             ("\\.inc$"      . f90-mode)
            )
           auto-mode-alist))

    ;; Make M-Backspace behave in Xemacs as it does in GNU Emacs. The default
    ;; behavior is apparently a long-known bug the fix for which wasn't
    ;; propagated from fortran.el to f90.el.
    ;; (http://list-archive.xemacs.org/xemacs-patches/200109/msg00026.html):
    (add-hook 'f90-mode-hook
              (function (lambda ()
                 (define-key f90-mode-map [(meta backspace)] 'backward-kill-word)
    )))

Other best practices
--------------------

When implementing :code:`IF` or :code:`SELECT` blocks always write
code for all cases -- including the default or else case.
This should be done even when that code is only a call to
raise an error that the case should not have been reached.
If you see a missing case anywhere then do add it. These
failsafes are essential in a large multi-purpose multi-user
code like the |PC|.

If a case is supposed to do nothing and it may be unclear that the
coder has recognized this fact then make it explicit by adding the
default case with a comment like

.. code:: fortran

    ! Do Nothing

The compiler will clean away any such empty blocks.


General changes to the code
---------------------------

It is sometimes necessary to do major changes to the code.
Since this may affect many people and may even be controversial
among the developers, such changes are restricted to the time
of the next |PC| User Meeting.
Such meetings are advertised on `<http://www.nordita.org/software/pencil-code/>`_
under the news section.
Notes about previous such meetings can be found under
`<http://www.nordita.org/software/pencil-code/UserMeetings/>`_.

Major changes can affect those developers who have not checked in
their latest changes for some time.
Before doing such changes it is therefore useful to contact the people
who have contributed to the latest developments on that module.
If it is not functional or otherwise in bad shape, it should be
moved to :file:`experimental`, i.e. one says:

.. code:: bash

    svn mv file.f90 experimental/file.f90

or:

.. code:: bash

    git mv file.f90 experimental/file.f90

However, any such directory change constitutes a major change in
itself and should be performed in agreement with those involved in
the development.
Otherwise any file that has been changed in the meantime will end up
being outside revision control, which is to be avoided at all cost.
