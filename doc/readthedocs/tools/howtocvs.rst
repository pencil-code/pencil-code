.. _howtocvs:
**********************
How to CVS
**********************

I know… I know…  
CVS? Really? Yes, really. Some venerable projects still swear by it, and if you find yourself here, congratulations — you are about to pilot a vintage time machine of version control.  
Think of it as a classic car: a bit old-school, occasionally cranky, but if you know what you’re doing, it still gets you where you need to go.



Some documentation
==================

* `Quick guide <https://www.cs.umb.edu/~srevilak/cvs.html>`_

* `Longer guide <https://www.yolinux.com/TUTORIALS/LinuxTutorialCVSintro.html>`_

* `Full CVS manual <https://www.gnu.org/software/trans-coord/manual/cvs/>`_

Configuration
=============

Setting up the environment
-------------------------- 

Before CVS can take off, you need to tell your shell **where the repository lives**:

.. code:: bash

    $ CVSROOT=/path/to/cvsroot
    $ export CVSROOT

Replace `/path/to/cvsroot` with the actual path or remote location.  
   Pro tip: put this in your `.bashrc` to save future headaches.

Set your preferred editor:



.. code:: bash

    $ CVSEDITOR=vi




General syntax 
==============

CVS commands follow a simple-but-old-school pattern:

.. code:: bash

    $ cvs cvs-options subcommand subcommand-options

Where subcommand is the thing you are asking cvs to do.

* To see a full list of subcommands:

.. code:: bash

    $ cvs --help

* To get detailed syntax for a specific subcommand:

.. code:: bash
    
    $ cvs -H subcommand


.. note::

   Think of CVS syntax like a riddle written in the 1980s — it works if you know the secret handshake.


Dowloading files
================

* Checkout an existing directory:

.. code:: bash
    
    $ cvs co dir

* Checkout without recursion:

.. code:: bash
    
    $ cvs co -l dir

* Update your working directory:

.. code:: bash
    
    $ cvs up -d

* List contents of the current directory:

.. code:: bash

    $ cvs ls

.. note::

   CVS doesn’t automatically tell you everything — it expects you to know what you’re doing. Like an old-school mentor.


Uploading changes
=================

After modifying your files:

* Commit changes:

.. code:: bash
    
    $ cvs ci -m 'commit message'

.. note::

   Remember: CVS commits are sacred. It trusts you to preserve history… or at least not destroy it entirely.


Checking differences
====================

See what has changed in your working copy:

.. code:: bash

    $ cvs status

.. note::

   Pro tip: `cvs diff` is your friend. Treat it like a telescope for observing the tiny ripples you’ve made in the timeline.



Adding files and directories
============================

* Add a new directoy:

.. code:: bash

    $ mkdir newdir
    $ cvs add newdir


* Add a new file:

.. code:: bash

    # create file
    $ cvs add newfile
    $ cvs ci -m 'message' 


Merging revisions
==================


If your file is out of date with the repository version, CVS will require you to **merge changes manually**. Think of it as carefully combining two timelines without causing a paradox.

.. note::

   CVS doesn’t do hand-holding. If two changes collide, **you** decide which timeline survives.


Resolving conflicts
===================

CVS inserts **conflict markers** when it detects overlapping edits:

.. code::

    <<<<<file 
    some output
    ==================
    some other output
    >>>>>

Edit the file, remove the conflict markers, and then commit.

.. note::

   Patience is a virtue. CVS expects you to be the arbiter of history.


Deleting files
===============

.. code:: bash

    $ rm filename           # remove local copy first
    $ cvs delete filename   # mark file for removal from repository
    $ cvs commit            # commfinalize deletionit

.. note::

   CVS likes a ceremonial approach — remove locally first, then let it know remotely, then commit.  


Other useful commands
======================

* Show differences between your local copy and the repository version:

.. code:: bash

    $ cvs diff filename
    $ cvs diff -r 1.2 filename           # Compare with version 1.2
    $ cvs diff -r 1.2 -r 1.3 filename    # Compare version 1.2 with 1.3

* Show commit log for a file:

.. code:: bash

    $ cvs log filename

* Annotate a file (see who changed each line):

.. code:: bash 

    $ cvs annotate filename

.. note::
   Useful for figuring out “who touched this timeline and why.”

