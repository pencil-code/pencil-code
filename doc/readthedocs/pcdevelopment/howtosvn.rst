.. _howtosvn:

*****************************************
svn: subversion, to who did what and when
*****************************************

svn is simple -- no need to make life more complicated.

Useful Documentation
====================

* `svnHowto <https://svnhowto.com/>`_ : step by step guide

* `svn tutorial <https://www.geeksforgeeks.org/svn/svn-tutorial/>`_

Configuration
=============

None (usually)

.. _howtosvn-svnhub:

svngit
------

`svngit <https://github.com>`_ since the svn bridge on github no longer works.


Cloning a Repository (checking out the code)
============================================

If the userID is AxelBrandenburg, then just write

.. code:: bash

    $ svn checkout --username=AxelBrandenburg https://pencil-code.org/svn/trunk pencil-code

You can also leave out "--username=AxelBrandenburg"

Checking Status & Differences
=============================

Checking the status of the svn repository:

.. code:: bash

    $ svn status

It will show a lot of files that are not under svn. To avoid this, just say

.. code:: bash

    $ svn status | grep -vE "^\?"

* To update the code, say

.. code:: bash

    $ svn up

You can either be directly under the pencil-code directory or under its src directory, for example.

* If you have your own changes, those file will appear with an M in front.
You can check them, or some of them, in like so

.. code:: bash

    $ svn ci src/magnetic.f90

If you are already in the src directory, you should omit src.

* To check the history of the file, say

.. code:: bash

    $ svn log src/magnetic.f90

It might be good to pipe the result into :command:`more` or into a file.

Checking Differences
--------------------

Checking for differences with local repository:

.. code:: bash

    $ svn diff -r41894 magnetic.f90

You will then see all changes since revision r41894.
If you want to see the difference between two subsequent changes, say

.. code:: bash

    $ svn diff -r41618:41638 magnetic.f90

.. note::

    Peek into the timeline before changing history


Pulling & Stashing
==================

In svn, you would usually say

.. tip::

    It is good to keep your version up-to-date.

* just update the code with

.. code:: bash
        
    $ svn up

As stated above, svn would update from the directory that you are in,
and other ones underneath.

* Before you update the code, it would be good to review your changes:

.. code:: bash
    
    $ svn diff magnetic.f90

would show you what has been changed.
If your changes are the ones you want to check in, just say

.. code:: bash
    
    $ svn ci magnetic.f90

In that case, your default editor will pop up, and there you would write your check-in message.

Keeping Uncommitted Changes
---------------------------

Your local changes will be preserved. But there can be conflicts, if somebody does changes in the same lines as your.
The conflicts are marked and explained.
You'd need to remove all conflict markers and make sure your version has been corrected before you can check them in.


Pushing Changes
================

.. attention::

    This is the same as checking in your changes.


Moving Files & Directories
==========================

This is straightforward.
Here an example:

.. code:: bash

    $ svn mv magnetic.f90 obsolete
    $ svn ci -m "Just checking whether I can more the file to the obsolete folder"

Branching
=========


