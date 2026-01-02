.. _howtogit:

***********************************
Git: The Time Machine for Your Code
***********************************




Git is like a time machine for your code. It remembers everything you’ve done, lets you undo your worst decisions, and even allows you to create alternate timelines without breaking reality (usually).

.. _howtogit-quick:

Super Quick workflow
=====================



If you are already familiar with Git and have everything configured for |PC|, this is the **minimal, paradox-resistant workflow** for committing and pushing changes.

Rule number one: **the timeline lies unless you check it**.  
Rule number two: check it with ``git status``.

.. code:: bash

    git status                          # Confirm branch, directory, and current state
    git add fileA fileB                 # Stage the files you want to commit
    git status                          # Did I stage what I meant to stage?
    
    git commit -m "Explanatory message" # Commit with a clear, descriptive message
    git status                          # Sanity check: working tree clean?

    git stash                           # Hide any remaining local changes (if any)
    git status                          # Nothing suspicious left behind?

    git pull --rebase                   # Sync with the main timeline (mandatory for |PC|)
    git status                          # Still on the right branch? Good.

    git push                            # Push your commits
    git status                          # Victory lap: confirm success

    git stash pop                       # Restore uncommitted local changes
    git status                          # Did everything come back alive?


.. note::

    If ``git status`` shows a clean working tree before pulling,  the ``git stash`` / ``git stash pop`` steps can be safely skipped.

.. note::

    ``git status`` is cheap, fast, and harmless. Running it too often has **no known side effects**.



For a full explanation of each step —  and how to resolve temporal paradoxes —
continue reading below.


History Rewriting Hazard (Git + SVN Server)
============================================

.. warning::

    **This server supports both Git and SVN.**
    This is not a normal Git setup.
    Certain Git workflows that are perfectly valid in pure Git
    **do not translate safely into SVN**.

    Problems have occurred in the past due to this mismatch.
    They are not theoretical.

The Problem
-----------

In a pure Git server, Git is very good at protecting history. Git history is represented as a **directed acyclic graph (DAG)**.
Multiple branches and merges are first-class concepts, and Git preserves all
ancestry information.

On our server, however, Git history is mirrored into **SVN**, which has a **strictly linear and immutable history**.

When these two models collide, **valid Git operations can be flattened into SVN
in ways that lose structural information**, especially around merges.

The issue is not that Git rewrites history by itself, but that **SVN cannot
represent Git’s merge topology**.


How Git Merges Actually Work
----------------------------

A Git merge can occur explicitly (``git merge``) or implicitly
(e.g. via ``git pull`` when local commits exist).

Consider the following situation:

* You start from commit ``A``
* You commit ``B`` locally (no push yet)
* While you work locally, someone pushes commit ``C`` to ``master``
* You then pull and push

Before the merge, the server history is:

.. code-block:: text

    C
    |
    A

Your local history is:

.. code-block:: text

    B
    |
    A

After a **Git merge**, the history becomes:

.. code-block:: text

          M   (merge commit)
         / \
        B   C
        |
        A

Here, **nothing is rewritten**:
all commits (``A``, ``B``, ``C``) remain intact, and the merge commit ``M``
records both parent histories.

This is correct and expected Git behavior.

Why This Fails with SVN
----------------------

SVN does not store merge topology.
It expects a **single immutable line of commits**.

When the Git–SVN bridge linearizes the above history, it must choose one path.
In doing so, **merge relationships are discarded**, and the resulting SVN history
may appear reordered or incomplete.

This is where confusion and apparent “history rewriting” arise —
not in Git, but during **SVN linearization**.

Rebase vs Merge in This Context
-------------------------------

A rebase produces a linear history by construction.

Using the same example, rebasing ``B`` onto the updated ``master`` results in:

.. code-block:: text

    B
    |
    C
    |
    A

This history is linear and can be represented safely in SVN.
For this reason, **rebasing avoids Git–SVN translation issues**, even though it
rewrites `local`` Git history (in this example, commit B is reapplied after C,
even though it was originally created before it). 

This reordering is explicit, predictable, and compatible with SVN’s linear
model, whereas merge topology is not.

Merges, while correct in Git, **preserve topology that SVN cannot encode**.


What Is Actually Dangerous
--------------------------

On this server, the following actions are problematic **on the main working
branch**:

* ``git pull`` without ``--rebase`` when local commits exist
* Merge commits on ``master`` that introduce multiple parents
* Any operation that rewrites *already published* history
  (force-push, amend, rebase after pushing)

These are not *bad Git practices*.
They are simply **incompatible with a Git repository that is mirrored into SVN**.

When these situations occur, the result is **history corruption**:
commits may disappear from the visible history, timelines can appear to fork
incorrectly, and parts of the repository history may become effectively
unreachable from SVN.

Once this happens, recovery is painful and error-prone, and in some cases
not possible without manual intervention or loss of information.


Practical Rule
--------------

.. important::

    On ``master``:
        **Avoid merge commits. Prefer rebasing.**

    On feature branches:
        **Use normal Git workflows.**

In practice:

* Use ``git pull --rebase`` on ``master``
* Integrate branches carefully, being aware of SVN limitations
* If unsure, inspect your history before pushing


If in doubt:

.. code:: bash

    git status

If you are confused:

.. code:: bash

    git status

If something looks strange:

.. code:: bash

    git status

You can also try: 

.. code:: bash

    git status
    git log --graph --oneline --decorate



Useful Documentation
====================


* `GitHowto <https://githowto.com/>`_ : step by step guide

* `Git tutorial <https://www.geeksforgeeks.org/git/git-tutorial/>`_

Configuration
=============

Before your computer can communicate with a remote Git host, it needs some form of identity badge — otherwise the server just sees an unknown lifeform trying to break in. Depending on the host, this badge might be an SSH key, a personal access token, or the old-fashioned username and password combo.

For GitHub, SSH keys are the preferred psychic paper: once set up, you can clone, pull, and push without reintroducing yourself every five minutes. Other servers have their own customs, but the idea is the same — prove who you are so your time-traveling code changes don’t get rejected on arrival.

.. _howtogit-github:

GitHub
------

`GitHub <https://github.com>`_ is one of the most widely used Git hosting platforms, offering remote repositories, collaboration tools, and SSH-based access.


SSH keys
^^^^^^^^

* Generate a key if you do not already have one on your computer:

  .. code:: bash

    $ ssh-keygen -t ed25519 -C "email@email"    


* Add the key to GitHub. Better to check GitHub documentation: `How to add SSH key <https://docs.github.com/en/authentication/connecting-to-github-with-ssh/adding-a-new-ssh-key-to-your-github-account>`_

But here is the current procedure:

1. Copy the SSH public key to your clipboard:

    .. code:: bash

        $ cat ~/.ssh/id_ed25519.pub

2. Go to GitHub -> Settings -> SSH and GPG keys -> New SSH key

3. Add a title and paste the public key content.

Cache your credentials
^^^^^^^^^^^^^^^^^^^^^^

To cache your password for #### seconds, you can configure: 

.. code:: bash

    git config --global credential.helper 'cache --timeout=####'

Pencil Code Credentials
------------------------

For |PC|, there’s a small twist: the gate looks like GitHub’s, but the key is different.
Even though the repository lives on |github|, you can’t use the usual GitHub credentials for write access.
To obtain permission to push changes, follow the steps in :ref:`download-writeaccess`.

Cloning a Repository
====================

This section covers how to clone a generic repository — we’ll call it :file:`project`.
If you want to clone the |PC| repository specifically, see :ref:`download`.


If the repository is **public** (for example, a public project on GitHub), you can clone it without any credentials or access configuration:

.. code:: bash

    $ git clone https://github.com/project/project.git 

This gives you your own local copy of the code universe, ready for exploration, but you **won’t be able to push any changes** unless you have the appropriate access (SSH keys, token, or username/password, depending on the service).


.. note::

    Congratulations, you are now in **Look-but-don’t-touch mode**. Feel free to explore the code, but the timeline is locked until you have write access!


If you do have access configured (for example, via SSH keys), you can clone using the SSH URL:

.. code:: bash

    $ cd my-dir
    $ git clone git@github.com:project/project.git

   
If using HTTPS instead of SSH:

.. code:: bash

    $ git clone https://github.com/project/project.git


The repository folder will now appear in your directory.
Think of it as opening a portal into the repository’s timeline.


Checking Status & Differences
=============================

Checking the status of the git repository:

.. code:: bash

    $ git status

It will show:

* The branch you are working on

* Untracked files (in red if any): files not tracked by Git

* Changes not staged for commit (in red if any): files modified but not added to the commit list

* Changes to be committed (in green if any): files staged but not yet committed


If you have a bunch of untracked files that you don’t want cluttering your status output,  
you can ignore them by using:


.. code:: bash

    git status -uno


Checking Differences
--------------------

Checking for differences with local repository:

.. code:: bash

    $ git diff

Check staged differences:

.. code:: bash

    $ git diff --staged

.. note::

    Peek into the timeline before changing history

Check differences between your local branch and the remote branch:

.. code:: bash

    $ git diff origin


Pulling & Stashing
==================


.. tip::

    Keep your timeline up to date—pull before working, stash experiments if needed.

* *Do not use* the basic pull if you have local changes:

.. code:: bash
        
    $ cd my/git/dir
    $ git pull

* Recommended pull (especially if you have unsynchronized changes):

.. code:: bash
    
    $ git pull --rebase

Working with multiple contributors may result in overlapping changes.  
The ``rebase`` option reapplies your commits on top of the latest changes from the remote branch, keeping a linear history without unnecessary merge commits.  

It works smoothly if changes do not overlap. Otherwise, **don't panic!** Everything has a solution.


Keeping Uncommitted Changes
---------------------------


Sometimes you’re working on something experimental, but suddenly you need to pull updates from the remote or switch branches. You don’t want to commit half-baked changes, and you don’t want to lose your work. Enter `git stash`—your own little **time-travel pocket dimension** for code.  


* Temporarily protect local changes before pulling or pushing:

.. code:: bash

    $ git stash        # hide your uncommitted changes

Now you can pull or push safely.

* List of your stashed experiments: 

.. code:: bash

    $ git stash list   # see all your stashed experiments

    $ git stash apply  # restore the latest stash without removing it


* Restore uncommitted changes:


.. code:: bash

    $ git stash pop   # restore the latest stash and remove it from the stash list

* Restore the latest stash:


.. code:: bash


    $ git stash apply  # restore the latest stash without removing it


.. note::

    Think of stash as hiding experiments in a TARDIS pocket dimension.


Advanced stash tips:
^^^^^^^^^^^^^^^^^^^^

* Name your stash to remember what’s inside:

.. code:: bash

    $ git stash push -m "experiment with time loops"

* Stash only specific files:

.. code:: bash

    $ git stash push path/to/file1 path/to/file2

* Drop a stash you no longer need:

.. code:: bash

    $ git stash drop stash@{0}

.. note::

    Use stashes wisely—too many, and your TARDIS starts to feel cluttered.




Staging Changes
===============

Before your changes can travel to the master timeline (the remote repository), Git requires a **pre-flight check**: this is the staging phase. Think of it as placing your edits into a sonic-proof capsule before sending them through the TARDIS.

.. code:: bash

    $ git add file_to_commit    # stage a single file
    $ git add .                 # stage all changes in current directory
    $ git add dir_to_add/       # stage all files in a specific folder

.. note::

    Staging lets you **choose exactly which changes** go into your next commit. You can have some edits ready for the next time jump while leaving experimental work behind.

.. changed

Pro tip: use `git status` after staging to double-check what’s staged and what’s still wandering in the timeline uncommitted:

.. code:: bash

    $ git status

.. note::

    This prevents “Oops! I committed that half-baked code” moments—every Time Lord needs a careful plan before hopping timelines.

Advanced tip: you can stage multiple sets of changes separately and then commit each with a different message. This lets you break your work into logical, focused commits instead of dumping everything into one messy time capsule.


.. code:: bash

    # Stage first set of changes (a file and a directory)
    $ git add file1.py
    $ git add big_dir/
    $ git commit -m "Implementing feature X"

    # Stage second set of changes (just a file)
    $ git add file2.py
    $ git commit -m "Fixing bug in feature Y"

    # Stage third set of changes (two files)
    $ git add file3.py
    $ git add file4.py
    $ git commit -m "Updating documentation"

.. note::

    Each `git add` is like sealing a small time capsule, and each `git commit -m` sends all the added files and directories safely into the master timeline. Your commit history will be clean, readable, and easy to navigate.


Interactive Staging with `git add -p`
-------------------------------------

Sometimes you’ve been tinkering in the same file and only part of your changes are ready for the next commit. Enter **interactive staging**:

.. code:: bash

    $ git add -p file_to_commit

This command will break your changes into **hunks** (chunks of modified lines) and ask you what to do with each:

* **y** – stage this hunk
* **n** – do not stage this hunk
* **s** – split the hunk into smaller pieces
* **q** – quit, do nothing
* **?** – show help

.. note::

    Think of `git add -p` as using a sonic screwdriver to precisely select which edits travel through time. You can send just the ready parts while leaving experimental changes safely behind.

.. changed

Pro tip: use this for clean, logical commits. You’ll thank yourself (and future developers) when browsing `git log`.



Pushing Changes
================

.. attention::

    Always pull (preferably with rebase) before pushing to avoid paradoxes.


Normal push sequence:

.. code:: bash

    $ git pull --rebase                       # update first!
    $ git add file_to_commit                  # stage the file 
    $ git commit -m "message of the commit"   # comment for the posterity
    $ git push                                # push to remote

and voilà!

.. admonition:: Don't panic!

    If this doesn't work... don't panic... check possible solutions in `Conflicts`_.

Discarding / Restoring / Canceling Changes
==========================================


Discarding Local Changes
------------------------

To discard local modifications:

.. code:: bash

    $ git restore working_on_it


Canceling Staged Changes
------------------------

Before committing staged changes:

.. code:: bash

    $ git restore --staged working_on_it

This will unstage changes without modifying the local file. To fully restore, refer to `Discarding Local Changes`_.


Canceling a Commit
------------------

This will undo the last commit (use with caution):

.. code:: bash

    $ git revert HEAD


Moving Files & Directories
==========================

Moving directories or file with git can be a bit tricky. The easiest way  (always check that your version is up to date beforehand!) is using Git itself:



.. code:: bash

    $ git mv <source> <destination>
    $ git commit -m "move directory/file to another location/name"
    $ git push


.. note::

    No ``git add`` needed. Teleport files like a sonic screwdriver.



Branching
=========


Branching is like opening an alternate timeline where you can experiment, build features, or break things gloriously *without* endangering the master universe (``master``). The idea is to keep these branches short-lived and focused—if your branch lasts longer than some house plants, you might actually be developing a completely different project.

When your work is done, you merge your branch back into ``master`` and pretend everything went according to plan.

Before you start, it's wise to check where you are:

.. code:: bash

    $ git status

The master branch is called ``master``. Feature branches can be named however you like—ideally something more helpful than ``new-stuff`` or ``pls-work``.

Basic commands:

* List all local branches:

    .. code:: bash

        $ git branch

* Create a new branch:

    .. code:: bash 

        $ git branch my-branch

* Switch to an existing branch:

    .. code:: bash

        $ git checkout my-branch

    or 

    .. code:: bash

        $ git switch my-branch

* Create and switch to a new branch:

    .. code:: bash

        $ git checkout -b my-branch

    or

    .. code:: bash

        $ git switch -c my-branch


* Rebase onto another branch:

    .. code:: bash

        $ git rebase my-branch

    Careful with this one. Can generate conflicts.

* Delete a branch, but only if it has been fully merged.

    .. code:: bash

        $ git branch -d my-branch

* Forcefully deletes a branch (use with care!)

    .. code:: bash

        $ git branch -D my-branch


* Merge into ``master``:

    .. code:: bash

        $ git switch master
        $ git merge my-branch

    .. attention:

        This merge will not work with the |PC|, please check sec :ref:`merge_pencil`


.. important::

    Always ensure you know which branch you are on before committing, pulling, or pushing.

Tips for working with Branches
------------------------------

A classic branching horror story goes like this: you create your branch, happily work on your changes for a while, and when you finally try to rebase onto ``master``, you discover that ``master`` has evolved into a completely different timeline. Now you’re staring at a kaiju-sized merge conflict wondering if you should fake your own death and start a new career.

To avoid this future therapy bill, the best practice is to regularly merge ``master`` into your branch:

.. code:: bash

    $ git switch documentation  # make sure your are on your branch
    $ git merge master          # merge master into your branch

By doing this often, any conflicts you hit will be smaller, friendlier, and less likely to question your life choices.

If you keep merging as you work, merging your branch later will feel less like boss-level combat and more like a polite handshake.



Pushing branches
----------------

Most of the time, you’ll work on your feature branch locally and then merge it into ``master`` when everything is ready. However, sometimes you need to **share your branch with others**, create a **pull request**, or simply **back it up to the remote repository**.

When you push a branch to the server **for the first time**, Git doesn’t know where to send it yet. So you must explicitly set the upstream:


.. code:: bash

    $ git push --set-upstream origin documentation

From that moment on, Git will remember the connection between your local ``documentation`` branch and the remote one, so you can simply:

.. code:: bash

    $ git push

.. note::

    The first push is like introducing your branch to the server: *"Hello, I exist now!"* — after that, Git will remember the relationship and stop asking awkward questions.


.. _merge_pencil:

How to merge your branch with the |PC| master
----------------------------------------------



Merging in the |PC| universe isn’t your regular “two lines diverged in a repo” situation.  
Because |PC| exists in a peculiar hybrid space-time where both ``svn`` and ``git`` coexist (through the miracle—or curse—of SubGit), every interaction with the repository must go through the central server at `<https://pencil-code.org>`_.  

This means that a normal merge won’t work. You need to follow the proper temporal protocols.

To keep your branch from tearing a hole in the space–code continuum, proceed as follows:



1. **Synchronize your branch with master — align your timelines**

    .. code:: bash

        $ git switch your-branch   # make sure you are on your branch
        $ git merge master         # merge latest timeline updates

    Congratulations, your branch is now aligned with the latest master timeline.
    Reality remains stable—for now.


2. **Merge into master — but not the fast-forward kind**

    A fast-forward merge may look tempting: quick, simple, elegant.  
    Unfortunately, in the |PC| multiverse, it’s also forbidden. SubGit guards the gate and will smite any attempt to rewrite the sacred SVN trunk.

    So instead, perform a :command:`non Fast-Forward merge` — the Git equivalent of gently folding timelines together rather than shoving one into the other.


    .. code:: bash

        $ git switch master             # make sure you are on master
        $ git merge your-branch --no-ff # no Fast forward, no paradoxes

    This will keep the history intact and prevent the repository from imploding into a causal loop.



3. **Push your changes to the central repository**

    .. code:: bash

        $ git push

    
    If everything worked, your branch is now part of master, history is safe, and you’ve successfully avoided the “Temporal Merge Conflict of Doom.”



The merge failed! (or, “I think we broke the timeline...”)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


If you ignored the “no Fast-Forward” prophecy and pushed anyway,  
Git will retaliate with an ancient curse that looks like this:


.. code:: 

    remote: 
    remote: SubGit ERROR REPORT (SubGit version 3.3.17 ('Bobique') build #4463):
    remote: 
    remote: You've received this message because SubGit (http://subgit.com/) is installed in your repository
    remote: and an error that needs to be dealt with has occurred in SubGit translation engine.
    remote: 
    remote: The following ref update is disallowed:
    remote:   refs/heads/master: leads to replacement of SVN branch 'trunk'
    remote: 
    remote: If changes were forcefully pushed to Git repository, try to merge them with the upstream instead;
    remote: If changes were result of fast-forward merge, retry merge with --no-ff option.
    remote: 
    remote: You can allow branch replacements by adjusting SubGit configuration file as follows:
    remote:   'svn.allowBranchReplacement = true' in remote mirror mode;
    remote:   'git.<ID>.allowBranchReplacement = true' in local mirror mode.
    remote: 
    usage: git credential-cache [<options>] <action>

        --[no-]timeout <n>    number of seconds to cache credentials
        --[no-]socket <path>  path of cache-daemon socket

    git credential-cache --timeout=9999

     store: 3: store: not found
    To https://pencil-code.org/git/
     ! [remote rejected]     master -> master (pre-receive hook declined)
    error: failed to push some refs to 'https://pencil-code.org/git/'

Don’t panic. The timeline can be repaired.


**Steps to fix your mistake and restore the flow of time:**

1. **Rewind to before the paradox**

    First, make sure you’re standing on the ``master`` branch (``git status`` will confirm your position in time).

    .. code:: bash

        $ git reset --hard origin/master  # return to the moment before the merge

2. **Update master — in case someone else tinkered with the timeline**

    .. code:: bash

        $ git pull

3. **Merge again, correctly this time**

    .. code:: bash

        $ git merge your-branch --no-ff

4. **Push, and watch as the timelines gracefully align**

    .. code:: bash

        $ git push

If you followed these steps, the merge should succeed and the repository will continue to exist in a stable reality.  

.. admonition:: Remember: 
    
    *merging with care is cheaper than rebuilding the universe.*  
    And whatever you do—never fast-forward past a fixed point in time.



History / Log
=============


Think of ``git log`` as the journal of your time-travel adventures: every change, every experiment, every “oops” that you later rewrote into a perfectly reasonable commit message. It lets you see what happened, when it happened, and who to glare at (even if it's just past-you).


Get a list of changes:

.. code:: bash

    $ git log

Some options:

* One line history and some options:

.. code:: bash

    $ git log --oneline         
    $ git log --oneline --max-count=2
    $ git log --oneline --since="5 minutes ago"
    $ git log --oneline --until="5 minutes ago"
    $ git log --oneline --author="Your Name"
    $ git log --oneline --all
    $ git log --pretty=format:"%h %ad | %s%d [%an]" --date=short





Conflicts
=========

Ah, Git conflicts—the stuff of nightmares that makes seasoned developers break out in cold sweats. Don’t worry, you’re not alone; I panic too.  

The good news is that most conflicts are avoidable if you follow a few simple rules of time-travel hygiene:  

* Always check your ``git status`` to know exactly which branch you’re meddling in.  
* Pull the latest changes before making your own edits.  
* Prefer ``rebase`` over messy merges whenever possible.  
* Read Git’s error messages carefully—they are surprisingly good at telling you exactly what to do (and they won’t judge you for your past mistakes).  

Follow these, and you’ll face fewer conflicts, less panic, and a lot more sanity.

Ignore these rules at your own peril: suddenly you’re in a parallel universe of code, facing monstrous conflicts that make you question every life choice, swear at your computer, and consider rewriting the project in interpretive dance instead of text.

Common Git Conflicts
--------------------

1. **Simple line conflicts**  
   Two changes on the same line. Resolve manually, then `git add` and continue.  
   .. note:: Imagine your past self arguing with your present self.

2. **File deleted vs. modified**  
   One deleted a file, another changed it. Decide if the file lives or dies.  
   .. note:: Like erasing a timeline — TARDIS advised.

3. **Directory vs. file**  
   A folder appears where a file existed. Rename or move one to resolve.  
   .. note:: Parallel universe tried to overwrite your living room with a closet.

4. **Multiple commits changing same lines**  
   Happens when rebasing long-lived branches. Resolve incrementally.  
   .. note:: Untangle the time knots carefully, one thread at a time.

5. **Binary files**  
   Git cannot merge them. Pick one version manually.  
   .. note:: Binary files are like Daleks — they don’t negotiate.


Step-by-Step Conflict Resolution
--------------------------------


``Push`` did not work
^^^^^^^^^^^^^^^^^^^^^^^^^^

After adding and committing your files, you tried to push your changes
and got the dreaded error:

.. code:: bash

    $ git push
        To https://pencil-code.org/git/
         ! [rejected]            master -> master (fetch first)
        error: failed to push some refs to 'https://pencil-code.org/git/'
        hint: Updates were rejected because the remote contains work that you do not
        hint: have locally. This is usually caused by another repository pushing to
        hint: the same ref. If you want to integrate the remote changes, use
        hint: 'git pull' before pushing again.
        hint: See the 'Note about fast-forwards' in 'git push --help' for details.

This happens when someone else has updated the remote branch since you last
pulled. Git is politely asking you to reconcile timelines before pushing
your changes — basically, don’t try to overwrite someone else’s work
with a vortex manipulator.
To fix it, pull the remote changes and rebase your commits on top:

.. code:: bash

    $ git pull --rebase
    $ git push

Example scenario:

* You added a new function ``compute_flux()`` in ``hydro.f90``.

* Meanwhile, a colleague added ``update_boundary()`` to the same file
  and pushed it.

* ``git push`` will be rejected until you ``git pull --rebase`` and
  integrate your function with theirs.

* If both edits touch the same lines, Git will pause and ask you to
  resolve conflicts manually — the next bullet points will guide you
  through that process.
This method works perfectly if your changes don’t overlap with the
remote edits. Otherwise, brace yourself for some conflict resolution fun.



Rebase paused due to conflicts (same lines touched)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If your edits overlap with the remote changes — for example, both you
and a colleague modified the same line in ``hydro.f90`` — Git will
pause the rebase and flag a conflict:

.. code:: bash

    $ git status
    # both modified: hydro.f90

Git inserts conflict markers in the file, like this:

.. code:: text

    <<<<<<< HEAD
    your change here
    =======
    colleague's change here
    >>>>>>> branch-to-rebase

At this point, you have to decide how to merge the two edits. Options:
* Keep your change, discard theirs.

* Keep theirs, discard yours.

* Combine both changes intelligently.

Once resolved, mark the file as resolved and continue the rebase:

.. code:: bash

    $ git add hydro.f90
    $ git rebase --continue

Then verify your changes:

.. code:: bash

    $ git log --oneline

And finally, push the integrated timeline:

.. code:: bash

    $ git push --force-with-lease origin master

Resolving conflicts when merging branches
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
* Abort the merge:

    .. code:: bash 

        $ git merge --abort
        $ git status
* Resolve the conflict by editing files and committing:

    .. code:: bash

        $ git add resolved_file
        $ git commit





.. admonition:: Remember

    Remember: conflicts may feel terrifying, but with careful time-travel hygiene, they are just minor bumps in the TARDIS ride of development.



Advanced topics
================

Pro Tips
-----------


A few extra moves that make you feel like a Git Time Lord:

* **.gitignore** – prevent unwanted files from sneaking into your timeline:

.. code:: bash

    # Example .gitignore
    *.log
    *.tmp
    

.. note::

    Think of it as shielding Daleks and temporary logs from your timeline.

* **Undo a commit** (`git reset`) – sometimes past-you made a mistake:

.. code:: bash

    $ git reset HEAD~1  # undo last commit but keep changes
    $ git reset --hard HEAD~1  # undo last commit and discard changes

.. note::

    Like a mini TARDIS to erase recent misadventures.

* **Check remotes** (`git remote -v`) – know which time portals your repo talks to:

.. code:: bash

    $ git remote -v

.. note::

    Useful before pushing to avoid accidentally sending code to a parallel universe.

.. _git-svn-postmortem:


Show git branch in the bash prompt
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

When working with Git, it is extremely useful to know **which branch you are on at all times**.
Displaying the current branch directly in your shell prompt helps avoid accidental commits
to the wrong timeline.

The configuration below has been tested on **GNU/Linux** systems using Bash.

To enable this feature, add the following lines to your :file:`.bashrc` file:



.. code:: bash

    #################
    # Show git branch in the bash prompt when in git-controlled directory   
    #
    source /etc/bash_completion.d/git-prompt

    # Colors
    RED="\[\033[0;31m\]"
    GREEN="\[\033[0;32m\]"
    YELLOW="\[\033[0;33m\]"
    BLUE="\[\033[0;34m\]"
    RESET="\[\033[0m\]"

    # Prompt with color and git branch
    export PS1="\u@\h${RESET}:${YELLOW}\w${RESET}${RED}\$(__git_ps1 ' (%s)')${RESET}\$ "
    ##########


Git on an SVN-Backed Server: A Post-Mortem
--------------------------------------------

This section is for advanced users who want to understand **how Git history can be
rewritten or lost on this server**, even without force-pushing or obvious errors.

No commands shown here are *wrong Git*.
They are simply incompatible with the Git–SVN bridge used by this repository.

This is a post-mortem, not a how-to.

What Actually Happened
^^^^^^^^^^^^^^^^^^^^^^^

In multiple incidents, parts of the repository history were lost or rewritten after
apparently normal Git operations.

The common pattern was:

* A developer pulled remote changes using ``git pull`` (without ``--rebase``)
* Git created a merge commit locally
* That merge commit was later pushed to the server
* The Git–SVN bridge attempted to linearize the resulting DAG
* Commits that could not be represented in SVN were silently dropped

No force-push was involved.
No error was reported.
From Git’s point of view, everything succeeded.

From SVN’s point of view, history was rewritten.

Why Git Could Not Protect You
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Git assumes that the remote endpoint understands Git’s commit graph.
SVN does not.

SVN cannot represent:

* Multiple parents
* Non-linear history
* Rewritten commit ancestry

When Git history containing merges is pushed through the bridge,
the bridge must **choose a linearization strategy**.

That strategy is not guaranteed to preserve all commits.

In other words:

Git did exactly what it was told to do.  
SVN accepted what it could understand.  
Everything else was discarded.

Why This Is Worse Than a Force-Push
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

A force-push is loud.
It rewrites references explicitly and is usually blocked or noticed.

This failure mode is silent.

* The push succeeds
* No warnings are emitted
* The repository appears healthy
* Missing commits are only noticed much later

By the time the problem is discovered, local clones may already have diverged.

Lessons Learned
^^^^^^^^^^^^^^^^^


* Git’s safety guarantees stop at the Git boundary
* A Git–SVN bridge is a translation layer, not a time machine
* Merge commits are fundamentally incompatible with SVN history
* Linear history is not a preference here — it is a requirement

This is why the documentation insists on:

* ``git pull --rebase``
* No merges on the main branch
* Frequent use of ``git status``

These are not stylistic choices.
They are **damage control**.

Final Note
^^^^^^^^^^^

If you are accustomed to advanced Git workflows (feature branches, merge commits,
interactive rebases), be aware that this server **cannot support them safely**.

The Prime Timeline here is linear.

Alternate timelines belong in local clones only.