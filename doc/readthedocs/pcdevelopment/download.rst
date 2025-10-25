.. _download:

************************************
Download the Pencil Code
************************************


The Pencil Code **lives** on GitHub — but unlike most repositories, it does **not socialize** using normal GitHub authentication. You can access it using **SVN** or **Git**. Yes, both. Because history is messy and developers have opinions. 

Depending on your intentions, you can access the code in two distinct modes:

* **Read-only mode** — *observe the code like a respectful museum visitor.*
* **Write access** — *join the crew and start rearranging the furniture.*


For general time-travel concepts like branches, commits, and rebasing paradoxes, see :ref:`howtogit`.  
This document focuses specifically on obtaining the Pencil Code using your preferred weapon of choice.

.. note::
   If you're wondering *“SVN or Git — which one should I use?”*  
   The answer is simple:

   **SVN** — stable, linear, sensible. Like eating with a fork. 
    
   **Git** — interdimensional branching chaos. Like eating with chopsticks while skateboarding.

   *(Both approaches are valid lifestyles.)*


Below is the **Pencil-Code-specific access protocol** for both read-only and write access.



Read Only (a.k.a. “Look but Don’t Touch” Mode)
==============================================

If you just want to explore the code without leaving a trace, you can clone or checkout anonymously. No identity badge required.

*SVN (for historians and minimalists):*

.. code:: bash

    svn checkout https://pencil-code.org/svn/trunk pencil-code

*Git (for people who enjoy rewriting history even when they don't have write access):*

.. code:: bash

    git clone https://pencil-code.org/git/ pencil-code

Congratulations! You now possess a **spectator copy**. Browse freely.  
Modify locally if you like — but **pushing changes will be politely rejected**.

.. tip::
   Think of this as **read-only Jedi training mode**. You may swing the lightsaber, but it won’t cut anything yet.


Write Access (a.k.a. “Let Me In, I Want to Contribute!”)
========================================================

So you want to **actually contribute**? Excellent — but that means the system must **know who you are**.

Although the repository *lives* on GitHub, the usual GitHub SSH key setup (see :ref:`howtogit-github`) **will NOT work here**.

Instead, follow the official Pencil Code procedure:

1. Register at `<https://account.pencil-code.org>`_
2. Apply for write access to the ``main`` repository
3. Ideally, choose the **same USERNAME as your GitHub account** — unless you enjoy future identity crises

Once approved, use the appropriate access method:

*SVN (steady and predictable):*

.. code:: bash

    svn checkout --username=USERNAME https://pencil-code.org/svn/trunk pencil-code

*Git (branch enthusiast edition):*

.. code:: bash

    git clone https://USERNAME@pencil-code.org/git/ pencil-code

.. warning::
   After write access is enabled, **use your powers wisely**.  
   SVN users will quietly carry on as before.  
   Git users will immediately create three branches and rewrite history twice before lunch.


You are now equipped to enter the Pencil multiverse — choose your timeline wisely.


What Next?
==========

You now have the code.  
To actually *work* with it, hop over to :ref:`howtogit` for branching, staging, rebasing, and mild emotional breakdowns.

Happy hacking — and may your merges be ever conflict-free.
Remember: with great power comes great responsibility… and occasional merge conflicts.
