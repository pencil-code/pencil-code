.. _python-pastmeetings:

Past Meetings
=============

Meeting archive for the Python postprocessing development group.

Each entry includes date, participants, and a summary of decisions.



.. rubric:: 2021 Sept 17

**Participants:** 

**Minutes:**

1. Rules for adding new routines and extending existing routines.

In the long term we should have test functions. (WD, SC)
Build on the existing test in pc/python/pencil/tests
Send regular emails to PC Python group to stick to coding rules. (IRL, SC)


2. Unified format for doc string.

Have all sections for every function:
General definition, Signature, Parameters, Returns, Return Type, Examples, (Notes)
Add section to the style guide about the documentation. (IRL)
Send email about unified doc string. (IRL, SC)


3. Usage of Black.

Add to the style guide that everyone should run Black when
making changes and additions. (IRL)
Send an email to the PC Python discussion group. (WD)


4. pv plots: should we extend?

Email authors of the visu/pv_* routines. (SC)



.. rubric:: 2021 Apr 14

**Participants:** Illa, Fred, Simon, Wolfgang,

**Minutes:**

* Wolfgang: updated us on the following topics:

1. import cycle: 

- he develop a script to check this isue 'python/utils/cyclic-imports'
- he did a lot of clean-up
- extracted 'is_sim_dir()' to a new package 'pencil/util' that must not import from other packages.
- removed 'lapack/pidly/pidly.py'

2. started to write tests

3. tried to fix 'pencil/calc/example_shocktube.py'

4. reformatted 'pencil/calc/streamlines', 'pencil/calc/tensors' and 'pencil/read/power' with *black*

We discussed the usage of *black* to format the python files, and agreed on using it as a default.
*black* can be used by default in several text editors.

5. Imports

We discussed the replacement of
'from ..sim import simulations' with 'from pencil.sim import simulations'
the only problem this might have is any change of jerarchy on the different modules.

6. Backpack:

the current definition of this library is "Check for standard library"
we discussed to rename it to 'third_party'

7. Problems of 'import pencil as pc'

Currently it takes ~4 seg to import the pencil package, and this might be a problem if one just
want to run a test.
We discuss the idea of having a separate package structrue for an interactive 
use.

8. Description of submodules:

We expend quite some time discussing if we should move the description of submodules into respective `__init__.py` (from `pencil/__init__.py`)

- use of black
- rename backpack to third party?
- Replace 'from ..sim import simulations' with 'from pencil.sim import simulations'
- Distinguish between interactive and library usage ?

* Illa

9. Mailing list: 

we encourage everyone interested to join out mailing list

10. Newsletter:

We discussed the contents of the newsletter and the information we want to write

11. Files Documentation:

We discussed a general idea for the formatting of a python file documentation.
I suggested to follow the outline of 'python/pencil/read/power.py'


* Simon

12. He remove the file npfile.py, as discussed in previous meetings.

13. He changed read.var


.. rubric:: 2021 Feb 24

**Participants:** 

**Minutes:**

Topics:
* Python <-> IDL descriptions

We added this to the wiki: 
https://github.com/pencil-code/pencil-code/wiki/Pencil-IDL------Python-guide
Please contribute!

* Documentation: sphinx creation of the documentation

We discussed where to post the documentation and agreed on using readthedocs.org integrated with github.

* npfile.py is used in animate_slices and appears to be deprecated in scipy, Simon will check this issue.

* get_sim: Wolfgang is trying to get rid of it





.. rubric:: 2021 Jan 14

**Participants:** 

**Minutes:**


Topics:

* Replace package names to non-conflicting  ones?

* Documentation: use sphinx

* Imports: import in the functions or in the header

* Cyclic dependence structure in different imports.

* Newsletter



.. rubric:: 2020 Nov 23

**Participants:** Fred, Simon, Wladimir, Wolfgang, Alberto, Illa

**Minutes:**


We agree on going over the issues on the files: 
-pencil-code/python/meetings/Nov20_topics
-pencil-code/python/pencil/TO_DO_LIST.txt

And discuss the issues:
 
* Where to import modules: 

Right now the modules are imported inside the functions, but according to the 
PEP8 guidelines they should be imported at the beginning of the file.

Simon proposed to leave it as it is, since it was agreed upon on a Pencil Code 
meeting 2 years ago, and is cleaner. Plus it is less messy when using 
autocomplition in ipython.

Wolfang and Illa proposed to move to the PEP8 convention and find a solution 
for ipython.

* Reload:

Functions won't reload in an ipython session.
Fred pointed out that the file names are the same as functions names, making 
it imposible to reload, and suggested to change file names.
We proposed different ways of reloading, using importlib, but it didn's seem to 
work.

* Standarize formatting:

Wolfgang suggested to use "black" for formatting the code. Also as a way to 
standarize the formatting. 
We agreed to give it a try, and he would give us more info on the next meeting.

* Automatic documentation:

Illa suggested to use either Epydoc or sphinx to create automatic documentation 
of the python code.
We agreed on the idea, and Illa will research the different options, keeping in 
mind doing the minimal work on the already existing documentation inside the 
functions.

* Tutorials:

We discuss the use of Jupyter versus plain scripts to add into the tutorials, 
and agreed on having both version.
The Jupyter notebooks can be translated into scripts using pandoc or jupyter nboconvert --to script

Alberto will try to add new tutorials.
We also agree on changing the dir names inside the tutorials to make them more 
relevant to the content. Creating a "Getting started" dir and examples of 
different topics.
Simon will try to change the look of the tutorials dir.

* Functions:

We will try to add the relevant IDL functionality.
We agree on adding the missing IDL functions to the TO_DO_LIST and keep a check 
there.

