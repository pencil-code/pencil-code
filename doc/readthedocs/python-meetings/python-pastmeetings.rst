.. _python-pastmeetings:

Past Meetings
=============

Meeting archive for the Python postprocessing development group.

Each entry includes date, participants, and a summary of decisions.

2026 Jan 29
------------

**Participants:** 

Alberto Roper Pol, Illa R. Losada, Simo Candelaresi, Kishore, Reihaneh

**Minutes:**

1. Points in the agenda:

* Development of Python tests on Pencil Code and exploring code coverage tools.

* Presentation of tutorials for particular samples developed by each person on their favourite sample/exercise with Pencil Code.

* Updates on the documentation (tutorials and updates on readthedocs)

* Inconsistencies in the array formatting with different I/O options.

* Check :file:`pencil-code/python/pencil/TO_DO_LIST.txt`. The last update is from Jun 29, 2024, is there any further update?

* Changing behaviour of .keys() methods (of Averages, Timeseries etc.). Currently it just prints the keys, but it would be more useful to return a list

* Decide on a convention for extra debug output (controlled by the ``quiet`` keyword in functions like ``pc.read.var`` and ``pc.read.grid``). Currently, whether it defaults to ``False`` or ``True`` varies from function to function, which is confusing and annoying. Ideally we would control such output by a module-wide flag (such that something like ``pc.shut_up = True`` at the beginning of a script would have the effect of setting ``quiet=True`` in all the ``pc`` functions).

- Changing the behaviour of ``pc.io.get_value_from_file``. Currently, whether it has been successful of not is indicated by a boolean return value (as if it is a C program). I think it makes more sense to explicitly raise an error message. I presume most users would want to know that the function is failing, rather than dealing with crytic errors later on.
- Rename ``pc.sim.__Simulation__`` to ``pc.sim.Simulation``. Leading underscores usually signify that something is not meant to be publicly accessed; why do we not want users to directly use this class?
- Supported Python versions: is it enough to go by the versions currently supported by upstream (<https://devguide.python.org/versions/>), or do we want to support Python versions older than 5 years?
- The item **Interactive vs. library mode** in :file:`pencil-code/doc/readthedocs/python-meetings/python-decisions.rst` is now obsolete due to the implementation of lazy loading.

- Folders in `python/` (like `meetings` and `docs`) are detected as Python modules if `$PENCIL_HOME/python` is added to `$PYTHONPATH`. Options to avoid this:




2025 Nov 26
------------

**Participants:** 
Alberto Roper Pol, Carolina Ortiz, Oindrila Ghosh, Kishore, Reihaneh, Illa Losada, Simon Candelaresi, Alireza, Mgurgeni, Nordita Participants, Isak Stomberg. 

**Minutes:**

1. Meeting Format.

a) Frequency of the Pencil Code Python Meetings should be every two months and should be kept up.
b) Rotate organizers for the meeting. Next meeting organized by Alberto.
c) Topics to be covered during the meetings: Tests, Tutorials, libraries with pip, formatting.

2. Tests.

We currently have a series of Python tests that can be run locally.
- Task: Further development would include code coverage and format checking done automatically on a testing machine, maybe github: Kishore and Simon.

3. Pip: issues with library installation in Ubuntu.

- Task: Add to the documentation how to create your own Python environment. (Simon)
- Task: Upload PC Python as separate pip package. (Alberto)

4. Documentation.

- Task: Update the tutorials. (Simon, Isak)
- Task: Move the tutorials into Read The Docs. (Illa)

5. Real life simulation purely in Jupyter with life visualization.

- Task: Add it to read the docs. (Isak)
- Task: Can we use simulation objects? (Simon)

6. Tutorials.

We currently have information dispersed.
- Task: Coding style should be put last. (Illa)
- Task: Link tutorials in the Pencil Code from Read The Docs. (Simon)

7. Inconsistency in array format when reading with io_dist or hdf5.
- Task: Make sure that the h5 file load follows same. (Kishore)






2021 Sept 17
------------

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


2021 Apr 14
------------

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


2021 Feb 24
------------

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





2021 Jan 14
------------

**Participants:** 

**Minutes:**


Topics:

* Replace package names to non-conflicting  ones?

* Documentation: use sphinx

* Imports: import in the functions or in the header

* Cyclic dependence structure in different imports.

* Newsletter



2020 Nov 23
------------

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

