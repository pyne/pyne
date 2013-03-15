FluDAG: DAGMC Interface for Fluka
==========================================================

The DAGMC interface for Fluka will be either a C++ main that calls
flukam_ with a geant4 flag, or will make use of the main routine
internal to flukam_ and the linking capabilities of the lfluka script
provided in the  FLUKA installation's flutil subdirectory.

The call to flukam_ is usually #define'd to flukam (no trailing underscore)
and is instantiated through the fluka library.

C++ main
--------
For some testing purposes it is convenient to use the C++ main.  In 
this case mainFluDag.cpp is the C++ file that contains a main() method
and performs DagMC initialization prior to calling flukam_.  If flukam_ 
is called with its integer argument flag set to the value of 3 it will
bypass some of its internal calls and instead make calls into
linked in object code.  The object code can be compiled from C++ source
(or Fortran).

Advantages to this method are:
	* A cmake build system is easily implemented
        * The fluka library can be linked using g++, with
          no overt reference to fortran  
        * The source code is all in one language, and the passing of
          arguments is more straightforward and flexible

Disadvantages
	* It is not possible to use FLUKA's method of linking in user
          routines that can pass information to flukam_
	* Workflow to prepare the correctly defined geometry cards for
          FLUKA is multi-step


FLUKA main
-----------
The FLUKA library has a main() method defined.  FLUKA provides a script to 
link the fluka library as a main, along with predefined FLUKA user object 
codes that can be compiled from user-written Fortran source.  The link command
looks like this:

${FLUPRO}/flutil/lfluka -o $@ -a MOAB -a dagmc -m fluka obj/userini.o  -L ${FLUPRO} -L ${MOAB_LIBDIR} obj/fluka_funcs.o obj/WrapInit.o obj/createRegionAndMatFiles.o


This is how FLUKA is intended to be linked in with other codes, namely GEANT4.  
However testing of this methodology showed that the input geometry file was not
read in early enough for material assignments to be set up programmatically
in a user-defined function.

Another issue with the lfluka link method is that creating a cmake system for it
will require expertise, assuming it's possible.

Planning
--------
This work is reviewed and updated weekly or biweekly.  The following sections 
provide an overview of how the code project is to be planned, tracked, tested,
and documented.

Code Project
~~~~~~~~~~~~
* Document User feature table
* Document Workflow

* Document Developer feature table
 
* Building
  - Currently two build systems are in use:  We started with GNU make, but are moving to cmake where possible.
  - A README file in the source directory contains build notes and instructions.

* Testing
  - The choice of cmake as a build system also permits implementation of ctest.

Source Control
~~~~~~~~~~~~~~
The FluDAG project is under Git source control, under svalinn/DAGMC as a publicly
visible repository.  It has been cloned and placed uder the julry repository, as 
julry/DAGMC/FluDAG, for local collaboration.

Git
___
- FluDAG Source Code
- Doxygen project file Doxyfile
- CMakeLists.txt 
- Makefile
- README files
- Sample input files
- Sample geometry files
Collaboration and Visibility
____________________________
- Publicly visible svalinn repository, pushed weekly
- Local collaboration branch in julry, pushed as often as needed

The code is self-documenting, with Doxygen providing the ability to collect and view
the hyperlinked code documentation in a browser window in file list or graphic form.

Besides the files above, project planning and status documentation will be placed in 
the CNERG pages, which is part of the same git repository.

Issue Tracking
~~~~~~~~~~~~~~
Github, Huboard
Closed issues will be reviewed for inclusion in the documentation.  They may be 
documented within code files, in separate documentation files, or be in Sphinx-compiled
doc-src files, depending on the issue.

Documentation
~~~~~~~~~~~~~
Code
____
Code comments will be formatted and tagged to be used by Doxygen in order to 
document functions and classes as we go.
We will use doxywizerd to run Doxygen 1.7.1 on the source file tree.
Class relationship diagrams and tables can be extracted via doxygen analysis.
The file named "Doxyfile" has been placed in the git repository.   This file  
contains the setting used by doxywizard to run doxygen and will have to be
updated by individuals cloning or forking this repository.

The html directory containing the graphical lists, trees, and diagrams produced 
by doxygen will not be placed under git, since it is reproducable from the Doxyfile.

Project
_______
The User Features and Workflow should be documented via Sphinx in the doc-src directory.
The Developer Feature status, and its map to the User features should remain in the 
git repo.

