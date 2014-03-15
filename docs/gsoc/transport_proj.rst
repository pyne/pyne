===================================
Transport Solver Components
===================================
The Transport Solver Components project will create the structure and first pieces to allow PyNE to have plug-and-play deterministic transport solver functionality. 

Brief explanation:
------------------
Any deterministic transport code needs components to deal with the spatial variable, the angular variable, and the energy variable. Each of these are typically solved with different numerical methods. This project is to implement and test at least one of these components. The goal is to build up a suite of solver options from which a user can select to build their code. 

Implementing a solver in PyNE will involve method selection, code implementation, testing, and documentation. Determining a structure for solver addition will pave the way for additional solvers to be added in the future. 

Expected results:
------------------
The ability to access one solver component for a determinisitic transport method through PyNE and creation of the associated documentation.

Knowledge Prerequisites:
------------------------

Necessary:

*  Intermediate level familiarity with at least one of Python, C++, or fortran. 
*  Beginner knowledge of version control.

Helpful:

*  Familiarity with numerical solvers.

Mentor:
-------
Professor Rachel Slaybaugh and Postdoctoral Scholar Kathryn Huff, University of California, Berkeley.

More Resources:
---------------
*  Source code : github.com/pyne/pyne
*  Example gallery : http://pyne.io/gallery/index.html

