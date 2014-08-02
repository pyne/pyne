===================
Integrate AHOTN
===================
The AHOTN project will allow PyNE users to access a spatial solver that is used in the computational transport field.

Brief explanation:
------------------
The Arbitrarily High Order Method of the Nodal type (AHOTN) is a nodal method used for spatial discretization of the SN deterministic transport equation. We have an existing, standalone implementation of AHOTN that is written in fortran. This project would incorporate that implemation such that it would be accessible from PyNE.

Bringing the AHOTN code into PyNE will likely involve segmenting it into reasonable modules from which we can build an application (the stand-alone CLI), possibly mixing in some existing PyNE notions. Determining how to do this will set an important precedent for incorporating other codes into PyNE. 

Expected results:
------------------
The ability to use the AHOTN method through PyNE and creation of the associated documentation.

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
*  The associated PyNE ticket: https://github.com/PyNE/PyNE/issues/219
*  Source code : github.com/pyne/pyne
*  Example gallery : http://pyne.io/gallery/index.html

