=======================================
Hexahedral Mesh Sampling
=======================================
The Hexahedral (hex) Mesh Sampling project will provide a tool for generating source distributions in cartesian coordiantes to first order accuracy given values on hexahedral meshes. This provides source definition flexibility to PyNE users.

Brief explanation:
------------------
Plasma power density is defined in a "plasma" coordinate system that is related to real space through Fourier coefficients. This defines the source on a set of known flux surfaces. We would like to be able to take this source data and transform it to cartesian space such that it can be used by solvers more easily. The transformation involves some coordinate mapping and numerical integration. The result can be formatted as a cumulative distribution function for Monte Carlo sampling, or left as a meshed distribution of values. 

This general idea was originally implemented in matlab. This project is to convert the original code to Python, improve the implementation, and perform some experiments to ensure it functions correctly. A stretch goal and/or continuing work would be to implement adaptive mesh selection. 

Expected results:
------------------
The ability to use the hex mesh interface through PyNE and creation of the associated documentation.

Knowledge Prerequisites:
------------------------

Necessary:

*  Intermediate level familiarity with Python.
*  Beginner knowledge of version control.

Helpful:

*  Familiarity with numerical intergration techniques.

Mentor:
-------
Professor Rachel Slaybaugh and Postdoctoral Scholar Kathryn Huff, University of California, Berkeley.

More Resources:
---------------
*  Source code : github.com/pyne/pyne
*  Example gallery : http://pyne.io/gallery/index.html
*  R.N. Slaybaugh, P.P.H. Wilson, L.A. El-Guebaly, E.P. Marriott. “Three-Dimensional Neutron Source Models for Toroidal Fusion Energy Systems.” *Fusion Engineering and Design* **84** (2009) 1774-1778.

