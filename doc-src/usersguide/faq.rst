Frequently Asked Questions
===============================

* Q1_ How does the implicit complement work?
* Q2_ Can I automate the CUBIT conversion process?

.. _Q1:

Q1: How does the implicit complement work?
----------------------------------------

**A1** Since your geometry has been imprinted and merged, all surfaces fall into one of two categories:

1. surfaces that are shared by two volumes
2. surfaces that are used in only one volume

All surfaces in this second group are, by definition, on the boundary
of the complement region since they have an explicitly defined region
on one side and nothing on the other.  The implicit complement is
formed automatically from this list of surfaces.

From the point of view of MCNP, and extra cell is created to represent
the implicit complement.  Your output file will include one entry for
each volume in your geometry, including the graveyard, plus one extra
entry representing the complement.

.. _Q2:

Q2: The CUBIT conversion process is tedious and frustrating. Is there a way to avoid or automate most of this work?
----------------------------------------

**A2** Yes, a script has been written that will automate all the
necessary steps to perform the CUBIT conversion process with limited
user input. Information about this script can be found on the
AutomatedCubitConversion website.

