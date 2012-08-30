.. _usersguide_simplesim:

============================================================
Simple Simulation Input Definitions -- :mod:`pyne.simplesim`
============================================================

not intended to replace an understanding of the manual for any code you want to
run a simulation in.

Examples that are too complicated are useless, but all features still need to
be shown. It's better to have many small examples than to have a few big ones.

========
Overview
========

.. image:: simplesimdiagram.svg


=====
Usage
=====

Card storage :
    Cards are all stored in ordered dictionaries. To see the names and cards
    that have been added to a definition, one can examine these dictionaries.
    For example, the user can get the names of all cells via the following::

        print sys.cells.keys()

    The dictionaries are read-only, though the cards within the dictionary are
    not. That is, a user cannot remove an element of the dictionary, but can
    modify any card contained within the dictionary. The dictionaries
    themselves can only be modified via add/remove methods in the definition.


========
Examples
========


******
Godiva
******
Features: surfaces, cells, ...

***********************
Simple infinite lattice
***********************

--DigitalWorkshop


******
Burnup
******

Explain a possible use scenario for manipulating a reactor: have a script that
creates the reactor, and have scripts that analyze results, etc...
