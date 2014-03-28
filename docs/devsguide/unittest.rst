.. _devsguide_unittest:

================
How to Unit Test  
================

First, install nose:
http://nose.readthedocs.org/en/latest/

To perform all tests:

    cd tests/
    nosetests

This will recursively look through the currently directory, open up every file
named test_* and run every function (or method) named test_*.

Nosetests can also take file(s) as an argument. For example, to just run the
mcnp and material module tests:

    nosetests test_mcnp.py test_material.py

A clean build/nucdatamake should yield a version of PyNE with all tests
passing.
 
