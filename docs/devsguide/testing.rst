.. _devsguide_unittest:

================
How to Test PyNE
================
This guide will teach you the basics of how to test PyNE code.

------------
Unit Testing
------------

First, install nose:
http://nose.readthedocs.org/en/latest/

To perform all unit tests::

    $ cd tests/
    $ nosetests

This will recursively look through the currently directory, open up every file
named test_* and run every function (or method) named test_*.

Nosetests can also take file(s) as an argument. For example, to just run the
mcnp and material module tests:

    nosetests test_mcnp.py test_material.py

A clean build/nucdatamake should yield a version of PyNE with all tests
passing.
 
---------------
Example Testing
---------------
The examples directory should also be kept up-to-date as much as possible.
PyNE examples are either in Python files or IPython notebooks. This means that
to test the examples requires a recent version of IPython.  Furthermore, the 
examples themseleves may have many other optional dependencies.  Don't be alarmed
if testing the examples fails do to a lack of having a dependency installed.
For this reason, testing the examples is not as important as unit tests, but still
should be done occassionally.

To run the examples automatically, go to the examples directory and run the 
``execer.py`` file from the root pyne dir.

.. code-block:: bash

    $ cd examples
    $ ../execer.py

----------------
Tutorial Testing
----------------
Tutorial testing is very similar to example testing except that all of the 
tutorials are IPython notebooks.

To run the tutorials automatically, go to the tutorial directory and run the 
``execer.py`` file from the root pyne dir.

.. code-block:: bash

    $ cd tutorial
    $ ../execer.py

-----------------------
Putting It All Together
-----------------------
If you'd like to run all of the tests automatically from the root pyne dir, 
you can chain the following BASH commands together::

    $ cd tests && nosetests && cd ../examples && ../execer.py || cd ../tutorial && \
      ../execer.py || cd ..

Happy testing!

