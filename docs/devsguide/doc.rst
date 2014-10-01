.. _devsguide_doc:

====================
How to Document PyNE
====================
Documentation takes many forms. This will guide you through the steps of 
successful documentation.

Docstrings
----------
No matter what language you are writing in, you should always have documentation
strings along with you code. This is so important that it is part of the 
style guide.  When writing in Python or Cython, your docstrings should be in 
reStructured Text using the numpydoc format. When writing in C, C++, or 
Fortran you should write your docstrings in Doxygen.

Auto-Documentation Hooks
------------------------
The docstrings that you have written will automatically be connected to the 
website, once the appropriate hooks have been setup.  At this stage, all 
documentation lives within pyne's top-level ``docs`` directory. 
PyNE uses the sphinx tool to manage and generate the documentation, which 
you can learn about from `the sphinx website <http://sphinx-doc.org/>`_.
If you want to generate the documentaion, first pyne itself must be installed 
and then you may run the following command from the ``docs`` dir:

.. code-block:: bash

    ~/pyne/docs $ make html

For each new 
module, you will have to supply the appropriate hooks. This should be done the 
first time that the module appears in a pull request.  From here, call the 
new module ``mymod``.  The following explains how to add hooks based on language:

Python & Cython Hooks
......................
Python and Cython documentation lives in the ``docs/pyapi`` directory.  
First create a file in this directory that represents the new module called
``mymod.rst``.  
The ``docs/pyapi`` directory matches the structure of the ``pyne`` directory.
So if your module is in a sub-package, you'll need to go into the sub-package's 
directory before creating ``mymod.rst``.
The contents of this file should be as follows:

**mymod.rst:**

.. code-block:: rst

    .. _pyne_mymod:

    ======================================
    My Awesome Module -- :mod:`pyne.mymod`
    ======================================

    .. currentmodule:: pyne.mymod

    .. automodule:: pyne.mymod
        :members:

This will discover all of the docstrings in ``mymod`` and create the appropriate 
webpage. Now, you need to hook this page up to the rest of the website.

Go into the ``index.rst`` file in ``docs/pyne`` or other subdirectory and add 
``mymod`` to the appropriate ``toctree`` (which stands for table-of-contents tree).
Note that every sub-package has its own ``index.rst`` file.

C++ et al. Hooks
................
This is very similar to the Python hooks except that all of the pages live in the 
``docs/cppapi`` directory.  Again, create a ``mymod.rst`` file in the appropriate 
place. This file will have contents as follows:

**mymod.rst:**

.. code-block:: rst

    My Awesome Module
    =====================================

    .. autodoxygenindex:: mymod.h
        :source: pyne_mymod

This must be performed for every header file.  Again, open up the ``index.rst`` file
edit the ``toctree`` to indclue ``mymod`` where appropriate.

User's Guide
----------------------
The user's guide is available for additions via the ``docs/usersguide`` directory.
This is a more free-form and high-level introduction to pyne topics than anywhere
else. Simply write rst files and add them to the ``index.rst``.

Gallery
-------
The gallery is a collection of IPython notebooks which provide simple examples of 
how to perform tasks with pyne.

Theory Manual
-------------
The theory manual is a required document for our quality assurance. These 
documents have a well defined structure that you may see in the 
``docs/theorymanual/theory_template.rst.t`` file.
