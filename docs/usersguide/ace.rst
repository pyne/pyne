.. _usersguide_ace:

==================
ACE Cross Sections
==================

.. currentmodule:: pyne.ace

.. automodule:: pyne.ace

For a complete specification for the classes in the ``ace`` module, please
refer to the Library Reference entry for :ref:`pyne_ace`.

****************************
Example Use of Library Class
****************************

To load data from an ACE file, one needs to simply initialize an instance of
the Library class specifying the path to the Library file.

.. code-block:: ipython

   In [1]: from pyne import ace

   In [2]: libFile = ace.Library('endf70a')

   In [3]: libFile.read()
