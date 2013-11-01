.. _pyne_simplesim:

================================================
Simple Simulation Input -- :mod:`pyne.simplesim`
================================================

.. currentmodule:: pyne.simplesim

.. automodule:: pyne.simplesim
  
The ``simplesim`` package provides modules for the definition of a system (its
geometry and material composition), the definition of simulation parameters,
and for the generation of input files to codes that can perform such
simulations (MCNP, Serpent, MCODE, etc.). The package is imported as such::

    from pyne import simplesim
    
The package is broken up into three modules: ``cards``, ``definition``, and
``inputfile``:

* :py:mod:`cards`: A ``card`` can be considered the basic unit of input; each
  piece of input is represented by a card. The word ``card`` is taken from its
  usage with MCNP.  The ``card`` classes are defined using the same
  fields/properties used in an MCNPX input, though the definition of a ``card``
  should not require knowledge about the format of MCNP cards (but perhaps
  their content). Each card has a method to generate a string formatted for
  MCNP input. The :py:mod:`nestedgeom` module exists to help with cell card
  input.
* :py:mod:`definition`: The ``definition`` module allows a user to create a
  system (reactor) definition through the combination of surface, region,
  material, and cell cards, and simulation definitions using a slew of other
  cards. Essentially, an object of a class in the ``definition`` module stores
  all the ``cards`` that are needed to define a system and simulation.
* :py:class:`inputfile`: Generates an input file for a neutronics code given
  a system definition and a options definition, both of which are classes in
  ``definition``.

The reference for each of the modules described above can be found using the
links below. Additionally, there is a page that shows an inheritance diagram
for all the classes in each module.

.. toctree::
    :maxdepth: 1

    cards 
    nestedgeom
    definition
    inputfile
    inheritance

For an explanation of how to use the module, with full examples, visit
:ref:`usersguide_simplesim`. If you are interested in developing this package
further or are looking for information about the state of development, visit
:ref:`devsguide_simplesim` and/or contact the developers.

The original author of this package is Chris Dembia, and it was developed in
consultation with Anthony Scopatz and Paul Wilson. See :ref:`dev_team` for
contact information.

.. moduleauthor:: Chris Dembia <cld72@cornell.edu>


