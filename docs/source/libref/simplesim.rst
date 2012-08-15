.. _pyne_simplesim:

===================================
Simple Simulation Input Definitions
===================================

.. currentmodule:: pyne.simplesim

.. automodule:: pyne.simplesim
  
The ``simplesim`` package provides modules for the definition of a system (its
geometry and material composition), the definition of simulation parameters,
and for the generation of input files to codes that can perform such
simulations (MCNP, Serpent, MCODE, etc.). The package is imported as such::

    from pyne import simplesim

The three objectives of this module, in increasing specificity, are:

1- To provide an easy to use object-oriented interface to commonly used codes
   in nuclear science and engineering.
2- To abstract the definition of a system and simulation parameters away from
   the input syntax of a specific code. It is truly valuable to be able to
   define a reactor in a way that can be used to generate input for two
   different codes. However, the definition is not necessarily generalizable
   across different codes.
3- To formulate the input to such codes in a modular way that is persisent and
   amenable to modification (e.g. for parameter space studies).

To achieve the second objective, the package has been broken up into three
modules: ``cards``, ``definition``, and ``inputfile``:

* ``cards``: A `card` can be considered the basic unit of input; each piece of
  input is represented by a card. The word `card` is taken from its usage with
  MCNP.  The ``card`` classes are defined using the same fields/properties used
  in an MCNPX input, though the ``card`` classes contain no information about the
  format of MCNP cards. The format of such cards is managed by the
  ``inputfile`` module described below.
* ``definition``: The ``definition`` module allows a user to create a system
  (reactor) definition through . Essentially, an object of a class in the
  ``definition`` module stores all the cards, from ``cards`` that are needed
  to define a system and simulation.
* ``inputfile``: Generates an input file for any neutronics code given a system
  definition and a options definition, both of which are classes in
  ``definition``.

The author believes that an additional module, ``run``, may provide the
much-desired functionality of executing the correct program (MCNPX, etc.) and
connecting the run to the appropriate input parsers provided in other modules
of `PyNE`. Furthermore, such a module could contain a method to perform a
parameter study simply with one method call (something like ``keff_vector =
run.parameter_study(param_name, param_values``). This module has not been
implemented.

The most important thing I have to say about this plan for the module is that
all the modules will be written for the abilities of MCNPX, though the exact
format of the MCNPX input will be properly modularized. Perhaps it is ideal to
create some very general definition that could be used for any code, but I
think that is a little ambitious and idealistic, and will ultimately result in
sloppy attempts to match the general definition to specific codes or will
hamper the functionality of the general definition. Thus, I think MCNPX should
be the "first child"/"first class citizen" (incorrect usage) of this package.
Classes can be written in ``inputfile`` allowing the system and option
definitions to be used to generate input files for other codes, but will
require special consideration. Two cases to consider are that densities and
reflecting boundary conditions are specified in different places between MCNPX
and Serpent. The code will be written so that these settings are placed in the
proper place for the MCNPX input, and the Serpent inputfile class will have to
pick that information from its place in the MCNPX cards and place it in the
right place for a Serpent input. The Serpent inputfile class will likely need
to contain many exceptions for inputs to MCNPX that Serpent cannot handle (e.g.
importances on cell cards). MCNPX is a widely used code, and a product taht can
cleanly generate MCNPX input but that cannot generate input to other codes as
well is more valuable, I think, than a product that falls short for all codes
(including MCNPX).

Though the ``cards`` module contains the most basic information needed to use
this class, its reference below is given after the references for
``definition`` and ``input``, since this module is by far the largest.

The original author of this package is Chris Dembia, and it was developed in
consultation with Anthony Scopatz and Paul Wilson. See :ref:`dev_team` for
contact information.

*****************
Definition Module
*****************

.. automodule:: pyne.simplesim.definition
    :member-order: 'bysource'
    :members:
    :inherited-members:

************
Input Module
************

.. automodule:: pyne.simplesim.inputfile
    :member-order: 'bysource'
    :members:
    :inherited-members:

************
Cards Module
************

.. automodule:: pyne.simplesim.cards
    :member-order: 'groupwise'
    :members:
    :inherited-members:
    :show-inheritance:

