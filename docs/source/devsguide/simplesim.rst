.. _devsguide_simplesim:

============================================================
Simple Simulation Input Definitions -- :mod:`pyne.simplesim`
============================================================

.. currentmodule:: pyne.simplesim


************
How it works
************
number referencing
    Similar to LaTeX TODO

the mcnp(float_format, sim) thing
    Passing sim and float_format

*************
Adding a card
*************

When adding a card to :py:mod:`pyne.simplesim.cards`, the following tasks need
to be completed:

* Choose the appropriate base class (most likely :py:class:`IMisc`), and
  subclass from this.
* If the card is `unique`, such as the Criticality card, then make sure to pass
  the unique name to the base class constructor, and to pass the
  ``unique=True`` keyword argument.
* Make sure that the base class constructor is called using :py:meth:`super`.
  If this is not done, future functionality may break the functionality you
  contribute.
* Be sure to give complete docstrings for at least the class and ``__init__``.
  Be sure to give at least one example for each method the user is expected to
  use.
* Place a ``.. inheritance-diagram::`` directive in the class docstring. See
  other classes for a guide.
* Place the card under the appropriate category in the index in the module's
  docstring.

Adding a new abstract base class
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
If introducing a new abstract base class or card category (i.e.
:py:class:`IMisc`): 

* Make sure to prepend the name of the class with `I`, a notation which stands
  for interface in Java.
* Following the model set by other abstract base class's docstrings, make it
  clear to the user that they do not interact with the class on their own.
* The constructor takes ``*args`` and ``**kwargs``, and passes ``*args`` and
  ``**kwargs`` to the next constructor in the MRO (method resolution order)
  through the call to :py:meth:`super`. One reason for this is that
  :py:class:`ICard` takes a ``unique`` keyword and it is not enjoyable to carry
  this keyword through all subclass constructors.


****************
Design decisions
****************

Abstraction
^^^^^^^^^^^
MCNP
   The most important thing to the initial developer was to have something that
   worked for MCNP, given his limited development time. Therefore, the content
   of the cards is very similar to that which is found in MCNP. Ideally, the
   cards would be created independently from a specific code. However, doing so
   would likely lead to making sacrifices for all codes, including MCNP. The
   way the code is written now, it surely works well for MCNP, though
   sacrifices may need to be made for other codes.
:py:meth:`ICard.mcnp`
   It was initially advised that the MCNP card functionality would be separated
   from the card classes. It seemed preferable to place all code-specific
   functionality in one class (i.e. :py:class:`MCNPInput` for MCNP). Given time
   constraints, it was not possible to do this.  Furthermore, it seems that the
   structure as it is now, where each card has a :py:meth:`mcnp` method, makes
   a lot of sense: it allows us to take advantange of inheritance and similar
   structure between cards in the same category.

Separation of :py:class:`SimulationDefinition` and :py:class:`InputFile`
   Since :py:class:`SimulationDefinition` is subclassed to create
   :py:class:`MCNPDefinition`, it seems to make senes that    TODO


Referencing other cards using strings or the objects themselves
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

float_format % uses dollar rather than str.format() because the % syntax is
more conventional outside of python.



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


The three objectives of this module, in increasing specificity, are:
TODO similar to latex, don't require the user to track numbers.
1- To provide an easy to use object-oriented interface to commonly used codes
   in nuclear science and engineering.
2- To abstract the definition of a system and simulation parameters away from
   the input syntax of a specific code. It is truly valuable to be able to
   define a reactor in a way that can be used to generate input for two or more
   different codes. However, the definition is not necessarily generalizable
   across different codes.
3- To formulate the input to such codes in a modular way that is persisent and
   amenable to modification (e.g. for parameter space studies).

To achieve the second objective, the package has been broken up into three
modules: ``cards``, ``definition``, and ``inputfile``:

* :py:mod:`cards`: A ``card`` can be considered the basic unit of input; each
  piece of input is represented by a card. The word ``card`` is taken from its
  usage with MCNP.  The ``card`` classes are defined using the same
  fields/properties used in an MCNPX input, though the ``card`` classes contain
  no information about the format of MCNP cards. The format of such cards is
  managed by the ``inputfile`` module described below. The :py:mod:`nestedgeom`
  module exists to help with cell card input.
* :py:mod:`definition`: The ``definition`` module allows a user to create a
  system (reactor) definition through . Essentially, an object of a class in
  the ``definition`` module stores all the cards, from ``cards`` that are
  needed to define a system and simulation.
* :py:class:`inputfile`: Generates an input file for any neutronics code given
  a system definition and a options definition, both of which are classes in
  ``definition``.


***********
Limitations
***********

Saving and loading definitions
    This is the biggest limitation, that definitions and cards cannot be saved
    in a persistent manner. Work has been done on this front, however. TODO
    JSON
    




*******************
Further Development
*******************
- Initial discussions about this package set forth a much larger task: creating
  an abstract syntax that  is truly abstracted from a specific code (i.e.
  MCNP). This is certainly possible for the cards associated with the system
  definition, but doing so with the simulation definition seems very difficult
  and likely not possible.
- An additional module, ``run``, may provide the
  much-desired functionality of executing the correct program (MCNP, etc.) and
  connecting the run to the appropriate input parsers provided in other modules
  of `PyNE`. Furthermore, such a module could contain a method to perform a
  parameter study simply with one method call (something like ``keff_vector =
  run.parameter_study(param_name, param_values``).
- **TODO** nodes for small items can be found scattered among the code.







