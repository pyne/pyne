.. _devsguide_simplesim:

================================================
Simple Simulation Input -- :mod:`pyne.simplesim`
================================================

.. currentmodule:: pyne.simplesim

.. automodule:: pyne.simplesim

This guide is intended to be the friend of anyone who wants to develop
:py:mod:`pyne.simplesim` further, or is disgusted by its currently limited
support. The guide is not comprehensive, and any prospective developer should
feel comfortable contacting Chris Dembia with any questions.

************
How it works
************
Referencing a card
    Nowhere does a user ever need to specify a card number, which at sight has
    little meaning, provides no insight to what it represents, and can make
    modification a real pain. So, as is done in LaTeX, cards are given names
    (labels), and when that card number is needed on another card, the user
    references the other card by its name (or with the card object itself).
    Now, the cards themselves live in the definition. The definition stores
    cards in ordered dictionaries, and right now the card numbers are the
    card's index in the dictionary (+ 1 because of zero-indexing). That means
    that at the time a card is written, it needs access to access the
    simulation, which knows the card numbers. Therefore, the
    :py:meth:`cards.ICard.mcnp` method takes in the simulation.
    These methods are called by the input file object, which has the
    simulation.

*************
Adding a card
*************

When adding a card to :py:mod:`pyne.simplesim.cards`, the following tasks need
to be completed:

* Choose the appropriate base class (most likely
  :py:class:`cards.IMisc`), and subclass from this.
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
* Each card has a :py:meth:`comment` and :py:meth:`mcnp` method. The comment
  should be defined for each card, but the mcnp card only needs to be defined
  for cards that MCNP supports. To support a different code, i.e. Serpent,
  define a :py:meth:`serpent` method.
* Place the card under the appropriate category in the index in the module's
  docstring.
* Add tests in `pyne/tests/test_simplesim.py`.

Adding a new abstract base class
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
If introducing a new abstract base class or card category (i.e.
:py:class:`cards.IMisc`): 

* Make sure to prepend the name of the class with `I`, a notation which stands
  for interface in Java.
* Following the model set by other abstract base class's docstrings, make it
  clear to the user that they do not interact with the class on their own.
* The constructor takes ``*args`` and ``**kwargs``, and passes ``*args`` and
  ``**kwargs`` to the next constructor in the MRO (method resolution order)
  through the call to :py:meth:`super`. One reason for this is that
  :py:class:`cards.ICard` takes a ``unique`` keyword and it is
  not enjoyable to carry this keyword through all subclass constructors.


****************
Design decisions
****************

MCNP
   The most important thing to the initial developer was to have something that
   worked for MCNP, given his limited development time. Therefore, the content
   of the cards is very similar to that which is found in MCNP. Ideally, the
   cards would be created independently from a specific code. However, doing so
   would likely lead to making sacrifices for all codes, including MCNP. The
   way the code is written now, it surely works well for MCNP, though
   sacrifices may need to be made for other codes.

   The following paragraph was written before development of the module:
   `Perhaps it is ideal to create some very general definition that could be
   used for any code, but I think that is a little ambitious and idealistic,
   and will ultimately result in sloppy attempts to match the general
   definition to specific codes or will hamper the functionality of the general
   definition. Thus, I think MCNP should be the "first child"/"first class
   citizen" (incorrect usage) of this package.  Classes can be written in`
   :py:mod:`pyne.simplesim.inputfile` `allowing the system and option
   definitions to be used to generate input files for other codes, but will
   require special consideration. Two cases to consider are that densities and
   reflecting boundary conditions are specified in different places between
   MCNP and Serpent. The code will be written so that these settings are placed
   in the proper place for the MCNP input, and the Serpent inputfile class will
   have to pick that information from its place in the MCNP cards and place it
   in the right place for a Serpent input. The Serpent inputfile class will
   likely need to contain many exceptions for inputs to MCNP that Serpent
   cannot handle (e.g.  importances on cell cards). MCNP is a widely used code,
   and a product taht can cleanly generate MCNP input but that cannot generate
   input to other codes as well is more valuable, I think, than a product that
   falls short for all codes (including MCNP).`


:py:meth:`cards.ICard.mcnp`
   It was initially advised that the MCNP card functionality would be separated
   from the card classes. It seemed preferable to place all code-specific
   functionality in one class (i.e.
   :py:class:`inputfile.MCNPInput` for MCNP). Given time
   constraints, it was not possible to do this.  Furthermore, it seems that the
   structure as it is now, where each card has a :py:meth:`mcnp` method, makes
   a lot of sense: it allows us to take advantange of inheritance and similar
   structure between cards in the same category.

Separation of :py:class:`definition.SimulationDefinition` and :py:class:`inputfile.IInputFile`
   At first glance it may not make sense to separate the defintion from the
   input file object.
   Since :py:class:`definition.SimulationDefinition` is
   subclassed to create :py:class:`definition.MCNPSimulation`,
   there is MCNP-related information in more than one place, and it may seem
   that all the code-specific information could be contained in one class.
   However, there should be code-independent simulation definition parameters
   that belong in `definition.SimulationDefintion`. It is
   possible that it becomes clear that the input file task can be roped into
   the simulation definition, but the additional modularization should be easy
   to 'undo'. Also, it was thought that the input file classes would eventually
   have the ability to parse in an input file, and so it might make sense to
   have a separate class.

Referencing other cards using strings or the objects themselves
    Many cards need information about other cards. For example, tallies need to
    know the card number of the surfaces or cells being tallied. There are at
    least two ways to provide this information to the tally card: give the
    actual card object, or just give the string name of the card to be tallied.
    In both cases, the card looks up the surface or cell number in the system
    defintion, using the ``sim`` input of the :py:meth:`mcnp` method. Right
    now, the cell card takes in the actual cards, not just card names, but all
    simulation-related cards take in card names, not hte actual cards. The
    reason for using card names only is that it makes it easier for the user to
    modify an input they've already created, after the user no longer has the
    original cards they created. Rather than doing something like
    ``Tally(sys.cells['cA'])`` to create a tally, the user can just do
    ``Tally('cA')``.
    
    This further has the advantage that the user can define
    the tally before the cell 'cA' is even created, putting less restrictions
    on how the user defines the system. The disadvantage here is that passing
    the actual object allows for type checking, and puts less restriction on
    the names used for cards. For example, some cards accept either a cell or
    universe. If providing a name only, then the code does not know whether the
    string is the name of a cell or universe. If providing the objects
    themselves, the object contains the knowledge of its type.

    This issue has been managed for tallies by creating the
    :py:mod:`pyne.simplesim.nestedgeom` module. In this module, there are
    classes such as :py:class:`nestedgeom.Cell`,
    :py:class:`nestedgeom.Surf`,
    :py:class:`nestedgeom.Univ` that only hold a string name of
    a cell/surface/universe (they are not cells/surfaces/universes themselves),
    but make it clear thta the string name represents a cell/surface/universe.
    The problem here is that the notation becomes quite verbose.

Float formatting
    The author did not want to presume how detailed the user wanted to be about
    the floats that appear in the input file. Therefore, the user can set the
    ``float_format`` in the input file. The syntax used for specifying the
    float format is what is typically found in C++ and MATLAB (i.e. '%.5g').
    Python now expects users to use their own syntax (i.e. '{:.5g}'. Here, the
    former has been used to make it easier for people who are new to Python to
    use this package. This does cause code to be more verbose, as for
    everything except these floats, string formatting is done with
    :py:meth:`str.format`.  It is possible that there are different categories
    of floats, and that the user would want each of these types of floats to be
    formatted differently.  This is not supported currently.

Default arguments
    Since codes like MCNP have default arguments themselves when arguments are
    left out, those arguments are left out the cards in this package. That is,
    the default values that a particular code uses are not stored in the
    package, as the defaults may change or may be different across different
    codes. Instead, default arguments are set to None, which usually means to
    omit that input from the input file. This may not be true for all cards.


***********
Limitations
***********

Saving and loading definitions
    This is the biggest limitation, that definitions and cards cannot be saved
    in a persistent manner. Work has been done on this front, however. A
    subclass of the JSONEncoder can be found in
    :py:mod:`pyne.simplesim.definition`. A crude way to allow for this saving
    and loading would be to define a ``__repr__`` method for all cards.
    
Surfaces
    Only a few surfaces are implemented right now; there needs to be many more.

MCNP card order
   In MCNP, data cards that refer to cell cards by the order of cell cards
   (i.e.  Volume) actually depend on the order in which the cell cards are
   listed, and not at all on the cell card numbers. In the current
   implementation, cells (well, all cards) are always printed in order of
   increasing card number, so this is not an issue. However, if the user
   supplies custom cards with card numbers out of order, then this would
   present an issue.

Card numbers
   It is assumed for the most part in the definitions that cards are numbered
   as they are in MCNP. This is fine for the simulations, since those are
   expected to be subclassed for a specific code. However, the system supposed
   to be completely code-independent.

MCNP Macrobody facets
   In MCNP, the user is able to specify a facet of a macrobody surface on a
   cell card (i.e. 3.2 is facet 2 of surface 3). Preliminary work has been done
   to implement this, but this is not functional.


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
- Universes, lattices, and fill are done either using IMisc cards or on the
  :py:class:`cards.Cell` card. However, it is preferable that
  there are actual Universe and Lattice cards/objects (this is not the same
  thing as the current Universes and Lattice cards, that are just an
  implementation of MCNP's U and LAT cards) that can be used to fill cells
  instead of regions.
- Geometry transformations. Right now, a crude form of geometry transformations
  are allowed; the user can shift and stretch surfaces or regions, so long as
  the surface supports it. Rotation is not yet supported. The surfaces are
  still tied to how they are defined in MCNP, and so for most surfaces rotation
  is not supported.
- Serpent provides a few macros, like the pin macro, that makes the
  specification of surfaces and cells a lot easier. These can be implemented
  here using static Region methods that create a region, and the required
  surfaces.
- It would be nice, for MCNP materials, if this package could automatically
  pick the correct nuclide library table identifiers for a given temperature
  (e.g. .71c for 600 K), though this would change over time. A manual override
  should also be allowed.
- MCNP cards are mostly printed in a sloppy way (except material cards); it
  would be nice to exploit the automated card string generation in order to
  have the strings print in a pretty way.
- **TODO** nodes for small items can be found scattered among the code.
- There is a `github wiki <https://github.com/fitze/pyne/wiki/simplesim>`_ page
  that has a list of some ideas for future development.


.. moduleauthor:: Chris Dembia <cld72@cornell.edu>




