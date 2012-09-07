"""The :py:mod:`cards` module can be imported as such::

    from pyne.simplesim import cards

The user creates cards, typically only working with the constructors, and then
adds the cards to a :py:class:`pyne.simplesim.definition.IDefinition`. The user
can modify cards later on by modifying attributes after the card has been
created.

The first section on this page, Index, list all the cards with their category.
The second section is Reference, which contains the reference for all the
classes. In the Reference, the classes are listed in alphabeticla order.

.. moduleauthor:: Chris Dembia <cld72@cornell.edu>


*****
Index
*****

This section presents a listing of the cell cards, with links to their
reference, based on the category to which the cell cards belong. The categories
are important because cards must have a unique name within their category.

An inheritance diagram of all the classes in this module can be found at
:ref:`pyne_simplesim_inheritance`. All cards are subclassed from
:py:class:`ICard`, which is not part of any of the categories below. The only
exception is the :py:class:`Particle` class, which is a helper class used
whenever a particle must be specified as part of a card.

------------
System cards
------------
These are the cards that are provided in a
:py:class:`pyne.simplesim.definition.SystemDefinition`. See further down for
more subcategories. Cell cards are added to a simulation using
:py:meth:`pyne.simplesim.definition.SystemDefinition.add_cell`, and cell card
names must be unique among cell cards.
Material cards are added using
:py:meth:`pyne.simplesim.definition.SystemDefinition.add_material`, and
material card names must be unique among material cards.

.. autosummary::
    :nosignatures:

    Cell
    CellMCNP
    Material


Surface cards
^^^^^^^^^^^^^
Surface cards are added to a system using 
:py:meth:`pyne.simplesim.definition.SystemDefinition.add_surface`.

.. autosummary::
    :nosignatures:

    ISurface
    IAxisSurface
    AxisCylinder
    AxisPlane


Macrobodies
^^^^^^^^^^^
These are surfaces built from simpler surfaces, and accordingly have (will
have) facet functionality. These cards are still added via
:py:meth:`pyne.simplesim.definition.SystemDefinition.add_surface`, and so their
names must be unique among the surfaces.

.. autosummary::
    :nosignatures:

    IMacrobody
    Facet
    Parallelepiped
    Cuboid


Region
^^^^^^
Cells require a region, and regions are constructed using surfaces. These
classes are used to form a binary tree. The & and | operators are overloaded to
facilitate the construction of solid regions from surfaces. These cards can be
named if so desired, but otherwise are not really part of a system like cell,
surface, and material cards are.

.. autosummary::
    :nosignatures:

    IRegion
    IRegionBool
    RegionAnd
    RegionOr
    RegionLeaf


----------------
Simulation cards
----------------
These are the cards that are provided in a
:py:class:`pyne.simplesim.definition.SimulationDefinition`. Some of the cards
can only be used in subclasses of :py:class:`SimulationDefinition`, such as
:py:class:`Transformation`.

The following two cars form their own categories, and are the only card in
their respective category. Distributions are added to simulations using
:py:meth:`pyne.simplesim.definition.SimulationDefinition.add_dist`, and
transformations are added using
:py:meth:`pyne.simplesim.definition.MCNPSimulation.add_transformation` (they
are only available for MCNP simulations). Card names must be unique within
their own category.

.. autosummary::
    :nosignatures:

    Distribution
    Transformation


Source
^^^^^^
These cards are used to define particle sources.
Add these cards to a simulation using
:py:meth:`pyne.simplesim.definition.SimulationDefinition.add_source`. Card
names must be unique among source cards.

.. autosummary::
    :nosignatures:

    ISource
    GeneralSource
    Criticality
    CriticalityPoints


Tally
^^^^^
Particle tallies/detectors.
Add these cards to a simulation using
:py:meth:`pyne.simplesim.definition.SimulationDefinition.add_tally`. Card names
must be unique among tally cards.

.. autosummary::
    :nosignatures:

    ITally
    ICellSurfTally
    SurfaceCurrent
    CellFlux
    CellEnergyDeposition
    CellFissionEnergyDeposition
    CellPulseHeight
    CellChargeDeposition
    IDetector
    PointDetector
    RingDetector


Miscellaneous
^^^^^^^^^^^^^
These cars do not have card numbers and do not fit into any other category.
Add these cards to a simulation using
:py:meth:`pyne.simplesim.definition.SimulationDefinition.add_misc`. Card names
must be unique among misc scards.

.. autosummary::
    :nosignatures:

    IMisc
    ScatteringLaw
    EnergyGrid
    ICellMod
    Universes
    Fill
    Lattice
    Volume
    Area
    FissionTurnoff
    PhotonWeight
    DetectorContribution
    IUniqueParticle
    Temperature
    TemperatureTimes
    ICellModParticle
    Importance
    ExponentialTransform
    ForcedCollision
    WeightWindowBound
    WeightWindowEnergies
    WeightWindowTimes
    DXTRANContribution
    DXTRANSpheres
    Vector


------------
Custom cards
------------
A custom card is used when other cards, within a given category, do not satisfy
the needs of the user. They are added to the system or simulation using the
``add_*`` method for their category, and the names must be unique within their category.

.. autosummary::
    :nosignatures:

    Custom
    CellCustom
    SurfaceCustom
    MaterialCustom
    SourceCustom
    DistributionCustom
    TallyCustom
    MiscCustom
    TransformationCustom


-------------
Unimplemented
-------------
The following cards are defined but are unimplemented.

.. autosummary::
    :nosignatures:

    MaterialMCNP
    Comment
    Mode
    NeutronPhysics
    PhotonPhysics
    ElectronPhysics
    ProtonPhysics
    Burnup


Nested Geometry
^^^^^^^^^^^^^^^
These cards would allow the preferred way of constructing universes and
lattices. Currently, universe and lattice functionality `is` implemented, but
not through these cards.

.. autosummary::
    :nosignatures:

    IUniverse
    UniverseByRegion
    UniverseByLattice
    Lattice
    LatticeByArray


*********
Reference
*********

All classes are listed in alphabetical order. Cards that are available in MCNP
have their MCNP card name in the reference; feel free to search for it to find
the related simplesim card.

"""

# TODO if i make a material card, make sure to revert back to the Material
# object that Anthony has once they implement a comment.
# TODO allow card suppression (don't delete but don't print it).
# TODO allow readfile commands.
# TODO emphasize how easy it is to modify the cards, just subclass it
# individually and use your own.
# Disadvantage of actually providing the cell card is now if you want to modify
# something later, you have to go grab the appropriate cell card, since you
# don't have the object anymore.
# TODO check WWGT default arg effect in WWN.
# TODO improve inheritance around ICellMod; lots of copied code from
# ICellModParticle.
# TODO big opportunity to pretty-up what the mcnp() output looks like.
# TODO consistent plural card names when appropriate.
# Refactor CellMCNP so all relevant classes have a method _mcnp_cell_comment
# and _mcnp_cell_card
# TODO in Cell, require that the material is not pyne.card.Material, but is a
# card here (error check).
# TODO names can be unique to tally type even.
# TODO for read/writing, using repr makes recreating really easy; otherwise
# decoding would take some effort. repr also means that i probably should use
# name referencing everywhere, not object passing. however using json() would
# create something more human-readable.

import abc
import collections
import copy
import warnings

import numpy as np

from pyne import nucname
from pyne import material
import pyne.simplesim.nestedgeom as ng

def _validate_name(value, isunique=False):
    if value.find(' ') != -1:
        raise ValueError("The property ``name`` cannot contain spaces. "
                         "User provided {0}.".format(value))
    if isunique:
        raise StandardError("This is a unique card, meaning only one card"
                            " of this type can be found in a ``definition``. "
                            "Accordingly, the name is read-only.")
    if value == '':
        raise ValueError("The ``name`` property of the cannot be empty.")


class ICard(object):
    """This class is not used by the user. Abstract base class for all cards.
    All cards have a name property and a comment() method.

    """
    # This line makes this class an abstract base class (ABC).
    __metaclass__ = abc.ABCMeta
    # Useful for cell cards and the TMP card in MCNP.
    # k is Boltzmann's constant in eV/K, and T is temperature in K
    kelvin2kT = 8.6173423e-11
    secs2shakes = 1e+8

    def __init__(self, name, unique=False, *args, **kwargs):
        """ 

        Parameters
        ----------
        name : str
            Name of this instance of this card, used for referencing this card
            from others and modifying this card after it has been added to a
            definition (see :py:mod:`pyne.simplesim.definition`). The only
            known restriction on the format of the name is that it cannot
            contain spaces. 

        """
        # Somebody on the internet said to use super() even when only
        # subclassing from object. I'm commenting this out because it prevents
        # MaterialCustom from working.
        #super(ICard, self).__init__()
        self.unique = unique
        if self.unique:
            self._name = name
        else:
            self.name = name

        # TODO Do we actually want to do this? Perhaps not; it may lead to
        # recursion.
    def __str__(self):
        return self.comment()

    # All subclasses must define a comment() method.
    @abc.abstractmethod
    def comment(self):
        """All cards define a comment describing the content of the card."""
        raise NotImplementedError

    def mcnp(self, float_format, sim):
        # sim is an instance of
        # :py:class:`pyne.simplesim.definition.SimulationDefinition`.
        raise NotImplementedError("Card {0}.".format(self.name))

    def json(self):
        """Returns something that can be used to encode this class in a JSON
        file.

        """
        # return self.__repr__()
        return self.__dict__

    @property
    def name(self):
        return self._name

    @name.setter
    def name(self, value):
        _validate_name(value, self.unique)
        self._name = value


class Cell(ICard):
    """A cell is a region of space filled with a material. If requesting a void
    cell, the ``material``, ``density``, and ``density_units`` attributes are
    all set to None (as by default).

    """
    def __init__(self, name, region, material=None, density=None,
                 density_units=None, universe=None, fill=None, lattice=None,
                 *args, **kwargs):
        """\

        Parameters
        ----------
        name : str
            See :py:class:`ICard`.
        region : :py:class:`Region` subclass
            Defines the region of space that this cell occupies (see
            :py:class:`Region`).
        material : :py:class:`pyne.simplesim.cards.Material`, None for void
            A material definition using the :py:mod:`pyne.material` module.
            For use here, the material's :py:attr:`name` property must be set
            to something other than '' and must be unique within the materials.
            See :py:class:`pyne.material.Material`.
        density : float, None for void
            Density for the material, in units of density_units.
        density_units : str, None for void
            Either 'g/cm^3', or 'atoms/b/cm'.
        universe : 2-element tuple, optional
            To make this cell part of a universe, provide the name of that
            universe. Universes are often composed of more than one cell, and
            so this name is likely to appear on other cell cars as well. Tuple
            contains:

            1. (str) name of universe
            2. (bool) True for truncating, False otherwise.

            Refer to :py:class:`Universes` for more information.
        fill : str, or 4-element tuple, optional
            The name of a universe with which to fill this cell. Typically,
            cells that are filled with a universe are void themselves. If
            tuple (only makes sense for lattice cells), it contains:

            1. 2-element list of x-index bounds ([], or [0, 0] for 1 row/col)
            2. 2-element list of y-index bounds ([], or [0, 0] for 1 row/col)
            3. 2-element list of z-index bounds ([], or [0, 0] for 1 row/col)
            4. 1-D, 2-D, or 3-D list of universe names, or None to use the
               real world universe. See examples.

        lattice : str, optional
            Makes this cell is a repeating lattice. Can be either
            'hexahedra' for 6-sided elements or 'hexagonal' for 8-sided
            elements.

        Examples
        --------
        Suppose we have surface cards ``surfA`` and ``surfB``, and material
        card ``matA``. The following creates a void cell on the
        negative side of ``surfA`` and the positive side of ``surfB`` (see
        :py:class:`Region` to learn how to create regions)::

            cellA = Cell('A', surfA.neg & surfB.pos)

        The following cell is filled with ``matA``::

            cellB = Cell('B', surfA.neg & surfB.pos, matA, 10.0, 'g/cm^3')
            cellB = Cell('B', surfA.neg & surfB.pos, matA, 0.5, 'atoms/b/cm')

        Note that if a material is specified, a density and density units must
        also be provided The following is not allowed::

            cellB = Cell('B', surfA.neg & surfB.pos, matA)
            cellB = Cell('B', surfA.neg & surfB.pos, matA, density=1)
            cellB = Cell('B', surfA.neg & surfB.pos, matA, 
                    density_units='g/cm^3')
            cellB = Cell('B', surfA.neg & surfB.pos, density=1)
            cellB = Cell('B', surfA.neg & surfB.pos, 
                    density_units='g/cm^3')

        The following creates a new universe, and places two cells in it::

            cellC = Cell('C', surfA.neg & surfB.pos, matA, density=1,
                    density_units='g/cm^3',
                    universe=('unitcell', True))
            cellD = Cell('D', surfC.neg, matC, density=2,
                    density_units='g/cm^3',
                    universe=('unitcell', False))

        Note that the universe isn't actually `created` until this card is
        added to a system::

            sys = definition.SystemDefinition()
            sys.add_cell(cellD)

        Now, we can create a cell filled with the 'unitcell' universe::

            cellE = Cell('E', fill='unitcell')

        If we had two universes 'A' and 'B', we could fill a finite lattice
        cell as follows::
       
            #                 element (0, 0, 0)   element (4, 0, 0)
            #                /                   /
            univ_names = [['A', 'B', 'A', 'B', 'A'],
                          ['B', 'A', 'B', 'A', 'B']]
            #               ^                   ^
            #               element (0, 1, 0)   element (4, 1, 0)
            cellF = Cell('F', surfA.neg, lattice='hexahedra',
                    fill=([0, 4], [0, 1], [], univ_names))

        Here is an example for a 3-D finite lattice::

            #                 (0, 0, -1)  (2, 0, -1)
            #                /           /
            univ_names = [[['A',  'B',  'A'], 
                           ['B',  'A',  'B']],
                          [['A',  None, 'A'],
                           ['B',  'A',  'B']]]
            #                ^           ^
            #                (0, 1, 0)  (2, 1, 0)
            #
            cellF = Cell('F', surfA.neg, lattice='hexahedra',
                    fill=([0, 2], [0, 1], [-1, 0], univ_names))

        """
        # TODO decide how I will do cross-referencing.
        super(Cell, self).__init__(name, *args, **kwargs)
        self.region = region
        self.material = material
        self.density = density
        self.universe = universe
        self.fill = fill
        self.lattice = lattice
        self.density_units = density_units
        if ((self.material and not (self.density and self.density_units)) or
                (self.density and not (self.density_units and self.material)) or
                (self.density_units and not (self.density and self.material))):
            raise ValueError("If specifying a material, ``material``, "
                    "``density``, and ``density_units`` must all be "
                    "specified.")

    def comment(self, period=True):
        string = "Cell {0!r}: region {1}, ".format(
                self.name, self.region.comment())
        if self.material and self.density and self.density_units:
            string += "material {0!r} density {1:g} {2}".format(
                    self.material.name, self.density, self.density_units)
        else:
            string += "void"
        # Univ, fill, lattice
        if self.universe:
            card = Universes(self.name, self.universe[0], self.universe[1])
            string += " univ{0}".format(card._comment_unit(self.name))
        if self.fill:
            if type(self.fill) is not tuple:
                # if IS tuple, mcnp_long_tuple is called.
                card = Fill(self.name, self.fill)
                string += " fill{0}".format(card._comment_unit(self.name))
        # Lattice
        if self.lattice:
            card = Lattice(self.name, self.lattice)
            string += " lattice{0}".format(card._comment_unit(self.name))
        if self.fill and type(self.fill) is tuple:
            string += " fill long format"
        if period: string += "."
        return string

    def mcnp(self, float_format, sim):
        # Card number.
        if self.name not in sim.sys.cells:
            raise StandardError("Cell {0!r} not in simulation.".format(
                self.name))
        formatstr = "{{: <{0}d}}".format(
                int(np.log10(len(sim.sys.cells))) + 1)
        string = formatstr.format(sim.sys.cell_num(self.name))
        if self.material and self.density and self.density_units:
            # Material number.
            string += " {0} ".format(sim.sys.material_num(self.material.name))
            # Density, with units prefix.
            if self.density_units == 'g/cm^3':
                string += float_format % -self.density
            else:
                string += float_format % self.density
        else:
            # Void.
            string += " 0"
        # Print surfaces.
        string += " {0}".format(self.region.mcnp(sim))
        # Univ, fill, lattice
        if self.universe:
            card = Universes(self.name, self.universe[0], self.universe[1])
            string += " U={0}".format(card._mcnp_unit(
                float_format, sim, self.name))
        if self.fill:
            if type(self.fill) is not tuple:
                # if IS tuple, mcnp_long_tuple is called.
                card = Fill(self.name, self.fill)
                string += " FILL={0}".format(card._mcnp_unit(
                    float_format, sim, self.name))
        # Lattice
        if self.lattice:
            card = Lattice(self.name, self.lattice)
            string += " LAT={0}".format(card._mcnp_unit(
                float_format, sim, self.name))
        if isinstance(self, CellMCNP):
            string += self._mcnp_keywords(float_format, sim)
        if self.fill and type(self.fill) is tuple:
            string += self._mcnp_long_fill(float_format, sim)
        return string

    def _mcnp_long_fill(self, float_format, sim):
        x_idxs = self.fill[0]
        y_idxs = self.fill[1]
        z_idxs = self.fill[2]
        univ_names = self.fill[3]
        string = "\n{0}FILL={1[0]}:{1[1]} {2[0]}:{2[1]} {3[0]}:{3[1]}".format(
                5 * " ",
                x_idxs if len(x_idxs) == 2 else [0, 0],
                y_idxs if len(y_idxs) == 2 else [0, 0],
                z_idxs if len(z_idxs) == 2 else [0, 0])
        n_dim = (len(x_idxs) == 2) + (len(y_idxs) == 2) + (len(z_idxs) == 2)
        if n_dim == 3: counter = z_idxs[0]
        for zdim in univ_names:
            # Allow easy looping if 1 dimension only.
            if type(zdim) is str: zdim = [zdim]
            if n_dim == 2: string += "\n{0}".format(4 * " ")
            if n_dim == 3:
                string += "\n{0}$ k = {1}".format(5 * " ", counter)
                counter += 1
            for ydim in zdim:
                # Allow easy looping if 2 dimensions only.
                if type(ydim) is str: ydim = [ydim]
                if n_dim == 3: string += "\n{0}".format(4 * " ")
                for univ_name in ydim:
                    string += " {0}".format(
                        sim.sys.universe_num(univ_name) if univ_name else 0)
        return string

    @property
    def region(self): return self._region

    @region.setter
    def region(self, obj): self._region = obj

    @property
    def material(self): return self._material

    @material.setter
    def material(self, obj):
        # This check is redundant.
        if obj and type(obj) is not Material:
            raise ValueError("The property ``material`` must be instance "
                    "of ``pyne.simplesim.cards.Material``. "
                    "User provided {0}.".format( obj))
        if obj and obj.name == '':
            raise ValueError("The ``name`` property of the material cannot "
                    "be empty.")
        self._material = obj

    @property
    def density(self): return self._density

    @density.setter
    def density(self, value): self._density = value

    @property
    def density_units(self): return self._density_units

    @density_units.setter
    def density_units(self, value):
        if (value and value != 'g/cm^3' and value != 'atoms/b/cm'):
            raise ValueError("The property ``density_units`` must be either "
                    "'g/cm^3' or 'atoms/b/cm'. User provided "
                    "{0!r}".format(value))
        self._density_units = value

    @property
    def universe(self): return self._universe

    @universe.setter
    def universe(self, value): self._universe = value

    @property
    def fill(self): return self._fill

    @fill.setter
    def fill(self, value): self._fill = value

    @property
    def lattice(self): return self._lattice

    @lattice.setter
    def lattice(self, value): self._lattice = value


class CellMCNP(Cell):
    """A cell card with keyword options that are available in MCNP.
    Thus, it
    only makes sense to use this card if writing an input for MCNP. A number of
    the keyword arguments are for a particular particle. The particles
    available are given in :py:class:`Particle`. The user provides the full
    name of the particle, as given as keys in :py:attr:`Particle.mcnp_abbrev`.
    The card will then use the appropriate particle designator when writing the
    card.
    
    See :py:class:`Cell` for universes, fill, and lattices.
    
    Note this card was written with MCNPX version 2.7 in mind.

    .. inheritance-diagram:: pyne.simplesim.cards.CellMCNP

    """
    # TODO Sphinx documentation should not list all keyword arguments.
    # TODO Transformation: allow default arguments.
    def __init__(self, name, region, material=None,
                 density=None, density_units=None,
                 temperature=None, volume=None,
                 importance=None,
                 exp_transform=None,
                 forced_coll=None,
                 weight_win_bound=None,
                 dxtran_contrib=None,
                 photon_weight=None,
                 fission_turnoff=None,
                 det_contrib=None,
                 transformation=None,
                 user_custom=None,
                 **kwargs
                 ):
        """\

        Parameters
        ----------
        name : str
            See :py:class:`ICard`.
        region : :py:class:`Region`
            See :py:class:`Cell`.
        material : :py:class:`pyne.simplesim.cards.Material`, None for void
            See :py:class:`Cell`.
        density : float, None for void
            See :py:class:`Cell`.
        density_units : str, None for void
            See :py:class:`Cell`.
        temperature : float, otional [Kelvin]
            Temperature of the cell for thermal treatment,
            :py:class:`Temperature`, **TMP**.
        volume : float, optional [cm^3]
            Volume of the cell, :py:class:`Volume`, **VOL**.
        importance : 2-element tuple, or list of tuples, optional
            Particle importance, :py:class:`Importance`, **IMP**. The tuple
            contains:
            
            1. (str) particle name (see :py:class:`Particle`)
            2. (int) the importance

            Refer to :py:class:`Importance` for more information.
            To specify this input for more
            than one particle, provide a list of these tuples.
        exp_transform : 4-element tuple, or list of tuples, optional
            An exponential transform, :py:class:`ExponentialTransform`,
            **EXT**. The tuple contains:

            1. (str) particle name (see :py:class:`Particle`)
            2. (str/float) stretch
            3. (str) direction
            4. (str) sign ('toward', 'away')
            
            Refer to :py:class:`ExponentialTransform` for the form of these
            inputs.
            To specify this input for more
            than one particle, provide a list of these tuples.
        forced_coll : 3-element tuple, or list of tuples, optional
            Forced collisions, :py:class:`ForcedCollision`, **FCL**. The tuple
            contains:

            1. (str) particle name (see :py:class:`Particle`)
            2. (float) probability
            3. (bool) only entering the cell triggers forced collision

            Refer to :py:class:`ForcedCollision` for the
            form of these inputs.
            To specify this input for more
            than one particle, provide a list of these tuples.
        weight_win_bound : 4-element tuple, or list of tuples, optional
            Weight window lower bound, :py:class:`WeightWindowBound`, **WWN**.
            The tuple contains:

            1. (str) particle name (see :py:class:`Particle`)
            2. (int) energy index
            3. (int) time index
            4. (float/str) lower bound, or 'killall'

            Refer to :py:class:`WeightWindowBound` for the form of these
            inputs.
            To specify this input for more than one particle, or energy/time
            index, provide a list of these tuples.
        dxtran_contrib : 3-element tuple, or list of tuples, optional
            Probability of contribution to a DXTRAN sphere,
            :py:class:`DXTRANContribution`, **DXC**.  The tuple contains:
            
            1. (str) particle name (see :py:class:`Particle`)
            2. (str/None) DXTRAN sphere name (None for all spheres)
            3. (float) probabilility of contribution
            
            Refer to :py:class:`DXTRANContribution` for the form of these
            inputs.
            To specify this input for more than one particle or DXTRAN sphere,
            provide a list of these tuples.
        photon_weight : 2-element tuple, or list of tuples, optional
            Threshold weight of photons that are produced at neutron
            collisions, :py:class:`PhotonWeight`, **PWT**. The tuple contains:

            1. (str/float) 'off', 'one', or weight threshold parameter
            2. (bool) ``pre_weight``; relevant is 1st element is float
               (optional)

            Refer to :py:meth:`PhotonWeight.set` for the form of these
            inputs.
        fission_turnoff : str, optional
            Fission turnoff, :py:class:`FissionTurnoff`, **NONU**. The allowed
            values are: 'capture-gamma', 'real-gamma', and 'capture-nogamma'.
            Refer to :py:meth:`FissionTurnoff` for more information.
        det_contrib : tuple of str and float, optional
            Detector contribution, :py:class:`DetectorContribution`, **PD**.
            The tuple contains:

            1. (str) name of :py:class:`IDetector` tally obtaining contribution 
               from this cell
            2. (float) probability of contribution
            
            Refer to :py:meth:`DetectorContribution` for the form of these
            inputs.
            To specify this input for more than one tally, provide a list of
            these tuples.
        transformation : str or 4-element tuple, optional
            Cell transformation, **TRCL**. If str, it is the name of a
            :py:class:`Transformation` card. If tuple, it contains:
            
            1. (list/:py:class:`np.array`) displacement vector
            2. (3 x 3 list/:py:class:`np.array`/:py:class:`np.matrix`) 
               rotation matrix
            3. (bool) If True, 'aux-in-main', else 'main-in-aux' (optional)
            4. (bool) Rotation in degrees if True, in cosines if False
               (optional)

             Refer to :py:class:`Transformation` for the form of these inputs
             and the default values for the optional arguments.
        user_custom : str, optional
            This string is appended to the end of the MCNP card. It is possible
            that the user will want to append a string to the end of the MCNP
            card, given the limited support of keyword arguments. This is
            perhaps most useful if the user wants to specify a keyword for a
            particle that is not supported by any of the keyword arguments
            above.

        Examples
        --------
        The following sets the temperature of the cell to 600 K, its volume
        to 1.5 cm^3, and the neutron impotartance to 1.::

            cellA = CellMCNP('A', surfA.neg & surfB.pos,
                    matA, 10.0, 'g/cm^3',
                    temperature=600, volume=1,
                    importance=('neutron', 1))

        The following sets the neutron importance to 1 and the photon
        importance to 0::

            cellA = CellMCNP('A', surfA.neg & surfB.pos,
                    matA, 10.0, 'g/cm^3',
                    importance=[('neutron', 1), ('photon', 0)])

        The following sets an exponential transform for neutrons with stretch
        'capture-to-total' toward the current direction of travel of the
        particle, and a differente transform for photons::

            cellA = CellMCNP(..., exp_transform=[
                    ('neutron', 'capture-to-total', 'currdir', 'toward'),
                    ('photon', 0.5, 'x', 'away')])

        The direction can be the name of a vector (e.g. `vec1`) on the
        :py:class:`Vector` card, as long as the :py:class:`Vector` has been
        added to the simulation before an input is generated::

            cellA = CellMCNP(..., 
                    exp_transform=('neutron', 0.5, 'vec1', 'away'))
            vec = Vector()
            vec.set('vec1', [0, 0, 0])

        We can request forced collision for neutrons that enter the cell and
        as part of weight games::

            cellA = CellMCNP(..., forced_coll=('neutron', 0.5, False))

        The following specifies the lower bound for neutrons at energy index,
        3, time index 1 (assuming that the :py:class:`WeightWindowEnergies`
        card is in the simulation)::

            cellA = CellMCNP(..., 
                    weight_win_bound=('neutron', 3, 1, 'killall'))

        The following specifies the probability that neutrons from this cell
        will tally in DXTRAN sphere 'sph1' (assuming a sphere with this name will be
        in the simulation, on the :py:class:`DXTRANSpheres` card)::
        
            cellA = CellMCNP(..., dxtran_contrib=('neutron', 'det1', 0.5))
            sim = definition.MCNPDefinition(...)
            det1 = PointDetector('det1', ...)
            sim.add_tally(det1)

        Here are two examples of specifying a photon weight threshold (the bool
        can be provided if the first element is a float, but defaults to False
        otherwise)::

            cellA = CellMCNP(..., photon_weight=('one',))
            cellA = CellMCNP(..., photon_weight=(0.5,))
            cellA = CellMCNP(..., photon_weight=(0.5, True))

        The following turns off fission in this cell, but still requests that
        gammas are generated by fission interactions::
        
            cellA = CellMCNP(..., fission_turnoff='capture-gamma')

        The following specifies the probability that this cell contributes to
        detector `det1`, which must be in the simulation::

            cellA = CellMCNP(..., det_contrib=('det1', 0.5))

        A transformation can be provided by referring to a
        :py:class:`Transformation` card::

            cellA = CellMCNP(..., transformation='transA')
            transA = Transformation('transA', ...)
            sim = definition.MCNPSimulation(...)
            sim.add_transformation(transA)

        Alternatively, a transformation can be specified right on the keyword.
        Here we have used the default values for the last two inputs::

            cellA = CellMCNP(..., transformation=([1, 0, 0], np.eye(3))

        Here, we want to change the 4th element from its default, value, so we
        must specify the 3rd element even though it has its default value::

            cellA = CellMCNP(..., transformation=([1, 0, 0], np.eye(3),
                    True, True))

        If the user wants to supply an exponential transform keyword, with a
        transform of '0.7V2', on their own, they can do the following::

            cellA = CellMCNP('A', surfA.neg & surfB.pos,
                    matA, 10.0, 'g/cm^3',
                    importance=[('neutron', 1), ('photon', 0)],
                    user_custom='EXT:N 0.7V2')

        and the user_custom string will be printed at the end of the cell card.

        See :py:class:`Cell` for more examples.

        """
        # TODO check that the importance isn't set for the same particle
        # multiple times.
        super(CellMCNP, self).__init__(name, region, material, density,
                                       density_units, **kwargs)
        # Assign keyword arguments.
        self.temperature = temperature
        self.volume = volume
        self.importance = importance
        self.exp_transform = exp_transform
        self.forced_coll = forced_coll
        self.weight_win_bound = weight_win_bound
        self.dxtran_contrib = dxtran_contrib
        self.photon_weight = photon_weight
        self.fission_turnoff = fission_turnoff
        self.det_contrib = det_contrib
        self.transformation = transformation
        self.user_custom = user_custom

    def comment(self):
        string = super(CellMCNP, self).comment(period=False)
        # temperature
        if self.temperature:
            card = Temperature(self.name, self.temperature)
            string += " TMP=" + card._comment_unit(self.name)
        # volume
        if self.volume:
            card = Volume(self.name, self.volume)
            string += " VOL=" + card._comment_unit(self.name)
        # importance
        if self.importance:
            for entry in self._make_list(self.importance):
                card = Importance(entry[0], self.name, entry[1])
                string += " IMP:{0}={1}".format(
                        Particle(entry[0]).mcnp(),
                        card._comment_unit(self.name))
        # exp_transform
        if self.exp_transform:
            for entry in self._make_list(self.exp_transform):
                card = ExponentialTransform(
                        entry[0], self.name, entry[1], entry[2], entry[3])
                string += " EXT:{0}={1}".format(
                        Particle(entry[0]).mcnp(),
                        card._comment_unit(self.name))
        # forced_coll
        if self.forced_coll:
            for entry in self._make_list(self.forced_coll):
                card = ForcedCollision(entry[0], self.name, entry[1], entry[2])
                string += " FCL:{0}={1}".format(
                        Particle(entry[0]).mcnp(),
                        card._comment_unit(self.name))
        # weight_win_bound
        if self.weight_win_bound:
            for entry in self._make_list(self.weight_win_bound):
                card = WeightWindowBound(
                        entry[0], entry[1], entry[2], self.name, entry[3])
                string += " WWN({0},{1}):{2}={3}".format(
                        entry[1], entry[2], Particle(entry[0]).mcnp(),
                        card._comment_unit(self.name, entry[1], entry[2]))
        # dxtran_contrib
        if self.dxtran_contrib:
            dxtran_contrib = self._make_list(self.dxtran_contrib)
            for entry in self._make_list(self.dxtran_contrib):
                card = DXTRANContribution(entry[0], entry[1], self.name, entry[2])
                string += " DXC{0!r}:{1}={2}".format( 
                        entry[1], Particle(entry[0]).mcnp(),
                        card._comment_unit(self.name))
        # photon_weight
        if self.photon_weight:
            # TODO this might cause issues.
            card = PhotonWeight()
            card.set(self.name, *list(self.photon_weight))
            string += " PWT={0}".format(card._comment_unit(self.name))
        # fission_turnoff
        if self.fission_turnoff:
            card = FissionTurnoff(self.name, self.fission_turnoff)
            string += " NONU={0}".format(card._comment_unit(self.name))
        # det_contrib
        if self.det_contrib:
            for entry in self._make_list(self.det_contrib):
                card = DetectorContribution(entry[0], self.name, entry[1])
                string += " PD for tally {0!r}={1}".format(
                        entry[0], card._comment_unit(self.name))
        # transform
        if self.transformation:
            # For brevity.
            entry = self.transformation
            string += " TRCL "
            if type(entry) is str:
                string += "{0!r}".format(entry)
            else:
                card = Transformation('temp', *list(entry))
                string += card._comment_unit()[:-1]
        # user_custom
        if self.user_custom: string += " and user's custom input"
        return string + "."

    def mcnp(self, float_format, sim):
        return super(CellMCNP, self).mcnp(float_format, sim)

    def _mcnp_keywords(self, float_format, sim):
        string = ""
        # temperature
        if self.temperature:
            card = Temperature(self.name, self.temperature)
            string += " TMP={0}".format(
                    card._mcnp_unit(float_format, sim, self.name))
        # volume
        if self.volume:
            card = Volume(self.name, self.volume)
            string += " VOL={0}".format(
                    card._mcnp_unit(float_format, sim, self.name))
        # importance
        if self.importance:
            for entry in self._make_list(self.importance):
                card = Importance(entry[0], self.name, entry[1])
                string += " IMP:{0}={1}".format(
                        Particle(entry[0]).mcnp(), 
                        card._mcnp_unit(float_format, sim, self.name))
        # exp_transform
        if self.exp_transform:
            for entry in self._make_list(self.exp_transform):
                card = ExponentialTransform(
                        entry[0], self.name, entry[1], entry[2], entry[3])
                string += " EXT:{0}={1}".format(
                        Particle(entry[0]).mcnp(), 
                        card._mcnp_unit(float_format, sim, self.name))
        # forced_coll
        if self.forced_coll:
            for entry in self._make_list(self.forced_coll):
                card = ForcedCollision(entry[0], self.name, entry[1], entry[2])
                string += " FCL:{0}={1}".format(
                        Particle(entry[0]).mcnp(), 
                        card._mcnp_unit(float_format, sim, self.name))
        # weight_win_bound
        if self.weight_win_bound:
            for entry in self._make_list(self.weight_win_bound):
                card = WeightWindowBound(
                        entry[0], entry[1], entry[2], self.name, entry[3])
                # Need the next line to obtain the linear index.
                card._find_n_energies(sim)
                string += " WWN{0}:{1}={2}".format( 
                        card.i_linear(entry[1], entry[2]),
                        Particle(entry[0]).mcnp(), 
                        card._mcnp_unit(float_format, self.name, entry[1],
                            entry[2]))
        # dxtran_contrib
        if self.dxtran_contrib:
            for entry in self._make_list(self.dxtran_contrib):
                card = DXTRANContribution(entry[0], entry[1], self.name,
                            entry[2])
                string += " DXC{0}:{1}={2}".format( 
                        card.sph_index(sim), Particle(entry[0]).mcnp(), 
                        card._mcnp_unit(float_format, sim, self.name))
        # photon_weight
        if self.photon_weight:
            card = PhotonWeight()
            card.set(self.name, *list(self.photon_weight))
            string += " PWT={0}".format(
                    card._mcnp_unit(float_format, sim, self.name))
        if self.fission_turnoff:
            card = FissionTurnoff(self.name, self.fission_turnoff)
            string += " NONU={0}".format(
                card._mcnp_unit(float_format, sim, self.name))
        # det_contrib
        if self.det_contrib:
            for entry in self._make_list(self.det_contrib):
                card = DetectorContribution(entry[0], self.name, entry[1])
                string += " PD{0}={1}".format(
                    sim.tally_num(entry[0]) * 10 + 5,
                    card._mcnp_unit(float_format, sim, self.name))
        # transform
        if self.transformation:
            # For brevity.
            entry = self.transformation
            string += " "
            if type(entry) is str:
                string += "TRCL={0}".format(sim.transformation_num(entry))
            else:
                if len(entry) == 4 and entry[3]:
                    string += "*"
                card = Transformation('temp', *list(entry))
                string += "TRCL ({0})".format(card._mcnp_unit(float_format))
        # user_custom
        if self.user_custom: string += " {0}".format(self.user_custom)
        return string
                
    def _make_list(self, arg):
        if type(arg) is list: return arg
        else: return [arg]

    @property
    def temperature(self): return self._temperature

    @temperature.setter
    def temperature(self, value):
        self._temperature = value

    @property
    def volume(self): return self._volume

    @volume.setter
    def volume(self, value): self._volume = value

    @property
    def importance(self): return self._importance

    @importance.setter
    def importance(self, value): self._importance = value

    @property
    def exp_transform(self): return self._exp_transform

    @exp_transform.setter
    def exp_transform(self, value): self._exp_transform = value

    @property
    def forced_coll(self): return self._forced_coll

    @forced_coll.setter
    def forced_coll(self, value): self._forced_coll = value

    @property
    def weight_win_bound(self): return self._weight_win_bound

    @weight_win_bound.setter
    def weight_win_bound(self, value): self._weight_win_bound = value

    @property
    def dxtran_contrib(self): return self._dxtran_contrib

    @dxtran_contrib.setter
    def dxtran_contrib(self, value): self._dxtran_contrib = value
    
    @property
    def photon_weight(self): return self._photon_weight
    
    @photon_weight.setter
    def photon_weight(self, value):
        if value and type(value) is not tuple:
            raise ValueError("Photon weight keyword argument must be a tuple.")
        self._photon_weight = value
    
    @property
    def fission_turnoff(self): return self._fission_turnoff
    
    @fission_turnoff.setter
    def fission_turnoff(self, value): self._fission_turnoff = value

    @property
    def det_contrib(self): return self._det_contrib
    
    @det_contrib.setter
    def det_contrib(self, value): self._det_contrib = value

    @property
    def transformation(self): return self._transformation
    
    @transformation.setter
    def transformation(self, value): self._transformation = value

    @property
    def user_custom(self): return self._user_custom
    
    @user_custom.setter
    def user_custom(self, value): self._user_custom = value
    

class IUniverse(ICard):
    """This class is not used by the user. Abstract base class for all
    universe cards.

    """
    # TODO
    # Ideally universes would be their own objects.
    __metaclass__ = abc.ABCMeta

    def __init__(self, name, *args, **kwargs):
        pass


class UniverseByRegion(IUniverse):
    """Unimplemented. Use cell card keywords."""
    # TODO
    def __init__(self, name, region):
        pass


class UniverseByLattice(IUniverse):
    """Unimplemented. Use cell card keywords."""
    # TODO
    def __init__(self, name, lattice):
        pass


class Lattice(ICard):
    """Unimplemented. Use cell card keywords."""
    # TODO
    def __init__(self, name, geom, universe):
        pass


class LatticeByArray(ICard):
    """Unimplemented. Use cell card keywords."""
    # TODO
    # TODO support of 3D arrays.
    def __init(self, name, geom, xindices, yindices, zindices,
               universe_array):
        pass
        

class Material(ICard):
    """In MCNP, this is the **M** card.
    Adds the attributes
    :py:attr:`description`, :py:attr:`tables` and the methods
    :py:meth:`comment` and :py:meth:`mcnp` to
    :py:class:`pyne.material.Material`. The :py:attr:`name` must be provided
    before the card is added to a system. The user can specify a description
    that is printed in an appropriate location in the input file.

    .. inheritance-diagram:: pyne.simplesim.cards.Material

    """
    def __init__(self, mat, name=None, description=None, tables=None, *args, **kwargs):
        """See :py:class:`pyne.material.Material` for superclass parameters.

        Parameters
        ----------
        mat : pyne.material.Material
            Material to base this card on.
        name : str, optional
            Must either be present in mat.attrs or supplied here.
        description : str, optional
            A description of this material that perhaps explains where the
            material came from (whether it's recycled, any references, etc.).
            If not supplied here, may be present in mat.attrs.
        tables : dict of (nucname, str) pairs
            Sometimes it is necessary to specify a library/table identifier for
            a given nuclide. These can be provided in this dictionary. Leave
            out the period. If not supplied here, may be present in mat.attrs.  
            See examples.

        Examples
        --------
        The usage of this card is meant to ensure that only valid instances of
        :py:class:`pyne.material.Material` are used as cards.  Additionally, 
        there are two new methods::

            originstory = "I found this water in a well a few years ago."
            h2o = Material(pyne.material.from_atom_frac({10010: 2.0, 'O16': 1.0}),
                           name='water', description=originstory)
            h2o.tables = {10010: '71c'}
            sys.add_material(h2o)

        Alternatively, the tables can be specified with the constructor::

            h2o = Material(pyne.material.from_atom_frac({10010: 2.0, 'O16': 1.0}),
                           name='water', tables={10010: '71c'})

        The ``nucname``'s used for ``tables`` can be different from those used
        for ``comp``::

            h2o = Material(pyne.material.from_atom_frac({10010: 2.0, 'O16': 1.0}),
                           name='water', tables={'H1': '71c'})

        """
        name = name if name is not None else  mat.attrs['name']
        self._mat = mat
        super(Material, self).__init__(name=name, description=description,
                                       tables=tables, *args, **kwargs)
        self.mat = mat
        if description is not None:
            self.description = description
        if tables is not None:
            self.tables = tables
        # Find longest table ID. Used in card printing for prettiness.

    def comment(self): 
        if self.name == '':
            raise ValueError("The ``name`` property of the material cannot be empty.")
        s = "Material {0!r}".format(self.name)
        if self.description:
            s += ": {0}".format(self.description)
        else:
            s += "."
        return s

    def mcnp(self, float_format, sim):
        s = "M{0}".format(sim.sys.material_num(self.name))
        for nuc, den in self._mat.to_atom_frac().items():
            # ZAID.
            s += "\n     {: 6d}".format(nucname.mcnp(nuc))
            # Table ID. Loop allows flexible keys for tables.
            flag = False 
            for key in self.tables:
                if nucname.mcnp(key) == nucname.mcnp(nuc):
                    flag = True
                    s += ".{0}".format(self.tables[key])
            if not flag:
                # +1 for he decimal point.
                s += (self._max_table_len + 1) * " "
            # Concentration/density.
            s += 2 * " " + float_format % den
            # Nuclide name.
            s += " $ {0}".format(nucname.name(nuc))
        return s

    @property
    def mat(self):
        return self._mat

    @mat.setter
    def mat(self, value):
        _validate_name(self.mat.attrs['name'], self.unique)
        self._mat = value
        if 'tables' in value.attrs:
            self._max_table_len = max([len(v) for v in value.attrs['tables'].values()])
        else:
            self._max_table_len = 0

    @property
    def name(self):
        return self._mat.attrs['name']

    @name.setter
    def name(self, value):
        _validate_name(value, self.unique)
        self._mat.attrs['name'] = value

    @property
    def description(self):
        if 'description' not in self._mat.attrs:
            self._mat.attrs['description'] = None
        return self._mat.attrs['description']

    @description.setter
    def description(self, value):
        self._mat.attrs['description'] = value

    @property
    def tables(self):
        if 'tables' not in self._mat.attrs:
            self._mat.attrs['tables'] = {}
        return self._mat.attrs['tables']

    @tables.setter
    def tables(self, value):
        self._mat.attrs['tables'] = value
        self._max_table_len = max([len(val) for val in value.values()])


class MaterialMCNP(Material):
    """Unimplemented; would automate the selection of table identifiers. """
    # TODO
    pass


class ISurface(ICard):
    """This class is not used by the user. Abstract base class for all
    surface cards.

    The Surface superclass contains properties to set the surface as reflecting
    or white. For codes other than MCNPX, reflecting or white surfaces may be
    specified on a separate boundary condition card (i.e. in Serpent) or may
    not even be available. For other codes, then, the appropriate
    :py:mod:`inputfile` class needs to pick up this information and print the
    appropriate string to the code's input file, or in the latter case return
    an exception.

    .. inheritance-diagram:: pyne.simplesim.cards.ISurface

    """
    # TODO support rotation.
    __metaclass__ = abc.ABCMeta
    
    def __init__(self, name, reflecting=False, white=False, **kwargs):
        """
        Parameters
        ----------
        name : str
            See :py:class:`ICard`.
        reflecting : bool, optional
            The surface has a reflective boundary condition.
        white : bool, optional
            The surface has a white boundary condition (reflection with a
            cosine distribution with respect to the surface normal).

        """
        super(ISurface, self).__init__(name, **kwargs)
        self.reflecting = reflecting
        self.white = white 
        if self.reflecting and self.white:
            raise ValueError("The user set the surface to be reflecting AND "
                    "white, but can only be neither or one of the two.")

    @abc.abstractmethod
    def comment(self, title):
        return "{0} {1!r}:{2}{3} ".format(title, self.name, 
                " reflecting." if self.reflecting else "",
                " white." if self.white else "")

    @abc.abstractmethod
    def mcnp(self, float_format, sim, keystring):
        if self.name not in sim.sys.surfaces:
            raise StandardError("Surface {0!r} not in simulation.".format(
                self.name))
        formatstr = "{0}{1}{{: <{2}d}} {3:<4}".format(
                "*" if self.reflecting else "",
                "+" if self.white else "",
                int(np.log10(len(sim.sys.surfaces))) + 2 - \
                    (self.reflecting or self.white),
                keystring)
        return formatstr.format(sim.sys.surface_num(self.name))

    def shift(self, vector):
        """Translates the surface. This is an abstract method, and must be
        defined by each surface card. Shifts may not be permitted for all
        surfaces in all directions, and in such cases an exception is raised.
        
        Parameters
        ----------
        vector : 3-element list or :py:class:`np.array`, float [centimeters]
            The elements specify translation in the x, y, and z directions, 
            in this order.

        Examples
        --------
        Both of the following lines shifts the surface along the x axis by 3 cm::

            surf.shift([3, 0, 0])
            surf.shift(np.array([3, 0, 0]))
        
        """
        raise NotImplementedError("Shifting not implemented for "
                "{0}.".format(self))

    def stretch(self, vector):
        """Stretches (scales) the surface from the origin. This is an abstract
        method, and must be defined by each surface card. Stretches may not be
        permitted for all surfaces in all directions, and in such cases an
        exception is raised.

        Parameters
        ----------
        vector : 3-element list or :py:class:`np.array`, float [unitless]
            The elements specify a stretch in the x, y, and z directions, in
            this order. A zero in any of the directions indicates that no
            stretch is done in that direction. Negative values are allowed, and
            represent reflections.

        Examples
        --------
        Both of the following lines stretch the surface along the y axis by a
        factor of 2. The x and z directions are unaffected::

            surf.stretch([0, 2, 0])
            surf.stretch(np.array([0, 2, 0]))
        
        """
        raise NotImplementedError("Stretching not implemented for "
                "{0}.".format(self))

    @property
    def neg(self):
        """A property that creates and returns a
        :py:class:`RegionLeaf` that can then be used in
        boolean arithmetic with subclasses of
        :py:class:`Region`. The region is define as the
        space on the side of the surface that has a negative sense.

        In the expected typical usage of the :py:mod:`pyne.simplesim` package,
        regions are constructed using these properties.

        For more information, see :py:class:`Region` and
        :ref:`usersguide_simplesim`.

        Examples
        --------
        The following shows a simple case of how a more complex region can be
        constructed from regions returned by this property::

            reg1 = surf1.neg
            reg2 = surf2.neg
            reg3 = reg1 & reg2
            reg4 = reg1 | reg2

        """
        return RegionLeaf(self, False)

    @property
    def pos(self):
        """Similar to :py:attr:`neg`, except the resulting
        :py:class:`RegionLeaf` is on the side of the surface with a positive
        sense.

        Examples
        --------
        The following shows a simple case of how a more complex region can be
        constructed from regions returned by this property::

            reg1 = surf1.pos
            reg2 = surf2.pos
            reg3 = reg1 & reg2
            reg4 = reg1 | reg2
        
        """
        return RegionLeaf(self, True)
    
    @property
    def reflecting(self): return self._reflecting

    @reflecting.setter
    def reflecting(self, value):
        if type(value) is not bool:
            raise TypeError("The property ``reflecting`` must be of "
                    "boolean type. User provided {0}.".format(value))
        self._reflecting = value

    @property
    def white(self): return self._white

    @white.setter
    def white(self, value):
        if type(value) is not bool:
            raise TypeError("The property ``white`` must be of "
                    "boolean type. User provided {0}.".format(value))
        self._white = value


class IAxisSurface(ISurface):
    """This class is not used by the user. Abstract base class for all simple
    axis-aligned surfaces. Accordingly, such classes share the cartesian_axis
    property.

    .. inheritance-diagram:: pyne.simplesim.cards.IAxisSurface
    
    """
    __metaclass__ = abc.ABCMeta

    def __init__(self, name, cartesian_axis, *args, **kwargs):
        """
        Parameters
        ----------
        name : str
            See :py:class:`ICard`.
        cartesian_axis : str
            Either 'x', 'X', 'y', 'Y', 'z', or 'Z'. Regardless of input, it is
            stored as lower-case. The meaning depends on the surface.
        reflecting : bool, optional
            See :py:class:`ISurface`
        white : bool, optional
            See :py:class:`ISurface`

        """
        super(IAxisSurface, self).__init__(name, *args, **kwargs)
        self.cartesian_axis = cartesian_axis

    @abc.abstractmethod
    def comment(self, *args):
        return super(IAxisSurface, self).comment(*args)

    @abc.abstractmethod
    def mcnp(self, *args):
        return super(IAxisSurface, self).mcnp(*args)
    
    @abc.abstractmethod
    def shift(self, vector):
        """See :py:meth:`ISurface.shift`."""
        raise NotImplementedError

    @abc.abstractmethod
    def stretch(self, vector):
        """See :py:meth:`ISurface.stretch`."""
        raise NotImplementedError

    @property
    def cartesian_axis(self): return self._cartesian_axis

    @cartesian_axis.setter
    def cartesian_axis(self, value):
        if type(value) is not str:
            raise ValueError("AxisCylinder's cartesian_axis property must be "
                    "a string. User provided {0}.".format(value))
        if (value.lower() != 'x' and value.lower() != 'y' and 
                value.lower() != 'z'):
            raise ValueError("AxisCylinder's cartesian_axis property must be "
                    "'x', 'X', 'y', 'Y', or 'z', 'Z'. "
                    "User provided '{0}'.".format(value))
        self._cartesian_axis = value.lower()


class AxisCylinder(IAxisSurface):
    """Cylinder aligned with and centered on one of the Cartesian axes.
    
    .. inheritance-diagram:: pyne.simplesim.cards.AxisCylinder
    
    """

    # TODO if this is shifted, then it becomes not an axis-cylinder.
    def __init__(self, name, cartesian_axis, radius, **kwargs):
        """
        Parameters
        ----------
        name : str
            See :py:class:`ICard`.
        cartesian_axis : str
            The axis with which the cylinder is aligned and centered.
            See :py:class:`IAxisSurface`.
        radius : float [centimeters]
            Radius of the cylinder.
        reflecting : bool, optional
            See :py:class:`ISurface`
        white : bool, optional
            See :py:class:`ISurface`

        Examples
        --------
        The following creates a cylinder aligned with the z axis, going through
        centered on the z axis, with a radius of 0.4 cm::

            cyl = AxisCylinder('mycyl', 'z', 0.4)

        """
        super(AxisCylinder, self).__init__(name, cartesian_axis, **kwargs)
        self.radius = radius

    def comment(self):
        string = super(AxisCylinder, self).comment("Axis cylinder")
        return string + ("aligned and centered on {0} axis, "
                "with radius {1:g} cm (diameter {2:g} cm).".format(
                    self.cartesian_axis, self.radius, 2 * self.radius))

    def mcnp(self, float_format, sim):
        string = super(AxisCylinder, self).mcnp(float_format, sim,
                "C{0}".format(self.cartesian_axis.upper()))
        string += float_format % self.radius
        return string
    
    def shift(self, vector):
        """See :py:meth:`ISurface.shift`. Axis cylinders can only be shifted along
        their axis, and even in such cases the shift has no effect. However,
        such a shift must be permitted in case this surface is part of a region
        that is being shifted.

        Examples
        --------
        The following is okay (where we have imported :py:mod:`numpy` as ``np``)::

            cyl = AxisCylinder('mycyl', 'z', 0.4)
            cyl.shift([0, 0, 3])
            cyl.shift(np.array([0, 0, 3]))

        The following do not work:

            cyl.shift([3, 0, 0])
            cyl.shift([0, 3, 3])
            
        
        """
        # Flag for exception.
        iserror = False
        if self.cartesian_axis == 'x' and (vector[1] != 0 or vector[2] != 0):
            iserror = True
            dirs = ('x', 'y', 'z')
        if self.cartesian_axis == 'y' and (vector[0] != 0 or vector[2] != 0):
            iserror = True
            dirs = ('y', 'x', 'z')
        if self.cartesian_axis == 'z' and (vector[0] != 0 or vector[1] != 0):
            iserror = True
            dirs = ('z', 'x', 'y')
        if iserror:
            raise ValueError("A cylinder aligned with the {0[0]} axis cannot "
                    "be shifted in the {0[1]} or {0[2]} directions.".format(
                        dirs))

    def stretch(self, vector):
        """See :py:meth:`ISurface:stretch`. Axis cylinders can be stretched in
        the direction of their axis, which has no effect (permitted in case
        this surface is part of a region that is being stretched), or can be
        stretched `uniformly` in the plane perpendicular to its axis.
        
        Examples
        --------
        The following stretches are okay for a cylinder aligned with the x axis
        (where we have imported :py:mod:`numpy` as ``np``)::
            
            cyl = AxisCylinder('mycyl', 'z', 0.4)
            cyl.stretch([0, 0, 2])
            cyl.stretch([3, 3, 0])
            cyl.stretch(np.array([3, 3, 2]))

        However, the following would cause the cylinder to lose its
        circular cross section, which cannot be accommodated::

            cyl.stretch([0, 3, 0])
            cyl.stretch([2, 3, 1])
        
        """
        # TODO allow some slop between the same two values for a uniform
        # perpendicular stretch.
        # Flag for exception.
        iserror = False
        # 'out' is used in the exception below.
        if self.cartesian_axis == 'x':
            if vector[1] != vector[2]:
                iserror = True
                out = ('y', vector[1], 'z', vector[2], 'x')
            elif vector[1] != 0:
                self.radius *= vector[1]
        if self.cartesian_axis == 'y':
            if vector[0] != vector[2]:
                iserror = True
                out = ('x', vector[0], 'z', vector[2], 'y')
            elif vector[0] != 0:
                self.radius *= vector[0]
        if self.cartesian_axis == 'z':
            if vector[0] != vector[1]:
                iserror = True
                out = ('x', vector[0], 'y', vector[1], 'z')
            elif vector[0] != 0:
                self.radius *= vector[0]
        if iserror:
            raise ValueError("Stretches perpendicular to the axis must be "
                    "uniform in the two perpendicular directions. User "
                    "provided {0[0]} stretch {0[1]:g} and {0[2]} stretch "
                    "{0[3]:g} for a {0[4]}-aligned cylinder.".format(out))

    @property
    def radius(self):
        return self._radius

    @radius.setter
    def radius(self, value):
        if value <= 0:
            raise ValueError("The ``radius`` property must be "
                    "positive. User provided {0:g}.".format(value))
        self._radius = value


class AxisPlane(IAxisSurface):
    """Plane perpendicular to one of the Cartesian axes.

    .. inheritance-diagram:: pyne.simplesim.cards.AxisPlane
    
    """

    def __init__(self, name, cartesian_axis, position, **kwargs):
        """
        Parameters
        ----------
        name : str
            See :py:class:`ICard`.
        cartesian_axis : str
            The axis to which the plane is perpendicular.
            See :py:class:`IAxisSurface`.
        position : float [centimeters]
            Position of the plane along :py:attr:`cartesian_axis`.
        reflecting : bool, optional
            See :py:class:`ISurface`
        white : bool, optional
            See :py:class:`ISurface`

        Examples
        --------
        The following creates a plane perpendicular to the x axis, 3 cm along
        the positive x axis with a reflecting boundary condition::

           plane = AxisPlane('myplane', 'x', 3, reflecting=True) 

        """
        super(AxisPlane, self).__init__(name, cartesian_axis, **kwargs)
        self.position = position
    
    def comment(self):
        string = super(AxisPlane, self).comment("Axis plane")
        return string + "{0} = {1:g} cm.".format(
                self.cartesian_axis, self.position)

    def mcnp(self, float_format, sim):
        string = super(AxisPlane, self).mcnp(float_format, sim,
                "P{0}".format(self.cartesian_axis.upper()))
        string += float_format % self.position
        return string

    def shift(self, vector):
        """See :py:meth:`ISurface.shift`. Axis planes can be shifted in any
        direction, but only shifts along their axis have an effect.

        Examples
        --------
        The following has the effect of shifting the plane's position to x = 6
        cm::

           plane = AxisPlane('myplane', 'x', 3)
           plane.shift([3, 0, 0])

        The following has no effect, but is allowed::

           plane.shift([0, 3, 2])

        """
        if self.cartesian_axis == 'x':   self.position += vector[0]
        elif self.cartesian_axis == 'y': self.position += vector[1]
        elif self.cartesian_axis == 'z': self.position += vector[2]

    def stretch(self, vector):
        """See :py:meth:`ISurface.stretch`. Axis planes can be stretched in any
        direction, but only stretches along their axis have an effect. The
        position of the plane is scaled by the stretch factor.

        Examples
        --------
        The following has the effect of moving the plane's position to x = 9 cm::

            plane = AxisPlane('myplane', 'x', 3)
            plane.stretch([3, 0, 0])
        
        The following has no effect, but is allowed::

            plane.stretch([0, 3, 2])

        """
        if self.cartesian_axis == 'x' and vector[0] != 0:
            self.position *= vector[0]
        elif self.cartesian_axis == 'y' and vector[1] != 0:
            self.position *= vector[1]
        elif self.cartesian_axis == 'z' and vector[2] != 0:
            self.position *= vector[2]

    @property
    def position(self): return self._position

    @position.setter
    def position(self, value): self._position = value


class IMacrobody(ISurface):
    """This class is not used by the user. Abstract base class for all
    macrobody cards. Macrobodies are an MCNP concept. The classes here account
    for the notion of facets.

    .. inheritance-diagram:: pyne.simplesim.cards.IMacrobody

    """
    __metaclass__ = abc.ABCMeta

    # TODO abstract method for obtaining "sub"-surfaces. This is one reason why
    # you'd pass objects rather than just object names around: I could pass
    # IMacrobody.top, etc.
    # TODO facets are not functional because they have the same name as their
    # macrobody, so when the System Definition looks it up, it only finds the
    # macrobody and not the facet. I don't know how to solve this yet.
    def __init__(self, name, *args, **kwargs):
        super(IMacrobody, self).__init__(name, *args, **kwargs)

    @abc.abstractmethod
    def comment(self, *args):
        return super(IMacrobody, self).comment(*args)

    @abc.abstractmethod
    def mcnp(self, *args):
        return super(IMacrobody, self).mcnp(*args)

    def facet(self, descriptor):
        """
        Returns
        -------
        facet : :py:class:`Facet`
            The facet can be used as a normal surface, though logically there
            is no MCNP surface card for just a facet. The facet can be used in
            constructing a :py:class:`IRegion`.

        """
        raise NotImplementedError("Facets not implemented for {0!r}.".format(
            self))


class Facet(ISurface):
    """Used in conjunction with macrobodies. See
    :py:meth:`pyne.simplesim.definition.SystemDefinition.surface_num` to see
    how the appropriate surface number is obtained. The name of the surface is
    the same as that of the macrobody.
    
    """
    # TODO ideally we would dynamically subclass this using whatever class
    # surface is.
    # TODO facets do not work right now because they have the same name as an
    # already-existant surface.
    def __init__(self, macrobody, descriptor, number):
        """The user does not use this constructor. The user creates facets
        using :py:meth:`IMacrobody.facet`.

        """
        raise Exception("Facets are not implemented properly yet.")
        super(Facet, self).__init__(macrobody.name)
        self.macrobody = macrobody
        self.descriptor = descriptor
        self.number = number

    def comment(self):
        return "{0}. facet {0!r}.".format(
                self.macrobody.comment(), self.descriptor)

    def mcnp(self):
        raise Exception("There is no notion of a facet surface card.")

    def shift(self, vector):
        self.macrobody.shift(vector)

    def stretch(self, vector):
        self.macrobody.stretch(vector)


class Parallelepiped(IMacrobody):
    """Rectangular parallelepiped in which all surfaces are parallel to the
    cartesian axes.

    .. inheritance-diagram::pyne.simplesim.cards.Parallelepiped

    """
    facets = {'east'  : 1,
              'west'  : 2,
              'north' : 3,
              'south' : 4,
              'top'   : 5,
              'bottom': 6}

    def __init__(self, name, xmin, xmax, ymin, ymax, zmin, zmax, **kwargs):
        """
        Parameters
        ----------
        name : str
            See :py:class:`ICard`.
        xmin, xmax, ymin, ymax, zmin, zmax : float [centimeters]
            Bounds of the parallelepiped in the given direction. The min value
            must be less than the max value. Setting both min and max in a
            given direction to 0 indicates the parallelepiped is infinite in
            that direction.
        reflecting : bool, optional
            See :py:class:`ISurface`.
        white : bool, optional
            See :py:class:`ISurface`.

        Examples
        --------
        The following creates a cube centered at the origin with 4 cm sides::

            pp = Parallelepiped('mypp', -2, 2, -2, 2, -2, 2)

        """
        super(Parallelepiped, self).__init__(name, **kwargs)
        self.xlims = np.array([xmin, xmax])
        self.ylims = np.array([ymin, ymax])
        self.zlims = np.array([zmin, zmax])

    def comment(self):
        string = super(Parallelepiped, self).comment("Parallelepiped")
        return string + ("[{0[0]:g}, {0[1]:g}] x "
                "[{1[0]:g}, {1[1]:g}] x [{2[0]:g}, {2[1]:g}] cm.".format(
                    self.xlims, self.ylims, self.zlims))

    def mcnp(self, float_format, sim):
        string = super(Parallelepiped, self).mcnp(float_format, sim, "RPP")
        formatstr =  "{0} {0}  {0} {0}  {0} {0}".format(float_format)
        return string + formatstr % (
                tuple(self.xlims) + tuple(self.ylims) + tuple(self.zlims))

    def facet(self, descriptor):
        """See :py:meth:`IMacrobody.facet`.
        
        Parameters
        ----------
        descriptor : str
            One of the following:

            - 'east': x = xmax plane
            - 'west': x = xmin plane
            - 'north': y = ymax plane
            - 'south': y = ymin plane
            - 'top': z = zmax plane
            - 'bottom': z = zmin plane

        Returns
        -------
        facet : :py:class:`Facet`
            The facet can be used as a normal surface, though logically there
            is no MCNP surface card for just a facet.

        """
        if descriptor not in self.facets:
            raise ValueError("Facet descriptor {0!r} unacceptable.".format(
                descriptor))
        return Facet(self, descriptor, self.facets[descriptor])

    def shift(self, vector):
        """See :py:meth:`ISurface.shift`.
        
        Examples
        --------
        The following::

            pp = Parallelepiped('mypp', -2, 2, -2, 2, -2, 2)
            pp.shift([2, 0, 0])

        creates a parallelepiped bounded by [0, 4] x [-2, 2] x [-2, 2].

        """
        self.xlims += vector[0]
        self.ylims += vector[1]
        self.zlims += vector[2]

    def stretch(self, vector):
        """See :py:meth:`ISurface.stretch`. Handling reflections (negative
        stretch factors) requires additional consideration for this surface,
        but is implemented.
        
        Examples
        --------
        The following::

            pp = Parallelepiped('mypp', 0, 4, -2, 2, -2, 2)
            pp.stretch([2, 0, 3])

        creates a parallelepiped bounded by [0, 8] x [-2, 2] x [-6, 6].
        Consider the reflection of the following parallelepiped
        about the z axis::
        
            pp = Parallelepiped('mypp', 0, 4, -2, 2, -3, 6)
            pp.stretch([0, 0, -1])

        This results in bounds of [0, 4] x [-2, 2] x [-6, 3].
        
        """
        if vector[0] != 0:
            if vector[0] > 0: self.xlims *= vector[0]
            else: self.xlims = vector[0] * self.xlims[::-1]
                # Stretch factor is negative, swap limits.
        if vector[1] != 0:
            if vector[1] > 0: self.ylims *= vector[1]
            else: self.ylims = vector[1] * self.ylims[::-1]
        if vector[2] != 0:
            if vector[2] > 0: self.zlims *= vector[2]
            else: self.zlims = vector[2] * self.zlims[::-1]

    @property
    def xlims(self): return self._xlims

    @xlims.setter
    def xlims(self, value):
        if value[0] > value[1]:
            raise ValueError("The value of xmin, {0:g}, is greater than "
                    "that of xmax, {1:g}.".format(value[0], value[1]))
        self._xlims = value

    @property
    def ylims(self): return self._ylims

    @ylims.setter
    def ylims(self, value):
        if value[0] > value[1]:
            raise ValueError("The value of ymin, {0:g}, is greater than "
                    "that of ymax, {1:g}.".format(value[0], value[1]))
        self._ylims = value

    @property
    def zlims(self): return self._zlims

    @zlims.setter
    def zlims(self, value):
        if value[0] > value[1]:
            raise ValueError("The value of zmin, {0:g}, is greater than "
                    "that of zmax, {1:g}.".format(value[0], value[1]))
        self._zlims = value


class Cuboid(Parallelepiped):
    """Same exact thing as a :py:class:`Parallelepiped`. This class is provided
    because the name is shorter, and thus may be preferred by those who fancy
    brevity.

    """
    def __init__(self, *args, **kwargs):
        super(Cuboid, self).__init__(*args, **kwargs)


class IRegion(ICard):
    """This class is not used by the user. Abstract base class for
    all regions. Represents a volume (space) confined by unions and
    intersections of surfaces.
    
    """
    # TODO transformation functions
    # TODO Complement functionality can be added by overloading the
    # __not__ operator and defining a complement boolean property that is set
    # by the __not__ operator.
    __metaclass__ = abc.ABCMeta

    def __init__(self, name, *args, **kwargs):
        super(IRegion, self).__init__(name, *args, **kwargs)

    @abc.abstractmethod
    def comment(self):
        raise NotImplementedError

    def set(self, region):
        # TODO deepcopy?
        if issubclass(region, RegionLeaf):
            self = RegionLeaf(region.surface, region.pos_sense, self.name)
        elif issubclass(region, RegionOr):
            self = RegionOr(region.left_child, region.right_child, self.name)
        elif issubclass(region, RegionAnd):
            self = RegionAnd(region.left_child, region.right_child, self.name)

    def __and__(self, arg):
        """Calls :py:meth:`intersect`."""
        return self.intersect(arg)

    def __or__(self, arg):
        """Calls :py:meth:`union`."""
        return self.union(arg)

    def intersect(self, arg):
        """
        Parameters
        ----------
        arg : :py:class:`IRegion` subclass
            The region with which to intersect this region.

        Returns
        -------
        reg : :py:class:`RegionAnd`
            A region that has been intersected with the input.

        """
        return RegionAnd(self, arg)

    def union(self, arg):
        """
        Parameters
        ----------
        arg : :py:class:`IRegion` subclass
            The region with which to union this region.

        Returns
        -------
        reg : :py:class:`RegionAnd`
            A region that has been unioned with the input.

        """
        return RegionOr(self, arg)

    def shift(self, vector):
        """Shifts all the surfaces that this region is composed of. Note that
        the surfaces themselves are modified, and are not copied, so the
        dimensions of the surfaces in other regions is also modified.

        Parameters
        ----------
        vector : 3-element list or :py:class:`np.array`
            See :py:class:`ISurface`.
        
        """
        # TODO make copies of the surfaces?
        self.left_child.shift(vector)
        self.right_child.shift(vector)


    def stretch(self, vector):
        """Stretches all the surfaces that this region is composed of. Note that
        the surfaces themselves are modified, and are not copied, so the
        dimensions of the surfaces in other regions is also modified.

        Parameters
        ----------
        vector : 3-element list or :py:class:`np.array`
            See :py:class:`ISurface`.

        """
        # TODO make copies of the surfaces?
        self.left_child.stretch(vector)
        self.right_child.stretch(vector)

    def walk(self, leaf_func, and_func=None, or_func=None):
        """Walks the Region tree, calling ``leaf_func`` if the Region is a
        ``RegionLeaf``, ``and_func`` if the Region is a ``RegionAnd``, and
        ``or_func`` if the Region is a ``RegionOr``. This is used in
        ``pyne.definition`` to discover all surfaces in the Region.

        """
        if isinstance(self, RegionLeaf):
            leaf_func.im_func(leaf_func.im_self, self)
        else:
            self.left_child.walk(leaf_func)
            if and_func and isinstance(self, and_func):
                and_func.im_func(and_func.im_self, self)
            if or_func and isinstance(self, or_func):
                or_func.im_func(or_func.im_func, self)
            self.right_child.walk(leaf_func)


class IRegionBool(IRegion):
    """This class is not used by the user. Abstract base class for
    :py:class:`RegionAnd` and :py:class:`RegionOr`.

    """
    __metaclass__ = abc.ABCMeta

    def __init__(self, left_child, right_child, name='<Empty>', *args, **kwargs):
        super(IRegionBool, self).__init__(name, *args, **kwargs)
        self.left_child = left_child
        self.right_child = right_child

    @abc.abstractmethod
    def comment(self, midchar):
        return ("(" + self.left_child.comment() + " " + midchar + " " +
                self.right_child.comment() + ")")

    @abc.abstractmethod
    def mcnp(self, sim, midchar):
        return ("(" + self.left_child.mcnp(sim) + midchar + 
                self.right_child.mcnp(sim) + ")")

    @property
    def left_child(self):
        """Left child"""
        return self._left_child

    @left_child.setter
    def left_child(self, value):
        self._left_child = value

    @property
    def right_child(self):
        """Right child."""
        return self._right_child

    @right_child.setter
    def right_child(self, value):
        self._right_child = value


class RegionAnd(IRegionBool):
    """The intersection of two regions, held in :py:attr:`left_child` and
    :py:attr:`right_child`.

    """
    def comment(self):
        return super(RegionAnd, self).comment('&')

    def mcnp(self, sim):
        return super(RegionAnd, self).mcnp(sim, ':')


class RegionOr(IRegionBool):
    """The union of two regions, held in :py:attr:`left_child` and
    :py:attr:`right_child`.

    """
    def comment(self):
        return super(RegionOr, self).comment('|')

    def mcnp(self, sim):
        return super(RegionOr, self).mcnp(sim, ' ')


class RegionLeaf(IRegion):
    """The region of space to one side of a surface (constructed by calling
    :py:attr:`ISurface.neg` or :py:attr:`ISurface.pos`).

    """
    def __init__(self, surface, pos_sense, name='<Empty>'):
        """
        Parameters
        ----------
        surface : :py:class:`ISurface` subclass
            The surface that defines the edge of this region.
        pos_sense : bool
            True if this region is on the side of the surface with a positive
            sense, False otherwise.
        name : str, optional
            The name of this surface. Currently, it is not really needed.

        """
        super(RegionLeaf, self).__init__(name)
        self.surface = surface
        self.pos_sense = pos_sense

    def comment(self):
        if self.pos_sense:  prefix = '+'
        else:               prefix = '-'
        return prefix + self.surface.name

    def mcnp(self, sim):
        if self.pos_sense:  prefix = ''
        else:               prefix = '-'
        return "{0}{1}".format(prefix, sim.sys.surface_num(self.surface.name))

    def shift(self, vector):
        """Calls the shift method of the surface."""
        self.surface.shift(vector)

    def stretch(self, vector):
        """Calls the stretch method of the surface."""
        self.surface.stretch(vector)

    @property
    def surface(self):
        return self._surface

    @surface.setter
    def surface(self, value):
        self._surface = value

    @property
    def pos_sense(self):
        return self._pos_sense

    @pos_sense.setter
    def pos_sense(self, value):
        # TODO this is probably not okay by the proponents of duck typing.
        if type(value) is not bool:
            raise TypeError("User provided a value for pos_sense that is "
                    "not of boolean type.")
        self._pos_sense = value


class IMisc(ICard):
    """This class is not used by the user. Abstract base class for unnumbered
    cards that do not fit in other categories.

    """
    __metaclass__ = abc.ABCMeta

    def __init__(self, name, *args, **kwargs):
        super(IMisc, self).__init__(name, *args, **kwargs)

    @abc.abstractmethod
    def comment(self):
        raise NotImplementedError


class ScatteringLaw(IMisc):
    """Scattering law for a material. Unique card for a given material, with
    name `scatlaw-<matname>`. In MCNP, this is the **MT** card.

    .. inheritance-diagram:: pyne.simplesim.cards.ScatteringLaw

    """
    def __init__(self, mat_name, libraries):
        """
        Parameters
        ----------
        mat_name : str
            Name of the material for which this card applies.
        libraries : dict of :py:class:`nucname`: str pairs
            The keys are the nuclides on the material for which a library is
            being provided, and the values are the appropriate library
            identifiers as strings.

        Examples
        --------
        This specifies hydrogen bound in water, in MCNP::

           h2o = Material(name='water')
           h2o.from_atom_frac({10010: 2.0, 'O16': 1.0})
           sys.add_material(h2o)
           sl = ScatteringLaw('water', {'H1': 'lwtr.16t'})

        """
        super(ScatteringLaw, self).__init__('scatlaw-{0}'.format(mat_name),
                unique=True)
        self.mat_name = mat_name
        self.libraries = libraries

    def comment(self):
        string = "Scattering law {0!r}:".format(self.name)
        for nuc, lib in self.libraries.items():
            string += " {0}: {1},".format(nucname.name(nuc), lib)
        return string[:-1] + "."

    def mcnp(self, float_format, sim):
        string = "MT{0}".format(sim.sys.material_num(self.mat_name))
        for nuc, lib in self.libraries.items():
            string += " {0}".format(lib)
        return string

    @property
    def libraries(self): return self._libraries

    @libraries.setter
    def libraries(self, value): self._libraries = value


class ISource(ICard):
    """This class is not used by the user. Abstract base class for all source
    cards.

    """
    __metaclass__ = abc.ABCMeta

    def __init__(self, name, *args, **kwargs):
        super(ISource, self).__init__(name, *args, **kwargs)

    @abc.abstractmethod
    def comment(self):
        raise NotImplementedError


class GeneralSource(ISource):
    """A general source. Unique with name `gensource`. In MCNP, this is the
    **SDEF** card. See the respective code (e.g. MCNP) for the default behavior
    of the arguments. Note that the x, y, and z coordinates of the
    source are separate inputs because they might each be on their own
    distributions. To learn about Python keyword arguments, visit
    `the Python docs
    <docs.python.org/tutorial/controlflow.html#keyword-arguments>`_. Each
    keyword argument is also an attribute of this class, and can be modified
    individually after the card is created::

        src = GeneralSource()
        src.cell = 'fuel'

    A source is added to a simulation as follows::
        
        sys = definition.SystemDefinition()
        ...
        sim = definition.MCNPSimulation(sys)
        sim.add_source(src)

    .. inheritance-diagram:: pyne.simplesim.cards.GeneralSource

    Limitations:

    - Keyword cannot depend on other keywords, as they can in MCNP.
    
    """
    # TODO non-ambiguous names required.
    # TODO cell card for repeated structures is not simple.
    def __init__(self, **keyword_args):
        # Using 'keyword_args' instead of 'kwargs' to help PyNovices.
        """
        Parameters
        ----------
        particle : str, int, :py:class:`Distribution` name, optional (default: None)
            Name of a particle, 'spont-fiss-by-neut', 'spont-fiss-by-hist', or
            'spont-phot', integer to specify a heavy ion ZAID, or the name of a
            :py:class:`Distribution` card. If
            the string is the name of a distribution in the simulation, then
            that distribution is used, and the other possible inputs are
            ignored. Therefore, be careful with how distribution cards are
            name.
        cell : str, :py:class:`Distribution` name, optional (default: None)
            Name of a cell in which the source exists, or the name of a
            :py:class:`Distribution`. If both a cell and distribution have the
            same card name (ambiguous input), an exception is raised. **CEL**
            in MCNP.
        surface : str, :py:class:`Distribution` name, optional (default: None)
            Name of a surface, or the name of a :py:class:`Surface`. If both a
            surface and a distribution have the same card name (ambiguous
            input), an exception is
            raised. **SUR** in MCNP.
        energy : float, :py:class:`Distribution` name, optional (default: None) [MeV] 
            Particle energy, or a distribution of energies. **ERG** in MCNP.
        time : float, :py:class:`Distribution` name, (default: None) [seconds]
            Time at which the source particles are created, or a distribution
            of times. **TME** in MCNP.
        ref_dir : 3-element list, optional (default: None)
            Reference vector for sampling particle's direction of travel.
            **VEC** in MCNP.
        cosine : float, :py:class:`Distribution` name, optional (default: None)
            Cosine of the angle between the reference vector
            :py:attr:`ref_dir` and the particle's direction of travel, or
            distribution. **DIR**
            in MCNP. **DIR** in MCNP.
        normal_sign : str, optional (default: None)
            'pos' for a positive surface normal, 'neg' for a negative surface
            normal. **NRM** in MCNP.
        ref_pos : 3-element list, :py:class:`Distribution` name, optional (default: None) [cm]
            Reference point for position sampling. **POS** in MCNP, or
            distribution. For
            cylinders, it is the point at the center of the base of the
            cylinder.
        offset : float, :py:class:`Distribution` name, optional (default: None) [cm]
            **EXT** in MCNP. Meaning varies, or distribution. For cylinders, it
            is the length of the cylinder.
        radius : float, :py:class:`Distribution` name, optional (default: None) [cm]
            Radial distance of sampled points from :py:attr:`ref_pos` or
            :py:attr:`axis`, or distribution. **RAD** in MCNP. For cylinders,
            it is the radius of the cylinder.
        axis : 3-element list, :py:class:`Distribution` name, optional (default: None)
            Reference vector or distribution for use in conjunction with
            :py:attr:`offset` and :py:attr:`radius`. **AXS** in MCNP. For
            cylinders, it is the axis of the cylinder.
        x : float, :py:class:`Distribution` name, optional (default: None) [cm]
            x-coordinate or distribution for source position sampling. **X** in
            MCNP.
        y : float, :py:class:`Distribution` name, optional (default: None) [cm]
            y-coordinate or distribution for source position sampling. **Y** in
            MCNP.
        z : float, :py:class:`Distribution` name, optional (default: None) [cm]
            z-coordinate or distribution for source position sampling. **Z** in
            MCNP.
        cookie_cutter : str, :py:class:`Distribution` name, optional (default: None)
            Name of a cell to use as cookie-cutter, or distribution Be wary of
            naming cells and distributions with the same name. **CCC** in MCNP.
        area : float, optional (default: None) [cm^2]
            Area of the surface, if this is a surface source. **ARA** in MCNP.
        weight : float, optional (default: None)
            Explicit particle weight. **WGT** in MCNP.
        transformation : str, optional (default: None)
            Name of a :py:class:`Transformation` in the simulation.
        eff : float, optional (default: None)
            **EFF** in MCNP. Explicit value. See MCNP manual.
        beam_emit : 3-element tuple, optional (default: None)
            **BEM** in MCNP. The tuple contains:

            1. parameter for x coordinate [pi-cm-rad]
            2. parameter for y coordinate [pi-cm-rad]
            3. distance [cm]

            See MCNP manual.
        beam_aper : 2-element tuple, optional
            **BAP** in MCNP. The tuple contains:

            1. parameter for x [cm]
            2. parameter for y [cm]

            See MCNP manual.
        user_custom : str, optional
            Appended to the end of :py:meth:`mcnp` output. The user can use
            this to account for functionality that is not otherwise provided by
            this class.

        Examples
        --------
        The following creates the default source::

            gs = GeneralSource()

        The following creates a photon source in cell 'cellA'::

            gs = GeneralSource(particle='photon', cell='cellA')

        The cookie-cutter input operates the same way as the cell input. The
        following uses a :py:class:`Distribution` for both the particle and the
        cell, assuming that the user has created distributions with the names
        'partdist' and 'celldist', and that there is no cell with the name
        'celldist'::

            gs = GeneralSource(particle='partdist', cell='celldist')

        Surfaces work the same way as cells. The following shows the
        specification of energies from a distribution, as well as a time::

            gs = GeneralSource(energy='energydist', time=2e10)

        The following shows the specification of a reference direction and the
        cosine from that reference vector::

            gs = GeneralSource(ref_dir=[1, 0, 0], cosine=0.5)

        The following shows how a distribution for ``x``, ``y``, and ``z`` can
        be requested (make sure to add the distributions to the simulation)::

            gs = GeneralSource(x='xdist', y='ydist', z='zdist')
            xd = Distribution('xdist', [0, 10], [0, 1])
            yd = Distribution('ydist', [-5, 5], [0, 1])
            zd = Distribution('zdist', [0, 1], [0, 1])

        The following shows how a transformation can be applied to the source
        (assuming the transformation is added to the simulation)::

            tr = Transformation('source', ...)
            gs = GeneralSource(transformation='source')


        The following shows how ``beam_emit`` and ``beam_aper`` are used::

            gs = GeneralSource(beam_emit=(1, 2, 3), beam_aper=(4, 5))

        All variables for which an example was not given function similarly to
        the other variables.

        """
        super(GeneralSource, self).__init__('gensource', unique=True)
        self.particle = keyword_args.get('particle', None)
        self.cell = keyword_args.get('cell', None)
        self.surface = keyword_args.get('surface', None)
        self.energy = keyword_args.get('energy', None)
        self.time = keyword_args.get('time', None)
        self.ref_dir = keyword_args.get('ref_dir', None)
        self.cosine = keyword_args.get('cosine', None)
        self.normal_sign = keyword_args.get('normal_sign', None)
        self.ref_pos = keyword_args.get('ref_pos', None)
        self.offset = keyword_args.get('offset', None)
        self.radius = keyword_args.get('radius', None)
        self.axis = keyword_args.get('axis', None)
        self.x = keyword_args.get('x', None)
        self.y = keyword_args.get('y', None)
        self.z = keyword_args.get('z', None)
        self.cookie_cutter = keyword_args.get('cookie_cutter', None)
        self.area = keyword_args.get('area', None)
        self.weight = keyword_args.get('weight', None)
        self.transformation = keyword_args.get('transformation', None)
        self.eff = keyword_args.get('eff', None)
        self.beam_emit = keyword_args.get('beam_emit', None)
        self.beam_aper = keyword_args.get('beam_aper', None)
        self.user_custom = keyword_args.get('user_custom', None)
        # Old constructor:
        #def __init__(self, cell=None, surface=None, energy=None, time=None, 
        #             ref_dir=None, cosine=None, normal_sign=None, ref_pos=None,
        #             offset=None, radius=None, axis=None, x=None, y=None,
        #             z=None, cookie-cutter=None, area=None, weight=None,
        #             transformation=None, ec*, direction*, normal_sign, 

    def comment(self):
        string = "General source {0!r}:".format(self.name)
        if self.particle:
            string += " particle={0},".format(self.particle)
        if self.cell:
            string += " cell={0},".format(self.cell)
        if self.surface:
            string += " surface={0},".format(self.surface)
        if self.energy:
            string += " energy={0},".format(self.energy)
        if self.time:
            string += " time={0},".format(self.time)
        if self.ref_dir:
            string += " ref. dir={0},".format(self.ref_dir)
        if self.cosine:
            string += " cosine={0},".format(self.cosine)
        if self.normal_sign:
            string += " {0} surf. normal,".format(self.normal_sign)
        if self.ref_pos:
            string += " ref. pos={0},".format(self.ref_pos)
        if self.offset:
            string += " offset={0},".format(self.offset)
        if self.radius:
            string += " radius={0},".format(self.radius)
        if self.axis:
            string += " axis={0},".format(self.axis)
        if self.x:
            string += " x={0},".format(self.x)
        if self.y:
            string += " y={0},".format(self.y)
        if self.z:
            string += " z={0},".format(self.z)
        if self.cookie_cutter:
            string += " c. cutter cell={0},".format(self.cookie_cutter)
        if self.area:
            string += " area={0},".format(self.area)
        if self.weight:
            string += " weight={0},".format(self.weight)
        if self.transformation:
            string += " trans.={0},".format(self.transformation)
        if self.eff:
            string += " eff.={0},".format(self.eff)
        if self.beam_emit:
            string += " beam emit.={0},".format(self.beam_emit)
        if self.beam_aper:
            string += " beam aper.={0},".format(self.beam_aper)
        if self.user_custom:
            string += " and user input."
        return string[:-1] + "."

    def mcnp(self, float_format, sim):
        vecformat = "{0} {0} {0}".format(float_format)
        string = "SDEF"
        if self.particle:
            string += " PAR="
            if self.particle in sim.dists:
                string += "D{0}".format(sim.dist_num(self.particle))
            elif self.particle == 'spont-fiss-by-neut': string += "SF"
            elif self.particle == 'spint-fiss-by-hist': string += "-SF"
            elif self.particle == 'spont-phot':         string += "SP"
            elif type(self.particle) is str:
                string += Particle(self.particle).mcnp()
            else:
                string += "{0}".format(self.particle)
        if self.cell:
            string += " CEL="
            #if isinstance(self.cell, Cell): 
            #    string += "{0}".format(sim.sys.cell_num(self.cell.name))
            #elif isinstance(self.cell, Distribution):
            #    string += "{0}".format(sim.dist_num(self.cell.name))
            if self.cell in sim.sys.cells and self.cell in sim.dists:
                raise Exception("Name {0!r} is both a cell and dist.".format(
                    self.cell))
            elif self.cell in sim.sys.cells:
                string += "{0}".format(sim.sys.cell_num(self.cell))
            elif self.cell in sim.dists:
                string += "D{0}".format(sim.dist_num(self.cell))
            else:
                raise Exception("Name {0!r} is not a cell or dist.".format(
                    self.cell))
        if self.surface:
            string += " SUR"
            if self.surface in sim.sys.surfaces and self.surface in sim.dist:
                raise Exception("Name {0!r} is both surface and dist.".format(
                    self.surface))
            elif self.surface in sim.sys.surfaces:
                string += "{0}".format(sim.sys.surface_num(self.surface))
            elif self.surface in sim.dists:
                string += "D{0}".format(sim.dist_num(self.surface))
            else:
                raise Exception("Name {0!r} is not a surface or dist.".format(
                    self.surface))
        if self.energy:
            string += " ERG="
            if type(self.energy) is str:
                string += "D{0}".format(sim.dist_num(self.energy))
            else: string += float_format % self.energy
        if self.time:
            string += " TME="
            if type(self.time) is str:
                string += "D{0}".format(sim.dist_num(self.time))
            else: string += float_format % (self.time * self.secs2shakes)
        if self.ref_dir:
            string += " VEC=" + vecformat % tuple(self.ref_dir)
        if self.cosine:
            string += " DIR="
            if type(self.cosine) is str:
                string += "D{0}".format(sim.dist_num(self.cosine))
            else: string += float_format % self.cosine
        if self.normal_sign:
            string += " NRM="
            if   self.normal_sign == 'neg': string += "-1"
            elif self.normal_sign == 'pos': string += "+1"
            else:
                raise Exception("Input ``normal_sign`` cannot be {0}.".format(
                    self.normal_sign))
        if self.ref_pos:
            string += " POS="
            if type(self.ref_pos) is str:
                string += "D{0}".format(sim.dist_num(self.ref_pos))
            else: formatstr = vecformat % tuple(self.ref_pos)
        if self.offset:
            string += " EXT="
            if type(self.offset) is str:
                string += "D{0}".format(sim.dist_num(self.offset))
            else: string += float_format % self.offset
        if self.radius:
            string += " RAD="
            if type(self.radius) is str:
                string += "D{0}".format(sim.dist_num(self.radius))
            else: string += float_format % self.radius
        if self.axis:
            string += " AXS="
            if type(self.axis) is str:
                string += "D{0}".format(sim.dist_num(self.axis))
            else: string += vecformat % self.axis
        if self.x:
            string += " X="
            if type(self.x) is str:
                string += "D{0}".format(sim.dist_num(self.x))
            else: string += float_format % self.x
        if self.y:
            string += " Y="
            if type(self.y) is str:
                string += "D{0}".format(sim.dist_num(self.y))
            else: string += float_format % self.y
        if self.z:
            string += " Z="
            if type(self.z) is str:
                string += "D{0}".format(sim.dist_num(self.z))
            else: string += float_format % self.z
        if self.cookie_cutter:
            string += " CCC="
            if (self.cookie_cutter in sim.sys.cells and 
                    self.cookie_cutter in sim.dists):
                raise Exception("Name {0!r} is both a cell and dist.".format(
                    self.cookie_cutter))
            elif self.cookie_cutter in sim.sys.cells:
                string += "{0}".format(sim.sys.cell_num(self.cookie_cutter))
            elif self.cookie_cutter in sim.dists:
                string += "D{0}".format(sim.dist_num(self.cookie_cutter))
            else:
                raise Exception("Name {0!r} is not a cell or dist.".format(
                    self.cookie_cutter))
        if self.area:
            string += " ARA="
            #if type(self.area) is str:
            #    string += "D{0}".format(sim.dist_num(self.area))
            #else:
            string += float_format % self.area
        if self.weight:
            string += " WGT="
            string += float_format % self.weight
        if self.transformation:
            string += " TR={0}".format(
                    sim.transformation_num(self.transformation))
        if self.eff:
            string += " EFF="
            string += float_format % self.eff
        if self.beam_emit:
            formatstr = " BEM={0} {0} {0}".format(float_format)
            string += formatstr % self.beam_emit
        if self.beam_aper:
            formatstr = " BAP={0} {0}".format(float_format)
            string += formatstr % self.beam_aper
        if self.user_custom:
            string += " {0}".format(self.user_custom)
        return string

    @property
    def particle(self): return self._particle

    @particle.setter
    def particle(self, value): self._particle = value

    @property
    def cell(self): return self._cell

    @cell.setter
    def cell(self, value): self._cell = value

    @property
    def surface(self): return self._surface

    @surface.setter
    def surface(self, value): self._surface = value

    @property
    def energy(self): return self._energy

    @energy.setter
    def energy(self, value): self._energy = value

    @property
    def time(self): return self._time

    @time.setter
    def time(self, value): self._time = value

    @property
    def ref_dir(self): return self._ref_dir

    @ref_dir.setter
    def ref_dir(self, value): self._ref_dir = value

    @property
    def cosine(self): return self._cosine

    @cosine.setter
    def cosine(self, value): self._cosine = value

    @property
    def normal_sign(self): return self._normal_sign

    @normal_sign.setter
    def normal_sign(self, value): self._normal_sign = value

    @property
    def ref_pos(self): return self._ref_pos

    @ref_pos.setter
    def ref_pos(self, value): self._ref_pos = value

    @property
    def offset(self): return self._offset

    @offset.setter
    def offset(self, value): self._offset = value

    @property
    def radius(self): return self._radius

    @radius.setter
    def radius(self, value): self._radius = value

    @property
    def axis(self): return self._axis

    @axis.setter
    def axis(self, value): self._axis = value

    @property
    def x(self): return self._x

    @x.setter
    def x(self, value): self._x = value

    @property
    def y(self): return self._y

    @y.setter
    def y(self, value): self._y = value

    @property
    def z(self): return self._z

    @z.setter
    def z(self, value): self._z = value

    @property
    def cookie_cutter(self): return self._cookie_cutter

    @cookie_cutter.setter
    def cookie_cutter(self, value): self._cookie_cutter = value

    @property
    def area(self): return self._area

    @area.setter
    def area(self, value): self._area = value

    @property
    def weight(self): return self._weight

    @weight.setter
    def weight(self, value): self._weight = value

    @property
    def transformation(self): return self._transformation

    @transformation.setter
    def transformation(self, value): self._transformation = value

    @property
    def eff(self): return self._eff
    
    @eff.setter
    def eff(self, value): self._eff = value

    @property
    def beam_emit(self): return self._beam_emit
    
    @beam_emit.setter
    def beam_emit(self, value): self._beam_emit = value

    @property
    def beam_aper(self): return self._beam_aper
    
    @beam_aper.setter
    def beam_aper(self, value): self._beam_aper = value

    @property
    def user_custom(self): return self._user_custom

    @user_custom.setter
    def user_custom(self, value): self._user_custom = value


class Distribution(ICard):
    """A distribution used in conjunction with source variables on the
    :py:class:`GeneralSource` card. The value of a
    source variable can be sampled from a distribution, described by this
    class. Distributions are added to a simulation using
    :py:attr:`definition.MCNPSimulation.add_dist`. In MCNP, this card is
    implemented with a pair of **SI** and **SP** cards.

    Limititations:

    - Embedded sources are not supported; to use this the user must provide an
      explicit string for ``keys``.
    - In MCNP, particles also have a number associated with them. It is
      necessary to use those numbers, instead of the abbreviations like 'N' or
      'H' (see :py:attr:`Particle.mcnp_abbrev`), to avoid conflict with the key
      settings. Specifying a particle name like 'proton' will not cause this
      code to find the particle number, but in such cases the user can provide
      the particle number on their own. Alternatively, the user can provide
      their ``keys`` input explicitly as a string .
    - Particles that are heavy ions are not supported through
      ``pyne.nucname``, but the user can provide the ZAID integer as part of the
      ``keys`` list

    """
    mcnp_keys = {'histogram': 'H',
                 'discrete' : 'L',
                 'pdf_indep': 'A',
                 'dist'     : 'S'}
    mcnp_vals = {'prob'     : 'D',
                 'cumprob'  : 'C',
                 'propvol'  : 'V',
                 'partint'  : 'W'}
    analytic = {'maxwell'   :  2,
                'watt'      :  3,  
                'gauss'     :  4,  
                'evap'      :  5,  
                'muir'      :  6,  
                'exp_decay' :  7,  
                'power'     : 21, 
                'exp'       : 31, 
                'gauss_beam': 41} 

    def __init__(self, name, keys, vals, key_setting=None, val_setting=None,
                 **kwargs):
        """The ``keys`` are the independent variable in the distribution, and
        the ``vals`` are probability values for the independent variable.

        Parameters
        ----------
        name : str
            See :py:class:`ICard`. Used to reference this card from the
            :py:class:`GeneralSource` card.
        key_setting : str, optional
            Specifies what the keys represent. The following settings are
            allowed:

            - 'histogram': bin upper boundaries. **H** in MCNP.
            - 'discrete': source variable values. **L** in MCNP.
            - 'pdf_indep': independent variable for probability density
              function. **A** in MCNP.
            - 'dist': names of other distributions. **S** in MCNP.
            - 'analytic': ``vals`` contains the parameters of a built-in
              analytic probability density function, and the function is
              determined by the ``key_setting``. In this case, ``keys`` must be
              an empty list ``[]`` or ``list()``, or the bounds of the domain
              for the function.

        val_setting : str, optional
            Specifies what the values represent. The following settings are
            allowed:

            - 'prob': bin probabilities. ** D** in MCNP.
            - 'cumprob': cumulative bin probabilities. **C** in MCNP.
            - 'propvol': probability proportional to cell vol. **V** in MCNP.
            - 'partint': intensities of particle sources. Values can be the
              name of a cell for spontaneous fission or spontaneous photon
              sources. **W** in MCNP.
            - if ``key_setting`` == 'analytic', this value is one of the
              following, with parameters specified in a list on ``vals``:
                 
                - 'maxwell': 1 parameter
                - 'watt': 2 parameters
                - 'gauss': 2 parameters
                - 'evap': 1 parameter
                - 'muir': 2 parameters
                - 'exp_decay': 1 parameter
                - 'power': 1 parameter
                - 'exp': 1 parameter
                - 'gauss_beam': 2 parameters

            All the analytic parameters have defaults, and the defaults can be
            invoked if ``vals`` is an empty list. See MCNP manual for more
            information.
        keys : list, str
            If a string, the string is printed after the ``key_setting``, if
            provided. A string is used if the class does not provide the
            desired functionality, and the user wants to provide their own
            input. If a list, elements can be:
            
            - particle names (keys of :py:attr:`Particle.mcnp_abbrev`),
            - 'spont-fiss-by-neut', 'spont-fiss-by-hist', 'spont-phot',
            - a 3-element list or :py:class:`np.array`,
            - the name of another :py:class:`Distribution`, or
            - anything else.

            The list should be as long as ``vals``, unless ``key_setting`` is
            ``analytic``, in which case it is an empty list or the bounds of
            the domain of the analytic function.
        vals : list, str
            If a string, the string is printed after the ``val_setting``, if
            provided. A string is used if the class does not provide the
            desired functionality, and the user wants to provide their own
            input. If a list, elements can be:
            
            - cell name,
            - analytic function parameters, or
            - anything else (a normal bin probability).

            The list should be as long as ``keys``, unless ``key_setting`` is
            'analytic', in which case the list contains the 1 or 2 parameters
            for the analytic function, or can be an empty list to use the
            default parameters.

        Examples
        --------
        The following might be used for the ``x`` source variable::

            sd = Distribution('distA', [-2, 2], [0, 1])
            sd = Distribution('distA', [-2, 2], [0, 1], 'histogram', 'prob')

        The following specifies multiple particle types and their intensities,
        also showing how a cell name can be provided as an intensity value, and
        a heavy ion can be specified::

            sd = Distribution('distA', 
                    ['neutron', 'photon', 'spont-fiss-by-neut', 92238],
                    [1, 1, 'cellA', 2],
                    key_setting='discrete',
                    val_setting='partint')

        The names of other distributions can be provided as well. It is assumed
        in this example that 'distB' and 'distC' are the names of
        :py:class:`Distribution`'s that will be added to the simulation before
        an input is created::

            sd = Distribution('distA', ['distB', 'distC'], [0.3, 0.7], 'dist')

        The following two commands specify an unbounded analytic distribution
        using the default parameters::

            sd = Distribution('distA', [], [], 'analytic', 'maxwell')
            sd = Distribution('distA', [], [], key_setting='analytic',
                    val_setting='maxwell')

        The following is a bounded analytic distribution, with parameters
        specified::

            sd = Distribution('distA', [0, 10], [1, 3], 'analytic', 'watt')

        """
        super(Distribution, self).__init__(name, **kwargs)
        self.key_setting = key_setting
        self.val_setting = val_setting
        self.keys = keys
        self.vals = vals

    def comment(self):
        vecformat = " {0} {0} {0}"
        string = "Source distribution {0!r}:".format(self.name)
        string += " key setting {0},".format(
                self.key_setting if self.key_setting else "default")
        string += " val setting {0},".format(
                self.val_setting if self.val_setting else "default")
        string += " KEYS:"
        # Independent variable.
        if type(self.keys) is str:
            string += " {0}".format(self.keys)
        elif type(self.keys) is list and len(self.keys) == 0:
            string += " default,"
        else:
            for key in self.keys:
                string += " {0},".format(key)
        # Probability values.
        string += " VALS:"
        if type(self.vals) is str:
            string += " {0}".format(self.vals)
        elif type(self.vals) is list and len(self.vals) == 0:
            string += " default,"
        else:
            for val in self.vals:
                string += " {0},".format(val)
        return string[:-1] + "."

    def mcnp(self, float_format, sim):
        vecformat = " {0} {0} {0}"
        # Independent variable.
        info = "SI{0}".format(sim.dist_num(self.name))
        if self.key_setting:
            if self.key_setting != 'analytic':
                info += " {0}".format(self.mcnp_keys[self.key_setting])
        if type(self.keys) is str:
            info += " {0}".format(self.keys)
        else:
            for key in self.keys:
                if key in Particle.mcnp_abbrev:
                    info += " {0}".format(Particle(key).mcnp())
                elif key == 'spont-fiss-by-neut': info += " SF"
                elif key == 'spont-fiss-by-hist': info += " -SF"
                elif key == 'spont-phot':         info += " SP"
                elif type(key) is list or isinstance(key, np.ndarray):
                    info += vecformat % tuple(key)
                elif key in sim.dists:
                    info += " {0}".format(sim.dist_num(key))
                else:
                    info += " {0}".format(key)
        # Probability values.
        prob = "SP{0}".format(sim.dist_num(self.name))
        if self.val_setting:
            if self.key_setting == 'analytic':
                prob += " {0}".format(-self.analytic[self.val_setting])
            else:
                prob += " {0}".format(self.mcnp_vals[self.val_setting])
        if type(self.vals) is str:
            prob += " {0}".format(self.vals)
        else:
            for val in self.vals:
                if val in sim.sys.cells:
                    prob += " {0}".format(-sim.sys.cell_num(val))
                elif type(val) is float:
                    prob += " " + float_format % val
                else:
                    prob += " {0}".format(val)
        # Don't print an SI card if ``keys`` is an empty list.
        if type(self.keys) is list and len(self.keys) == 0:
            string = prob
        else:
            string = "{0}\n{1}".format(info, prob)
        return string

    @property
    def key_setting(self): return self._key_setting

    @key_setting.setter
    def key_setting(self, value):
        acceptable = self.mcnp_keys.keys() + ['analytic']
        if value and value not in acceptable:
            raise ValueError("The ``key_setting`` {0!r} is "
                "unacceptable.".format(value))
        self._key_setting = value

    @property
    def val_setting(self): return self._val_setting

    @val_setting.setter
    def val_setting(self, value):
        acceptable = self.mcnp_vals.keys() + self.analytic.keys()
        if value and value not in acceptable:
            raise ValueError("The ``val_setting`` {0!r} is "
                "unacceptable.".format(value))
        self._val_setting = value

    @property
    def keys(self): return self._keys

    @keys.setter
    def keys(self, value): self._keys = value

    @property
    def vals(self): return self._vals

    @vals.setter
    def vals(self, value): self._vals = value


class Criticality(ISource):
    """A criticality source of neutrons. Unique card with name `criticality`.
    In MCNP, this is the **KCODE** card.
    
    .. inheritance-diagram:: pyne.simplesim.cards.Criticality

    """
    def __init__(self, n_histories=1000, keff_guess=1.0,
            n_skip_cycles=30, n_cycles=130):
        """
        Parameters
        ----------
        n_histories : int, optional
            Number of particle histories to run in each cycle.
        keff_guess : float, optional
            Initial guess for the effective multiplication constant of the
            system.
        n_skip_cycles : int, optional
            The number of cycles to skip.
        n_cycles : int, optional
            The total number of cycles to simulate (skipped + active).

        Examples
        --------
        The following::

            critsrc = Criticality(2000, 1.5, 30, 300)

        creates a criticality source with 2000 histories per cycle, an initial
        k_eff guess of 1.5, 30 skipped cyles, and 300 total cycles.
        
        """
        super(Criticality, self).__init__('criticality', unique=True) 
        self.n_histories = n_histories
        self.keff_guess = keff_guess
        self.n_skip_cycles = n_skip_cycles
        self.n_cycles = n_cycles

    def comment(self):
        return ("Criticality source {0!r}: n_histories: {1}, keff_guess: {2:g}"
                ", n_skip_cycles: {3}, n_cycles: {4}.".format(self.name,
                    self.n_histories, self.keff_guess, self.n_skip_cycles,
                    self.n_cycles))

    def mcnp(self, float_format, sim):
        formatstr = "KCODE {0} {1} {2} {3}".format(self.n_histories,
                float_format, self.n_skip_cycles, self.n_cycles)
        return formatstr % self.keff_guess

    @property
    def n_histories(self): return self._n_histories

    @n_histories.setter
    def n_histories(self, value):
        if type(value) is not int:
            raise ValueError("The property ``n_histories`` must be an "
                    "integer. User provided {0}.".format(value))
        if value <= 0:
            raise ValueError("The property ``n_histories`` must be positive. "
                "User provided {0}.".format(value))
        self._n_histories = value

    @property
    def keff_guess(self): return self._keff_guess
 
    @keff_guess.setter
    def keff_guess(self, value):
        if value < 0:
            raise ValueError("The property ``keff_guess`` must be "
                    "non-negative. User provided {0:g}.".format(value))
        self._keff_guess = value

    @property
    def n_skip_cycles(self): return self._n_skip_cycles

    @n_skip_cycles.setter
    def n_skip_cycles(self, value):
        if type(value) is not int:
            raise ValueError("The property ``n_skip_cycles`` must be an "
                    "integer. User provided {0}.".format(value))
        if value <= 0:
            raise ValueError("The property ``n_skip_cycles`` must be positive. "
                "User provided {0}.".format(value))
        self._n_skip_cycles = value

    @property
    def n_cycles(self): return self._n_cycles

    @n_cycles.setter
    def n_cycles(self, value):
        if type(value) is not int:
            raise ValueError("The property ``n_cycles`` must be an "
                    "integer. User provided {0}.".format(value))
        if value < self.n_skip_cycles:
            raise ValueError("The property ``n_cycles`` must be equal to or "
                    "greater than ``n_skip_cycles``. "
                    "User provided {0}.".format(value))
        # If n_cycles is greater or equal to n_skip_cycles, it is positive.
        self._n_cycles = value


class CriticalityPoints(ISource):
    """Initial source points for neutrons generated by a criticality source.
    Unique card with name 'criticalitypoints'. In MCNP, this is the **KSRC**
    card.

    .. inheritance-diagram:: pyne.simplesim.cards.CriticalityPoints
    
    """
    def __init__(self, points=[[0, 0, 0]]):
        """
        Parameters
        ----------
        points : list of 3-element lists, optional [centimeters]
            A list of 3-element lists (or numpy arrays) specifying
            initial neutron source points in 3-D space.

        Examples
        --------
        The following specifies two initial source points at (1, 2, 3) cm and
        at (3.141..., 2.718..., 0) cm, where we have imported :py:mod:`numpy`
        as ``np``::

            critpts = CriticalityPoints([[1, 2, 3],
                                        np.array([np.pi, np.e, 0])])

        """
        super(CriticalityPoints, self).__init__('criticalitypoints',
                                                unique=True)
        self.points = points

    def comment(self):
        string = "Criticality points {0!r}:".format(self.name)
        counter = 0
        for point in self.points:
            counter += 1
            string += " ({0[0]:g}, {0[1]:g}, {0[2]:g}) cm{1}".format(point,
                "," if counter < len(self.points) else ".")
        return string

    def mcnp(self, float_format, sim):
        # Will be editing the points list through popping, so make a copy.
        points = copy.deepcopy(self.points)
        formatstr = "KSRC {0} {0} {0}".format(float_format)
        string = formatstr % tuple(points.pop(0))
        formatstr = "\n{0}{1} {1} {1}".format(5 * " ", float_format)
        for point in points:
            string += formatstr % tuple(point)
        return string

    @property
    def points(self): return self._points

    @points.setter
    def points(self, value):
        for point in value:
            if len(point) is not 3:
                raise ValueError("Length of all point lists/arrays must be 3."
                        " User provided a point {0}.".format(point))
        self._points = value


class ITally(ICard):
    """This class is not used by the user. Abstract base class for
    tallies (MCNP) or detectors (Serpent).

    """
    __metaclass__ = abc.ABCMeta

    def __init__(self, name, particle, alt_units=False, *args, **kwargs):
        """
        Parameters
        ----------
        name : str
            See :py:class:`ICard`. Used for, e.g., tally multiplier cards.
        particle : str, list of str
            Name of the particle to tally, e.g. 'neutron', 'photon', electron',
            or 'proton'. For some tallies, may be a list of particle names.
        alt_units : bool, optional
            If set to True and the tally can use alternative units, alternative
            units are used for the tally. See subclasses.

        """
        super(ITally, self).__init__(name, *args, **kwargs)
        if type(particle) is list:
            self.particle = [Particle(el) for el in particle]
        else:
            self.particle = Particle(particle)
        self.alt_units = alt_units

    @abc.abstractmethod
    def comment(self, title):
        string = "{0} tally {1!r} of ".format(title, self.name)
        if type(self.particle) is not list:
            string += self.particle.name
            # 'all' is only for CellEnergyDeposition.
            if self.particle.name != 'all': string += "s"
        else:
            pcounter = 0
            for par in self.particle:
                pcounter += 1
                string += "{0}s".format(par)
                if pcounter < len(self.particle): string += ", "
        string += "{0}:".format(
                " (in alt. units)" if self.alt_units else "")
        return string
    
    @abc.abstractmethod
    def mcnp(self, float_format, sim, num, pre='', post=''):
        if pre != '' and self.alt_units:
            raise Exception("Can't have * and + as tally prefix.")
        if self.alt_units: pre = '*'
        string = "{0}F{1}{2}".format(
                pre, 10 * sim.tally_num(self.name) + num, post)
        if type(self.particle) is not list and self.particle.name == 'all':
            pass
        else:
            # 'all' is only for CellEnergyDeposition.
            string += ":"
            if type(self.particle) is list: plist = self.particle
            else:                           plist = [self.particle]
            pcounter = 0
            for part in plist:
                pcounter += 1
                string += part.mcnp()
                if pcounter < len(plist): string += ","
        return string + " "

    @property
    def particle(self): return self._particle

    @particle.setter
    def particle(self, value):
        self._particle = value

    @property
    def alt_units(self): return self._alt_units

    @alt_units.setter
    def alt_units(self, value): self._alt_units = value


class ICellSurfTally(ITally):
    """This class is not used by the user. Abstract base class for
    tallies over cells and surfaces, as opposed to detector tallies. In MCNP,
    these are the **F1**, **F2**, **F4**, **F6**, **F7** and **F8** tallies.
    Some of these are for cells, and some are for surfaces.

    """
    __metaclass__ = abc.ABCMeta

    def __init__(self, name, particle, cards, *args, **kwargs):
        """
        Parameters
        ----------
        name : str
            See :py:class:`ITally`.
        particle : str
            See :py:class:`ITally`.
        cards : str name of :py:class:`Cell` or :py:class:`ISurface`, OR :py:class:`pyne.simplesim.nestedgeom.IUnit`, list, list of lists
            **Basic** If tallying 1 cell/surface, the input is that
            cell/surface card. If tallying multiple cells/surfaces, the
            individual cell/surface cards are provided in a list. To obtain the
            average tally across multiple cells/surfaces, these cell/surface
            cards are provided in their own list, within the outer list. To
            avoid ambiguity, if only one set of averages is desired, then this
            set must be nested in two lists. See the examples.

            **Unit** For a more complex tally, specifically for nested geometry
            or what is called `repeated structures` in MCNP, the user can
            supply an instance of a subclass of
            :py:class:`pyne.simplesim.nestedgeom.IUnit`, or such an instance as
            part of the list of card names. However, a unit cannot be averaged
            with other units (e.g. nested within another list); there is
            functionality within the units for averaging/unioning. There is a
            good amount of documentation in that module, and there are examples
            below.
        alt_units : bool, optional
            See :py:class:`ITally`.

        Examples
        --------
        The following gives the tally in cell 'cA'::

            tally = CellFlux('fuel', 'neutron', 'cA')

        The following two cards give the tally in surface 'sA' and 'sB', and
        the average across surfaces 'sB' and 'sC'::

            tally = SurfaceFlux('fuel', 'photon', ['sA', 'sB', ['sB', 'sC']],
                    average=True)

        To obtain the average across surfaces 'sA' and 'sB' only, a nested list
        is required to distinguish the case of tallying on 'sA' and 'sB'
        individually::

            tally = SurfaceFlux('fuel', 'neutron', [['sA', 'sB']])

        Alternatively, the cards can be specified by something like the
        following (see :py:mod:`pyne.simplesim.nestedgeom`)::

            import pyne.simplesim.nestedgeom as ng
            unit = ng.Surf('sA') < ng.FCell('cA')
            tally = SurfaceFlux('fuel', 'neutron', unit)

        The basic input and the unit input can be mixed, and two different
        units can be requested::

            tally = SurfaceFlux('fuel', 'neutron', ['sA', unit])
            unit2 = ng.Surf('sA') < ng.FCell('cA').lat(ng.Lin(2))
            tally = SurfaceFlux('fuel', 'neutron', [unit, unit2])

        """
        super(ICellSurfTally, self).__init__(name, particle, *args, **kwargs)
        self.cards = cards
        # Affects behavior in comment() and mcnp(), etc.
        self.card_type = kwargs.get('card_type', None)
        if self.card_type != 'cell' and self.card_type != 'surface':
            raise Exception("Card type must be 'cell' or 'surface. "
                    "Got {0}.".format(self.card_type))

    @abc.abstractmethod
    def comment(self, title, union_type):
        string = super(ICellSurfTally, self).comment(title) + " "
        if self.card_type == 'cell':      classcheck = Cell
        elif self.card_type == 'surface': classcheck = ISurface
        if type(self.cards) is not list: 
            if isinstance(self.cards, ng.IUnit):
                string += self.cards.comment()
            elif type(self.cards) is str:
                string += "{0} {1!r}".format(self.card_type, self.cards)
            else:
                raise ValueError("Expected list, string, or IUnit; "
                        "got {0}".format(self.cards))
        elif type(self.cards) is list:
            if type(self.cards[0]) is not list:
                string += "{0}s ".format(self.card_type)
            outcounter = 0
            for obj in self.cards:
                outcounter += 1
                if type(obj) is not list:
                    if isinstance(obj, ng.IUnit):
                        string += obj.comment()
                    elif type(obj) is str:
                        string += "{0!r}".format(obj)
                    else:
                        raise ValueError("Expected list, string, or IUnit; "
                                "got {0}".format(obj))
                elif type(obj) is list:
                    # Not a unit.
                    string += "{0} in ".format(union_type)
                    incounter = 0
                    for avgobj in obj:
                        if type(avgobj) is not str:
                            raise("Expected string; got {0}".format(avgobj))
                        incounter += 1
                        # Must be a cell/surface card name.
                        string += "{0!r}".format(avgobj)
                        if incounter < len(obj): string += ", "
                if outcounter < len(self.cards): string += "; "
        return string

    @abc.abstractmethod
    def mcnp(self, float_format, sim, num, **kwargs):
        string = super(ICellSurfTally, self).mcnp(float_format, sim, num,
                                                  **kwargs)
        if self.card_type == 'cell':      getname = sim.sys.cell_num
        elif self.card_type == 'surface': getname = sim.sys.surface_num
        if type(self.cards) is list: clist = self.cards
        else:                        clist = [self.cards]
        outcounter = 0
        for obj in clist:
            outcounter += 1
            if type(obj) is not list:
                if isinstance(obj, ng.IUnit):
                    string += obj.mcnp(float_format, sim)
                elif type(obj) is str:
                    string += " {0:d}".format(getname(obj))
                else:
                    raise ValueError("Expected string or IUnit; got "
                            "{0}".format(obj))
            elif type(obj) is list:
                string += " ("
                for avgobj in obj:
                    # Must be a cell/surface card name.
                    if type(avgobj) is not str:
                        raise ValueError("Expected string; got "
                                "{0}".format(avgobj))
                    string += " {0:d}".format(getname(avgobj))
                string += ")"
        return string

    def _unique_card_list(self):
        # TODO Broken after implementing IUnit (does not expect card to contain
        # any IUnits) but seems unneeded.
        # Returns a unique list of all the card names provided in self.cards.
        # This method, in the future, may be called by
        # :py:class:`pyne.simplesim.SimulationDefinition` for error-checking.
        # TODO this is ideally recursive, and maybe can be implemented in a
        # cleaner way.
        if type(self.cards) is not list:
            return [self.cards]
        elif type(self.cards) is list:
            card_list = list()
            for obj in self.cards:
                if type(obj) is list:
                    for avgobj in obj:
                        if (type(avgobj) is not list and
                                avgobj not in card_list):
                            card_list += [avgobj]
                elif obj not in card_list:
                    # Not a list, must be a cell/surface name.
                    card_list += [obj]
            return card_list
        else:
            raise ValueError("Expected cell or surface name, or list,"
                " got {0}.".format(type(self.cards)))

    @property
    def cards(self): return self._cards

    @cards.setter
    def cards(self, value): self._cards = value

    @property
    def alt_units(self): return self._alt_units

    @alt_units.setter
    def alt_units(self, value): self._alt_units = value


class SurfaceCurrent(ICellSurfTally):
    """Surface current tally. In MCNP, this is the **F1** card.

    .. inheritance-diagram:: pyne.simplesim.cards.SurfaceCurrent

    """
    mcnp_num = 1

    def __init__(self, name, particle, cards, total=False, **kwargs):
        """
        Parameters
        ----------
        name : str
            See :py:class:`ITally`.
        particle : str
            See :py:class:`ITally`.
        cards : str name of :py:class:`ISurface`, OR :py:class:`pyne.simplesim.nestedgeom.IUnit`, list, list of lists
            See :py:class:`ICellSurfTally`.
        total : bool, optional
            Include a tally for the total current across all surfaces
            specified on this card (note: NOT across all surfaces in the
            `problem`).
        alt_units : bool, optional
            If True, Tally is additionally weighted by particle energy.
            See :py:class:`ITally`.

        Examples
        --------
        The following requests the tally in surface 'sA', surface 'sB', as well
        as the total across 'sA' and 'sB'::

            tally = SurfaceCurrent('fuel', 'electron', ['sA', 'sB'],
                    total=True)

        In the following, the tally is also weighted by particle energy::

            tally = SurfaceCurrent('fuel', 'photon', [['sA', 'sB']],
                    alt_units=True)

        """
        super(SurfaceCurrent, self).__init__(name, particle, cards,
                                             card_type='surface', **kwargs)
        self.total = total

    def comment(self):
        string = super(SurfaceCurrent, self).comment("Surface current",
                                                     'total')
        if self.total: string += "; and total of all provided."
        else:          string += "."
        return string

    def mcnp(self, float_format, sim):
        string = super(SurfaceCurrent, self).mcnp(float_format, sim,
                self.mcnp_num)
        if self.total: string += " T"
        return string


class IAverageTally(ICellSurfTally):
    """This class is not used by the user. Abstract base class for
    tallies of averaged quantities. In MCNP, these are the **F2**, **F4**,
    **F6**, **F7** and **F8** tallies. Some of these are for cells, and some
    are for surfaces.

    """
    __metaclass__ = abc.ABCMeta

    def __init__(self, name, particle, cards, average=False, **kwargs):
        """
        Parameters
        ----------
        name : str
            See :py:class:`ITally`.
        particle : str
            See :py:class:`ITally`.
        cards : str name of :py:class:`Cell` or :py:class:`ISurface`, OR :py:class:`pyne.simplesim.nestedgeom.IUnit`, list, list of lists
            See :py:class:`ICellSurfTally`.
        average : bool, optional
            Include a tally for the average flux across all cells/surfaces
            specified on this card (note: NOT across all cells in the
            `problem`).
        alt_units : bool, optional
            If set to True and the tally can use alternative units, alternative
            units are used for the tally. See subclasses, and :py:class:`ITally`.

        Examples
        --------
        The following requests the tally in cell 'cA', cell 'cB', as well as
        the average across 'cA' and 'cB'::

            tally = CellEnergyDeposition('fuel', 'neutron', ['cA', 'cB'],
                    average=True)

        """
        super(IAverageTally, self).__init__(name, particle, cards, **kwargs)
        self.average = average

    @abc.abstractmethod
    def comment(self, title):
        avgstr = 'avg.'
        string = super(IAverageTally, self).comment(title, avgstr)
        if self.average: string += "; and {0} of all provided.".format(avgstr)
        else:            string += "."
        return string

    def mcnp(self, float_format, sim, num, **kwargs):
        string = super(IAverageTally, self).mcnp(float_format, sim, num,
                                                 **kwargs)
        if self.average: string += " T"
        return string

    @property
    def average(self): return self._average

    @average.setter
    def average(self, value): self._average = value


class SurfaceFlux(IAverageTally):
    """Surface flux tally. In MCNP, this is the **F2** card.

    .. inheritance-diagram:: pyne.simplesim.cards.SurfaceFlux

    """
    mcnp_num = 2

    def __init__(self, name, particle, cards, **kwargs):
        """
        Parameters
        ----------
        name : str
            See :py:class:`ITally`.
        particle : str
            See :py:class:`ITally`.
        cards : str name of :py:class:`ISurface`, OR :py:class:`pyne.simplesim.nestedgeom.IUnit`, list, list of lists
            See :py:class:`IAverageTally`.
        average : bool, optional
            See :py:class:`IAverageTally`.
        alt_units : bool, optional
            If True, Tally is additionally weighted by particle energy.
            See :py:class:`ITally`.

        Examples
        --------
        The following requests the tally in surface 'sA', surface 'sB', as well
        as the average across 'sA' and 'sB'::

            tally = SurfaceFlux('fuel', 'electron', ['sA', 'sB'],
                    average=True)

        In the following, the tally is also weighted by particle energy::

            tally = SurfaceFlux('fuel', 'proton', [['sA', 'sB']],
                    alt_units=True)
        
        See base classes for more examples.

        """
        super(SurfaceFlux, self).__init__(name, particle, cards,
                                          card_type='surface', **kwargs)

    def comment(self):
        return super(SurfaceFlux, self).comment("Surface flux")

    def mcnp(self, float_format, sim):
        return super(SurfaceFlux, self).mcnp(float_format, sim, self.mcnp_num)

    @property
    def total(self): return self._total

    @total.setter
    def total(self, value): self._total = value


class CellFlux(IAverageTally):
    """Cell flux tally. In MCNP, this is the **F4** card.

    .. inheritance-diagram:: pyne.simplesim.cards.CellFlux

    """
    mcnp_num = 4

    def __init__(self, name, particle, cards, **kwargs):
        """
        Parameters
        ----------
        name : str
            See :py:class:`ITally`.
        particle : str
            See :py:class:`ITally`.
        cards : str name of :py:class:`Cell`, OR :py:class:`pyne.simplesim.nestedgeom.IUnit`, list, list of lists
            See :py:class:`IAverageTally`.
        average : bool, optional
            See :py:class:`IAverageTally`.
        alt_units : bool, optional
            If True, Tally is additionally weighted by particle energy.
            See :py:class:`ITally`.

        Examples
        --------
        The following requests the tally in cell 'cA', cell 'cB', as well as
        the average across 'cA' and 'cB'::

            tally = CellFlux('fuel', 'electron', ['cA', 'cB'],
                    average=True)

        In the following, the tally is also weighted by particle energy::

            tally = CellFlux('fuel', 'proton', [['cA', 'cB']],
                    alt_units=True)
        
        See base classes for more examples.

        """
        super(CellFlux, self).__init__(name, particle, cards, card_type='cell',
                                       **kwargs)

    def comment(self):
        return super(CellFlux, self).comment("Cell flux")

    def mcnp(self, float_format, sim):
        return super(CellFlux, self).mcnp(float_format, sim, self.mcnp_num)


class CellEnergyDeposition(IAverageTally):
    """Energy deposition tally. In MCNP, this is the **F6** card. In MCNP, it
    is not permitted to use a particle 'all' and also to use alternative units.

    .. inheritance-diagram:: pyne.simplesim.cards.CellEnergyDeposition

    """
    mcnp_num = 6

    # TODO in mcnp input, prevent particle all and alt_units
    def __init__(self, name, particles, cards, **kwargs):
        """
        Parameters
        ----------
        name : str
            See :py:class:`ITally`.
        particles : str, list of str
            See :py:class:`ITally`. For this tally, the user can specify the
            particle type as a list of strs to tally more than one type of
            particle. Also, the additional value of 'all' is allowed, and
            specifies collision heating. As may be expected, 'all' cannot be
            provided as part of a list.
        cards : str name of :py:class:`Cell`, OR :py:class:`pyne.simplesim.nestedgeom.IUnit`, list, list of lists
            See :py:class:`IAverageTally`.
        average : bool, optional
            See :py:class:`IAverageTally`.
        alt_units : bool, optional
            If True, alternative units are used for the tally. In MCNP, the
            default units are MeV/g and the alternative units are jerks/g.
            See :py:class:`ITally`.

        Examples
        --------
        The following requests the energy deposited by neutrons in cell 'cA'::

            tally = CellEnergyDeposition('energy', 'neutron', 'cA')


        The following requests the energy deposited by neutrons and protons in
        cell 'cA'::

            tally = CellEnergyDeposition('energy', ['neutron', 'proton'],
                    'cA')
        
        The following requests the energy deposited by all particles in cell
        'cA'::

            tally = CellEnergyDeposition('energy', 'all', 'cA')

        The following are not allowed in MCNP::
            
            tally = CellEnergyDeposition('energy', ['neutron', 'all'], 'cA')
            tally = CellEnergyDeposition('energy', 'all', 'cA', alt_units=True)

        See base classes for more examples.

        """
        super(CellEnergyDeposition, self).__init__(name, particles, cards,
                                                   card_type='cell', **kwargs)
        # TODO move this error check to the MCNP method.
        if (type(self.particle) is not list and self.particle.name == 'all' and
                self.alt_units):
            raise ValueError("The particle cannot be 'all' if alt_units is "
                    "True.")

    def comment(self):
        return super(CellEnergyDeposition, self).comment("Energy deposition")

    def mcnp(self, float_format, sim):
        if type(self.particle) is not list and self.particle.name == 'all':
            pre = '+'
        else:
            pre = ''
        return super(CellEnergyDeposition, self).mcnp(float_format, sim,
                                                      self.mcnp_num,
                                                      pre=pre)

    @property
    def particle(self): return self._particle

    @particle.setter
    def particle(self, value):
        # Make sure that if there is more than one particle, that none of them
        # has name 'all'.
        if type(value) is list:
            for par in value:
                if par.name == 'all':
                    raise ValueError("The particle 'all' cannot be given with"
                            " other particles.")
        self._particle = value


class CellFissionEnergyDeposition(IAverageTally):        
    """Fission energy deposition tally. In MCNP, this is the **F7** card. The
    particle is necessarily neutron.

    .. inheritance-diagram:: pyne.simplesim.cards.CellFissionEnergyDeposition

    """
    mcnp_num = 7
    # TODO prevent user from specifying a different particle.
    def __init__(self, name, cards, **kwargs):
        """
        Parameters
        ----------
        name : str
            See :py:class:`ITally`.
        cards : str name of :py:class:`Cell`, OR :py:class:`pyne.simplesim.nestedgeom.IUnit`, list, list of lists
            See :py:class:`IAverageTally`.
        average : bool, optional
            See :py:class:`IAverageTally`.
        alt_units : bool, optional
            If True, alternative units are used for the tally. In MCNP, the
            default units are MeV/g and the alternative units are jerks/g.
            See :py:class:`ITally`.

        Examples
        --------
        The following requests the tally in cell 'cA', cell 'cB', as well as
        the average across 'cA' and 'cB'::

            tally = CellFissionEnergyDeposition('fuel', ['cA',
                    'cB'], average=True)

        In the following, the alternate units are used::

            tally = CellFissionEnergyDeposition('fuel', [['cA',
                    'cB']], alt_units=True)
        
        See base classes for more examples.

        """
        super(CellFissionEnergyDeposition, self).__init__(name, 'neutron',
                cards, card_type='cell', **kwargs)

    def comment(self):
        return super(CellFissionEnergyDeposition, self).comment(
                "Fission energy deposition")

    def mcnp(self, float_format, sim):
        return super(CellFissionEnergyDeposition, self).mcnp(
                float_format, sim, self.mcnp_num)


class CellPulseHeight(IAverageTally):
    """Pulse height tally in cells. In MCNP, this is the **F8** card. For a
    charge deposition tally, see :py:class:`CellChargeDeposition`.

    .. inheritance-diagram:: pyne.simplesim.cards.CellPulseHeight

    """
    mcnp_num = 8

    def __init__(self, name, particles, cards, **kwargs):
        """
        Parameters
        ----------
        name : str
            See :py:class:`ITally`.
        particles : str, list of str
            See :py:class:`ITally`. Multiple particles can be provided in a
            list of str. In MCNP, if only 'proton', or 'electron' is
            specified, both are automatically included.
        cards : str name of :py:class:`Cell`, OR :py:class:`pyne.simplesim.nestedgeom.IUnit`, list, list of lists
            See :py:class:`IAverageTally`.
        average : bool, optional
            See :py:class:`IAverageTally`.
        alt_units : bool, optional
            If True, alternative units are used for the tally. In MCNP, the
            default units are pulses and the alternative units are MeV.
            See :py:class:`ITally`.

        Examples
        --------
        The following requests the tally in cell 'cA' and cell 'cB' for both
        protons and electrons, and requests units of MeV::

            tally = CellPulseHeight('fuel', ['proton', 'electron'], ['cA',
                    'cB'], alt_units=True)

        See base classes for more examples.

        """
        super(CellPulseHeight, self).__init__(name, particles, cards,
                                              card_type='cell', **kwargs)

    def comment(self):
        return super(CellPulseHeight, self).comment("Pulse height")

    def mcnp(self, float_format, sim, **kwargs):
        return super(CellPulseHeight, self).mcnp(float_format, sim,
                                                 self.mcnp_num,
                                                 **kwargs)

    @property
    def particle(self): return self._particle

    @particle.setter
    def particle(self, value):
        # Copied from CellEnergyDeposition.particle
        self._particle = value


class CellChargeDeposition(CellPulseHeight):
    """Charge deposition tally in cells. In MCNP, this is the **+F8** card. No
    alternative units are available.

    .. inheritance-diagram:: pyne.simplesim.cards.CellChargeDeposition

    """
    # TODO it doesn't make sense that the user can provide the particle type
    # here, at least for MCNP.
    def __init__(self, name, particles, cards, **kwargs):
        """
        Parameters
        ----------
        name : str
            See :py:class:`ITally`.
        particles : str, list of str
            See :py:class:`ITally`. Multiple particles can be provided in a
            list of str. In MCNP, if only 'proton', or 'electron' is
            specified, both are automatically included.
        cards : str name of :py:class:`Cell`, OR :py:class:`pyne.simplesim.nestedgeom.IUnit`, list, list of lists
            See :py:class:`IAverageTally`.
        average : bool, optional
            See :py:class:`IAverageTally`.

        Examples
        --------
        The following requests the tally in cell 'cA' and cell 'cB' for both
        protons and electrons::

            tally = CellChargeDeposition('fuel', ['proton', 'electron'], 
                    ['cA', 'cB'])

        See base classes for more examples.

        """
        super(CellChargeDeposition, self).__init__(name, particles, cards,
                                                   **kwargs)

    def comment(self):
        return super(CellChargeDeposition, self).comment("Charge deposition")

    def mcnp(self, float_format, sim):
        return super(CellChargeDeposition, self).mcnp(float_format, sim,
                                                      pre='+')


class IDetector(ITally):
    """This class is not used by the detector. Abstract base class for detector
    tallies.
    
    """
    mcnp_num = 5

    def __init__(self, name, particle, args_per_set, *args, **kwargs):
        """
        Parameters
        ----------
        particles : str, list of str
            See :py:class:`ITally`. In MCNP, this tally is only for neutrons and
            photons.
        sep_direct : bool, optional (default: True)
            In MCNP, the direct contribution to the tally is printed
            separately. Set to False to disable the separate printing. 

        """
        super(IDetector, self).__init__(name, particle, **kwargs)
        self.detectors = []
        self.sep_direct = kwargs.get('sep_direct', True)
        self._args_per_set = args_per_set
        for i_set in range(len(args) / self._args_per_set):
            i_start = self._args_per_set * i_set
            self.add(*args[i_start:i_start+self._args_per_set])

    @abc.abstractmethod
    def comment(self, name, post=None):
        string = super(IDetector, self).comment(name)
        if post: string += "{0}".format(post)
        counter = 0
        for det in self.detectors:
            counter += 1
            string += self._comment_unit(det)
            if counter < len(self.detectors): string += ";"
        string += "; direct contrib is {0}separate.".format(
                '' if self.sep_direct else 'not ')
        return string
    
    @abc.abstractmethod
    def _comment_unit(self):
        raise NotImplementedError

    @abc.abstractmethod
    def mcnp(self, float_format, sim, num, **kwargs):
        string = super(IDetector, self).mcnp(float_format, sim, num, **kwargs)
                
        for det in self.detectors:
            if det is self.detectors[0]: string += 1 * " "
            else:                        string += "\n{0}".format(5 * " ")
            string += self._mcnp_unit(float_format, sim, det)
        if not self.sep_direct: string += " ND"
        return string
    
    def _mcnp_unit(self):
        raise NotImplementedError

    @property
    def detectors(self): return self._detectors

    @detectors.setter
    def detectors(self, value): self._detectors = value

    @property
    def sep_direct(self): return self._sep_direct

    @sep_direct.setter
    def sep_direct(self, value): self._sep_direct = value


class PointDetector(IDetector):
    """A point detector tally. In MCNP, this is the **F5** card. This is not to
    be confused with the more general use of the term `Detector` in Serpent.

    .. inheritance-diagram:: pyne.simplesim.cards.PointDetector

    """
    Point = collections.namedtuple('Point', ['pos', 'soe_rad', 'soe_units'])
    """Named tuple with which input entries are stored."""

    def __init__(self, name, particles, *args, **kwargs):
        """\

        Parameters
        ----------
        name : str
            See :py:class:`ITally`.
        particles : str, list of str
            See :py:class:`ITally`. In MCNP, this tally is only for neutrons and
            photons.
        position : 3-element list, or :py:class:`np.array`
            Location of the point detector.
        soe_radius : float
            Radius of the sphere of exclusion.
        soe_units : str
            Units for the radius of the sphere of exclusion. Either 'cm' or
            'mfp' for mean free paths.
        *args : position, soe_radius, soe_units, ...
            To request more than one detector, supply the last 3 arguments for
            the additional detectors. See examples. This can also be done using
            :py:meth:`add`.
        sep_direct : bool, optional
            See :py:class:`IDetector`.
        alt_units : bool, optional
            If True, Tally is additionally weighted by particle energy.
            See :py:class:`ITally`.

        Examples
        --------
        The following creates a single point detector at the origin, without a
        sphere of exclusion::

            det = PointDetector('point', 'neutron', [0, 0, 0], 0, 'cm')

        The following creates a detector at (1, 1, 1) cm with a sphere of
        exclusion with a radius of 1 cm, where we have imported :py:mod:`numpy`
        as ``np``)::

            det = PointDetector('point', 'neutron', np.array([1, 1, 1]), 1,
                    'cm')

        The radius for the sphere of exclusion here is 3 mfp::

            det = PointDetector('point', 'neutron', [1, 0, 0], 3, 'mfp')

        This is an example of requesting two point detectors::
        
            det = PointDetector('point', 'photon', [0, 0, 0],  0, 'cm',
                                                   [1, 0, 0], 3, 'mfp')

        Here, it is requested that the direct contribution is not tallied
        separately and that alternative units are used for the output::

            det = PointDetector('point', 'photon', [0, 0, 0], 0, 'cm',
                    sep_direct=False, alt_units=True)

        """
        super(PointDetector, self).__init__(name, particles, 3, *args,
                                            **kwargs)

    def add(self, position, soe_radius, soe_units):
        """
        Parameters
        ----------
        position : 3-element list, or :py:class:`np.array`
        soe_radius : float
        soe_units : str

        Examples
        --------
        The last example above can also be achieved by the following::

            det = PointDetector('point', 'photon', sep_direct=False,
                    alt_units=True)
            det.add([0, 0, 0], 0, 'cm')

        """
        if soe_units != 'cm' and soe_units != 'mfp':
            return ValueError("The input ``soe_units`` must be 'cm' or "
                    "'mfp'. User provided {0}.".format(soe_units))
        self.detectors += [self.Point(position, soe_radius, soe_units)]
        
    def comment(self):
        return super(PointDetector, self).comment("Point detector")

    def _comment_unit(self, det):
        return (" point ({0[0]:g}, {0[1]:g}, {0[2]:g}) cm, radius "
                "{1:g} {2}".format(det.pos, det.soe_rad, det.soe_units))

    def mcnp(self, float_format, sim):
        return super(PointDetector, self).mcnp(float_format, sim,
                                               self.mcnp_num)

    def _mcnp_unit(self, float_format, sim, det):
        formatstr = "{0} {0} {0} {0}".format(float_format)
        return formatstr % (tuple(det.pos) + (det.soe_rad * (
            -1 if det.soe_units == 'mfp' else 1),))


class RingDetector(IDetector):
    """A ring detector tally. In MCNP, this is the **F5a** card. This is not to
    be confused with the more general use of the term `Detector` in Serpent.

    .. inheritance-diagram:: pyne.simplesim.cards.RingDetector

    """
    Ring = collections.namedtuple('Ring',
                ['pos', 'rad', 'soe_rad', 'soe_units'])
    """Named tuple with which input entries are stored."""
    
    def __init__(self, name, particles, cartesian_axis, *args, **kwargs):
        """\

        Parameters
        ----------
        name : str
            See :py:class:`ITally`.
        particles : str, list of str
            See :py:class:`ITally`. In MCNP for this tally, only neutrons and
            photons are allowed.
        cartesian_axis : str
            Cartesian axis string ('x', 'X', 'y', 'Y', 'z', or 'Z').
        position : float
            Position along the cartesian axis.
        radius : float
            Radius of the ring.
        soe_radius : float
            Radius of the sphere of exclusion.
        soe_units : float
            Units for the radius of the sphere of exclusion.
        *args : position, radius, soe_radius, soe_units, ...
            To request more than one detector, supply the last 4 arguments for
            the additional detectors. See examples. This can also be done using
            :py:meth:`add`.
        sep_direct : bool, optional
            See :py:class:`IDetector`.
        alt_units : bool, optional
            If True, Tally is additionally weighted by particle energy.
            See :py:class:`ITally`.

        Examples
        --------
        The following creates a single ring detector at x = 10.0 cm, with a
        2.0 cm radius, and a 1.0 cm radius sphere of exclusion::

            det = RingDetector('ring', 'neutron', 'x', 10.0, 2.0,  1.0, 'cm')

        In the following, the sphere of exclusion has a radius of 1.0 mfp::

            det = RingDetector('ring', 'neutron', 'x', 10.0, 2.0, 1.0, 'mfp')

        This is an example of requesting two ring detectors::

            det = RingDetector('ring', 'neutron', 'y', 10.0, 2.0, 1.0, 'mfp',
                                                       20.0, 3.0, 1.0, 'cm')

        Here it is requested that the direct contribution is not tallied
        separately, and alternative units are used for the output::
        
            det = RingDetector('ring', 'neutron', 'x', 10.0, 2.0, 1.0, 'cm'
                    sep_direct=False, alt_units=True)

        """
        super(RingDetector, self).__init__(name, particles, 4, *args, 
                                           **kwargs)
        self.cartesian_axis = cartesian_axis

    def add(self, position, radius, soe_radius, soe_units):
        """\

        Parameters
        ----------
        position : float
        radius : float
        soe_radius : float
        soe_units : float

        Examples
        --------
        The last example above can be achieved with the following::

            det = RingDetector('ring', 'neutron', 'x', sep_direct=False,
                    alt_units=True)
            det.add(10.0, 2.0, 1.0, 'cm')

        """
        if soe_units != 'cm' and soe_units != 'mfp':
            return ValueError("The input ``soe_units`` must be 'cm' or "
                    "'mfp'. User provided {0}.".format(soe_units))
        self.detectors += [self.Ring(position, radius, soe_radius, soe_units)]

    def comment(self):
        return super(RingDetector, self).comment("Ring detector",
                " along {0} axis.".format(self.cartesian_axis))

    def _comment_unit(self, det):
        return (" ring {0} = {1:g} cm, radius {2:g} cm, s.o.e. "
                "radius {3:g} {4}".format(self.cartesian_axis, det.pos,
                    det.rad, det.soe_rad, det.soe_units))

    def mcnp(self, float_format, sim):
        return super(RingDetector, self).mcnp(float_format, sim,
                self.mcnp_num,
                post=self.cartesian_axis.upper())

    def _mcnp_unit(self, float_format, sim, det):
        formatstr = "{0} {0} {0}".format(float_format)
        return formatstr % (det.pos, det.rad, det.soe_rad * (
            -1 if det.soe_units == 'mfp' else 1))

    @property
    def cartesian_axis(self): return self._cartesian_axis

    @cartesian_axis.setter
    def cartesian_axis(self, value): self._cartesian_axis = value


class EnergyGrid(IMisc):
    """Energy grid for tallies. In MCNP, this is the **E** card.

    .. inheritance-diagram:: pyne.simplesim.cards.EnergyGrid

    """
    # TODO make this a unique-by-tally card.
    def __init__(self, name, tally, energies):
        """\

        Parameters
        ----------
        name : str
            See :py:class:`ICard`.
        tally : str name of :py:class:`ITally`, None
            The name of the tally for which this energy grid should apply. If
            requesting for this grid to apply to all tallies, then this is
            None.
        energies : list, :py:class:`np.array`
            The upper bounds of the energy groups.

        Examples
        --------
        The following energy grid applies to all tallies::

            egrid = EnergyGrid('grid0', None, [1e-4, 1, 100e3, 10e6])

        The following applies to tally 'tallyA'::

            egrid = EnergyGrid('grid0', 'tallyA',
                    np.array([1e-4, 1, 100e3, 10e6]))

        """
        super(EnergyGrid, self).__init__(name)
        self.tally = tally
        self.energies = energies

    def comment(self):
        string = "Energy grid {0!r} for ".format(self.name)
        if self.tally is None: string += "all tallies"
        else:                  string += "tally {0}".format(self.tally)
        return string + ": {0} groups.".format(len(self.energies))

    def mcnp(self, float_format, sim):
        string = "E{0}".format(
                sim.tally_num(self.tally) if self.tally else "0")
        for val in self.energies:
            string += (" " + float_format) % val
        return string

    @property
    def energies(self): return self._energies

    @energies.setter
    def energies(self, value): self._energies = value


class Comment(IMisc):
    """Unimplemented."""
    pass


class Mode(IMisc):
    """Unimplemented."""
    pass


class NeutronPhysics(IMisc):
    """Unimplemented."""
    pass


class PhotonPhysics(IMisc):
    """Unimplemented."""
    pass


class ElectronPhysics(IMisc):
    """Unimplemented."""
    pass


class ProtonPhysics(IMisc):
    """Unimplemented."""
    pass


class Transformation(ICard):
    """A coordinate transformation. In MCNP, this is the **TR** card. MCNP
    allows the specification of less than 9 values in the rotation matrix; this
    is not supported here. Note that if generating an input for MCNP, this card
    needs to be added to the simulation by
    :py:meth:`pyne.simplesim.definition.MCNPSimulation.add_transformation`.

    .. inheritance-diagram:: pyne.simplesim.cards.Transformation

    """
    # TODO support for MCNP's less-than-9-element transformation matrices.
    def __init__(self, name, displacement, rotation, aux_in_main=True,
                 degrees=False, **kwargs):
        """\

        Parameters
        ----------
        name : str
            See :py:class:`ICard`.
        displacement : 3-element list, :py:class:`np.array`, float [cm]
            Displacement vector.
        rotation : 3 x 3 list, :py:class:`np.array`, :py:class:`np.matrix`
            Rotation matrix.
        aux_in_main : bool, optional
            The displacement vector is in the main coordinate system and points
            to the origin of the transformed/auxiliary coordinate system. If
            False, then the displacement vector is in the
            transformed/auxiliary coordiante system and points to the origin of
            the main coordinate system.
        degrees : bool, optional
            If True, then the elements of the rotation matrix are in degrees.

        Examples
        --------
        Here's a number of examples::

            trans = Transformation('trans1', [1, 2, 3],
                    [[4, 5, 6], [7, 8, 9], [10, 11, 12]])
            trans = Transformation('trans1', [1, 2, 3],
                    np.array([[4, 5, 6], [7, 8, 9], [10, 11, 12]]))
            trans = Transformation('trans1', [1, 2, 3],
                    np.matrix([[4, 5, 6], [7, 8, 9], [10, 11, 12]]))
            trans = Transformation('trans1', [1, 2, 3],
                    [[4, 5, 6], [7, 8, 9], [10, 11, 12]],
                    aux_in_main=False)
            trans = Transformation('trans1', [1, 2, 3],
                    [[4, 5, 6], [7, 8, 9], [10, 11, 12]],
                    degrees=True)

        """
        super(Transformation, self).__init__(name, **kwargs)
        self.displacement = displacement
        # For numpy matrices.
        if hasattr(rotation, 'tolist'): self.rotation = rotation.tolist()
        else:                           self.rotation = rotation
        self.aux_in_main = aux_in_main
        self.degrees = degrees

    def comment(self):
        return "Transformation {0!r}: pos. of {1}".format(self.name,
                self._comment_unit())

    def _comment_unit(self):
        if self.aux_in_main: string = "aux origin in main"
        else:                string = "main origin in aux"
        string += " ({0[0]:g}, {0[1]:g}, {0[2]:g}) cm,".format(
                self.displacement)
        dirs = ['x', 'y', 'z']
        for idx in range(3):
            string += " {0}' <{1:g}, {2:g}, {3:g}>{4}".format(dirs[idx],
                    self.rotation[0][idx],
                    self.rotation[1][idx], 
                    self.rotation[2][idx],
                    " deg" if self.degrees else "")
            if idx != 2: string += ","
        return string + "."

    def mcnp(self, float_format, sim):
        string = "{0}TR{1}".format("*" if self.degrees else "",
                sim.transformation_num(self.name))
        string += self._mcnp_unit(float_format)
        return string

    def _mcnp_unit(self, float_format): 
        # Needed by CellMCNP.
        formatstr = " {0} {0} {0}".format(float_format)
        string = formatstr % tuple(self.displacement)
        for i_row in range(3):
            for i_col in range(3):
                string += " "
                string += float_format % self.rotation[i_row][i_col]
        return string + (" 1" if self.aux_in_main else " -1")

    @property
    def displacement(self): return self._displacement
    
    @displacement.setter
    def displacement(self, value): self._displacement = value

    @property
    def rotation(self): return self._rotation
    
    @rotation.setter
    def rotation(self, value): self._rotation = value

    @property
    def aux_in_main(self): return self._aux_in_main
    
    @aux_in_main.setter
    def aux_in_main(self, value): self._aux_in_main = value


class ICellMod(IMisc):
    """This class is not used by the user. Abstract base class for cards that
    can be specified in MCNP on both the cell card or in the data block.  All
    subclasses are unique, have a ``cell`` property, and have a similar form.
    Entries for a given cell can be modified by providing an input for the same
    cell.

    """
    # The mcnp() method of this class implements the jump feature, which
    # ICellModParticle does not.
    __metaclass__ = abc.ABCMeta

    def __init__(self, name, n_args_per_cell, *args, **kwargs):
        """\

        Parameters
        ----------
        name : str
            See :py:class:`ICard`.
        n_args_per_cell : int
            The number of arguments the subclass expects per cell.

        """
        super(ICellMod, self).__init__(name, unique=True)
        self.cells = []
        self._n_args_per_cell = n_args_per_cell
        if len(args) % n_args_per_cell != 0:
            raise StandardError("The length of ``*args`` must be a multiple "
                    "of {0}. Length is {1}.".format(n_args_per_cell, len(args)))

    def _process_varargs(self, args):
        for i_cell in range(len(args) / self._n_args_per_cell):
            i_start = self._n_args_per_cell * i_cell
            self.set(*args[i_start:i_start+self._n_args_per_cell])

    @abc.abstractmethod
    def set(self, cell):
        """Add an entry or modify an existing one. See subclasses."""

        if cell not in self.cells: self.cells += [cell]

    @abc.abstractmethod
    def comment(self, title, poststr=""):
        string = "{0} {1!r}:{2}".format(title, self.name, poststr)
        counter = 0
        for cell in self.cells:
            counter += 1
            string += " cell {0!r}".format(cell)
            string += self._comment_unit(cell)
            if counter < len(self.cells): string += ","
        return string + "."

    @abc.abstractmethod
    def _comment_unit(self):
        raise NotImplementedError

    @abc.abstractmethod
    def mcnp(self, float_format, sim, keystring):
        # TODO this ordering might not be correct, particularly once we add
        # support for universes, etc.
        string = "{0}".format(keystring)
        # TODO this should loop through in the print order.
        #
        # `gathering` empty cells (cells for which no value is given).
        empty_count = 0
        for key in sim.sys.cells: 
            if key in self.cells:
                if empty_count > 0:
                    string += " {0}J".format(empty_count)
                    empty_count = 0
                string += " " + self._mcnp_unit(float_format, sim, key)
            else:
                empty_count += 1
                # Account for running out of cells with empty baggage.
                if key is sim.sys.cells.keys()[-1]:
                    string += " {0}J".format(empty_count)
        return string

    def _mcnp_unit(self):
        raise NotImplementedError

    @property
    def cells(self): return self._cells

    @cells.setter
    def cells(self, value):
        for val in value:
            if type(val) is not str:
                raise Exception("Expected cell name, got {0}".format(val))
        self._cells = value


class Universes(ICellMod):
    """Universe names for each cell. Unique card with name `universes`. In
    MCNP, this is the **U** card. The user can initialize this card without
    providing any universe names.

    .. inheritance-diagram:: pyne.simplesim.cards.Universes

    """
    def __init__(self, *args):
        """
        Parameters
        ----------
        cell : str name of :py:class:`Cell` or subclass
            The name of the cell for which a universe name is being provided.
        univ_name : str
            The name of the universe that the cell is a part of.
        truncate : bool
            By default, the cell is truncated by the boundary of higher-level
            cells. Not an optional argument, but is typically True.
        *args : cell, univ_name, truncate, ...
            To provide universe names for more than one cell, supply the last
            three arguments for the other cells. See example. This can also be
            done using :py:meth:`set`.

        Examples
        --------
        The following adds the cell cards with names 'pincell' and
        'coolantcell' to the universe 'unitcell'. All other cell cards, as of
        yet, are then part of the real world universe::

            uni = Universes('pincell', 'unitcell', True,
                            'coolantcell', 'unitcell', True)

        """
        super(Universes, self).__init__('universes', 3, *args)
        self.univ_names = collections.OrderedDict()
        self.truncates = collections.OrderedDict()
        self._process_varargs(args)

    def set(self, cell, univ_name, truncate):
        """Add an entry or modify an existing one.

        Parameters
        ----------
        cell : str name of :py:class:`Cell` or subclass
        univ_name : str
        truncate : bool

        Examples
        --------
        The example above can be achieved by the following::

            uni = Universes()
            uni.set('pincell', 'unitcell', True)
            uni.set('coolantcell', 'unitcell', True)

        """
        super(Universes, self).set(cell)
        self.univ_names[cell] = univ_name
        self.truncates[cell] = truncate

    def comment(self):
        return super(Universes, self).comment("Universes")

    def _comment_unit(self, cell):
        return " {0} ({1}truncated)".format(self.univ_names[cell],
                "" if self.truncates[cell] else "not ")

    def mcnp(self, float_format, sim):
        return super(Universes, self).mcnp(float_format, sim, "U")

    def _mcnp_unit(self, float_format, sim, cell):
        if self.univ_names[cell] not in sim.sys.universes:
            raise Exception("This card has not been added to the simulation.")
        return "{0}".format(sim.sys.universe_num(self.univ_names[cell]) * (
            1 if self.truncates[cell] else -1))

    @property
    def univ_names(self): return self._univ_names

    @univ_names.setter
    def univ_names(self, value): self._univ_names = value

    @property
    def truncates(self): return self._truncates

    @truncates.setter
    def truncates(self, value): self._truncates = value


class Fill(ICellMod):
    """Allows the user to fill certain cells with a universe or lattice. Unique
    card with name `fill`. In MCNP, this is the **FILL** card. The user can
    initialize this card without providing any inputs. 
    

    .. inheritance-diagram:: pyne.simplesim.cards.Fill

    Limitations:

    - Combining transformations with the fill is not supported.

    """
    # The following is no longer true:
    #The :py:class:`Universe`
    #card, or keywords on the :py:class:`Cell` card, must be added to the
    #simulation before this card. Otherwise, the code will not be able to find
    #the universes specified on this card.
    def __init__(self, *args):
        """
        Parameters
        ----------
        cell : str name of :py:class:`Cell` or subclass
            The name of the cell being filled by a universe.
        univ_name : str
            The name of the universe filling this cell.
        *args : cell, univ_name, ...
            To fill more than one cell, supply the last 2 arguments for the
            other cells. See examples. This can also be done using
            :py:meth:`set`.

        Examples
        --------
        If 'unit' and 'cellB' are names of cells and 'unitcell' and 'otheruniv`
        are the names of universes, the user can do the following::

            fill = Fill('unit', 'unitcell', 'cellB', 'otheruniv')

        """
        super(Fill, self).__init__('fill', 2, *args)
        self.univ_names = collections.OrderedDict()
        self._process_varargs(args)

    def set(self, cell, univ_name):
        """Add an entry or modify an existing one.

        Parameters
        ----------
        cell : str name of :py:class:`Cell` or subclass
        univ_name : str

        Examples
        --------
        The example above can be achieved by the following::

            fill = Fill()
            fill.set('unit', 'unitcell')
            fill.set('cellB', 'otheruniv')

        """
        super(Fill, self).set(cell)
        self.univ_names[cell] = univ_name

    def comment(self):
        return super(Fill, self).comment("Fill")

    def _comment_unit(self, cell):
        return " {0}".format(self.univ_names[cell])

    def mcnp(self, float_format, sim):
        return super(Fill, self).mcnp(float_format, sim, "FILL")

    def _mcnp_unit(self, float_format, sim, cell):
        return "{0}".format(sim.sys.universe_num(self.univ_names[cell]))

    @property
    def univ_names(self): return self._univ_names

    @univ_names.setter
    def univ_names(self, value): self._univ_names = value


class Lattice(ICellMod):
    """Sets cells to be lattice cells. Unique card with name `lattice`. In MCNP,
    this is the **LAT** card. The user can initialize this card without
    providing any inputs.

    .. inheritance-diagram:: pyne.simplesim.cards.Lattice

    """
    acceptable = ['hexahedra', 'hexagonal']

    def __init__(self, *args):
        """
        Parameters
        ----------
        cell : str name of :py:class:`Cell` or subclass
            The name of the cell for which the lattice setting is being
            provided.
        setting : str
            'hexahedra', or 'hexagonal'.
        *args : cell, setting, ...
            To provide lattice settings for more than one cell, supply the last
            2 arguments for the other cells. This can also be done using
            :py:meth:`set`.

        Examples
        --------
        The following shows a typical usage, where 'pincell' is a cell
        representing a fuel pin, 'coolantcell' is a cell representing the
        coolant flowing around the pin, 'unit' is a cell filled with the
        'unitcell' universe, and is made into a hexahedra lattice::
            
            uni = Universes()
            uni.set('pincell', 'unitcell', True)
            uni.set('coolantcell', 'unitcell', True)
            fill = Fill('unit', 'unitcell')
            lat = Lattice('unit', 'hexahedra')

        Inputs for multiple cells can be provided::

            lat = Lattice('cellA', 'hexahedra', 'cellB', 'hexagonal')

        """
        super(Lattice, self).__init__('lattice', 2, *args)
        self.settings = collections.OrderedDict()
        self._process_varargs(args)

    def set(self, cell, setting):
        """Add an entry or modify an existing one.

        Parameters
        ----------
        cell : str name of :py:class:`Cell` or subclass
        setting : str

        Examples
        --------
        The example above can be achieved by the following::

            lat = Lattice()
            lat.set('cellA', 'hexahedra')
            lat.set('cellB', 'hexagonal')

        """
        super(Lattice, self).set(cell)
        self.settings[cell] = setting

    def comment(self):
        return super(Lattice, self).comment("Lattice")

    def _comment_unit(self, cell):
        return " {0}".format(self.settings[cell])

    def mcnp(self, float_format, sim):
        return super(Lattice, self).mcnp(float_format, sim, "LAT")

    def _mcnp_unit(self, float_format, sim, cell):
        return "{0:d}".format(1 if self.settings[cell] == 'hexahedra' else 2)

    @property
    def settings(self): return self._settings

    @settings.setter
    def settings(self, value):
        for val in value:
            if val not in self.acceptable:
                raise ValueError("``setting`` {0!r} not acceptable.".format(
                        val))
        self._settings = value


class Volume(ICellMod):
    """Cell volumes. Unique card with name `volume`. In MCNP, this is the
    **VOL** card. The user can initialize this card without providing any cell
    volumes.

    .. inheritance-diagram:: pyne.simplesim.cards.Volume

    """
    def __init__(self, *args, **kwargs):
        """
        Parameters
        ----------
        cell : str name of :py:class:`Cell` or subclass
            The name of the cell for which the volume is being provided.
        vol : float [centimeters^3]
            Volume of the cell.
        *args : cell, vol, ...
            To provide volumes for more than one cell, supply the last two
            arguments for the other cells. See examples. This can also be done
            using :py:meth:`set`.
        manual : bool, optional (default: False)
            If set to True, the code will not attempt to calculate cell volumes
            on its own.

        Examples
        --------
        The following are examples of the usage of this card::

            vol = Volume('cellA', 1)
            vol = Volume('cellA', 1, 'cellB', 2, manual=True)

        """
        super(Volume, self).__init__('volume', 2, *args)
        self.manual = kwargs.get('manual', False)
        self.vols = dict()
        self._process_varargs(args)

    def set(self, cell, vol):
        """Add an entry or modify an existing one.

        Parameters
        ----------
        cell : str name of :py:class:`Cell` or subclass
        vol : float

        Examples
        --------
        The example above can be achieved by the following::

            vol = Volume(manual=True)
            vol.set('cellA', 1)
            vol.set('cellB', 2)

        Previously provided values can be modified later on::

            vol.set('cellB', 3)

        """
        super(Volume, self).set(cell)
        self.vols[cell] = vol

    def comment(self):
        if self.manual: poststr = " (all manual)"
        else: poststr = ""
        return super(Volume, self).comment("Volume", poststr)
 
    def _comment_unit(self, cell):
        return " {0:g} cm^3".format(self.vols[cell])

    def mcnp(self, float_format, sim):
        string = "VOL"
        if self.manual: string += " NO"
        return super(Volume, self).mcnp(float_format, sim, string)

    def _mcnp_unit(self, float_format, sim, cell):
        return float_format % self.vols[cell]

    @property
    def manual(self): return self._manual

    @manual.setter
    def manual(self, value): self._manual = value

    @property
    def vols(self): return self._vols

    @vols.setter
    def vols(self, value): self._vols = value


class Area(ICellMod):
    """Cell surface areas. Unique card with name `area`. In MCNP, this is the
    **AREA** card. The user can initialize this card without providing any cell
    areas.

    .. inheritance-diagram:: pyne.simplesim.cards.Area

    """
    def __init__(self, *args):
        """
        Parameters
        ----------
        cell : str name of :py:class:`Cell` or subclass
            The name of the cell for which the area is being provided.
        area : float [centimeters^2]
            Surface area of the cell.
        *args : cell, area, ...
            To provide areas for more than one cell, suppy the last two
            arguments for the other cells. See examples. This can also be
            done using :py:meth:`set`.

        Examples
        --------
        The follow are examples of the usage of this card::

            area = Area('cellA', 10)
            area = Area('cellA', 10, 'cellB', 20)

        """
        super(Area, self).__init__('area', 2, *args)
        self.areas = dict()
        self._process_varargs(args)

    def set(self, cell, area):
        """Add an entry or modify an existing one.

        Parameters
        ----------
        cell : str name of :py:class:`Cell` or subclass
        area : float [centimeters^2]

        Examples
        --------
        The example above can be achieved by the following::

            area = Area()
            area.set('cellA', 10)
            area.set('cellB', 20)

        Previously provided values can be modified later on::

            area.set('cellB', 30)

        """
        super(Area, self).set(cell)
        self.areas[cell] = area
 
    def comment(self):
        return super(Area, self).comment("Area")

    def _comment_unit(self, cell):
        return " {0:g} cm^2".format(self.areas[cell])

    def mcnp(self, float_format, sim):
        return super(Area, self).mcnp(float_format, sim, "AREA")

    def _mcnp_unit(self, float_format, sim, cell):
        return float_format % self.areas[cell]

    @property
    def areas(self): return self._areas

    @areas.setter
    def areas(self, value): self._areas = value


class FissionTurnoff(ICellMod):
    """Fission turnoff options. Unique card with name `fissionturnoff`. In
    MCNP, this is the **NONU** card. If no arguments are provided, then the
    setting of 'capture-gamma' is applied to all cells.

    .. inheritance-diagram:: pyne.simplesim.cards.FissionTurnoff

    """
    def __init__(self, *args):
        """
        Parameters
        ----------
        cell : str name of :py:class:`Cell` or subclass
            The name of the cell for which input is provided.
        setting : str
            One of the following 3 strings: 'capture-gamma' to request fissions
            to count as capture, but still cause gamma radiation; 'real-gamma'
            for fission to actually occur and to also produce gamma radiation;
            and 'capture-nogamma' for fission to be counted as capture, and to
            not release gammas.
        *args : cell, setting, ...
            To provide a setting for more than one cell, supply the last two
            arguments for the other cells. See examples. This can also be done
            using :py:meth:`set`.

        Examples
        --------
        The following sets all cells to the 'capture-gamma' setting::

            fto = FissionTurnoff()

        The followings specifies settings for more than one cell::

            fto = FissionTurnoff('cellA', 'real-gamma', 'cellB',
                    'capture-nogamma')

        """
        super(FissionTurnoff, self).__init__('fissionturnoff', 2, *args)
        self.settings = dict()
        self._process_varargs(args)

    def set(self, cell, setting):
        """Add an entry or modify an existing one.

        Parameters
        ----------
        cell : str name of :py:class:`Cell` or subclass
        setting : str

        Examples
        --------
        The example above can be achieved by the following::
        
            fto = FissionTurnoff()
            fto.set('cellA', 'real-gamma')
            fto.set('cellB', 'capture-nogamma')

        Previously provided values can be modified later on::

            fto.set('cellB', 'real-gamma')

        """
        super(FissionTurnoff, self).set(cell)
        self.settings[cell] = setting
 
    def comment(self):
        string = "Fission turnoff"
        if len(self.cells) == 0:
            return string + " 'fissionturnoff': capture-gamma for all cells."
        else:
            return super(FissionTurnoff, self).comment("Fission turnoff")

    def _comment_unit(self, cell):
        return " {0}".format(self.settings[cell])

    def mcnp(self, float_format, sim):
        string = "NONU"
        if len(self.cells) == 0: return string
        return super(FissionTurnoff, self).mcnp(float_format, sim, string)

    def _mcnp_unit(self, float_format, sim, cell):
        if self.settings[cell] == 'capture-gamma':     return "0"
        elif self.settings[cell] == 'real-gamma':      return "1"
        elif self.settings[cell] == 'capture-nogamma': return "2"
        else:
            raise ValueError("Expected 'capture-gamma', 'real-gamma', "
                    "'capture-nogamma'. User provided {0!r}.".format(
                    self.settings[cell]))

    @property
    def settings(self): return self._settings

    @settings.setter
    def settings(self, value): self._settings = value


class PhotonWeight(ICellMod):
    """Controls production of neutron-induced photons through particle weight.
    Unique card with name `photonweight`. In MCNP, this is the **PWT** card.
    Unlike with other subclasses of :py:class:`ICellMod` and
    :py:class:`ICellModParticle`, this class cannot be formed via only the
    constructor; the :py:meth:`set` method must be used.

    .. inheritance-diagram:: pyne.simplesim.cards.PhotonWeight

    """
    Weight = collections.namedtuple('Weight', ['setting', 'pre_weight'])
    """Named tuple with which input entries are stored."""

    def __init__(self):
        # The reason for the basic constructor is that it is important that I
        # use an optional keyword argument here, and so a variable argument
        # list is not so pleasant to work with.
        super(PhotonWeight, self).__init__('photonweight', 2)
        self.weights = dict()

    def set(self, cell, setting, pre_weight=False):
        """Add an entry or modify an existing one.

        Parameters
        ----------
        cell : str name of :py:class:`Cell` or subclass
            The name of the cell for which the photon weight is being provided.
        setting : str/float
            The setting 'off' turns of photon production, the setting 'one'
            produces one photon per neutron collision if photon production is
            enabled, and if a float is provided, it is used to compute a weight
            threshold for the creation of photons.
        pre_weight : bool, optional
            If True, then the weight threshold incorporates the starting weight
            of the neutron that induced the photon. This argument is ignored if
            ``setting`` is of type ``str``.

        Examples
        --------
        The card is first initialized, after which weights can be added using
        this method. The following turns off neutron-induced photon
        production for 'cellA'::

            pw = PhotonWeight()
            pw.set('cellA', 'off')

        The following requests that one photon is produced for every neutron
        collision in 'cellB'::

            pw.set('cellB', 'one')

        An explicit weight can be given, as well::

            pw.set('cellC', 0.5)

        The starting neutorn weight can be incorporated in the weight
        threshold::

            pw.set('cellD', 0.5, pre_weight=True)

        Previously provided values may be modified later on::

            pw.set('cellD', 0.7, pre_weight=True)

        """
        super(PhotonWeight, self).set(cell)
        if type(setting) is float and setting < 0:
            raise ValueError("The weight must be non-negative.")
        self.weights[cell] = self.Weight(setting, pre_weight)

    def comment(self):
        return super(PhotonWeight, self).comment("Photon weight thresholds")

    def _comment_unit(self, cell):
        setting = self.weights[cell].setting
        if type(setting) is str: return " " + setting
        else: return " {0:g}{1}".format(setting, 
                " (pre-weighted)" if self.weights[cell].pre_weight else "")

    def mcnp(self, float_format, sim):
        if len(self.cells) == 0:
            raise Exception("No entries provided.")
        string = "PWT"
        return super(PhotonWeight, self).mcnp(float_format, sim, string)

    def _mcnp_unit(self, float_format, sim, cell):
        setting = self.weights[cell].setting
        # The following if-statement is for if setting is a float.
        if self.weights[cell].pre_weight: setting = -setting
        if setting == 'off':   return "-1.0E6"
        elif setting == 'one': return "0"
        else:                  return float_format % setting

    @property
    def weights(self): return self._weights

    @weights.setter
    def weights(self, value):
        self._weights = value


class DetectorContribution(ICellMod):
    """Probability of contribution of cell to a detector. Unique card with name
    `detcontrib-<detname>`. In MCNP, this is the **PD** card.

    .. inheritance-diagram:: pyne.simplesim.cards.DetectorContribution

    """
    def __init__(self, det_name, *args):
        """
        Parameters
        ----------
        det_name : str
            Name of the detector (subclass of :py:class:`IDetector`) for which
            this card stores probabilities of contribution.
        cell : str name of :py:class:`Cell` or subclass
            The name of the cell for which the probability of contribution is
            being provided.
        prob : float
            Probability that this cell contributes to this detector.
        *args : cell, prob, ...
            To provide probabilities for more than one cell, supply the last
            two arguments for the other cells. See examples. This can also be
            done using :py:meth:`set`.

        Examples
        --------
        Suppose a cell 'cellA' and a detector 'det1' have been
        defined and are in the simulation. Then, we can do the following::

            dc = DetectorContribution('det1', 'cellA', 0.5)

        or, if we also have a 'cellB':

            dc = DetectorContribution('det1', 'cellA', 0.5, 'cellB', 0.6)


        """
        self.det_name = det_name
        super(DetectorContribution, self).__init__(
                'detcontrib-{0}'.format(det_name), 2, *args)
        self.probs = dict()
        self._process_varargs(args)

    def set(self, cell, prob):
        """Add an entry or modify an existing one.

        Parameters
        ----------
        cell : str name of :py:class:`Cell` or subclass
        prob : float

        Examples
        --------
        The last example above can be achieved by the following::

            dc = DetectorContribution('det1')
            dc.set(cellA, 0.5)
            dc.set('cellB', 0.6)

        Previously provided values can be modified later on::

            dc.set('cellB', 0.7)

        """
        super(DetectorContribution, self).set(cell)
        self.probs[cell] = prob
 
    def comment(self):
        string = "Detector contribution" # to {0!r}".format(self.det_name)
        return super(DetectorContribution, self).comment(string)

    def _comment_unit(self, cell):
        return " {0}".format(self.probs[cell])

    def mcnp(self, float_format, sim):
        # TODO exception if the detector doesn't exist.
        string = "PD{0}".format(sim.tally_num(self.det_name) * 10 + 5)
        return super(DetectorContribution, self).mcnp(
                float_format, sim, string)

    def _mcnp_unit(self, float_format, sim, cell):
        return float_format % self.probs[cell]

    @property
    def probs(self): return self._probs

    @probs.setter
    def probs(self, value): self._probs = value


class IUniqueParticle(IMisc):
    """This class is not used by the user. Abstract base class for cards that
    are unique for a given particle.

    """
    # TODO not passig any *args or **kwargs.
    __metaclass__ = abc.ABCMeta

    def __init__(self, pre_name, particle, *args, **kwargs):
        """
        Parameters
        ----------
        pre_name : str
            First part of card :py:attr:`name`. Names are of the form
            ``<pre_name>-<particle>``.
        particle : str
            A particle string, taken from the keys of
            :py:attr:`Particle.mcnp_abbrev`.

        """
        super(IUniqueParticle, self).__init__(
                "{0}-{1}".format(pre_name, particle), unique=True)
        self.particle = Particle(particle)

    @property
    def particle(self): return self._particle

    @particle.setter
    def particle(self, value): self._particle = value


class Temperature(ICellMod):
    """Temperature of cells. Unique card with name `temperature`. If a time
    index is provided, the unique name is `temperature-idx<index>`. In MCNP, this
    is the **TMP** card and is used for the free gas thermal treatment. This
    card may require that :py:class:`TemperatureTimes` card is in the
    simulation.

    .. inheritance-diagram:: pyne.simplesim.cards.Temperature

    """
    # TODO check if THTME card exists.
    def __init__(self, *args, **kwargs):
        """
        Parameters
        ----------
        cell : str name of :py:class:`Cell` or subclass
            Name of the cell for which the temperature is being specified.
        temp : float [Kelvin]
            Temperature of the cell. For MCNP, this is converted into MeV.
        *args : cell, temp, ...
            To provide the temperature for more than one cell, supply the last
            three arguments for the other cells. See examples. This can also be
            done with :py:meth:`set`.
        index : int, optional
            The index (starting at 1) of a time on the
            :py:class:`TemperatureTimes` card for which the temperatures are
            provided. If no :py:class:`TemperatureTimes` card is used, the user
            does not provide this input, and it defaults to ``None``.

        Examples
        --------
        The following are examples of the usage of this card::

            temp = Temperature('cellA', 600)
            temp = Temperature('cellA', 600, 'cellB', 900, index=2)

        """
        self.index = kwargs.get('index', None)
        name = 'temperature'
        if self.index: name += '-idx{0}'.format(self.index)
        super(Temperature, self).__init__(name, 2, *args)
        self.temps = dict()
        self._process_varargs(args)

    def set(self, cell, temp):
        """Add an entry or modify an existing one.

        Parameters
        ----------
        cell : str name of :py:class:`Cell` or subclass
        temp : float [Kelvin]

        Examples
        --------
        The following does the same as the example above::

            temp = Temperature(index=2)
            temp.set('cellA', 600)
            temp.set('cellB', 900)

        Previously provided values can be modified later on::

            temp.set('cellB', 950)

        """
        super(Temperature, self).set(cell)
        if temp < 200:
            warnings.warn("Temperature set as less than 200 K. "
                    "Are you trying to specify temperature in degrees "
                    "Celcius, etc.? User provided {0:g}.".format(temp))
        if temp < 1:
            warnings.warn("Temperature set as less than 1 K. "
                    "Are you trying to specify temperature as 'kT'? "
                    "User provided {0:g}.".format(temp))
        self.temps[cell] = temp

    def comment(self):
        # TODO get access to sim and show the actual time, not just the index.
        string = "Temperatures"
        if self.index: string += " for time index {0}".format(self.index)
        return super(Temperature, self).comment(string)

    def _comment_unit(self, cell):
        return " {0:g} K".format(self.temps[cell])

    def mcnp(self, float_format, sim):
        string = "TMP"
        if self.index: string += "{0}".format(self.index)
        return super(Temperature, self).mcnp(float_format, sim, string)

    def _mcnp_unit(self, float_format, sim, cell):
        return float_format % (self.temps[cell] * self.kelvin2kT)

    @property
    def index(self): return self._index

    @index.setter
    def index(self, value): self._index = value

    @property
    def temps(self): return self._temps

    @temps.setter
    def temps(self, value): 
        self._temps = value

class TemperatureTimes(IMisc):
    """Times at which temperatuers are specified on the :py:class:`Temperature`
    card. Unique card with name `temptimes`. In MCNP, this is the **THTME**
    card.

    .. inheritance-diagram:: pyne.simplesim.cards.TemperatureTimes

    """
    # TODO ideally this would not be a separate card, and would be generated by
    # only using a Temperature card.
    def __init__(self, times):
        """
        Parameters
        ----------
        times : float [seconds]
            The times at which temperature cards are provided.

        Examples
        --------
        Here is an example of how this card is used::

            thtme = TemperatureTimes([1e10, 2e10])

        """
        super(TemperatureTimes, self).__init__('temptimes', unique=True)
        self.times = times

    def comment(self):
        string = "Temperature times {0!r} (in MeV):".format(self.name)
        for val in self.times: string += " {0:g}".format(val)
        return string + "."

    def mcnp(self, float_format, sim):
        string = "THTME"
        for val in self.times:
            string += " "
            string += float_format % (val * self.secs2shakes)
        return string

    @property
    def times(self): return self._times

    @times.setter
    def times(self, value): self._times = value


class ICellModParticle(IUniqueParticle):
    """This class is not used by the user. Abstract base class for cards that
    can be specified in MCNP on both the cell card or in the data block, and
    are unique by particle.  All subclasses have a ``particle`` and ``cell``
    property, and similar form.  All subclasses are unique for a given particle
    type. Entries for a given cell can be modified by providing an input for
    the same cell.

    """
    __metaclass__ = abc.ABCMeta

    def __init__(self, pre_name, particle, n_args_per_cell, *args,
                 **kwargs):
        """
        Parameters
        ----------
        pre_name : str
            See :py:class:`IUniqueParticle`.
        particle : str
            See :py:class:`IUniqueParticle`.
        cell : str name of :py:class:`Cell` or subclass.
            The name of the cell for which the card applies.
        n_args_per_cell : int
            The number of arguments the subclass expects per cell.

        """
        super(ICellModParticle, self).__init__(pre_name, particle)
        self.cells = []
        self._n_args_per_cell = n_args_per_cell
        if len(args) % n_args_per_cell != 0:
            raise StandardError("The length of ``*args`` must be a multiple "
                    "of {0}. Length is {1}.".format(n_args_per_cell, len(args)))

    def _process_varargs(self, args):
        for i_cell in range(len(args) / self._n_args_per_cell):
            i_start = self._n_args_per_cell * i_cell
            self.set(*args[i_start:i_start+self._n_args_per_cell])

    @abc.abstractmethod
    def set(self, cell):
        """Add an entry or modify an existing one. See subclasses."""
        if cell not in self.cells: self.cells += [cell]

    @abc.abstractmethod
    def comment(self, title):
        string = "{0} {1!r}:".format(title, self.name)
        counter = 0
        for cell in self.cells:
            counter += 1
            string += " cell {0!r}".format(cell)
            string += self._comment_unit(cell)
            if counter < len(self.cells): string += ";"
        return string + "."

    @abc.abstractmethod
    def _comment_unit(self):
        raise NotImplementedError

    @abc.abstractmethod
    def mcnp(self, float_format, sim, keystring):
        # TODO this ordering might not be correct, particularly once we add
        # support for universes, etc.
        string = "{0}:{1!s}".format(keystring, self.particle.mcnp())
        # TODO this should loop through in the print order.
        for key in sim.sys.cells: 
            if key in self.cells:
                string += " " + self._mcnp_unit(float_format, sim, key)
            else: 
                string += " 0"
        return string

    @abc.abstractmethod
    def _mcnp_unit(self):
        raise NotImplementedError

    @property
    def cells(self): return self._cells

    @cells.setter
    def cells(self, value):
        for val in value:
            if type(val) is not str:
                raise Exception("Expected cell name, got {0}".format(val))
        self._cells = value


class Importance(ICellModParticle):
    """Particle importance. Unique card for a given particle type, with name
    `importance-<particle>`. In MCNP, this is the **IMP** card. The user can
    initialize this card without providing any cell importances.  The typical
    usage is to provide the importance for all cells, though the card does not
    require this.
    .. inheritance-diagram:: pyne.simplesim.cards.Importance

    """
    def __init__(self, particle, *args):
        """
        Parameters
        ----------
        particle : str
            See :py:class:`ICellMod`.
        cell : str name of :py:class:`Cell` or subclass.
            See :py:class:`ICellMod`.
        imp : int
            Cell importance.
        *args : cell, imp, ...
            To provide importances for more than one cell, supply the last two
            arguments for the other cells. See examples. This can also be done
            using :py:meth:`set`.

        Examples
        --------
        The following specifies the neutron importance in cells 'cellA' and
        'cellB'::

            impn = Importance('neutron', 'cellA', 1, 'cellB', 2)

        If providing importances for many cells it may be easier to do the
        following::

            args = ['cellA', 1, 'cellB', 1, 'cellC', 0]
            impn = cards.Importance('neutron', *args)

        """
        super(Importance, self).__init__("importance", particle, 2, *args)
        self.imps = dict()
        self._process_varargs(args)

    def set(self, cell, imp):
        """Add an entry or modify an existing one.

        Parameters
        ----------
        cell : str name of :py:class:`Cell` or subclass.
        imp : int

        Examples
        --------
        The example above can be achieved by the following::

            impn = Importance('neutron')
            impn.set('cellA', 1)
            impn.set('cellB', 2)

        Previously provided values can be modified later on::

            impn.set('cellB', 1)

        """
        super(Importance, self).set(cell)
        self.imps[cell] = imp

    def comment(self):
        return super(Importance, self).comment("Importance")

    def _comment_unit(self, cell):
        return " " + str(self.imps[cell])

    def mcnp(self, float_format, sim):
        return super(Importance, self).mcnp(float_format, sim, "IMP")

    def _mcnp_unit(self, float_format, sim, cell):
        return "{0}".format(self.imps[cell])

    @property
    def imps(self): return self._imps

    @imps.setter
    def imps(self, value): self._imps = value


class ExponentialTransform(ICellModParticle):
    """An exponential transform that adjusts the total cross section by a given
    factor in a given direction. Unique card for a given particle type, with
    name `exptransform-<particle>`. In MCNP, this is the **EXT** card.  The
    user can initialize this card without requesting any exponential
    transforms.

    .. inheritance-diagram:: pyne.simplesim.cards.ExponentialTransform

    """
    def __init__(self, particle, *args):
        """
        Parameters
        ----------
        particle : str
            See :py:class:`ICellMod`.
        cell : str name of :py:class:`Cell` or subclass
            See :py:class:`ICellMod`.
        stretch : str or float
            The stretch factor. If 'capture-to-total', the stretch factor is
            the ratio of the capture cross section to the total cross section
            (referred to as Sigma_a in the MCNP manual). Otherwise, the factor
            is a number between 0 and 1.
        direction : str 
            If 'currdir' then the stretching is done in the particle's
            direction of travel. If 'x', 'X', 'y', 'Y', or 'z', 'Z', the
            stretching is with respect to the requested axis. Otherwise, it is
            the name of a vector on the :py:class:`Vector` card.
        sign : str
            If 'toward', the stretching is done toward the direction requested.
            If 'away', the stretching is done away from the direction requested.
        *args : cell, stretch, direction, sign...
            To request an exponential transform for more than one cell, supply
            the last four arguments for the additional cells. See examples.
            This can also be done using :py:meth:`set`.
 
        Examples
        --------
        Consider cell 'cellA'. The following requests a
        transformation that stretches by the ratio of the capture cross section
        to the total cross section, in the direction of the particle's travel.
        The name of the card will be `exptransform-neutron`::

            extn = ExponentialTransform('neutron', 'cellA', 'capture-to-total',
                    'currdir', 'toward')
            assert extn.name == 'exptransform-neutron'
        
        The following requests a transformation that stretches by a factor of
        0.5 in the direction of the particle's travel::

            extn = ExponentialTransform('neutron', 'cellA', 0.5, 'currdir',
                    'toward')

        The following requests a transformation with respect to the x axis::

            extn = ExponentialTransform('neutron', 'cellA', 0.5, 'x', 'toward')

        The following requests a transformation away from the origin::

            extn = ExponentialTransform('neutron', 'cellA', 0.5, 'vec1',
                    'away')

        where somewhere else the user has added::

            vec = Vector()
            vec.set('vec1', [0, 0, 0])

        to the simulation. If the user wants to request an exponential
        transform for cells 'cellB' and 'cellC' as well, they can do the
        following::

            extn = ExponentialTransform('neutron', 
                    'cellA', 'capture-to-total', 'currdir', 'toward', 
                    'cellB', 0.5, 'currdir', 'toward',
                    'cellC', 0.5, 'vec1', 'away')

        """
        super(ExponentialTransform, self).__init__('exptransform', particle,
                                                   4, *args)
        # Initialize properties for the subclass.
        self.stretchs = dict()
        self.directions = dict()
        self.signs = dict()
        # If information for multiple cells has been provided...
        self._process_varargs(args)

    def set(self, cell, stretch, direction, sign):
        """The user can add additional transforms, for additional cells, or
        modify existing values using this method. See above for a description
        of the input.

        Parameters
        ----------
        cell : str name of :py:class:`Cell` or subclass
        stretch : str or float
        direction : str 
        sign : str

        Examples
        --------
        The last example above can also be achieved by the following, assuming
        cell cards 'cellA', 'cellB' and 'cellC' have been created::

            extn = ExponentialTransform('neutron') 
            extn.set('cellA', 'capture-to-total', 'currdir', 'toward')
            extn.set('cellB', 0.5, 'currdir', 'toward')
            extn.set('cellC', 0.5, 'vec1', 'away')

        Previously provided values can be modified later on::

            extn.set('cellC', 0.7, 'vec1', 'away')

        """
        super(ExponentialTransform, self).set(cell)
        self.stretchs[cell] = stretch
        self.directions[cell] = direction
        self.signs[cell] = sign

    def comment(self):
        return super(ExponentialTransform, self).comment(
                "Exponential transform")

    def _comment_unit(self, cell):
        string = " stretch by "
        if type(self.stretchs[cell]) is str: string += self.stretchs[cell]
        else: string += "{0:g}".format(self.stretchs[cell])
        string += " "
        string += self.signs[cell] + " "
        if self.signs[cell] == 'away': string += "from "
        string += self.directions[cell]
        return string

    def mcnp(self, float_format, sim):
        return super(ExponentialTransform, self).mcnp(float_format, sim,
                "EXT")

    def _mcnp_unit(self, float_format, sim, cell):
        string = ""
        if self.signs[cell] == 'away': string += "-"
        if self.stretchs[cell] == 'capture-to-total': string += "S"
        else: string += float_format % self.stretchs[cell]
        if self.directions[cell] == 'currdir': pass
        elif (self.directions[cell].upper() == 'X' or
              self.directions[cell].upper() == 'Y' or
              self.directions[cell].upper() == 'Z'):
            string += self.directions[cell].upper()
        else:
            if 'vector' not in sim.misc:
                raise Exception("Vector card is needed for the Exponential "
                        "transform card.")
            vecname = self.directions[cell]
            string += "V{0}".format(sim.misc['vector'].index(vecname))
        return string

    # I'm aware that this spelling of stretchs is incorrect; it's for
    # consistency.
    @property
    def stretchs(self): return self._stretchs

    @stretchs.setter
    def stretchs(self, value): self._stretchs = value

    @property
    def directions(self): return self._directions

    @directions.setter
    def directions(self, value): self._directions = value

    @property
    def signs(self): return self._signs

    @signs.setter
    def signs(self, value):
        for arg in value:
            if arg != 'toward' and arg != 'away':
                raise ValueError("The value of ``sign`` must be 'toward' or "
                        "'away'. User provided {0!r}.".format(arg))
        self._signs = value


class ForcedCollision(ICellModParticle):
    """A forced collision setting. Unique card for a given particle type, with
    name `forcedcoll-<particle>`. In MCNP, this is the **FCL** card. The user
    can initialize this card without requesting forced collisions.

    .. inheritance-diagram:: pyne.simplesim.cards.ForcedCollision

    """
    def __init__(self, particle, *args):
        """
        Parameters
        ----------
        particle : str
            See :py:class:`ICellMod`.
        cell : str name of :py:class:`Cell` or subclass
            See :py:class:`ICellMod`.
        prob : float
            If 'none', there is no forced collision for this cell.
            If float (between 0 and 1), it is the survival probability as
            described in the MCNP manual.
        only_entering : bool
            If True, the card applies only to particles entering the cell. If
            False, the card applies to particles entering as well as those
            surviving weight games in the cell.
        *args : cell, prob, only_entering, ...
            To request forced collision for more than one cell, supply
            the last three arguments for the additional cells. See examples.
            This can also be done using :py:meth:`set`.

        Examples
        --------
        The following cards request forced collisions with neutrons in the cell
        defined by cell card 'cellA'::  

            fcl = ForcedCollision('neutron', 'cellA', 0.5, True)
            fcl = ForcedCollision('neutron', 'cellA', 0.5, False)

        If the user wants to request forced collisions for 'cellB' as well,
        they can do the following::

            fcl = ForcedCollision('neutron', 'cellA', 0.5, True,
                                             'cellB', 0.5, False)

        """
        super(ForcedCollision, self).__init__('forcedcoll', particle, 3, *args)
        self.probs = dict()
        self.only_enterings = dict()
        # If information for multiple cells has been provided...
        self._process_varargs(args)

    def set(self, cell, prob, only_entering):
        """The user can add additional forced collision probabilities, for
        additional cells, or modify existing data, using this method. See above
        for a description of the input.

        Parameters
        ----------
        cell : str name of :py:class:`Cell` or subclass
        prob : float
        only_entering : bool

        Examples
        --------
        The last example above can also be achieved by the following::

            fcl = ForcedCollision('neutron')
            fcl.set('cellA', 0.5, True)
            fcl.set('cellB', 0.5, False)

        Previously provided values can be modified later on::

            fcl.set('cellB', 0.7, False)

        """
        super(ForcedCollision, self).set(cell)
        self.probs[cell] = prob
        self.only_enterings[cell] = only_entering

    def comment(self):
        return super(ForcedCollision, self).comment("Forced collision")

    def _comment_unit(self, cell):
        if self.only_enterings[cell]: oestr = "entering only"
        else: oestr = "entering and weight games"
        return " prob {0:g} for {1}".format(self.probs[cell], oestr)

    def mcnp(self, float_format, sim):
        return super(ForcedCollision, self).mcnp(float_format, sim, "FCL")

    def _mcnp_unit(self, float_format, sim, cell):
        string = "-" if self.only_enterings[cell] else ""
        return string + float_format % self.probs[cell]

    @property
    def probs(self): return self._probs

    @probs.setter
    def probs(self, value): self._probs = value

    @property
    def only_enterings(self): return self._only_enterings

    @only_enterings.setter
    def only_enterings(self, value): self._only_enterings = value


class WeightWindowBound(ICellModParticle):
    """Cell-based weight window lower bounds. Unique card for a given particle
    type, with name `weightwinbound-<particle>`. In MCNP, this is the **WWN**
    card. In MCNP, these are typically generated automatically using its weight
    window generator.

    .. inheritance-diagram:: pyne.simplesim.cards.WeightWindowBound

    """
    def __init__(self, particle, *args):
        """
        Parameters
        ----------
        particle : str
            See :py:class:`ICellMod`.
        idx_energy : int
            Index of the energy, on a :py:class:`WeightWindowEnergies` card,
            for which a bound is being specified. If no energy bins are
            specified (e.g. 1 energy bin for weight windows), provide 1. The
            index starts at 1.
        idx_time : int
            Index of the time, on a :py:class:`WeightWindowTimes` card,
            for which a bound is being specified. If no time intervals are
            specified in the problem (e.g. 1 time interval for weight windows),
            provide 1. The index starts at 1.
        cell : str name of :py:class:`Cell` or subclass
            See :py:class:`ICellMod`.
        bound : str or float
            Lower bound for the weight to cause rouletting, or 'killall' to
            kill all particles entering the cell.
        *args : cell, idx_energy, idx_time, bound, ...
            Any number of sets of the previous four arguments can be provided
            as additional arguments. See examples. This can also be done using
            :py:meth:`set`.

        Examples
        --------
        The following specifies the weight window bounds for 'cellB' if the
        :py:class:`WeightWindowEnergies` and :py:class:`WeightWindowTimes`
        cards are not used::

            wwn = WeightWindowBound('neutron', 1, 1, 'cellB', 0.01)

        The value of ``bound`` can be 'killall'::

            wwn = WeightWindowBound('neutron', 1, 1, 'cellB', 'killall')

        The following::

            wwn = WeightWindowBound('neutron', 1, 1, 'cellB', 0.01,
                                               2, 3, 'cellB', 0.01)

        assumes that cards such as the following have been added to the
        simulation::

            wwe = WeightWindowEnergies('neutron', [1, 10])
            wwt = WeightWindowTimes('neutron', [1, 1e12, 2e12])

        """
        super(WeightWindowBound, self).__init__('weightwinbound', particle,
                                                4, *args)
        # Check for existence of weight window cards.
        self.idx_energys = []
        self.idx_times = []
        # _multi_dict is a method in this class.
        self.bounds = self._multi_dict(4)
        self._process_varargs(args)

    def set(self, idx_energy, idx_time, cell, bound):
        """Add an entry or modify an existing one.

        Parameters
        ----------
        idx_energy : int
        idx_time : int
        cell : str name of :py:class:`Cell` or subclass
        bound : str or float

        Examples
        --------
        The following does the same as the second example above::

            wwn = WeightWindowBound('neutron')
            wwn.set(1, 1, 'cellB', 0.01)
            wwn.set(2, 3, 'cellB', 0.01)

        Previously provided values can be modified later on::

            wwn.set(1, 1, 'cellB', 0.02)

        """
        super(WeightWindowBound, self).set(cell)
        if idx_energy not in self.idx_energys: self.idx_energys += [idx_energy]
        if idx_time not in self.idx_times:     self.idx_times += [idx_time]
        self.bounds[cell][idx_energy][idx_time] = bound

    def comment(self):
        string = "Weight window bounds {0!r} for {1}s:".format(
                self.name, self.particle)
        for i_e in self.idx_energys:
            string += " energy idx {0}:".format(i_e)
            for i_t in self.idx_times:
                string += " time idx {0}:".format(i_t)
                for cell in self.cells:
                    if self.bounds[cell][i_e][i_t]:
                        # The 3-dim dict's entries are initialized to None
                        string += " cell {0!r}:".format(cell)
                        string += self._comment_unit(cell, i_e, i_t)
                        string += ","
        # Change last character from a comma to a period.
        return string[:-1] + "."

    def _comment_unit(self, cell, i_e, i_t):
        bound = self.bounds[cell][i_e][i_t]
        if type(bound) is str: return " {0}".format(bound)
        else:                  return " {0:g}".format(bound)

    def mcnp(self, float_format, sim):
        # Prepare to obtain linear index.
        self._find_n_energies(sim)
        string = ""
        # Finally, create all necessary cards (one per linear index).
        counter = 0
        for i_t in self.idx_times:
            for i_e in self.idx_energys:
                # Start card, but only if any values are assigned for this idx.
                if self._n_vals_for(i_e, i_t) > 0:
                    counter += 1
                    string += "{0}WWN{1}:{2}".format(
                            "\n" if counter > 1 else "",
                            self.i_linear(i_e, i_t),
                            self.particle.mcnp())
                    # Check all cells in the system.
                    for key in sim.sys.cells:
                        # Is this cell on this card, and is there a bound
                        # defined for it, for this energy and time?
                        if key in self.cells and self.bounds[key][i_e][i_t]:
                            string += " " + self._mcnp_unit(
                                    float_format, key, i_e, i_t)
                        else:
                            string += " 0"
        return string

    def _mcnp_unit(self, float_format, cell, i_e, i_t):
        this_bound = self.bounds[cell][i_e][i_t]
        if type(this_bound) is float: return float_format % this_bound
        elif this_bound == 'killall': return "-1"
        else:                         raise ValueError("Unexpected input.")
 
    def _n_vals_for(self, i_e, i_t):
        # TODO may not need this check.
        n_vals = 0
        for cell in self.cells:
            if self.bounds[cell][i_e][i_t]: n_vals += 1
        return n_vals

    def _find_n_energies(self, sim):
        wwge_name = 'weightwingenenergy-{0}'.format(self.particle)
        wwe_name = 'weightwinenergy-{0}'.format(self.particle)
        if wwge_name in sim.misc and wwe_name in sim.misc:
            raise UserWarning("Both a WWGE and a WWE card have been added; "
                    "using WWGE, ignoring WWE.")
        if wwge_name in sim.misc:
            # Deal with MCNP default indices that are unhandled by the WWGT
            # card here.
            if (sim.misc[wwge_name].for_gen and
                    len(sim.misc[wwge_name].bounds) == 0):
                self._n_energies = 10
            else:
                self._n_energies = sim.misc[wwge_name].n_bounds
        elif wwe_name in sim.misc:
            self._n_energies = sim.misc[wwe_name].n_bounds
        else:
            raise Exception("No WWGT:{0} or WWT:{0} card found in the "
                    "simulation.".format(self.particle.mcnp()))

    def i_linear(self, i_e, i_t):
        return (i_t - 1) * self._n_energies + i_e
       
    def _multi_dict(self, n_dims):
        if n_dims <= 1:
            return None
        return collections.defaultdict(lambda: self._multi_dict(n_dims - 1))

    @property
    def idx_energys(self): return self._idx_energys

    @idx_energys.setter
    def idx_energys(self, value): self._idx_energys = value

    @property
    def idx_times(self): return self._idx_times

    @idx_times.setter
    def idx_times(self, value): self._idx_times = value

    @property
    def bounds(self): return self._bounds

    @bounds.setter
    def bounds(self, value): self._bounds = value


class WeightWindowEnergies(IUniqueParticle):
    """Upper energy bounds for weight windows. Unique card for a given particle
    type, with name `weightwinenergy-<particle>` or
    `weightwingenenergy-<particle>`. In MCNP, this is the **WWE**
    or **WWGE** card.

    .. inheritance-diagram:: pyne.simplesim.cards.WeightWindowEnergies

    """
    # The choice to combine this card with the WWGE card means that we can't
    # read in an MCNP input quite as easily? Raise an excpetion if a regular
    # and a generator card are provided by the user; or maybe just a warning,
    # and go with the indices provided by the former. it's possible that the
    # indices are the same. Why wouldn't they be?
    # TODO bounds [] can be provided if for_gen=False, in which case MCNP
    # generates a WWE card with 10 indices on it, that could then be used on a
    # WWN card, though if the user is using WWG then they shouldn't be using
    # WWN directly, but in any case the WWN might say that the WWGE card does
    # not have 10 indices, since only MCNP is handling that. Need to show this
    # to the WWN card. EDIT this is managed in WWN.mcnp()
    def __init__(self, particle, bounds, for_gen=False):
        """
        Parameters
        ----------
        particle : str
            See :py:class:`IUniqueParticle`.
        bounds : list, :py:class:`np.array` [MeV]
            Upper energy bounds of the weight windows.
        for_gen : bool, optional
            If set to True, this becomes the energies for the weight window
            generator energies card, **WWGE**.
            To request the default functionality that MCNP provides on the
            **WWGE** card when no bounds are specified, provide an empty list
            for ``bounds``.

        Examples
        --------
        The following specifies energy bins [0, 1] and [1, 10] MeV for photon
        weight windows::

            wwe = WeightWindowEnergies('photon', [1, 10])            

        The following specifies the energy bins on a weight window generator
        times card::

            wwt = WeightWindowEnergies('photon', [1, 1e12], for_gen=True)            

        The following, in MCNP, uses the default bins if using the generator
        option::
        
            wwt = WeightWindowEnergies('photon', [], for_gen=True)            

        """
        super(WeightWindowEnergies, self).__init__(
                'weightwin{0}energy'.format('gen' if for_gen else ''),
                particle)
        self.bounds = bounds
        self.for_gen=for_gen

    def comment(self):
        if self.for_gen and len(self.bounds) == 0:
            n_bins = "default"
        else:
            n_bins = len(self.bounds)
        return ("Weight window {0}energies {1!r} for {2}s: {3} bins.".format(
            "generator " if self.for_gen else "", self.name, self.particle,
            n_bins))

    def mcnp(self, float_format, sim):
        string = "WW{0}E:{1}".format("G" if self.for_gen else "",
                self.particle.mcnp())
        float_format = " " + float_format
        for bound in self.bounds: string += float_format % bound
        return string
        
    @property
    def bounds(self): return self._bounds

    @bounds.setter
    def bounds(self, value):
        self._bounds = value
        self.n_bounds = len(self.bounds)

    @property
    def for_gen(self): return self._for_gen

    @for_gen.setter
    def for_gen(self, value): self._for_gen = value


class WeightWindowTimes(IUniqueParticle):
    """Upper time bounds for weight windows. Unique card for a given particle
    type, with name `weightwintime-<particle>` or `weightwingentime-<particle`.
    In MCNP this is the **WWT** or **WWGT** card.

    .. inheritance-diagram:: pyne.simplesim.cards.WeightWindowTimes

    """
    def __init__(self, particle, bounds, for_gen=False):
        """
        Parameters
        ----------
        particle : str
            See :py:class:`IUniqueParticle`.
        bounds : list, :py:class:`np.array` [seconds]
            Upper time bounds of the weight windows.
        for_gen : bool, optional
            If set to True, this becomes the energies for the weight window
            generator energies card, **WWGE**.
            To request the default functionality that MCNP provides on the
            **WWGE** card when no bounds are specified, provide an empty list
            for ``bounds``.

        Examples
        --------
        The following specifies time intervals [-inf, 1] and [1, 1e12] sec for
        photon weight windows::
            
            wwt = WeightWindowTimes('photon', [1, 1e12])            

        The following specifies the time intervals on a weight window generator
        times card::

            wwt = WeightWindowTimes('photon', [1, 1e12], for_gen=True)            

        The following, in MCNP, uses the default bins if using the generator
        option::

            wwt = WeightWindowTimes('photon', [], for_gen=True)            

        """
        super(WeightWindowTimes, self).__init__(
                'weightwin{0}time'.format('gen' if for_gen else ''), particle)
        self.bounds = bounds
        self.for_gen = for_gen

    def comment(self):
        # TODO don't tie the comment to any certain units; units are
        # code-dependent.
        if self.for_gen and len(self.bounds) == 0:
            n_bins = "default"
        else:
            n_bins = len(self.bounds)
        return ("Weight window {0}times {1!r} for {2}s: "
                "{3} intervals.".format(
                    "generator " if self.for_gen else "",
                    self.name, self.particle, n_bins))

    def mcnp(self, float_format, sim):
        string = "WW{0}T:{1}".format("G" if self.for_gen else "",
                self.particle.mcnp())
        float_format = " " + float_format
        for bound in self.bounds: 
            string += float_format % (bound * self.secs2shakes)
        return string
        
    @property
    def bounds(self): return self._bounds

    @bounds.setter
    def bounds(self, value): 
        self._bounds = value
        self.n_bounds = len(self.bounds)

    @property
    def for_gen(self): return self._for_gen

    @for_gen.setter
    def for_gen(self, value): self._for_gen = value


class DXTRANContribution(ICellMod):
    """Contribution probabilities to DXTRAN sphere by cells. Unique card for a
    given particle type and DXTRAN sphere, with name
    `dxtrancont-<sph_name>-<particle>` (or more simply `dxtrancont-<particle>`
    if input applies for all spheres). In MCNP, this is the **DXC** card.

    .. inheritance-diagram:: pyne.simplesim.cards.DXTRANContribution

    """
    def __init__(self, particle, sph_name, *args):
        """
        Parameters
        ----------
        particle : str
            See :py:class:`ICellMod`.
        sph_name : str or None
            Name of the DXTRAN sphere, on the :py:class:`DXTRANSpheres` card,
            for which this card applies. To apply this card to all spheres, set
            to None.
        cell : str name of :py:class:`Cell` or subclass
            See :py:class:`ICellMod`.
        prob : float
            The probability of contribution to the sphere.

        Examples
        --------
        The following shows how this card is used to specify contribution
        probabilities for all DXTRAN spheres::

            dxc = DXTRANContribution('neutron', None, 'cellA', 0.5)

        The following applies for only 'sph1' on the :py:class:`DXTRANSpheres`
        card::

            dxc = DXTRANContribution('neutron', 'sph1', 'cellA', 0.5)

        The following provides probabilities for two different cells::

            dxc = DXTRANContribution('neutron', 'sph1', 'cellA', 0.5,
                                                        'cellB', 0.75)

        """
        # DXTRANSphere must be added first before this can be used?
        super(DXTRANContribution, self).__init__(
                'dxtrancont{0}-{1}'.format(
                    "-" + sph_name if sph_name else "", particle),
                2, *args)
        self.particle = Particle(particle)
        self.sph_name = sph_name
        self.probs = dict()
        self._process_varargs(args)

    def set(self, cell, prob):
        """Add an entry or modify an existing one.

        Parameters
        ----------
        cell : str name of :py:class:`Cell` or subclass
        prob : float

        Examples
        --------
        The following shows how this card can be constructed using this
        method::

            dxc = DXTRANContribution('neutron', 'sph1')
            dxc.set('cellA', 0.5)
            dxc.set('cellB', 0.75)

        Previously provided values can be modified later on::

            dxc.set('cellB', 0.8)

        """
        super(DXTRANContribution, self).set(cell)
        self.probs[cell] = prob

    def comment(self):
        string = "DXTRAN contribution"
        if self.sph_name: string += " for sphere {0!r}".format(self.sph_name)
        else:             string += " for all spheres"
        return super(DXTRANContribution, self).comment(string)

    def _comment_unit(self, cell):
        return " {0:g}".format(self.probs[cell])

    def mcnp(self, float_format, sim):
        string = "DXC{0}:{1}".format(
                self.sph_index(sim), self.particle.mcnp())
        return super(DXTRANContribution, self).mcnp(float_format, sim, string)

    def sph_index(self, sim):
        if self.sph_name:
            # Get sphere index.
            dxt_name = 'dxtranspheres-{0}'.format(self.particle)
            if dxt_name not in sim.misc:
                raise StandardError("To specify DXTRAN contributions for "
                        "{0}s, a {1} misc card must be in the "
                        "simulation.".format(self.particle, dxt_name))
            index = sim.misc[dxt_name].index(self.sph_name)
        else:
            # Applies to all spheres.
            index = ''
        return index

    def _mcnp_unit(self, float_format, sim, cell):
        return float_format % self.probs[cell]

    @property
    def particle(self): return self._particle

    @particle.setter
    def particle(self, value): self._particle = value

    @property
    def probs(self): return self._probs

    @probs.setter
    def probs(self, value): self._probs = value


class DXTRANSpheres(IUniqueParticle):
    """DXTRAN spheres. Unique card for a given particle type, with name
    `dixtranspheres-<particle>`. In MCNP, this is the **DXT** card. See the
    code's (e.g. MCNP's) manual for default values. All the spheres added to
    this card are named so that they can be referenced by other cards.

    .. inheritance-diagram:: pyne.simplesim.cards.DXTRANSpheres

    """
    Sphere = collections.namedtuple('Sphere',
            ['name', 'center', 'inrad', 'outrad'])
    """Named tuple with which input entries are stored."""

    def __init__(self, particle, *args, **kwargs):
        """
        Parameters
        ----------
        particle : str
            See :py:class:`IUniqueParticle`.
        sph_name : str
            Name of the sphere.
        center : 3-element list or :py:class:`np.array` [centimeters]
            Center of this sphere.
        inner_radius : float [centimeters]
            Inner radius of this sphere.
        outer_radius : float [centimeters]
            Outer radius of this sphere.
        *args : sph_name, point, inner_radius, outer_radius, ...
            To request more than one sphere, supply the last 4 arguments for
            the additional spheres. See examples. This can also be done using
            :py:meth:`set`.
        upper_cutoff : float, optional
            Upper weight cutoff, for all spheres.
        lower_cutoff : float, optional
            Lower weight cutoff, for all spheres.
        min_photon_weight : float, optional
            Only relevant if this card is for neutrons.

        Examples
        --------
        The following requests a DXTRAN sphere `sph1` at (1, 2, 3) cm, with
        inner radius 4 cm and outer radius 5 cm::

            dsph = DXTRANSpheres('neutron', 'sph1', [1, 2, 3], 4, 5)

        It is possible to define multiple spheres in the constructor; keyword
        arguments come after the sphere parameters::

            dsph = DXTRANSpheres('neutron', 'sph1', [1, 2, 3], 4, 5,
                                            'sph2', [4, 5, 6], 7, 8,
                                 upper_cutoff=0.1, lower_cutoff=0.05,
                                 min_photon_weight=0.5)

        """
        super(DXTRANSpheres, self).__init__('dxtranspheres', particle)
        self.upper_cutoff = kwargs.get('upper_cutoff', None)
        self.lower_cutoff = kwargs.get('lower_cutoff', None)
        self.min_photon_weight = kwargs.get('min_photon_weight', None)
        self.spheres = collections.OrderedDict()
        self._n_args_per_set = 4
        # TODO copied from ICellMod, ICellModParticle.
        for i_set in range(len(args) / self._n_args_per_set):
            i_start = self._n_args_per_set * i_set
            self.set(*args[i_start:i_start+self._n_args_per_set])

    def set(self, sph_name, center, inner_radius, outer_radius):
        """Add an entry or modify an existing one.

        Parameters
        ----------
        sph_name : str
            Name of the sphere.
        center : 3-element list or :py:class:`np.array` [centimeters]
            Center of this sphere.
        inner_radius : float [centimeters]
            Inner radius of this sphere.
        outer_radius : float [centimeters]

        Examples
        --------
        The last example above can be achieved as follows::

            dsph = DXTRANSpheres('neutron', upper_cutoff=0.1,
                                 lower_cutoff=0.05,
                                 min_photon_weight=0.5)
            dsph.set('sph1', [1, 2, 3], 4, 5)
            dsph.set('sph2', [4, 5, 6], 7, 8)

        Previously provided values can be modified later on::

            dsph.set('sph2', [4.5, 5.5, 6.5], 8.5, 9.5)

        """
        self.spheres[sph_name] = self.Sphere(
                sph_name, center, inner_radius, outer_radius)

    def comment(self):
        upcut = self.upper_cutoff
        lowcut = self.lower_cutoff
        mpw = self.min_photon_weight
        string = ("DXTRAN spheres {0!r}: "
                  "up. cut. {1}, low cut. {2}, min photon wgt. {3}.".format(
                      self.name,
                      upcut  if upcut  else 'default',
                      lowcut if lowcut else 'default',
                      mpw    if mpw    else 'default'))
        counter = 0
        for name in self.spheres:
            counter += 1
            string += " sphere " + self._comment_unit(name)
            if counter < len(self.spheres): string += ";"
        return string + "."

    def _comment_unit(self, name):
        sph = self.spheres[name]
        return ("{0!r} at ({1[0]:g}, {1[1]:g}, {1[2]:g}) cm, "
                "in. rad. {2:g} cm, out. rad. {3:g} cm".format(
                    name, sph.center, sph.inrad, sph.outrad))

    def mcnp(self, float_format, sim):
        mpw = self.min_photon_weight
        string = "DXT:{0}".format(self.particle.mcnp())
        for name in self.spheres:
            string += "  " + self._mcnp_unit(float_format, sim, name)
        string += 2 * " "
        if self.upper_cutoff: string += float_format % self.upper_cutoff
        else:                 string += "0"
        string += " "
        if self.lower_cutoff: string += float_format % self.lower_cutoff
        else:                 string += "0"
        string += " "
        if mpw:               string += float_format % mpw
        else:                 string += "0"
        return string

    def _mcnp_unit(self, float_format, sim, name):
        sph = self.spheres[name]
        formatstr = "{0} {0} {0} {0} {0}".format(float_format)
        return formatstr % (tuple(sph.center) + (sph.inrad, sph.outrad))

    def index(self, sph_name):
        return self.spheres.keys().index(sph_name) + 1

    @property
    def upper_cutoff(self): return self._upper_cutoff

    @upper_cutoff.setter
    def upper_cutoff(self, value): self._upper_cutoff = value

    @property
    def lower_cutoff(self): return self._lower_cutoff

    @lower_cutoff.setter
    def lower_cutoff(self, value): self._lower_cutoff = value

    @property
    def min_photon_weight(self): return self._min_photon_weight

    @min_photon_weight.setter
    def min_photon_weight(self, value): self._min_photon_weight = value

    @property
    def spheres(self): return self._spheres

    @spheres.setter
    def spheres(self, value): self._spheres = value


class Vector(IMisc):
    """Position vector. In MCNP, this is the **VECT** card. Unique card with name
    `vector`. Vector are added to the card via the :py:meth:`add`. The card is
    used with the :py:class:`ExponentialTransform` card in MCNP.

    .. inheritance-diagram:: pyne.simplesim.cards.Vector

    """
    # TODO don't require the user to make this card. The reason we're allowing
    # it is vector re-use.
    def __init__(self):
        super(Vector, self).__init__('vector', unique=True)
        self.vectors = collections.OrderedDict()

    def set(self, vecname, vector):
        """Adds a vector with name ``vecname`` to the card. If a vector with that
        vecname has already been added, that vector is then changed to the new
        value. Returns the index
        of the vector.

        Parameters
        ----------
        vecname : str
            Name of the vector. The names 'currdir', 'x', 'X', 'y', 'Y', 'z',
            and 'Z' will cause a conflict with the
            :py:class:`ExponentialTransform`, as they have other meanings on
            that card.
        vector : 3-element list, :py:class:`np.array`, etc. [centimeters]

        Returns
        -------
        index : int
            Index of the vector that ends up on the **VECT** card.

        Examples
        --------
        The following adds two vectors, where we have imported :py:mod:`numpy`
        as ``np``::

            vec = Vector()
            vec.set('origin', [0, 0, 0])
            vec.set('x-axis', np.array([1, 0, 0]))

        The vector 'origin' has index 1 and the vector 'x-axis' has index 2.

        """
        self.vectors[vecname] = vector
        return self.index(vecname)

    def comment(self):
        string = "Vector {0!r}:".format(self.name)
        counter = 0
        for key, val in self.vectors.iteritems():
            counter += 1
            string += " {0}: ({1[0]:g}, {1[1]:g}, {1[2]:g}) cm".format(
                    key, val)
            if counter < len(self.vectors): string += ","
        return string + "."

    def mcnp(self, float_format, sim):
        if len(self.vectors) == 0:
            raise StandardError("No vectors added.")
        string = "VECT"
        counter = 0
        for key, val in self.vectors.iteritems():
            counter += 1
            index = self.index(key)
            formatstr = " V{0} {1} {1} {1}".format(index, float_format)
            string += formatstr % tuple(val)
        return string
        
    def index(self, vecname):
        # MCNP is okay with index 0.
        return self.vectors.keys().index(vecname)


class Burnup(IMisc):
    """Unimplemented."""

    def __init__(self,
                 times=None,
                 power_fracs=None,
                 power=None,
                 materials=None,
                 matl_omit=None,
                 min_frac=None,
                 chain_convergence=None,
                 fp_tier=None,
                 out_order=None,
                 out_per_step=False,
                 model_opt=None,
                 volume=None,
                 conc_change=None):
        # TODO
        pass
    #def __init__(self, times=[1], power_fracs, power=1, materials=None,
    #        matl_omit=None, min_frac=1e-10, chain_convergence=1e-10, fp_tier=1,
    #        out_order='mass', out_per_step=False, model_opt='fatal',
    #        volume=None, conc_change=None)


class Custom(ICard):
    """Allows the user to specify a custom card. It is logical for this class,
    in the future, to have attributes like ``serpent``.
    
    """
    def __init__(self, *args, **kwargs):
        """ 
        Parameters
        ----------
        name : str
            Name of the card. Must be unique within the type of card.
        comment : str, optional
            A comment to be placed in the input file, if comments are
            requested.
        mcnp : str, optional
            The card to be printed in an MCNP input. Must be given for any card
            that will be used in generating an MCNP input. This excludes the
            card number, except for one case with tallies.

        Examples
        --------
        Except for materials, custom cards look like the following::

            cust = SurfaceCustom('customsurf', comment='My Surface', 
                                               mcnp='RPP ...')

        """
        # I am intentionally avoiding using super().
        ICard.__init__(self, *args, **kwargs)
        self.comment_string = kwargs.get('comment', None)
        self.mcnp_string = kwargs.get('mcnp', None)

    def comment(self):
        if not self.comment_string:
            return "Custom card."
        return self.comment_string

    def mcnp(self, float_format, sim, number=None, pre=None):
        """Raises an exception if the ``mcnp`` attribute is not defined."""
        if not self.mcnp_string:
            raise Exception("mcnp not defined on custom {0}.".format(
                self))
        return "{0}{1}{2}".format(
                pre if pre else "",
                "{0} ".format(number) if number else "",
                self.mcnp_string)


class CellCustom(Custom, Cell):
    """Custom :py:class:`Cell` card.
    
    .. inheritance-diagram:: CellCustom

    """
    # Provide number automatically.
    # TODO correct spacing given the number.
    def mcnp(self, float_format, sim):
        return super(CellCustom, self).mcnp(float_format, sim,
                sim.sys.cell_num(self.name))


class SurfaceCustom(Custom, ISurface):
    """Custom :py:class:`ISurface` card.
    
    .. inheritance-diagram:: SurfaceCustom

    """
    # Provide number automatically.
    # TODO correct spacing given the number.
    def mcnp(self, float_format, sim):
        return super(SurfaceCustom, self).mcnp(float_format, sim,
                sim.sys.surface_num(self.name))


class MaterialCustom(Custom, Material):
#class MaterialCustom(Material, Custom):
    """Custom :py:class:`Material` card.

    .. inheritance-diagram:: MaterialCustom

    NOTE The name of the MaterialCustom
    card must be given as a keyword argument, otherwise an error will arise::

        mc = MaterialCustom(name='matc', comment='Made in USA', mcnp='...')

    """
    def __init__(self, mat, name=None, *args, **kwargs):
        name = name if name is not None else  mat.attrs['name']
        self._mat = mat
        super(MaterialCustom, self).__init__(name=name, *args, **kwargs)

    # Provide number automatically.
    # TODO correct spacing given the number.
    def mcnp(self, float_format, sim):
        return super(MaterialCustom, self).mcnp(float_format, sim,
                sim.sys.material_num(self.name), pre="M")


class SourceCustom(Custom, ISource):
    """Custom :py:class:`ISource` card.
    
    .. inheritance-diagram:: SourceCustom
    
    """
    pass


class DistributionCustom(Custom, Distribution):
    """Custom :py:class:`Distribution` card.
    
    .. inheritance-diagram:: DistributionCustom

    """
    def __init__(self, *args, **kwargs):
        """
        Parameters
        ----------
        name : str
            See :py:class:`Custom`.
        comment : str, optional
            See :py:class:`Custom`.
        mcnp : list of 2-element tuples, optional
            Unlike with :py:class:`Custom`, this is a list of 2-element tuples,
            since a distribution in MCNP is done with a number of cards (i.e.
            SI, SP, SB).

        Examples
        --------
        The following is expected::

            dc = DistributionCustom('disti', comment='My dist.',
                    mcnp=[('SP', 'prob string')])
            dc = DistributionCustom('disti', comment='My dist.',
                    mcnp=[('SI', 'info string'), ('SP', 'prob string')])

        """
        super(DistributionCustom, self).__init__(*args, **kwargs)

    # Provide number automatically.
    # TODO correct spacing given the number.
    def mcnp(self, float_format, sim):
        # Use superclass's exception, but not its output.
        super(DistributionCustom, self).mcnp(float_format, sim)
        string = ""
        for pair in self.mcnp_string:
            prestr = pair[0]
            content = pair[1]
            string += "{0}{1} {2}\n".format(
                    prestr,
                    sim.dist_num(self.name),
                    content)

        return string


class TallyCustom(Custom, ITally):
    """Custom :py:class:`ITally` card.
    The card is printed in the order that it
    is added, but unlike with the other custom cards the user must provide
    their own tally number.

    .. inheritance-diagram:: TallyCustom

    """
    def __init__(self, *args, **kwargs):
        """
        Parameters
        ----------
        name : str
            See :py:class:`Custom`.
        comment : str, optional
            See :py:class:`Custom`.
        mcnp : str, optional
            See :py:class:`Custom`.
        mcnp_pre : str, optional
            A string (i.e. '+' or '*') to print before the tally number (i.e.
            'F14').
        tallyclass : :py:class:`ITally` class object, optional
            Since each type of tally is numbered separately, the user may want
            to specify which tally class to use so that the custom tally is
            added in the right place. If this is not provided, then the user
            must provide the tally number in the output string (e.g. in
            ``mcnp``).

        Examples
        --------
        The tallyclass does not need to be provided::

            tal1 = TallyCustom('tal1', comment="", mcnp="F14:N 1 2 3")

        But if it is provided, the class takes care of printing the begining of
        the card::

            tal1 = TallyCustom('tal1', comment="", mcnp=":N 1 2 3",
                    tallyclass=CellFlux)
            tal1 = TallyCustom('tal1', comment="", mcnp=":N 1 2 3",
                    tallyclass=CellFlux)
            tal1 = TallyCustom('tal1', comment="", mcnp=":N 1 2 3",
                    mcnp_pre='*',
                    tallyclass=CellFlux)

        """
        # I am intentionally avoiding using super().
        super(TallyCustom, self).__init__(*args, **kwargs)
        self.mcnp_pre = kwargs.get('mcnp_pre', "")
        self.tallyclass = kwargs.get('tallyclass', None)

    def mcnp(self, float_format, sim):
        if self.tallyclass:
            return super(TallyCustom, self).mcnp(float_format, sim,
                    10 * sim.tally_num(self.name) + self.tallyclass.mcnp_num,
                    pre=(self.mcnp_pre + "F"))
        else:
            return super(TallyCustom, self).mcnp(float_format, sim,
                    pre=self.mcnp_pre)


class MiscCustom(Custom, IMisc):
    """Custom :py:class:`IMisc` card.
    
    .. inheritance-diagram:: MiscCustom

    """
    pass


class TransformationCustom(Custom, Transformation):
    """Custom :py:class:`Transformation` card.

    .. inheritance-diagram:: TransformationCustom

    """
    # Provide number automatically.
    # TODO correct spacing given the number.
    def mcnp(self, float_format, sim):
        return super(TransformationCustom, self).mcnp(float_format, sim,
                sim.transformation_num(self.name), pre="TR")


class Particle(object):
    """Facilitates the use of particle abbreviations, and is an attribute in
    various cards (:py:class:`ITally` cards, and other :py:class:`IMisc` cards.)
    
    """
    # These are provided in the order in which they are listed in the MCNP
    # manual. TODO may need to escape some of these characters?
    mcnp_abbrev = {'neutron': 'N',
                   'anti-neutron': '-N',
                   'photon': 'P',
                   'electron': 'E',
                   'positron': '-E',
                   'muon': '|',
                   'anti-muon': '-|',
                   'tau': '*',
                   'electron_neutrino': 'U',
                   'anti-electron_neutrino': '-U',
                   'muon_neutrino': 'V',
                   'tau_neutrino': 'W',
                   'proton': 'H',
                   'anti-proton': '-H',
                   'lambda': 'L',
                   'sigma+': '+',
                   'sigma-': '-',
                   'cascade0': 'X',
                   'cascade-': 'Y',
                   'omega': 'O',
                   'lambdac+': 'C',
                   'cascadec+': '!',
                   'cascadec0': '?',
                   'lambdab0': '<',
                   'pion+': '/',
                   'pion-': '-/',
                   'neutral_pion': 'Z',
                   'kaon+': 'K',
                   'kaon-': '-K',
                   'K0-short': '%',
                   'K0-long': '^',
                   'D+': 'G',
                   'D0': '@',
                   'Ds+': 'f',
                   'B+': '>',
                   'B0': 'B',
                   'Bs0': 'Q',
                   'deuteron': 'D',
                   'triton': 'T',
                   'helium-3': 'S',
                   'helium-4': 'A',
                   'heavy_ions': '#'
                   }
    """Map of particle names to their abbreviations in MCNP."""

    def __init__(self, name):
        """ 
        Parameters
        ----------
        name : str
            Name of the particle. Acceptable names depend on the use. For
            example, the particles acceptable for MCNP are in
            :py:attr:`mcnp_abbrev`.
            
        """
        self.name = name

    def __str__(self):
        return self.name

    def mcnp(self):
        """ 
        Returns
        -------
        abbrev : str
            MCNP's one or two character abbreviation for a particle.

        """
        if self.name not in self.mcnp_abbrev:
            raise ValueError("Particle {0!r} not recognized.".format(
                    self.name))
        return self.mcnp_abbrev[self.name]

