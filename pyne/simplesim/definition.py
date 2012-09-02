#!/usr/bin/env python

"""The ``definition`` module can be imported as such::

    from pyne.simplesim import definition

********
Overview
********

There are primarily two types of
definitions: system definitions and simulation definitions:

    - System definitions describe the geometry and material composition of a
      system, and nothing else. It is intended that the definition is
      completely code-dependent (i.e. independent of codes such as MCNP),
      though the way in which geometry is specified is very similar to how it
      is done for MCNP.
    - Simulation definitions, which require a system definition, contain
      information about how to simulate the system: where to place particle
      sources, what information is desired from the simulation, etc. The
      simulation definition may have different information depending on the
      code being used to perform the simulation (e.g. MCNP, etc.), though it is
      likely that these codes do have some common elements. The common elements
      are to be found/placed in :py:class:`SimulationDefinition`, while the
      elements relevant to the abilities of specific codes would be found in
      classes such as :py:class:`MCNPSimulation`.


***********
Limitations
***********
- Saving and loading simulations is not implemented, though some preliminary
  work with JSON encoding has been performed. At the bottm of this module is a
  subclass of :py:class:`json.JSONEncoder` that is intended to aid the
  reading/writing of simulations.
- MCNP is the only code for which simulations have been explored.


An inheritance diagram of all the classes in this module can be found at
:ref:`pyne_simplesim_inheritance`.

.. moduleauthor:: Chris Dembia <cld72@cornell.edu>


*********
Reference
*********

"""
# TODO check overwriting warning.
# TODO test the exceptions for getting *_num() not in list/dict.
import abc
import collections
import pickle
import json
import warnings

import numpy as np

from pyne import material
from pyne.simplesim import cards

class IDefinition(object):
    """This class is not used by the user. Abstract base class for
    system and simulation definitions.

    """
    __metaclass__ = abc.ABCMeta

    def __init__(self, fname=None, verbose=True):
        """Creates a new definition or loads one from ``fname``.

        Parameters
        ----------
        fname : str, optional
            If not provided, then a new definition is created. Otherwise, an
            existing definition is loaded from file ``fname``.
        verbose : bool, optional
            Prints updates about what the class is doing.

        """
        self.verbose = verbose
        if fname is not None:
            if self.verbose:
                print "Opening definition stored in {0!r}.".format(fname)
            self._open(fname)
        else:
            if self.verbose:
                print "Creating a new definition."
            self._create_new()

    @abc.abstractmethod
    def _create_new(self):
        """Definition started from scratch. Initialize all fields. """
        raise NotImplementedError

    @abc.abstractmethod
    def save(self, fname):
        """Save object data to a JSON file. Unimplemented."""
        fid = open(fname, 'w')
        #pickle.dump(self, fid)
        fid.write(DefinitionEncoder(indent=3).encode(self))
        fid.close()

    @abc.abstractmethod
    def _open(self, fname):
        """Open object data from a JSON file. Unimplemented."""
        # TODO Ensure the file exists.
        fid = open(fname)
        # Process.
        #self = pickle.load(fname)
        fid.close()

    def _assert_unique(self, card_type, card):
        """Checks if the name on a card has already been used for another
        card. If so, a warning is issued.

        Parameters
        ----------
        card_type : str
            String identifying a type/category of card.
        card : :py:class:`pyne.simplesim.cards.ICard` or subclass

        Returns
        -------
        unique : bool
            True if card with the name is not in the simulation, otherwise
            False.

        """
        # TODO Don't require the user to pass the card_type string; just get
        # the class of the object passed in.
        if card_type == 'cell':
            dict_to_check = self.cells
        elif card_type == 'surface':
            dict_to_check = self.surfaces
        elif card_type == 'material':
            dict_to_check = self.materials
        elif card_type == 'source':
            dict_to_check = self.source
        elif card_type == 'tally':
            dict_to_check = self.tally
        elif card_type == 'misc':
            dict_to_check = self.misc
        elif card_type == 'transformation':
            dict_to_check = self.transformations
        elif card_type == 'dist':
            dict_to_check = self.dists
        else:
            raise ValueError("The input ``card_type`` must be either "
                    "'cell', 'surface', 'material', 'source', "
                    "'tally', or 'misc'.")
        if card.name in dict_to_check and self.verbose:
            warnings.warn("Card {0!r}, type {1!r} is already part of the "
                    "definition; overwriting.".format(card.name, card_type))
            #raise UserWarning("Card {0!r}, type {0!r} is already part of the "
            #        "definition; overwriting.".format(card.name, card_type))
            #raise Exception("The %s name %s has already been used for "
            #        "another %s" % (card_type, card.name, card_type))
            return False
        return True

    @property
    def verbose(self): return self._verbose

    @verbose.setter
    def verbose(self, value): self._verbose = value


class SystemDefinition(IDefinition):
    """This class creates a system definition as is done in MCNP: homogeneous
    regions in space in the reactor, called cells, are defined through the
    intersection, union, etc of surfaces and are filled by materials.

    """
    def __init__(self, **kwargs):
        """Creates a new reactor definition or loads one from a JSON file. See
        :py:class:`IDefinition`.
        
        """
        super(SystemDefinition, self).__init__(**kwargs)

    def _create_new(self):
        self._surfaces = collections.OrderedDict()
        self._materials = collections.OrderedDict()
        self._cells = collections.OrderedDict()
        self._universes = list()

    def add_cell(self, *args):
        """Adds a cell to the system, and also adds the materials and surfaces
        in those cells to the system. This method can also be used to replace a
        cell that has already been given to the system. If this cell is part of
        a universe, that universe is also registered with the system.

        Parameters
        ----------
        cell : :py:class:`pyne.simplesim.cards.Cell` or subclass
        *args : cell, ...
            Any number of cells.

        Returns
        -------
        num : int
            If only one card is added, the method returns the card number
            of that card. Be careful; the card number may change if other cards
            are deleted. May be useful when defining unique cards.

        Examples
        --------
        This card is used as follows::

            sys = SystemDefinition()
            sys.add_cell(cellA)
            sys.add_cell(cellB, cellC)

        """
        for cell in args:
            # Warn the user if a cell with this name is already in the system.
            # already.
            self._assert_unique('cell', cell)
            # Custom card; just add and move on.
            if isinstance(cell, cards.Custom):
                self._cells[cell.name] = cell
                return
            if self.verbose:
                print "Adding cell %s." % cell.name
            # Add all surfaces that aren't already added. Do this by walking the
            # region tree and calling _add_unique_surfaces() at the leaves.
            cell.region.walk(self._add_unique_surfaces)
            # Only add the material if this is not a void cell and if it doesn't
            # already exist.
            if (cell.material and cell.material.name not in self.materials):
                self.add_material(cell.material)
            # Constituent surfaces and material have been added, so we can added
            # the cell.
            # Add universes to system.
            if cell.universe:
                self._register_universe(cell.universe[0])
            self._cells[cell.name] = cell
        if len(args) == 1:
            return self.cell_num(args[0].name)

    def add_surface(self, *args):
        """This method is only used by the user for surfaces that are not on a
        cell card. Surfaces on a cell card are added automatically. Duplication
        is checked by card name, therefore surface names should be unique
        unless the user intends to replace a surface.

        Parameters
        ----------
        surface : :py:class:`pyne.simplesim.cards.ISurface` subclass
        *args : surface, ...
            Any number of surfaces.

        Returns
        -------
        num : int
            If only one card is added, the method returns the card number
            of that card. Be careful; the card number may change if other cards
            are deleted. May be useful when defining unique cards.

        Examples
        --------
        This card is used as follows::

            sys = SystemDefinition()
            sys.add_surface(surfA)
            sys.add_surface(surfB, surfC)

        """
        for surface in args:
            self._assert_unique('surface', surface)
            self._surfaces[surface.name] = surface
        if len(args) == 1:
            return self.surface_num(args[0].name)

    def _add_unique_surfaces(self, regionleaf):
        # Used by add_cell().
        name = regionleaf.surface.name
        if self.verbose:
            print "Trying to add surface {0!r}...".format(name)
        if name not in self.surfaces:
            self.add_surface(regionleaf.surface)
            if self.verbose:
                print "  Surface {0!r} added successfully.".format(name)
        else:
            if self.verbose:
                print ("  Surface {0!r} already exists in the "
                        "definition.".format(name))

    def add_material(self, *args):
        """This method is only used by the user for materials that are not on a
        cell card.  Materials on cell cards are added automatically.
        Alternatively, this method can be used to replace an existing material.
        Note that this must be a :py:class:`pyne.simplesim.cards.Material`, not
        just a :py:class:`pyne.material.Material`.

        Parameters
        ----------
        material : :py:class:`pyne.simplesim.cards.Material` subclass
        *args : material, ...
            Any number of materials.

        Returns
        -------
        num : int
            If only one card is added, the method returns the card number
            of that card. Be careful; the card number may change if other cards
            are deleted. May be useful when defining unique cards.

        Examples
        --------
        This card is used as follows::

            sys = SystemDefinition()
            sys.add_material(matA)
            sys.add_material(matB, matC)

        """
        for material in args:
            if material.name == None or material.name == '':
                raise ValueError("The ``name`` property of the material cannot "
                        "be empty.")
            self._assert_unique('material', material)
            self._materials[material.name] = material
        if len(args) == 1:
            return self.material_num(args[0].name)

    def _register_universe(self, univ_name):
        """This method is not used by the user. It is called by
        :py:meth:`cards.Cell.mcnp` and :py:meth:`cards.Universes._mcnp_unit` to
        register universe names with the system definition.

        """
        #self._assert_unique('universe', univ_name)
        if univ_name not in self.universes:
            self._universes += [univ_name]

    def cell_num(self, name):
        """Returns the cell number for the cell name provided. If a cell with
        that name is not in the system, an exception is raised.

        """
        # TODO removing cards.
        # Must add one because indices start at 0 but card numbers at 1.
        if name not in self.cells:
            raise StandardError("Cell {0!r} is not in the system.".format(
                    name))
        return self.cells.keys().index(name) + 1

    def surface_num(self, name):
        """Returns the surface number for the surface name provided. If a
        surface with that name is not in the system, an exception is raised.

        """
        # Must add one because indices start at 0 but card numbers at 1.
        if name not in self.surfaces:
            raise StandardError("Surface {0!r} is not in the system.".format(
                    name))
        num = self.surfaces.keys().index(name) + 1
        if isinstance(self.surfaces[name], cards.Facet):
            # This does not work, as self.surfaces is the macrobody with the
            # same name, not the facet.
            return num + self.surfaces[name].number / 10
        return num

    def material_num(self, name):
        """Returns the material number for the material name provided. If a
        material with that name is not in the system, an exception is raised.

        """
        # Must add one because indices start at 0 but card numbers at 1.
        if name not in self.materials:
            raise StandardError("Material {0!r} is not in the system.".format(
                    name))
        return self.materials.keys().index(name) + 1

    def universe_num(self, name):
        """Returns the universe number for the universe name provided. If a
        universe with that name is not in the system, an exception is raised.

        """
        if name not in self.universes:
            raise StandardError("Universe {0!r} is not in the system.".format(
                    name))
        return self.universes.index(name) + 1

    def remove_cell(self, name):
        """Removes a cell from the system. No dependency check is performed.

        Parameters
        ----------
        name : str
            Name of the cell to remove.

        Returns
        -------
        cell : :py:class:`pyne.simplesim.cards.Cell`
            The cell that has been removed from the system.

        """
        # TODO Perform a dependency check.
        return self._cells.pop(name)

    def remove_surface(self, name):
        """Removes a surface from the system. No dependency check is performed.

        Parameters
        ----------
        name : str
            Name of the surface to remove.

        Returns
        -------
        surface : :py:class:`pyne.simplesim.cards.ISurface`
            The surface that has been removed from the system.

        """
        # TODO Perform a dependency check.
        return self._surfaces.pop(name)

    def remove_material(self, name):
        """Removes a material from the system. No dependency check is performed.

        Parameters
        ----------
        name : str
            Name of the material to remove.

        Returns
        -------
        material : :py:class:`pyne.simplesim.cards.Material`
            The material that has been removed from the system.

        """
        # TODO Perform a dependency check.
        return self._materials.pop(name)

    def remove_universe(self, name):
        """Removes a universe from the system. No dependency check is performed.

        Parameters
        ----------
        name : str
            Name of the universe to remove.

        Returns
        -------
        universe : str
            The name of the universe that has been removed from the system.

        """
        # TODO Perform a dependency check.
        return self._universes.pop(self.universes.index(name))

    def save(self, fname):
        super(SystemDefinition, self).save(fname)

    def _open(self, fname):
        super(SystemDefinition, self)._open(fname)

    @property
    def surfaces(self):
        """Ordered dictionary of surfaces (:py:class:`cards.ISurface`)."""
        return self._surfaces

    @property
    def materials(self):
        """Ordered dictionary of materials
        (:py:class:`pyne.simplesim.cards.Material`).
        
        """
        return self._materials

    @property
    def cells(self):
        """Ordered dictionary of cells (:py:class:`pyne.simplesim.cards.Cell`).
        
        """
        return self._cells

    @property
    def universes(self):
        """List of universe names (see :py:class:`cards.Universes`)."""
        # TODO would become an ordered dictionary when Universes become their
        # own cards.
        return self._universes


class SimulationDefinition(IDefinition):
    """This definition contains source, source distribution, tally/detector,
    and miscellaneous cards. Code-specific functionality/cards is found in
    subclasses of this class.

    """
    def __init__(self, systemdef, **kwargs):
        """Creates a new definition or loads one from a file.
        
        Parameters
        ----------
        systemdef: :py:class:`SystemDefinition`
            A system definition for which this simulation is defined. This
            system definition is stored in :py:attr:`sys`, from where it can be
            modified.
        fname : str, optional
            See :py:class:`IDefinition`.
        verbose : bool, optional
            See :py:class:`IDefinition`.
        
        """

        super(SimulationDefinition, self).__init__(**kwargs)
        # TODO when saving a simulation def, the the system def should also be
        # saved, and so this assignment that happens on this next line
        # shouldn't happen this way.
        self.sys = systemdef

    def _create_new(self):
        """Initializes attributes/properties."""
        self._source = collections.OrderedDict()
        self._tally = collections.OrderedDict()
        self._misc = collections.OrderedDict()
        self._dists = collections.OrderedDict()

    def add_source(self, *args):
        """Adds source cards to the simulation, or replaces existing ones if
        the added card has the same name as a source card in the simulation.
        
        Parameters
        ----------
        source : :py:class:`pyne.simplesim.cards.ISource` or subclass
        *args : source, ...
            Any number of source cards.

        Examples
        --------
        This card is used as follows::

            sys = SystemDefinition()
            sim = SimulationDefinition(sys)
            sim.add_source(sourceA)
            sim.add_source(sourceB, sourceC)

        """
        for card in args:
            if not isinstance(card, cards.ISource):
                raise ValueError("Only cards subclassed from ``ISource`` can be "
                        "added by this method. User provided {0}.".format(card))
            self._assert_unique('source', card)
            self._source[card.name] = card

    def add_tally(self, *args):
        """Adds tally cards to the simulation, or replaces existing ones if
        the added card has the same name as a tally card in the simulation.
        

        Parameters
        ----------
        tally : :py:class:`pyne.simplesim.cards.ITally` or subclass
        *args : tally, ...
            Any number of tally cards.

        Examples
        --------
        This card is used as follows::

            sys = SystemDefinition()
            sim = SimulationDefinition(sys)
            sim.add_tally(tallyA)
            sim.add_tally(tallyB, tallyC)

        """
        for card in args:
            # TODO check cells and surfaces? only if moving to string refs.
            if not isinstance(card, cards.ITally):
                raise ValueError("Only cards subclassed from ``ITally`` can be "
                        "added by this method. User provided {0}.".format(card))
            self._assert_unique('tally', card)
            self._tally[card.name] = card

    def add_misc(self, *args):
        """Adds miscellaneous cards to the simulation, or replaces existing
        ones if the added card has the same name as a misc card in the
        simulation.

        Parameters
        ----------
        misc : :py:class:`pyne.simplesim.cards.IMisc` or subclass
        *args : misc, ...
            Any number of misc cards.

        Examples
        --------
        This card is used as follows::

            sys = SystemDefinition()
            sim = SimulationDefinition(sys)
            sim.add_misc(miscA)
            sim.add_misc(miscB, miscC)

        """
        for card in args:
            # TODO check references to tallies? only if moving to string refs.
            if not isinstance(card, cards.IMisc):
                raise ValueError("Only cards subclassed from ``IMisc`` can be "
                        "added by this method. User provided {0}.".format(card))
            # Add universes to system.
            if isinstance(card, cards.Universes):
                for uni in card.univ_names.values():
                    self.sys._register_universe(uni)
            self._assert_unique('misc', card)
            self._misc[card.name] = card

    def add_dist(self, *args):
        """Adds distribution cards to the simulation, or replaces existing ones if the added
        card has the same name as a distribution card in the simulation.

        Parameters
        ----------
        dist : :py:class:`cards.Distribution` or subclass
            The card to be added to the simulation.
        *args : dist, ...
            Any number of such cards.

        Returns
        -------
        num : int
            If only one card is added, the method returns the card number
            of that card. Be careful; the card number may change if other cards
            are deleted. May be useful when defining unique cards.

        Examples
        --------
        This card is used as follows::

            sys = SystemDefinition()
            sim = SimulationDefinition(sys)
            sim.add_dist(distA)
            sim.add_dist(distB, distC)

        """
        for card in args:
            if not isinstance(card, cards.Distribution):
                raise ValueError("Only ``Distribution`` cards can be "
                        "added by this method. User provided {0}.".format(card))
            self._assert_unique('dist', card)
            self._dists[card.name] = card
        if len(args) == 1:
            return self.dist_num(args[0].name)

    def dist_num(self, name):
        """Retrieve the number of a :py:class:`cards.Distribution` card in
        the MCNP input file.

        Parameters
        ----------
        name : str
            Name of the :py:class:`cards.Distribution`. Names must be unique
            across all distributions, and perhaps across cell and surface
            cards.

        Returns
        -------
        dist_num : int
            The distribution number.

        """
        if name not in self.dists:
            raise StandardError("Distribution {0!r} is not in "
                    "the simulation.".format(name))
        return self.dists.keys().index(name) + 1

    def remove_source(self, name):
        """Removes a source card from the system. No dependency check is
        performed.

        Parameters
        ----------
        name : str
            Name of the source card to remove.

        Returns
        -------
        source : :py:class:`pyne.simplesim.cards.ISource`
            The source card that has been removed from the system.

        """
        # TODO Perform a dependency check.
        return self._source.pop(name)

    def remove_tally(self, name):
        """Removes a tally card from the system. No dependency check is
        performed.

        Parameters
        ----------
        name : str
            Name of the tally card to remove.

        Returns
        -------
        tally : :py:class:`pyne.simplesim.cards.ITally`
            The tally card that has been removed from the system.

        """
        # TODO Perform a dependency check.
        return self._tally.pop(name)

    def remove_misc(self, name):
        """Removes a misc card from the system. No dependency check is
        performed.

        Parameters
        ----------
        name : str
            Name of the misc card to remove.

        Returns
        -------
        misc : :py:class:`pyne.simplesim.cards.IMisc`
            The misc card that has been removed from the system.

        """
        # TODO Perform a dependency check.
        return self._misc.pop(name)

    def remove_dist(self, name):
        """Removes a distribution card from the system. No dependency check is
        performed.

        Parameters
        ----------
        name : str
            Name of the distribution card to remove.

        Returns
        -------
        dist : :py:class:`pyne.simplesim.cards.Distribution`
            The distribution card that has been removed from the system.

        """
        # TODO Perform a dependency check.
        return self._dists.pop(name)

    def save(self, fname):
        super(SimulationDefinition, self).save(fname)

    def _open(self, fname):
        super(SimulationDefinition, self)._save(fname)

    @property
    def sys(self): return self._sys

    @sys.setter
    def sys(self, value): self._sys = value

    @property
    def source(self):
        """Ordered dictionary of source cards (from :py:class:`cards.ISource`)

        """
        return self._source

    @property
    def tally(self):
        """Ordered dictionary of tallies (from :py:class:`cards.ITally`)."""
        return self._tally

    @property
    def misc(self):
        """Ordered dictionary of misc. cards (from :py:class:`cards.IMisc`)."""
        return self._misc

    @property
    def dists(self):
        """Ordered dictionary of distribution cards (from
        :py:class:`cards.Distribution`).

        """
        return self._dists


class MCNPSimulation(SimulationDefinition):
    """Adds transformation and more tally functionality to
    :py:class:`SimulationDefinition`. See the constructor there.

    """
    def _create_new(self):
        """Initializes any attributes/properties."""
        super(MCNPSimulation, self)._create_new()
        self._transformations = collections.OrderedDict()
        self._tally_surfacecurrent = collections.OrderedDict()
        self._tally_surfaceflux = collections.OrderedDict()
        self._tally_cellflux = collections.OrderedDict()
        self._tally_cellenergydep = collections.OrderedDict()
        self._tally_cellfissiondep = collections.OrderedDict()
        self._tally_pulseheight = collections.OrderedDict()
        self._tally_detector = collections.OrderedDict()

    def add_tally(self, *args):
        """Adds tally cards to the simulation.

        Parameters
        ----------
        card : :py:class:`cards.ICard` or subclass
            The tally card to be added to the simulation.
        *args :
            Any number of such cards.

        """
        for card in args:
            # Can move this line out of the loop.
            super(MCNPSimulation, self).add_tally(card)
            if isinstance(card, cards.SurfaceCurrent):
                self._tally_surfacecurrent[card.name] = card
            elif isinstance(card, cards.SurfaceFlux):
                self._tally_surfaceflux[card.name] = card
            elif isinstance(card, cards.CellFlux):
                self._tally_cellflux[card.name] = card
            elif isinstance(card, cards.CellEnergyDeposition):
                self._tally_cellenergydep[card.name] = card
            elif isinstance(card, cards.CellFissionEnergyDeposition):
                self._tally_cellfissiondep[card.name] = card
            elif isinstance(card, (cards.CellPulseHeight,
                    cards.CellChargeDeposition)):
                self._tally_pulseheight[card.name] = card
            elif isinstance(card, cards.IDetector):
                self._tally_detector[card.name] = card
            elif isinstance(card, cards.Custom):
                # See if class attribute is not None.
                if card.tallyclass:
                    if issubclass(card.tallyclass, cards.SurfaceCurrent):
                        self._tally_surfacecurrent[card.name] = card
                    elif issubclass(card.tallyclass, cards.SurfaceFlux):
                        self._tally_surfaceflux[card.name] = card
                    elif issubclass(card.tallyclass, cards.CellFlux):
                        self._tally_cellflux[card.name] = card
                    elif issubclass(card.tallyclass, cards.CellEnergyDeposition):
                        self._tally_cellenergydep[card.name] = card
                    elif issubclass(card.tallyclass, cards.CellFissionEnergyDeposition):
                        self._tally_cellfissiondep[card.name] = card
                    elif issubclass(card.tallyclass, (cards.CellPulseHeight,
                            cards.CellChargeDeposition)):
                        self._tally_pulseheight[card.name] = card
                    elif issubclass(card.tallyclass, cards.IDetector):
                        self._tally_detector[card.name] = card
                    else:
                        raise Exception("Unrecognized tally card {0}.".format(card))
            else:
                raise Exception("Unrecognized tally card {0}.".format(card))

    def tally_num(self, name):
        """Retrieve the number of a :py:class:`cards.ITally` card in the MCNP
        input file. Note that :py:class:`pyne.simplesim.cards.TallyCustom`
        cards with :py:attr:`tallyclass` as None do not have a tally number.

        Parameters
        ----------
        name : str
            Name of the tally. Names must be unique across all tally types.

        Returns
        -------
        tally_num : int
            The tally number for a given tally type. For the first cell flux
            tally, a value of 1 is returned, not 14.

        """
        card = self.tally[name]
        # Must add one because indices start at 0 but card numbers at 1.
        if isinstance(card, cards.SurfaceCurrent):
            return self._tally_surfacecurrent.keys().index(name) + 1
        elif isinstance(card, cards.SurfaceFlux):
            return self._tally_surfaceflux.keys().index(name) + 1
        elif isinstance(card, cards.CellFlux):
            return self._tally_cellflux.keys().index(name) + 1
        elif isinstance(card, cards.CellEnergyDeposition):
            return self._tally_cellenergydep.keys().index(name) + 1
        elif isinstance(card, cards.CellFissionEnergyDeposition):
            return self._tally_cellfissiondep.keys().index(name) + 1
        elif isinstance(card, (cards.CellPulseHeight,
                cards.CellChargeDeposition)):
            return self._tally_pulseheight.keys().index(name) + 1
        elif isinstance(card, cards.IDetector):
            return self._tally_detector.keys().index(name) + 1
        elif isinstance(card, cards.Custom):
            # See if class attribute is not None.
            if card.tallyclass:
                if issubclass(card.tallyclass, cards.SurfaceCurrent):
                    return self._tally_surfacecurrent.keys().index(name) + 1
                elif issubclass(card.tallyclass, cards.SurfaceFlux):
                    return self._tally_surfaceflux.keys().index(name) + 1
                elif issubclass(card.tallyclass, cards.CellFlux):
                    return self._tally_cellflux.keys().index(name) + 1
                elif issubclass(card.tallyclass, cards.CellEnergyDeposition):
                    return self._tally_cellenergydep.keys().index(name) + 1
                elif issubclass(card.tallyclass, cards.CellFissionEnergyDeposition):
                    return self._tally_cellfissiondep.keys().index(name) + 1
                elif issubclass(card.tallyclass, (cards.CellPulseHeight,
                        cards.CellChargeDeposition)):
                    return self._tally_pulseheight.keys().index(name) + 1
                elif issubclass(card.tallyclass, cards.IDetector):
                    return self._tally_detector.keys().index(name) + 1
                else:
                    raise Exception("Unrecognized tally card {0}.".format(card))
        else:
            raise Exception("Unrecognized tally type. name {0}.".format(name))

    def remove_tally(self, name):
        """See :py:meth:`SimulationDefinition.remove_tally`."""
        card = self.tally[name]
        # Must add one because indices start at 0 but card numbers at 1.
        if isinstance(card, cards.SurfaceCurrent):
            self._tally_surfacecurrent.pop(name)
        elif isinstance(card, cards.SurfaceFlux):
            self._tally_surfaceflux.pop(name)
        elif isinstance(card, cards.CellFlux):
            self._tally_cellflux.pop(name)
        elif isinstance(card, cards.CellEnergyDeposition):
            self._tally_cellenergydep.pop(name)
        elif isinstance(card, cards.CellFissionEnergyDeposition):
            self._tally_cellfissiondep.pop(name)
        elif isinstance(card, (cards.CellPulseHeight,
                cards.CellChargeDeposition)):
            self._tally_pulseheight.pop(name)
        elif isinstance(card, cards.IDetector):
            self._tally_detector.pop(name)
        elif isinstance(card, cards.Custom):
            # See if class attribute is not None.
            if card.tallyclass:
                if issubclass(card.tallyclass, cards.SurfaceCurrent):
                    self._tally_surfacecurrent.pop(name)
                elif issubclass(card.tallyclass, cards.SurfaceFlux):
                    self._tally_surfaceflux.pop(name)
                elif issubclass(card.tallyclass, cards.CellFlux):
                    self._tally_cellflux.pop(name)
                elif issubclass(card.tallyclass, cards.CellEnergyDeposition):
                    self._tally_cellenergydep.pop(name)
                elif issubclass(card.tallyclass, cards.CellFissionEnergyDeposition):
                    self._tally_cellfissiondep.pop(name)
                elif issubclass(card.tallyclass, (cards.CellPulseHeight,
                        cards.CellChargeDeposition)):
                    self._tally_pulseheight.pop(name)
                elif issubclass(card.tallyclass, cards.IDetector):
                    self._tally_detector.pop(name)
        return super(MCNPSimulation, self).remove_tally(name)

    def add_transformation(self, *args):
        """Adds transformation cards to the simulation.

        Parameters
        ----------
        card : :py:class:`cards.Transformation` or subclass
            The card to be added to the simulation.
        *args : card, ...
            Any number of such cards.

        Returns
        -------
        num : int
            If only one card is added, the method returns the card number
            of that card. Be careful; the card number may change if other cards
            are deleted. May be useful when defining unique cards.

        """
        for card in args:
            if not isinstance(card, cards.Transformation):
                raise ValueError("Only ``Transformation``s can be "
                        "added by this method. User provided {0}.".format(card))
            self._assert_unique('transformation', card)
            self._transformations[card.name] = card
        if len(args) == 1:
            return self.transformation_num(args[0].name)

    def transformation_num(self, name):
        """Retrieve the number of a :py:class:`cards.Transformation` card in
        the MCNP input file.

        Parameters
        ----------
        name : str
            Name of the :py:class:`cards.Transformation`. Names must be unique
            across all transformations.

        Returns
        -------
        transformation_num : int
            The transformation number.

        """
        if name not in self.transformations:
            raise StandardError("Transformation {0!r} is not in "
                    "the simulation.".format(name))
        return self.transformations.keys().index(name) + 1

    def remove_transformation(self, name):
        """Removes a transformation card from the system. No dependency check is
        performed.

        Parameters
        ----------
        name : str
            Name of the transformation card to remove.

        Returns
        -------
        transformation : :py:class:`pyne.simplesim.cards.Transformation`
            The transformation card that has been removed from the system.

        """
        # TODO Perform a dependency check.
        return self._transformations.pop(name)

    @property
    def transformations(self):
        """Ordered dictionary of transformation cards (from
        :py:class:`cards.Transformation`).

        """
        return self._transformations


class DefinitionEncoder(json.JSONEncoder):
    """Used to encode definitions into JSON files. Unimplemented."""
    # TODO circular reference issue with Regions.

    def __init__(self, *args, **kwargs):
        """ """
        super(DefinitionEncoder, self).__init__(*args, **kwargs)

    def default(self, obj):
        try:
            if isinstance(obj, cards.ISurface):
                # TODO something like this to get around the circular
                # reference.
                return obj.name
            if isinstance(obj, material.Material):
                # Issue with stlconverters.
                return repr(obj)
            if hasattr(obj, '__dict__'):
                mydict = obj.__dict__
                for key, val in mydict.items():
                    if isinstance(val, np.ndarray):
                        mydict[key] = val.tolist()
                return mydict
            if isinstance(obj, np.ndarray):
                return obj.tolist()
            #if isinstance(obj, cards.Cell) or issubclass(obj, cards.Cell):
            #    print "cell"
            #    return ''
            return json.JSONEncoder.default(self, obj)
        except:
            print "exception: "
            print type(obj)
            print obj
            raise

