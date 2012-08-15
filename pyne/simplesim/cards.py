#!/usr/bin/env python

"""The ``cards`` module can be imported as such::

    from pyne.simulation import cards

The user does not interact with this module. Rather, the user interacts with
the ``definition`` module. Instances of classes in the ``definition`` module
contains instances of the classes in this module.  Below is the reference for
the module.

"""

import abc

from pyne import material

class ICard(object):
    """Abstract superclass for all cards. All cards have a name and number
    property, as well as a comment method.
    
    """
    __metaclass__ = abc.ABCMeta
    
    def __init__(self, name):
        """

        """
        self.name = name

    @abc.abstractmethod
    def comment(self):
        return

    @property
    def name(self):
        return self._name

    @name.setter
    def name(self, value):
        self._name = value

class CellVoid(ICard):
    """
    """
    def __init__(self, name, region):
        """
        """
        self.region(region)

   def comment(self):
       return

   @property
   def region(self):
       return self._region

   
class Cell(CellVoid):
    """
    """

    def __init__(self, name, region, material, density, density_units):
        """
        """
        super(Cell, self).__init__(name, region)
        self.material = material
        self.density = density
        self.density_units = density_units

    def comment(self):
        return

    @property
    def material(self):
        return self._material

    @property
    def density(self):
        """Density for the material, in units of self.density_units."""
        return self._density

    @property
    def density_units(self):
        """Either 'g/cm^3', or 'atoms/b/cm'."""
        return self._density_units


class CellSimpleVoid(CellVoid):
    """
    """

    def __init__(self, name, neg_surfs, pos_surfs):
        """There must be thought given to how default arguments are treated. It
        seems to make more sense that there are no kwargs in this module, that
        they are all in ``definition``, since ``definition`` contains the
        methods that the user calls.

        """
        # TODO this could be an overloaded constructor to CellVoid().
        # Create a region using the negative and positive surfaces.
        super(CellSimpleVoid, self).__init__(name, region)
        return

    def comment(self):
        return

#    @property
#    def neg_surface_cards(self):
#        """Surfaces with a negative sense for this cell."""
#        return self._neg_surface_cards
#
#    @property
#    def pos_surface_cards(self):
#        """Surfaces with a positive sense for this cell."""
#        return self._pos_surface_cards

class CellSimple(Cell):
    """
    """
    def __init__(self, name, neg_surfs, pos_surfs, material, density,
            density_units):
        # Create region from neg and pos surfaces.
        super(CellSimple, self).__init__(name, region,
                                         material, density, density_units)
        return


class CellSimpleVoidMCNP(CellSimpleVoid):
    """
    """

    def __init__(self, name, neg_surfs, pos_surfs,
                 temperature=None, volume=None,
                 neutron_imp=None, photon_imp=None, electron_imp=None,
                 proton_imp=None):
        super(CellSimpleVoidMCNP, self).__init__(name, neg_surfs,
                pos_surfs)
        # Assign keyword arguments.
        self.temperature = temperature
        self.volume = volume
        self.neutron_imp = neutron_imp
        self.photon_imp = photon_imp
        self.electron_imp = electron_imp
        self.proton_imp = proton_imp

    @property
    def temperature(self):
        """Temperature of the cell, in Kelvin."""
        return self._temperature

    @property
    def volume(self):
        """Volume of the cell, in cm^3. Optional."""
        return self._volume

    @property
    def neutron_imp(self):
        return self._neutron_imp

    @property
    def photon_imp(self):
        return self._photon_imp

    @property
    def electron_imp(self):
        return self._electron_imp

    @property
    def proton_imp(self):
        return self._proton_imp


class CellSimpleMCNP(CellSimpleVoidMCNP, CellSimple):
    """
    """

    def __init__(self, name, neg_surfs, pos_surfs, material,
                 temperature=None, volume=None,
                 neutron_imp=None, photon_imp=None, electron_imp=None,
                 proton_imp=None):
        """
        """
        # Based on Python's Method Resolution Order (MRO), the constructor for
        # CellSimpleVoidMCNP is called because it is listed first above.
        super(CellSimpleMCNP, self).__init__(name, neg_surfs,
                                             pos_surfs,
                                             temperature=temperature,
                                             volume=volume,
                                             neutron_imp=neutron_imp,
                                             photon_imp=photon_imp,
                                             electron_imp=electron_imp,
                                             proton_imp=proton_imp)
        # The following fields are not initialized via the constructor above.
        self.material = material
        self.density = density
        self.density_units = density_units
        return

    def comment(self):
        return


class IUniverse(Card):
    """

    """
    def __init__(self, name):
        pass


class IUniverseByRegion(IUniverse):
    """

    """
    def __init__(self, name, region):
        pass


class IUniverseByLattice(IUniverse):
    """

    """
    def __init__(self, name, lattice):
        pass


class Lattice(Card):
    """

    """
    def __init__(self, name, geom, universe):
        pass


class LatticeByArray(Card):
    """

    """
    # TODO support of 3D arrays.
    def __init(self, name, geom, xindices, yindices, zindices,
               universe_array):
        pass
        

class ISurface(ICard):
    """
    """
    # TODO support rotation.
    __metaclass__ = abc.ABCMeta
    
    def __init__(self, name, reflecting, white):
        """TODO The Surface superclass contains fields for reflecting and white
        surfaces. For codes other than MCNPX, reflecting or white surfaces may
        be specified on a separate boundary condition card (Serpent) or may not
        even be available. For other codes, then, the appropriate ``inputfile``
        class needs to pick up this information and print the appropriate
        string to the code's input file, or in the latter case return an
        exception.

        """
        super(ISurface, self).__init__(name)
        self.reflecting = reflecting
        self.white = white 

    @abc.abstractmethod
    def comment(self):
        return

    @property
    def neg(self):
        """
        Returns
        -------
        region : cards.Region
        """
        return RegionLeaf(self, False)

    @property
    def pos(self):
        """
        Returns
        -------
        region : cards.Region

        """
        return RegionLeaf(self, True)

    @abc.abstractmethod
    def shift(self, vector):
        """Shifts the surface."""
        return

    @abc.abstractmethod
    def stretch(self, vector):
        """Stretches the surface."""
        return
    
    @property
    def reflecting(self):
        return self._reflecting

    @reflecting.setter
    def reflecting(self, value):
        if value is not None and type(value) is not bool:
            raise TypeError("'reflecting' keyword argument must be "
                    "None or of boolean type.");
        self._reflecting = value

    @property
    def white(self):
        return self._white

    @white.setter
    def white(self, value):
        if value is not None and type(value) is not bool:
            raise TypeError("'white' keyword argument must be "
                    "None or of boolean type.");
        self._white = value


class IAxisSurface(ISurface):
    """Superclass for simple axis-aligned surfaces. Accordingly,
    such classes share the cartesian_axis property.
    
    """
    __metaclass__ = abc.ABCMeta
    def __init__(self, name, cartesian_axis, reflecting, white):
        """
        Parameters
        ----------

        """
        super(AxisSurface, self).__init__(name, reflecting, white)
        self.cartesian_axis = cartesian_axis

    @abc.abstractmethod
    def comment(self):
        return
    
    @abc.abstractmethod
    def shift(self, vector):
        """Shifts the surface."""
        return

    @abc.abstractmethod
    def stretch(self, vector):
        """Stretches the surface."""
        return

    @property
    def cartesian_axis(self):
        return self._cartesian_axis

    @cartesian_axis.setter
    def cartesian_axis(self, value):
        if type(value) is not str:
            raise ValueError("AxisCylinder's cartesian_axis property must be "
                    "a string.")
        if (value.lower() != 'x' and
                value.lower() != 'y' and 
                value.lower() != 'z')
            raise ValueError("AxisCylinder's cartesian_axis property must be "
                    "'x', 'X', 'y', 'Y', or 'z', 'Z'.")
        self._cartesian_axis = value


class AxisCylinder(IAxisSurface):
    """

    """
    # TODO if this is shifted, then it becomes not an axis-cylinder.
    def __init__(self, name, cartesian_axis, radius,
                 reflecting=None, white=None):
        """

        """
        super(AxisCylinder, self).__init__(name, cartesian_axis, 
                                           reflecting, white)
        self.radius = radius

    def comment(self):
        # TODO
        return
    
    def shift(self, vector):
        """Shifts the surface."""
        return

    def stretch(self, vector):
        """Stretches the surface."""
        return

    @property
    def radius(self):
        return self._radius

    @radius.setter
    def radius(self, value):
        if value <= 0
            raise ValueError("AxisCylinder's radius property must be "
                    "positive.")
        self._radius = value


class Plane(IAxisSurface):
    """

    """
    def __init__(self, name, cartesian_axis, position,
                 reflecting=None, white=None):
        """

        """
        super(Plane, self).__init__(name, cartesian_axis,
                                    reflecting, white)
        self.position = position

    @property
    def position(self):
        return self._position

    @position.setter
    def position(self, value):
        self._position = value


class Macrobody(ISurface):
    pass


class Parallelepiped(Macrobody):
    def __init__(self, name, xmin, xmax, ymin, ymax, zmin, zmax,
            reflecting=False, white=False):


class Cuboid(Parallelepiped):
    """
    """
    def __init__(self, name, xmin, xmax, ymin, ymax, zmin, zmax,
            reflecting=False, white=False):
        """
        """
        super(Cuboid, self).__init__(name, xmin, xmax, ymin, ymax, zmin, zmax,
                                     reflecting, white)

class Region(object):
    """Represents a volume (space) confined by unions and intersections of
    surfaces."""
    # TODO transformation functions
    # Cell cards are then formed by a region and a material.

    # TODO Complement functionality can be added by overloading the
    # __not__ operator and defining a complement boolean property that is set
    # by the __not__ operator.
    # TODO add transformation methods.
    def __and__(self, arg):
        return self.intersect(arg)

    def __or__(self, arg):
        return self.union(arg)

    def intersect(self, arg):
        return RegionAnd(self, arg)

    def union(self, arg):
        return RegionOr(self, arg)

class IRegionBool(Region):
    """Abstract class; should have no instances of this."""
    def __init__(self, left_region, right_region):
        self.left_region = left_region
        self.right_region = right_region

    @property
    def left_region(self):
        return self._left_region

    @left_region.setter
    def left_region(self, value):
        self._left_region = value

    @property
    def right_region(self):
        return self._right_region

    @right_region.setter
    def right_region(self, value):
        self._right_region = value


class RegionAnd(IRegionBool):
    pass


class RegionOr(IRegionBool):
    pass


class RegionLeaf(Region):

    def __init__(self, surface, pos_sense):
        self.surface = surface
        self.pos_sense = pos_sense

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
        if type(value) is not bool:
            raise TypeError("User specified a value for pos_sense that is "
                    "not of boolean type.")
        self._pos_sense = value



class Option(ICard):
    """All cards classified as data cards in MCNP, except for materials."""

    def __init__(self, name):
        """
        """
        return

    def comment(self):
        """This should be abstract as well."""
        return


class Source(ICard):
    """ """
    def __init__(self, name):
        pass




