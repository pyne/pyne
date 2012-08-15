#!/usr/bin/env python

"""The :py:mod:`cards` module can be imported as such::

    from pyne.simplesim import cards

The user does not interact with this module. Rather, the user interacts with
the :py:mod:`definition` module. Instances of classes in the :py:mod:`definition` module
contains instances of the classes in this module.  Below is the reference for
the module.

"""

# TODO in the Sphinx documentation, create a developers and a user's version of
# this documentation.
# TODO in the Sphinx documentation, provide a list of all classes without any
# docstrings. autosummary
# autosummary, :nosignatures: pyne.simplesim.cards
# TODO use sphinx domains where possible instead of double single quotes.
# TODO write a development guide next to the usersguide.
# TODO make error messages valuable: give back to the user their input.
# TODO sphinx inheritance diagrams.

import abc

import numpy as np

from pyne import material

class ICard(object):
    """This class is not used by the user. Abstract superclass for all cards.
    All cards have a name property and a comment() method.

    """

    # This line makes this class an abstract base class (ABC).
    __metaclass__ = abc.ABCMeta
    
    def __init__(self, name):
        """
        Parameters
        ----------
        name : str
            Name of this instance of this card, used for referencing this card
            from others and modifying this card after it has been added to a
            definition (see :py:mod:`pyne.simplesim.definition`). There is no
            known restriction on the characters that can be used in the name,
            but the name must be unique.

        """
        self.name = name

    # All subclasses must now define a comment() method.
    @abc.abstractmethod
    def comment(self):
        """All cards define a comment describing the content of the card."""
        raise NotImplementedError

    @property
    def name(self):
        return self._name

    @name.setter
    def name(self, value):
        if value == '':
            raise ValueError("The ``name`` property of the cell cannot "
                    "be empty.")
        self._name = value


class CellVoid(ICard):
    """An empty region of space; this cell does not contain a material."""

    def __init__(self, name, region):
        """
        Parameters
        ----------
        region : :py:class:`Region`
            Defines the region of space that this cell occupies (see
            :py:class:`Region`).

        Examples
        --------
        TODO

        """
        self.region = region

    def comment(self):
       # TODO Walk the region.
       return "Void cell. Comment unimplemented."

    @property
    def region(self):
        return self._region

    @region.setter
    def region(self, obj):
        self._region = obj

   
class Cell(CellVoid):
    """A cell is a region of space filled with a material.

    """
    def __init__(self, name, region, material, density, density_units):
        """
        Parameters
        ----------
        name : str
            See :py:class:`ICard`.
        region : :py:class:`Region`
            See :py:class:`CellVoid`
        material : :py:class:`pyne.material.Material`
            A material definition using the :py:mod:`pyne.material` module.
            For use here, the material's :py:attr:`name` property must be set to something
            other than '' and must be unique. See
            :py:class:`pyne.material.Material`.
        density : float
            Density for the material, in units of density_units.
        density_units : str
            Either 'g/cm^3', or 'atoms/b/cm'.

        Examples
        --------
        TODO
        
        """
        # TODO decide how I will do cross-referencing.

        super(Cell, self).__init__(name, region)
        self.material = material
        self.density = density
        self.density_units = density_units

    def comment(self):
        # TODO walk the region.
        # TODO print material description.
        return "Cell. comment unimplemented."

    @property
    def material(self):
        return self._material

    @material.setter
    def material(self, obj):
        if obj.name == '':
            raise ValueError("The ``name`` property of the material cannot "
                    "be empty.")
        self._material = obj

    @property
    def density(self):
        return self._density

    @density.setter
    def density(self, value):
        self._density = value

    @property
    def density_units(self):
        return self._density_units

    @density_units.setter
    def density_units(self, value):
        if density_units != 'g/cm^3' and density_units != 'atoms/b/cm':
            raise ValueError("The property ``density_units`` must be either "
                    "'g/cm^3' or 'atoms/b/cm'. User provided "
                    "'{0}'".format(value))
        self._density_units = value


class CellVoidMCNP(CellVoid):
    """A cell card with keyword options that are available in MCNP. Thus, it
    only makes sense to use this card if writing an input for MCNP. This is a
    void (no material) cell; see :py:class:`CellMCNP` for the corresponding
    cell filled with a material.
    
    The U, LAT, and FILL keywords are not available; as this functionality
        should be obtained by using Universe and Lattice cards.

    Note this card was written with MCNPX version 2.7 in mind.

    """
    # TODO Sphinx documentation should not list all keyword arguments.

    def __init__(self, name, region,
                 temperature=None, volume=None,
                 neutron_imp=None,
                 photon_imp=None,
                 electron_imp=None,
                 proton_imp=None,
                 proton_weight_lim=None,
                 neutron_exp_transform=None,
                 photon_exp_transform=None,
                 electron_exp_transform=None,
                 proton_exp_transform=None,
                 neutron_force_coll=None,
                 photon_force_coll=None,
                 electron_force_coll=None,
                 proton_force_coll=None,
                 neutron_weight_win_bound=None,
                 photon_weight_win_bound=None,
                 electron_weight_win_bound=None,
                 proton_weight_win_bound=None,
                 neutron_dxtran_contrib=None,
                 photon_dxtran_contrib=None,
                 electron_dxtran_contrib=None,
                 proton_dxtran_contrib=None,
                 fission_turnoff=None,
                 det_contrib=None,
                 transform=None
                 ):
        """
        Parameters
        ----------
        name : str
            See :py:class:`ICard`.
        region : :py:class:`Region`
            See :py:class:`CellVoid`
        temperature : float, otional [Kelvin]
            Temperature of the cell.
        volume : float, optional [cm^3]
            Volume of the cell.
        TODO

        Examples
        --------
        TODO

        """
        # TODO allow use of U, LAT, and FILL keywords?
        super(CellVoidMCNP, self).__init__(name, region)
        # Assign keyword arguments.
        self.temperature = temperature
        self.volume = volume
        self.neutron_imp = neutron_imp
        self.photon_imp = photon_imp
        self.electron_imp = electron_imp
        self.proton_imp = proton_imp

    def comment(self):
        # TODO walk the region.
        # TODO print material description.
        return "Cell. comment unimplemented."

    @property
    def temperature(self):
        return self._temperature

    @temperature.setter
    def temperature(self, value):
        if value < 200:
            raise UserWarning("Temperature set as less than 200 K. Are you "
                    "trying to specify temperature in degrees "
                    "Celcius, etc.? User provided %.4f." % value)
        if value < 1:
            raise UserWarning("Temperature set as less than 1 K. Are you "
                    "trying to specify temperature as 'kT'? "
                    "User provided %.4f." % value)
        self._temperature = value

    @property
    def volume(self):
        return self._volume

    @property
    def neutron_imp(self):
        return self._neutron_imp

    @neutron_imp.setter
    def neutron_imp(self, value):
        self._neutron_imp = value

    @property
    def photon_imp(self):
        return self._photon_imp

    @photon_imp.setter
    def photon_imp(self, value):
        self._photon_imp = value

    @property
    def electron_imp(self):
        return self._electron_imp

    @electron_imp.setter
    def electron_imp(self, value):
        self._electron_imp = value

    @property
    def proton_imp(self):
        return self._proton_imp

    @proton_imp.setter
    def proton_imp(self, value):
        self._proton_imp = value


class CellMCNP(CellVoidMCNP, Cell):
    """A cell card with keyword options that are available in MCNP. Thus, it
    only makes sense to use this card if writing an input for MCNP.    

    The U, LAT, and FILL keywords are not available; as this functionality
        should be obtained by using Universe and Lattice cards.

    Note this card was written with MCNPX version 2.7 in mind.

    """
    # TODO flesh out keyword arguments.
    def __init__(self, name, region, material,
                 temperature=None, volume=None,
                 neutron_imp=None,
                 photon_imp=None,
                 electron_imp=None,
                 proton_imp=None,
                 proton_weight_lim=None,
                 neutron_exp_transform=None,
                 photon_exp_transform=None,
                 electron_exp_transform=None,
                 proton_exp_transform=None,
                 neutron_force_coll=None,
                 photon_force_coll=None,
                 electron_force_coll=None,
                 proton_force_coll=None,
                 neutron_weight_win_bound=None,
                 photon_weight_win_bound=None,
                 electron_weight_win_bound=None,
                 proton_weight_win_bound=None,
                 neutron_dxtran_contrib=None,
                 photon_dxtran_contrib=None,
                 electron_dxtran_contrib=None,
                 proton_dxtran_contrib=None,
                 fission_turnoff=None,
                 det_contrib=None,
                 transform=None
                 ):
        """
        Parameters
        ----------
        name : str
            See :py:class:`ICard`.
        region : :py:class:`Region`
            See :py:class:`CellVoid`
        material : :py:class:`pyne.material.Material`
            See :py:class:`Cell`
        **kwargs : varies
            See :py:class:`CellVoidMCNP`.
            
        Examples
        --------
        TODO

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
        # The following fields are not initialized via the superclass
        # constructor above.
        self.material = material
        self.density = density
        self.density_units = density_units
        return

    def comment(self):
        # TODO walk the region.
        # TODO print material description.
        return "Cell. comment unimplemented."


class IUniverse(ICard):
    """This class is not used by the user. Abstract superclass for all
    universe cards.

    """
    def __init__(self, name):
        pass


class UniverseByRegion(IUniverse):
    """

    """
    # TODO
    def __init__(self, name, region):
        pass


class UniverseByLattice(IUniverse):
    """

    """
    # TODO
    def __init__(self, name, lattice):
        pass


class Lattice(ICard):
    """

    """
    # TODO 
    def __init__(self, name, geom, universe):
        pass


class LatticeByArray(ICard):
    """

    """
    # TODO support of 3D arrays.
    def __init(self, name, geom, xindices, yindices, zindices,
               universe_array):
        pass
        

class ISurface(ICard):
    """This class is not used by the user. Abstract superclass for all
    surface cards.

    The Surface superclass contains properties to set the surface as reflecting
    or white. For codes other than MCNPX, reflecting or white surfaces may be
    specified on a separate boundary condition card (i.e. in Serpent) or may
    not even be available. For other codes, then, the appropriate :py:mod:`inputfile`
    class needs to pick up this information and print the appropriate string to
    the code's input file, or in the latter case return an exception.

    """
    # TODO support rotation.
    __metaclass__ = abc.ABCMeta
    
    def __init__(self, name, reflecting, white):
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
        super(ISurface, self).__init__(name)
        self.reflecting = reflecting
        self.white = white 
        if self.reflecting and self.white:
            raise ValueError("The user set the surface to be reflecting AND "
                    "white, but can only be neither or one of the two.")

    @abc.abstractmethod
    def comment(self):
        raise NotImplementedError

    @abc.abstractmethod
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
        raise NotImplementedError

    @abc.abstractmethod
    def stretch(self, vector):
        """Stretches (scales) the surface from the origin. This is an abstract
        method, and must be defined by each surface card. Stretches may not be
        permitted for all surfaces in all directions, and in such cases an
        exception is raised.

        Parameters
        ----------
        vector : 3-element list or :py:class`np.array`, float [unitless]
            The elements specify a stretch in the x, y, and z directions, in
            this order. A zero in any of the directions indicates that no
            stretch is done in that direction.

        Examples
        --------
        Both of the following lines stretch the surface along the y axis by a
        factor of 2. The x and z directions are unaffected::

            surf.stretch([0, 2, 0])
            surf.stretch(np.array([0, 2, 0]))
        
        """
        raise NotImplementedError

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
        :py:mod:`RegionLeaf` is on the side of the surface with a positive
        sense.
        
        """
        return RegionLeaf(self, True)
    
    @property
    def reflecting(self):
        return self._reflecting

    @reflecting.setter
    def reflecting(self, value):
        if value is not None and type(value) is not bool:
            raise TypeError("The property ``reflecting`` must be "
                    "None or of boolean type. User provided "
                    "{0}.".format(value))
        self._reflecting = value

    @property
    def white(self):
        return self._white

    @white.setter
    def white(self, value):
        if value is not None and type(value) is not bool:
            raise TypeError("The property ``white`` must be "
                    "None or of boolean type. User provided "
                    "{0}.".format(value))
        self._white = value


class IAxisSurface(ISurface):
    """This class is not used by the user. Abstract superclass for all simple
    axis-aligned surfaces. Accordingly, such classes share the cartesian_axis
    property.
    
    """
    __metaclass__ = abc.ABCMeta
    def __init__(self, name, cartesian_axis, reflecting, white):
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
        super(IAxisSurface, self).__init__(name, reflecting, white)
        self.cartesian_axis = cartesian_axis

    @abc.abstractmethod
    def comment(self):
        raise NotImplementedError
    
    @abc.abstractmethod
    def shift(self, vector):
        """See :py:meth:`ISurface.shift`."""
        raise NotImplementedError

    @abc.abstractmethod
    def stretch(self, vector):
        """See :py:meth:`ISurface.stretch`."""
        raise NotImplementedError

    @property
    def cartesian_axis(self):
        return self._cartesian_axis

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
    """Cylinder aligned with and centered on one of the Cartesian axes."""

    # TODO if this is shifted, then it becomes not an axis-cylinder.
    def __init__(self, name, cartesian_axis, radius,
                 reflecting=None, white=None):
        """
        Parameters
        ----------
        name : str
            See :py:class:`ICard`.
        cartesian_axis : str
            The axis with which the cylinder is aligned and centered.
            See :py:class`IAxisSurface`.
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
        super(AxisCylinder, self).__init__(name, cartesian_axis, 
                                           reflecting, white)
        self.radius = radius

    def comment(self):
        return ("Axis cylinder %s: aligned and centered on %s axis, "
                "with radius %.4f cm (diameter %.4f cm)." %
                        (self.name, self.cartesian_axis,
                        self.radius, 2 * self.radius))
    
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
            raise ValueError("A cylinder aligned with the %s axis cannot "
                    "be shifted in the %s or %s directions." % dirs)

    def stretch(self, vector):
        """See :py:meth:`ISurface:stretch`. Axis cylinders can be stretched in
        the direction of their axis, which has no effect (permitted in case
        this surface is part of a region that is being stretched), or can be
        stretched `uniformly` in the plane perpendicular to its axis.
        
        Examples
        --------
        The following stretches are okay for a cylinder aligned with the x axis
        (where we have imported :py:mod`numpy` as ``np``)::
            
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
                    "provided %s stretch %.4f and %s stretch %.4f for a "
                    "%s-aligned cylinder." % out)

    @property
    def radius(self):
        return self._radius

    @radius.setter
    def radius(self, value):
        if value <= 0:
            raise ValueError("The ``radius`` property must be "
                    "positive. User provided %.4f." % value)
        self._radius = value


class AxisPlane(IAxisSurface):
    """
    .. inheritance-diagram:: pyne.simplesim.cards.AxisCylinder
    
    Plane perpendicular to one of the Cartesian axes.
    
    
    """

    def __init__(self, name, cartesian_axis, position,
                 reflecting=None, white=None):
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
        super(AxisPlane, self).__init__(name, cartesian_axis,
                                    reflecting, white)
        self.position = position
    
    def comment(self):
        return "Axis plane %s: %s = %.4f cm" % (self.name, self.cartesian_axis,
                self.position)
        pass

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
        if self.cartesian_axis == 'x':
            self.position += vector[0]
        elif self.cartesian_axis == 'y':
            self.position += vector[1]
        elif self.cartesian_axis == 'z':
            self.position += vector[2]

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
    def position(self):
        return self._position

    @position.setter
    def position(self, value):
        self._position = value


class IMacrobody(ISurface):
    """This class is not used by the user. Abstract superclass for all
    macrobody cards. Macrobodies are an MCNP concept.

    """
    def __init__(self, name, reflecting, white):
        """

        """
        super(IMacrobody, self).__init__(name, reflecting, white)

    @abc.abstractmethod
    def comment(self):
        raise NotImplementedError

class Parallelepiped(IMacrobody):
    """

    """
    def __init__(self, name, xmin, xmax, ymin, ymax, zmin, zmax,
                 reflecting=False, white=False):
        """

        Examples
        --------

        """
        super(Parallelepiped, self).__init__(name, reflecting, white)
        self.xlims = np.array([xmin, xmax])
        self.ylims = np.array([ymin, ymax])
        self.zlims = np.array([zmin, zmax])

    def comment(self):
        return

    @property
    def xlims(self):
        return self._xlims

    @xlims.setter
    def xlims(self, value):
        if value[0] > value[1]:
            raise ValueError("The value of xmin, %.4f, is greater than "
                    "that of xmax, %.4f." % (value[0], value[1]))
        self._xlims = value

    @property
    def ylims(self):
        return self._ylims

    @ylims.setter
    def ylims(self, value):
        if value[0] > value[1]:
            raise ValueError("The value of ymin, %.4f, is greater than "
                    "that of ymax, %.4f." % (value[0], value[1]))
        self._ylims = value

    @property
    def zlims(self):
        return self._zlims

    @zlims.setter
    def zlims(self, value):
        if value[0] > value[1]:
            raise ValueError("The value of zmin, %.4f, is greater than "
                    "that of zmax, %.4f." % (value[0], value[1]))
        self._zlims = value


class Cuboid(Parallelepiped):
    """Same exact thing as a Parallelepiped. This class is provided because the
    name is shorter, and thus may be preferred by those who fancy brevity.

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
    """
    """

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
        # TODO this is probably not okay by the proponents of duck typing.
        if type(value) is not bool:
            raise TypeError("User provided a value for pos_sense that is "
                    "not of boolean type.")
        self._pos_sense = value



class Option(ICard):
    """All cards classified as data cards in MCNP, except for materials."""

    def __init__(self, name):
        """
        """
        return

    @abc.abstractmethod
    def comment(self):
        raise NotImplementedError


class Source(ICard):
    """ """
    def __init__(self, name):
        pass




