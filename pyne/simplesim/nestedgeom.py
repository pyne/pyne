#!/usr/bin/env python

"""The classes in this module aid the specification of tally `units`  in
:py:class:`pyne.simplesim.cards.ICellSurfTally` (see examples below). The class
focuses on the ability to specify cells and surfaces that are nested within
universes or other cells in nearly arbitrarily complex ways. Apart from making
user input clearler, this class has methods that help with creating cell
comments and mcnp strings. 

The names of classes here are fairly short, so the user may not want to import
the entire namespace. A suggested way to import the module is::

    import pyne.simplesim.nestedgeom as ng


An inheritance diagram of all the classes in this module can be found at
:ref:`pyne_simplesim_inheritance`.


Usage Examples
--------------
The subclasses of :py:class:`ICellSurfTally` can take as input any subclass of
:py:class:`IUnit`. The following shows different ways in which complex units
can be built. Here it is assumed that there is a cell named '1' in the system,
a surface named '1' in the system, a universe named '1' in the system, etc. The
unions represent averaging or totaling the units, depending on the type of
tally.

+-----------+-------------------------------------------------------------+------------------------------------------------------------+-----------------------+
|category   |description                                                  |nestedgeom                                                  |mcnp                   |
+===========+=============================================================+============================================================+=======================+
|basic      |surface 1 (s1)                                               |``s1 = Surf('1')``                                          |1                      |
|           +-------------------------------------------------------------+------------------------------------------------------------+-----------------------+
|           |cell 1 (c1)                                                  |``c1 = Cell('1')``                                          |1                      |
|           +-------------------------------------------------------------+------------------------------------------------------------+-----------------------+
|           |universe 1 (u1)                                              |``u1 = Univ('1')``                                          |U=1                    |
+-----------+-------------------------------------------------------------+------------------------------------------------------------+-----------------------+
|union      |union of s1 and s2                                           |``Union(s1, s2)``                                           |(1 2)                  |
|           |                                                             +------------------------------------------------------------+                       |
|           |                                                             |``(s1 | s2)``                                               |                       |
|           +-------------------------------------------------------------+------------------------------------------------------------+-----------------------+
|           |union of c1, c2, and c3                                      |``Union(c1, c2, c3)``                                       |(1 2 3)                |
|           |                                                             +------------------------------------------------------------+                       |
|           |                                                             |``(c1 | c2 | c3)``                                          |                       |
|           +-------------------------------------------------------------+------------------------------------------------------------+-----------------------+
|           |union of cells in u1                                         |``Union(u1)``                                               |(U=1)                  |
+-----------+-------------------------------------------------------------+------------------------------------------------------------+-----------------------+
|nesting    |s1 in filled cell 2 (fc2)                                    |``fc2 = FCell('2'): s1 < fc2``                              |(1 < 1)                |
|           |                                                             +------------------------------------------------------------+                       |
|           |                                                             |``s1.of(fc2)``                                              |                       |
|           +-------------------------------------------------------------+------------------------------------------------------------+-----------------------+
|           |c1 in fc3 in fc3                                             |``c1 < fc2 < fc3``                                          |(1 < 2 < 3)            |
|           |                                                             +------------------------------------------------------------+                       |
|           |                                                             |``c1.of( fc2.of(fc3) )``                                    |                       |
|           +-------------------------------------------------------------+------------------------------------------------------------+-----------------------+
|           |s1 in union of fc2 and fc3                                   |``s1 < Union(fc2, fc3)``                                    |(1 < (2 3))            |
|           |                                                             +------------------------------------------------------------+                       |
|           |                                                             |``s1.of( Union(fc2, fc3) )``                                |                       |
|           +-------------------------------------------------------------+------------------------------------------------------------+-----------------------+
|           |s1 in (union of fc2 and fc3) in fc4                          |``s1 < Union(fc2, fc3) < fc4``                              |(1 < (2 3) < 4)        |
|           |                                                             +------------------------------------------------------------+                       |
|           |                                                             |``s1.of( Union(fc2, fc3).of(fc4) )``                        |                       |
|           +-------------------------------------------------------------+------------------------------------------------------------+-----------------------+
|           |s1 in u1 in fc4                                              |``s1 < u1 < fc4``                                           |(1 < U=1 < 4)          |
|           |                                                             +------------------------------------------------------------+                       |
|           |                                                             |``s1.of( u1.of(fc4) )``                                     |                       |
+-----------+-------------------------------------------------------------+------------------------------------------------------------+-----------------------+
|vectorized |(over s1 and s2) in fc2                                      |``ng.Vec(s1, s2) < fc2``                                    |(1 2 < 2)              |
|           |                                                             +------------------------------------------------------------+                       |
|           |                                                             |``ng.Vec(s1, s2).of(fc2)``                                  |                       |
|           |                                                             +------------------------------------------------------------+                       |
|           |                                                             |``(s1 & s2) < fc2``                                         |                       |
|           +-------------------------------------------------------------+------------------------------------------------------------+-----------------------+
|           |s1 in (over fc2 and fc3)                                     |``s1 < ng.Vec(fc2, fc3)``                                   |(1 < 2 3)              |
|           |                                                             +------------------------------------------------------------+                       |
|           |                                                             |``s1 < (fc2 & fc3)``                                        |                       |
+-----------+-------------------------------------------------------------+------------------------------------------------------------+-----------------------+
|lattice    |s1 in (union of fc2 and (fc3 in lattice elem 5)) in fc4      |``la = Lin(5): s1 < (fc2 | FCell('3', la)) < fc4``          |(1 < (2 3[5]) < 4      |
|           |                                                             +------------------------------------------------------------+                       |
|           |                                                             |``s1 < (fc2 | fc3.lat(la)) < fc4``                          |                       |
|           +-------------------------------------------------------------+------------------------------------------------------------+-----------------------+
|           |s1 in (fc2 in lattice x range 0:1, y range 0:2, z range 0:3) |``s1 < fc2.lat(Rng([0,1], [0,2], [0,3])``                   |(1 < 2[0:1 0:2 0:3])   |
|           +-------------------------------------------------------------+------------------------------------------------------------+-----------------------+
|           |s1 in (fc2 in lattice coord (0, 1, 2))                       |``s1 < fc2.lat(Cor([0, 1, 2]))``                            |(1 < 2[0 1 2])         |
|           +-------------------------------------------------------------+------------------------------------------------------------+-----------------------+
|           |s1 in (fc2 in lattice coords (0, 1, 2), (3, 2, 1))           |``s1 < fc2.lat( Cor([ [0, 1, 2], [3, 2, 1]]) )``            |(1 < 2[0 1 2, 3 2 1])  |
+-----------+-------------------------------------------------------------+------------------------------------------------------------+-----------------------+


"""

import numpy as np

class IUnit(object):
    """Abstract base class for tally units. The user does not use this class
    directly.
    
    """
    # TODO need to define a 'recursive' __repr__.
    def __init__(self, up=None, down=None):
        """Currently, the two keyword arguments are not actually used, but are
        sometimes set directly, after initialization.

        """
        self.up = up
        self.down = down

    def __lt__(self, next_level_up):
        """The user can use the < operator to perform the operation of
        :py:meth:`of`, e.g. (cells/surfs) < (cells/univ)
        
        """
        return self.of(next_level_up)

    def of(self, next_level_up):
        """Returns a unit where ``self`` must be in ``next_level_up`` to be
        tallied. ``self`` must be a cell or surface, or `vector` or `union`
        thereof (cell/surf), and
        ``next_level_up`` can be a cell or universe, or union thereof
        (cells/univ).

        """
        self.up = next_level_up
        self.up.down = self
        if self.down: return self.down
        else:         return self

    def __or__(self, right):
        """A convenient way to call :py:meth:`union`."""
        return self.union(right)

    def union(self, right):
        """Returns a :py:class:`Union` of ``self`` with ``right``."""
        if isinstance(self, Union):
            self.brothers += [right]
            return self
        else:
            return Union(self, right)

    def __and__(self, right):
        """A convenient way to call :py:meth:`vector`."""
        return self.vector(right)

    def vector(self, right):
        """Returns a :py:class:`Vec` of ``self`` with ``right``."""
        if isinstance(self, Vec):
            self.sisters += [right]
            return self
        else:
            return Vec(self, right)

    def comment(self, inner):
        return "{0}{1}{2}{3}".format(
                " (" if self.up and not self.down else "",
                inner,
                (" in" + self.up.comment()) if self.up else "",
                ")" if self.down and not self.up else "")

    def mcnp(self, float_format, sim, inner):
        return "{0}{1}{2}{3}".format(
                " (" if self.up and not self.down else "",
                inner,
                (" <" + self.up.mcnp(float_format, sim)) if self.up else "",
                ")" if self.down and not self.up else "")


class ICellSurf(IUnit):
    """Abstract base class for surfaces and cells in the lowest level of the
    nested geometry. The user directly uses the subclasses of this class. For
    cells in higher levels of nesting (e.g. closer to the real world), see
    :py:class:`FCell`.

    """
    def __init__(self, name):
        """
        Parameters
        ----------
        name : str
            Name of the surface or cell. Depending on the subclass, the name is
            looked up appropriately in the system definition to obtain the
            surface or cell number.

        """
        super(ICellSurf, self).__init__()
        self.name = name


class Surf(ICellSurf):
    """The user uses this class directly to reference a surface in the system
    definition for tallying.

    .. inheritance-diagram:: pyne.simplesim.nestedgeom.Surf

    """
    def comment(self):
        return super(Surf, self).comment(" surf {0!r}".format(self.name))

    def mcnp(self, float_format, sim):
        return super(Surf, self).mcnp(float_format, sim,
                " {0}".format(sim.sys.surface_num(self.name)))


class Cell(ICellSurf):
    """The user uses this class directly to reference a cell in the system
    definition for tallying.

    .. inheritance-diagram:: pyne.simplesim.nestedgeom.Cell

    """
    def comment(self, inner=None):
        return super(Cell, self).comment(" cell {0!r}{1}".format(self.name,
                inner if inner else ""))

    def mcnp(self, float_format, sim, inner=None):
        return super(Cell, self).mcnp(float_format, sim,
                " {0}{1}".format(sim.sys.cell_num(self.name),
                    inner if inner else ""))


class FCell(Cell):
    """This is subclassed from :py:class:`Cell`. Its name stands for filled
    cell. It is to be used for higher-level cells (closer to the real world),
    and has an additional attribute to specify specific lattice elements from
    this cell if it is a lattice cell.

    .. inheritance-diagram:: pyne.simplesim.nestedgeom.FCell

    """
    def __init__(self, name, lat_spec=None):
        """
        Parameters
        ----------
        name : str
            See :py:class:`ICellSurf`.
        lat_sec : subclass of :py:class:`LatticeSpec`

        """
        super(FCell, self).__init__(name)
        self.lat_spec = lat_spec

    def lat(self, lat_spec):
        self.lat_spec = lat_spec
        return self

    def comment(self):
        return super(FCell, self).comment(
                self.lat_spec.comment() if self.lat_spec else "")

    def mcnp(self, float_format, sim):
        return super(FCell, self).mcnp(float_format, sim,
                self.lat_spec.mcnp(float_format, sim) if self.lat_spec else "")


class Univ(IUnit):
    """A universe. It is used in higher levels of nesting (closer to the real
    world). The class, given the name of a universe in the system, looks up the
    appropriate universe number.

    .. inheritance-diagram:: pyne.simplesim.nestedgeom.Univ

    """
    def __init__(self, name):
        """
        Parameters
        ----------
        name : str
            Name of the universe in the system definition.  

        """
        super(Univ, self).__init__()
        self.name = name

    def comment(self):
        return super(Univ, self).comment(" univ {0!r}".format(self.name))

    def mcnp(self, float_format, sim):
        return super(Univ, self).mcnp(float_format, sim,
                " U={0}".format(sim.sys.universe_num(self.name)))


class Union(IUnit):
    """A union of surfaces, cells, or of a single universe.

    .. inheritance-diagram:: pyne.simplesim.nestedgeom.Union

    """
    def __init__(self, *args):
        """
        Parameters
        -----------
        *args : list of instances of :py:class:`IUnit` subclasses.

        """
        # all args must be a instance of Unit or its subclasses.
        self.brothers = args
        super(Union, self).__init__()

    def comment(self):
        string = ""
        counter = 0
        for bro in self.brothers:
            counter += 1
            string += bro.comment()
            if counter < len(self.brothers): string += ","
        return super(Union, self).comment(" union of ({0})".format(string))

    def mcnp(self, float_format, sim):
        string = ""
        for bro in self.brothers:
            string += bro.mcnp(float_format, sim)
        return super(Union, self).mcnp(float_format, sim, 
                " ({0})".format(string))


class Vec(IUnit):
    """A "vector" of surfaces or cells. This class named after the vectorized
    notation that can be used in MATLAB, as this class allows the specification
    of multiple units of input in a single unit. Typically, a :py:class:`Vec`
    is the first or last element in a nested unit. The following are two
    example usages of this class, and their equivalent in comment as two
    separate units::

        # Surf('A') < FCell('C')     Surf('B') < FCell('C')
        Vec(Surf('A'), Surf('B')) < FCell('C')
        # Surf('A') < FCell('C')     Surf('A') < FCell('D')
        Surf('A') < Vec(FCell('C'), FCell('D'))
    """
    # Named after matlab's vectorized notation
    def __init__(self, *args):
        # all args must be a instance of Unit or its subclasses.
        self.sisters = args
        super(Vec, self).__init__()

    def comment(self):
        string = ""
        counter = 0
        for sis in self.sisters:
            counter += 1
            string += sis.comment()
            if counter < len(self.sisters): string += ","
        return super(Vec, self).comment(" over ({0})".format(string))

    def mcnp(self, float_format, sim):
        string = ""
        for sis in self.sisters:
            string += sis.mcnp(float_format, sim)
        return super(Vec, self).mcnp(float_format, sim, string)


class ILatticeSpec(object):
    """Abstract base class for lattice element specifiers. The user does not
    use this class directly. There are 3 subclasses:

    - :py:class:`Lin` : a single linear index of a lattice element.
    - :py:class:`Rng` : x, y, and z index range of lattice elements.
    - :py:class:`Cor` : list coordinates of lattice elments.
    
    """
    def comment(self, inner):
        return "-lat {0}".format(inner)

    def mcnp(self, float_format, sim, inner):
        # Lattice specification goes in square brackets.
        return "[{0}]".format(inner)


class Lin(ILatticeSpec):
    """A single linear index of a lattice element.

    .. inheritance-diagram:: pyne.simplesim.nestedgeom.Lin

    """
    def __init__(self, linear_index):
        """
        Parameters
        ----------
        linear_index : int
            Linear (1-D) of a lattice element for the lattice cell that this
            specifier becomes a part of.

        """
        self.index = linear_index

    def comment(self):
        return super(Lin, self).comment("linear idx {0:d}".format(self.index))

    def mcnp(self, float_format, sim):
        return super(Lin, self).mcnp(float_format, sim,
                "{0:d}".format(self.index))


class Rng(ILatticeSpec):
    """A range of lattice elements.

    .. inheritance-diagram:: pyne.simplesim.nestedgeom.Rng

    """
    def __init__(self, x_bounds=[0, 0], y_bounds=[0, 0], z_bounds=[0, 0]):
        """
        Parameters
        ----------
        x_bounds : 2-element list of int, optional
            Something like [0,5] for lattice elements i=0 through i=5. First
            element is less than the second element. Don't
            specify to leave as [0, 0].
        y_bounds : 2-element list of int, optional
            First element is less than the second element. Don't specify to
            leave as [0, 0].
        z_bounds : 2-element list of int, optional
            First element is less than the second element. Don't specify to
            leave as [0, 0].

        Examples
        --------
        In the following, there is no y dimension::

            myrange = Rng([0, 5], z_bounds=[0, 2])

        """
        super(Rng, self).__init__()
        self.x_bounds = x_bounds
        self.y_bounds = y_bounds
        self.z_bounds = z_bounds

    def comment(self):
        return super(Rng, self).comment(
                "x range {0[0]:d}:{0[1]:d}, y range {1[0]:d}:{1[1]:d}, "
                "z range {2[0]:d}:{2[1]:d}".format(
                self.x_bounds, self.y_bounds, self.z_bounds))

    def mcnp(self, float_format, sim):
        return super(Rng, self).mcnp(float_format, sim,
                "{0[0]:d}:{0[1]:d} {1[0]:d}:{1[1]:d} {2[0]:d}:{2[1]:d}".format(
                self.x_bounds, self.y_bounds, self.z_bounds))


class Cor(ILatticeSpec):
    """A list of lattice element coordinates (in indices).

    .. inheritance-diagram:: pyne.simplesim.nestedgeom.Cor

    """
    def __init__(self, points=[0, 0, 0]):
        """
        Parameters
        ----------
        points : 3-element list of int, list of lists, optional
            Coordinates of a lattice element, or a list of coordinates.

        Examples
        --------
        The following work::

            latspec = Cor()
            latspec = Cor([1, 2, 3])
            latspec = Cor([ [1, 2, 3], [-1, 3, -2]])

        """
        super(Cor, self).__init__()
        # We want a nested list, even if the user doesn't provide it. If the
        # first element of the list is an int, then it's not a 3-element list
        # or numpy array, so we need to nest it in a loop for the methods to
        # work.
        if type(points[0]) is int: points = [points]
        self.points = points

    def comment(self):
        string = "coords"
        counter = 0
        for pt in self.points:
            counter += 1
            string += " ({0[0]:d}, {0[1]:d}, {0[2]:d})".format(pt)
            if counter < len(self.points): string += ","
        return super(Cor, self).comment(string)

    def mcnp(self, float_format, sim):
        string = ""
        counter = 0
        for pt in self.points:
            counter += 1
            string += " {0[0]:d} {0[1]:d} {0[2]:d}".format(pt)
            if counter < len(self.points): string += ","
        return super(Cor, self).mcnp(float_format, sim, string)



# # "raw" docstring for table above.
# ====
# MCNP
# ====
# 1
# 1
# U=1
# 
# (1 2)
# (1 2 3)
# (U=1)
# 
# (1 < 1)
# (1 < 2 < 3)
# (1 < (2 3))
# (1 < (2 3) < 4)
# 
# (1 2 < 2)
# (1 < 2 3)
# 
# 1 < (2 3[5]) < 4
# 1 < 2[0:1 0:2 0:3]
# 1 < 2[0 1 2]
# 1 < 2[0 1 2, 3 2 1]
# 
# ===========
# description
# ===========
# surface 1 (s1)
# cell 1 (c1)
# universe 1 (u1)
# 
# union of s1 and s2
# union of c1, c2, and c3
# union of cells in u1
# 
# s1 in filled cell 2 (fc2)
# c1 in fc3 in fc3
# s1 in union of fc2 and fc3
# s1 in (union of fc2 and fc3) in fc4
# s1 in u1 in fc4
# 
# (over s1 and s2) in fc2
# s1 in (over fc2 and fc3)
# 
# s1 in (union of fc2 and (fc3 in lattice elem 5)) in fc4
# s1 in (fc2 in lattice x range 0:1, y range 0:2, z range 0:3)
# s1 in (fc2 in lattice coord (0, 1, 2))
# s1 in (fc2 in lattice coords (0, 1, 2), (3, 2, 1))
# 
# ==========
# nestedgeom
# ==========
# s1 = Surf('1')
# c1 = Cell('1')
# u1 = Univ('1')
# 
# Union(s1, s2)  /  (s1 | s2)
# Union(c1, c2, c3)  /  (c1 | c2 | c3)
# Union(u1)
# 
# fc2 = FCell('2'): s1 < fc2  /  s1.of(fc2)
# c1 < fc2 < fc3  /  c1.of( fc2.of(fc3) )
# s1 < Union(fc2, fc3)  /  s1.of( Union(fc2, fc3) )
# s1 < Union(fc2, fc3) < fc4  /  s1.of( Union(fc2, fc3).of(fc4) )
# s1 < u1 < fc4  /  s1.of( u1.of(fc4) )
# 
# ng.Vec(s1, s2) < fc2  /  ng.Vec(s1, s2).of(fc2)  /  (s1 & s2) < fc2
# s1 < ng.Vec(fc2, fc3)  /  s1 < (fc2 & fc3)
# 
# la = Lin(5): s1 < (fc2 | FCell('3', la)) < fc4  /  s1 < (fc2 | fc3.lat(la)) < fc4
# s1 < fc2.lat(Rng([0,1], [0,2], [0,3])
# s1 < fc2.lat(Cor([0, 1, 2]))
# s1 < fc2.lat( Cor([ [0, 1, 2], [3, 2, 1]]) )







