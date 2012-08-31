import numpy as np

from pyne.simplesim import cards, definition

class IUnit(object):
    def __init__(self, up=None, down=None):
        self.up = up
        self.down = down

    def __lt__(self, next_level_up):
        return self.of(next_level_up)

    def of(self, next_level_up):
        self.up = next_level_up
        self.up.down = self
        if self.down: return self.down
        else:         return self

    def __and__(self, right):
        return self.union(right)

    def union(self, right):
        if isinstance(self, union):
            self.brothers += [right]
            return self
        else:
            return union(self, right)

    def comment(self):
        raise NotImplementedError

    def mcnp(self, float_format, sim, inner):
        return "{0}{1}{2}{3}".format(
                " (" if self.up and not self.down else "",
                inner,
                (" <" + self.up.mcnp(float_format, sim)) if self.up else "",
                ")" if self.down and not self.up else "")


class CellSurf(IUnit):
    def __init__(self, name):
        super(CellSurf, self).__init__()
        self.name = name

    def comment(self):
        return super(CellSurf, self).comment(" {0!r}".format(self.name))


class surf(CellSurf):
    def mcnp(self, float_format, sim):
        return super(surf, self).mcnp(float_format, sim,
                " {0}".format(sim.sys.surface_num(self.name)))


class cell(CellSurf):
    def mcnp(self, float_format, sim, inner=None):
        return super(cell, self).mcnp(float_format, sim,
                " {0}{1}".format(sim.sys.cell_num(self.name),
                    inner if inner else ""))


class ucell(cell):
    def __init__(self, name, lat_spec=None):
        super(ucell, self).__init__(name)
        self.lat_spec = lat_spec

    def mcnp(self, float_format, sim):
        return super(ucell, self).mcnp(float_format, sim,
                self.lat_spec.mcnp(float_format, sim))


class uni(IUnit):
    def __init__(self, name):
        super(uni, self).__init__()
        self.name = name

    def mcnp(self, float_format, sim):
        return super(uni, self).mcnp(float_format, sim,
                " U={0}".format(sim.sys.universe_num(self.name)))


class union(IUnit):
    def __init__(self, *args):
        # all args must be a instance of Unit or its subclasses.
        self.brothers = args
        super(union, self).__init__()

    def mcnp(self, float_format, sim):
        string = ""
        for bro in self.brothers:
            string += bro.mcnp(float_format, sim)
        return super(union, self).mcnp(float_format, sim, 
                " ({0})".format(string))


class vec(IUnit):
    def __init__(self, *args):
        self.sisters = args
        super(vec, self).__init__()

    def mcnp(self, float_format, sim):
        string = ""
        for sis in self.sisters:
            string += sis.mcnp(float_format, sim)
        return super(vec, self).mcnp(float_format, sim, string)


class LatticeSpec(object):
    def mcnp(self, float_format, sim, inner):
        # Lattice specification goes in square brackets.
        return "[{0}]".format(inner)


class lin(LatticeSpec):
    def __init__(self, linear_index):
        self.index = linear_index

    def mcnp(self, float_format, sim):
        return super(lin, self).mcnp(float_format, sim,
                "{0:d}".format(self.index))


class rng(LatticeSpec):
    def __init__(self, x_bounds, y_bounds, z_bounds):
        super(rng, self).__init__()
        self.x_bounds = x_bounds
        self.y_bounds = y_bounds
        self.z_bounds = z_bounds

    def mcnp(self, float_format, sim):
        return super(rng, self).mcnp(float_format, sim,
                "{0[0]:d}:{0[1]:d} {1[0]:d}:{1[1]:d} {2[0]:d}:{2[1]:d}".format(
                self.x_bounds, self.y_bounds, self.z_bounds))


class cor(LatticeSpec):
    def __init__(self, points):
        super(cor, self).__init__()
        # We want a nested list, even if the user doesn't provide it. If the
        # first element of the list is an int, then it's not a 3-element list
        # or numpy array, so we need to nest it in a loop for the methods to
        # work.
        if type(points[0]) is int: points = [points]
        self.points = points

    def mcnp(self, float_format, sim):
        string = ""
        counter = 0
        for pt in self.points:
            counter += 1
            string += " {0[0]:d} {0[1]:d} {0[2]:d}".format(pt)
            if counter < len(self.points): string += ","
        return super(cor, self).mcnp(float_format, sim, string)










