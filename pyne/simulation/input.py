#!/usr/bin/env python

"""This module employs the modules `reactordef` and `material` to generate
plaintext input files for a general code. Support is provided for MCNPX, and
support for Serpent is not complete but should be straightforward. The
extension to other codes may require more effort.

- Write out
- Read in a JSON file input def.
"""

import abc

class InputFile(object):
    """Create a more general name that encompasses reading an input file.
    """
    __metaclass__ = abc.ABCMeta

    def __init__(self):
        pass

    @abc.abstractmethod
    def _cell(self, cell):
        """Returns a cell card string."""
        return

    @abc.abstractmethod
    def _infinite_cylinder_surface(self, surface):
        """Returns an infinite cylinder surface card."""
        return

class JSONInput(InputFile):
    """TODO ideally we can exchange file formats through something like this
    intermediate."""

class MCNPInput(InputFile):
    """Contains a write method for each type of surface.
    """

    def __init__(self, reactordef, simulationdef):
        """

        Parameters
        ----------
        reactordef : pyne.ReactorDefinition

        """

    def _cell(self, cell):
        """Returns a cell card string."""
        return

    def _infinite_cylinder_surface(self, surface):
        """Returns an infinite cylinder surface card."""
        return




class SerpentInput(object):
    """Must find the cell used for a given material, and would need to create
    more than one material if necessary.

    """
    def __init__(self, reactordef, simulationdef):
        pass


