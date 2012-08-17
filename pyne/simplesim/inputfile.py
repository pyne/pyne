#!/usr/bin/env python

"""The ``inputfile`` module can be imported as such::

    from pyne.simplesim import inputfile

Below is the reference for this module.



This module employs the modules `reactordef` and `material` to generate
plaintext input files for a general code. Support is provided for MCNPX, and
support for Serpent is not complete but should be straightforward. The
extension to other codes may require more effort.

- Write out
- Read in a JSON file input def.
"""
# TODO need to be able to tell the user the numbers given to the different
# cards, for parsing.

import abc

class IInputFile(object):
    """Abstract base class for classes that take system and option definitions
    to create an input file for a certain code (e.g. MCNPX, Serpent, MCODE,
    etc.).

    """
    __metaclass__ = abc.ABCMeta

    def __init__(self, fname, simdef, comments=True):
        """

        Parameters
        ----------
        fname : str
            Filename/path at which to create the input file.
        simdef: :py:class:`SimulationDefinition` or subclass.
            TODO
        comments : bool, optional
            TODO

        """
        self.fname = fname
        self.sim = simdef

    def write(self):
        self.set_up()
        self._writesubclass()
        self.fid.close()

    def set_up(self):
        self.fid = open(self.fname, 'w')

    def clean_up(self):
        self.fid.close()
    
    @abc.abstractmethod
    def _write_subclass(self):
        return NotImplementedError

    @abc.abstractmethod
    def _cell(self, cell):
        """Returns a cell card string."""
        return


class MCNPInput(IInput):
    """Contains a write method for each type of surface.
    """

    def __init__(self, fname, simdef):
        """

        """
        super(MCNPInput, self).__init__(fname, simdef)

    def _write_subclass(self):
        # Write cell cards.
        for cell in self.sim.sys.cells:
            cell.write()
        # Write surface cards.
        for surf in self.sim.sys.surfaces:
            surf.write()
        # Write data cards.
        # Source cards.
        for src in self.sim.source:
            src.write()
        # Tally cards.
        for tall in self.sim.tally:
            tall.write()




    def _cell(self, cell):
        """Returns a cell card string given a Cell card."""
        return


class SerpentInput(IInput):
    """Must find the cell used for a given material, and would need to create
    more than one material if necessary.

    """
    pass


class JSONInput(IInput):
    """Ideally we can exchange file formats through something like this
    intermediate."""
    pass
