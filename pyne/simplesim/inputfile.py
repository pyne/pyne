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
# TODO check uniqueness
# TODO manage print order of cell cards, this may just be given by the ordered
# dict.
# TODO let users enter cell or surface cards to a certain location or number;
# describe how the remaining cards are modified.
# TODO let user specify number format for different types of floats
# TODO provision for printing materials...
# TODO filename should be with write
# TODO special treatment for Material card.
# TODO improve how the $ nucname.name is printed on material cards: alignment!
# TODO raise exception or warning whenever, if bypass_wrap, when the number of
# characters between newlines is larger than 80.

import abc
import datetime
import re
import textwrap
import warnings

import numpy as np

class IInputFile(object):
    """Abstract base class for classes that take system and option definitions
    to create an input file for a certain code (e.g. MCNPX, Serpent, MCODE,
    etc.).

    """
    __metaclass__ = abc.ABCMeta

    def __init__(self, simdef, comments=True, heading=None, plug=True,
                 float_format="% .5g"):
        """

        Parameters
        ----------
        simdef: :py:class:`SimulationDefinition` or subclass.
            TODO
        comments : bool, optional
            TODO

        """
        self.sim = simdef
        self.comments = comments
        self.heading = heading
        self.plug = plug
        self.float_format = float_format

    def write(self, fname):
        """
        fname : str
            Filename/path at which to create the input file.
        """
        self.fname = fname
        self.set_up()
        # Should write the plug in the appropriate place.
        self._write_subclass()
        self.fid.close()
        self.fname = None

    def set_up(self):
        self.fid = open(self.fname, 'w')

    def clean_up(self):
        self.fid.close()
    
    @abc.abstractmethod
    def _write_subclass(self):
        return NotImplementedError

    def _write_plug(self):
        if self.plug:
            self._write_plug_subclass(
                    "Generated with the Python package PyNE (pyne.github.com).")
    
    @abc.abstractmethod
    def _write_plug_subclass(self, string):
        return NotImplementedError

    @abc.abstractmethod
    def add_user_card(self, block, card, comment=None):
        return NotImplementedError

    @abc.abstractmethod
    def add_user_card_literal(self, block, string):
        return NotImplementedError


class MCNPInput(IInputFile):
    """Contains a write method for each type of surface.
    """
    # TODO user can overload commenting methods
    def __init__(self, simdef, comments=True, heading=None,
            description=None, plug=True, float_format="% .5g",
            cont_by_amp=False):
        """
        cont_by_amp : bool, optional


        """
        # TODO could cleanup keyword arguments wiht **kwarg.
        # Be careful with the order of inputs here for the kwargs.
        super(MCNPInput, self).__init__(simdef, comments, heading,
                plug, float_format)
        self.description = description
        self.cont_by_amp = cont_by_amp
        # Comment wrapping.
        self.commentwrap = textwrap.TextWrapper(
                width=80,
                initial_indent=   '  C ',
                subsequent_indent='  c     ',
                break_long_words=True)
        # Card wrapping.
        if self.cont_by_amp:
            self.cardwrap = textwrap.TextWrapper(
                    width=78,
                    #initial_indent=2 * ' ',
                    #subsequent_indent=2 * ' '
                    )
            self._card_end_line = ' &\n'
        else:
            self.cardwrap = textwrap.TextWrapper(
                    width=80,
                    #initial_indent=2 * ' ',
                    subsequent_indent=5 * ' ')
            self._card_end_line = '\n'

    def _write_subclass(self):
        # Header
        if self.heading:
            self._write_comment(self.heading)
        else:
            # MCNP files need a first 'comment' line.
            self._write_comment("datetime: %s" % str(datetime.datetime.now()))
        if self.description:
            self._write_comment(self.description)
        self._write_plug()

        # Write cell cards.
        self._write_deck_heading("Cell")
        self._write_dictionary(self.sim.sys.cells)

        # Write surface cards.
        self._new_line()
        self._write_deck_heading("Surface")
        self._write_dictionary(self.sim.sys.surfaces)

        # Write data cards.
        self._new_line()
        self._write_deck_heading("Data")
        # TODO WWN card will not work with the card wrapper as it is....
        # Material cards.
        self._write_data_heading("Material")
        self._write_dictionary(self.sim.sys.materials)
        # Source cards.
        self._write_data_heading("Source")
        self._write_dictionary(self.sim.source)
        # Tally cards.
        self._write_data_heading("Tally")
        self._write_dictionary(self.sim.tally)
        # Misc cards.
        self._write_data_heading("Miscellaneous")
        self._write_dictionary(self.sim.misc)

    def _write_dictionary(self, dictionary):
        for key, card in dictionary.iteritems():
            if self.comments:
                self._write_comment(card.comment())
            self._write_card(card)
            if self.comments:
                self._write_comment()

    def add_user_card(self, block, card, comment=None):
        # TODO
        # use textwrap
        pass


    def add_user_card_literal(self, block, string):
        # TODO
        pass

    def _write_plug_subclass(self, string):
        self._write_comment(string)

    def _write_deck_heading(self, heading):
        heading = "{0} Cards".format(heading)
        n_chars = len(heading)
        self._write_comment(n_chars * "=")
        self._write_comment(heading)
        self._write_comment(n_chars * "=")

    def _write_data_heading(self, heading):
        heading = "{0} Cards".format(heading)
        n_chars = len(heading)
        self._write_comment()
        self._write_comment(n_chars * "*")
        self._write_comment(heading)
        self._write_comment(n_chars * "*")

    def _write_comment(self, comment=""):
        if comment != "": strlist = self.commentwrap.wrap(comment)
        else:             strlist = ["  C"]
        for entry in strlist:
            self.fid.write(entry + '\n')

    def _write_card(self, card):
        string = card.mcnp(self.float_format, self.sim)
        if card.bypass_wrap:
            # Check that the number of characters between any two newlines is
            # not greater than 80. Got this line of code from activestate.com.
            starts = [0]
            starts += [match.start() for match in re.finditer('\n', string)]
            # Account for the fact that the string does not end in a newline.
            starts += [len(string) - 1]
            if (np.diff(np.array(starts)) > 79).any():
                warnings.warn("Card {0} contains lines longer than "
                        "79 columns.".format(card.name))
            self.fid.writelines(string)
        else:
            strlist = self.cardwrap.wrap(string)
            counter = 0
            for entry in strlist:
                counter += 1
                self.fid.write(entry)
                if counter < len(strlist):
                    self.fid.writelines(self._card_end_line)
        self._new_line()

    def _new_line(self):
        self.fid.write('\n')


class SerpentInput(IInputFile):
    """Must find the cell used for a given material, and would need to create
    more than one material if necessary.

    """
    pass

