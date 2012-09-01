#!/usr/bin/env python

"""The ``inputfile`` module can be imported as such::

    from pyne.simplesim import inputfile

Below is the reference for this module.

TODO
This module employs the modules `reactordef` and `material` to generate
plaintext input files for a general code. Support is provided for MCNPX, and
support for Serpent is not complete but should be straightforward. The
extension to other codes may require more effort.

An inheritance diagram of all the classes in this module can be found at
:ref:`pyne_simplesim_inheritance`.

.. moduleauthor:: Chris Dembia <cld72@cornell.edu>

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

import abc
import datetime
import re
import textwrap
import warnings

import numpy as np

from pyne.simplesim import definition, cards

class IInputFile(object):
    """Abstract base class for classes that take system and option definitions
    to create an input file for a certain code (e.g. MCNPX, Serpent, MCODE,
    etc.).

    """
    __metaclass__ = abc.ABCMeta

    def __init__(self, simdef, comments=True, title=None, plug=True,
                 float_format="% .5g", *args, **kwargs):
        """
        Parameters
        ----------
        simdef : :py:class:`definition.SimulationDefinition` or subclass.
            The simulation for which the user desires an input file.
        comments : bool, optional
            Display comments along with cards. The comments are rather long.
        title : str, optional
            A short (one line) title for the input.
        plug : bool, optional
            Prints a plug for PyNE somewhere in the input file =p. 
        float_format : str, optional
            The format with which floating point numbers should be printed.
            Uses familiar ``sprintf`` style format specifications that is
            common in many languages.

        """
        super(IInputFile, self).__init__()
        self.sim = simdef
        self.comments = comments
        self.title = title
        self.plug = plug
        self.float_format = float_format

    def write(self, fname):
        """
        Parameters
        ----------
        fname : str
            Filename/path at which to create the input file.

        """
        self.fname = fname
        self._set_up()
        # Should write the plug in the appropriate place.
        self._write_subclass()
        self._clean_up()
        self.fname = None

    def _set_up(self):
        self.fid = open(self.fname, 'w')

    def _clean_up(self):
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
    def add_user_card(self, *args):
        return NotImplementedError

    @abc.abstractmethod
    def add_user_literal(self, *args):
        return NotImplementedError


class MCNPInput(IInputFile):
    """An input file for MCNP. To write the file, call :py:meth:`write`. The
    user can modify the format of the input file by subclassing this class and
    overloading the _write_ methods.

    """
    def __init__(self, simdef, description=None, cont_by_amp=False, **kwargs):
        """
        Parameters
        ----------
        simdef : :py:class:`definition.MCNPSimulation`
            See :py:class:`IInputFile`. This cannot be a
            :py:class:`definition.SimulationDefinition`.
        description : str, optional
            An arbitrarily long commented description placed after the title
            but before any cards.
        cont_by_amp : bool, optional
            In MCNP, line continuation for lines with more than 80 characters
            can be achieved via 5 or more spaces on subsequent lines or by
            appending an ampersand (&) before the continued line.
        comments : bool, optional
            See :py:class:`IInputFile`.
        title : str, optional
            See :py:class:`IInputFile`. If not provided, the title defaults to
            a datetime string.
        plug : bool, optional
            See :py:class:`IInputFile`.
        float_format : str, optional
            See :py:class:`IInputFile`.

        """
        super(MCNPInput, self).__init__(simdef, **kwargs)
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
        # User cards. comments:
        self.user_cell_comm = []
        self.user_surface_comm = []
        self.user_data_comm = []
        # cards:
        self.user_cell_cards = []
        self.user_surface_cards = []
        self.user_data_cards = []
        # literal:
        self.user_cell_literal = []
        self.user_surface_literal = []
        self.user_data_literal = []

    def add_user_card(self, block, card, comment=None):
        """It is possible that the user wants to insert a card for which there
        is no :py:class:`ICard`. The user can supply such cards with this
        method or with :py:meth:`add_user_literal`. This method allows the user
        to supply a comment along with the card, and automatically wraps the
        card if it is longer than 80 characters. Any number of cards can be
        added to a block, and they are placed at the end of the block, but
        before text added by :py:class:`add_user_literal`. Cards are printed in
        the order they are added.

        Parameters
        ----------
        block : str
            The block to which the user desires to add a card. Either 'cell',
            'surface', or 'data'.
        card : str
            The card, including the card number.
        comment : str, optional
            A comment placed before the card.

        Examples
        --------
        The following shows how the user can provide a cell card::

            inp.add_user_card('cell', '11 0 -5 IMP:N=0', comment='Graveyard.')

        """
        if block == 'cell':
            self.user_cell_comm += [comment]
            self.user_cell_cards += [card]
        elif block == 'surface':
            self.user_surface_comm += [comment]
            self.user_surface_cards += [card]
        elif block == 'data':
            self.user_data_comm += [comment]
            self.user_data_cards += [card]
        else:
            raise ValueError("The input ``block`` must be 'cell', 'surface', "
                    "or 'data'. User provided {0!r}.".format(block))

    def add_user_literal(self, block, string):
        """This method is used instead of :py:meth:`add_user_card` if line
        wrapping is not desired, such as if the user wants to specify a
        multi-line card on their own. This can also be used by the user to
        provide comments. A newline is appended to the end of the
        string. Any number of cards can be added to a block, and they are
        placed after :py:class:`add_user_card` cards. Cards are printed in
        the order they are added.

        Parameters
        ----------
        block : str
            The block to which the user desires to add a card. Either 'cell',
            'surface', or 'data'.
        string : str
            The card, possibly including newlines.

        Examples
        --------
        In this case, line wrapping is `not` automatically performed::

            inp.add_user_literal('data', 'M1 1001 1\\n     8016 2')

        """
        if block == 'cell':
            self.user_cell_literal += [string]
        elif block == 'surface':
            self.user_surface_literal += [string]
        elif block == 'data':
            self.user_data_literal += [string]
        else:
            raise ValueError("The input ``block`` must be 'cell', 'surface', "
                    "or 'data'. User provided {0!r}.".format(block))

    def _write_subclass(self):
        # Header
        if self.title:
            if len(self.title) > 80:
                raise ValueError("MCNP titles must be less than 80 chars.")
            self.fid.write(self.title + '\n')
        else:
            # MCNP files need a first 'comment' line.
            self.fid.write("datetime: %s\n" % str(datetime.datetime.now()))
        if self.description:
            self._write_comment(self.description)
        self._write_plug()

        # Write cell cards.
        self._write_deck_heading("Cell")
        self._write_dictionary(self.sim.sys.cells)
        # User cards.
        self._write_user("cell", self.user_cell_comm, self.user_cell_cards)
        self._write_user_literal("cell", self.user_cell_literal)

        # Write surface cards.
        self._new_line()
        self._write_deck_heading("Surface")
        self._write_dictionary(self.sim.sys.surfaces)
        # User cards.
        self._write_user("surface", self.user_surface_comm,
                self.user_surface_cards)
        self._write_user_literal("surface", self.user_surface_literal)

        # Write data cards.
        self._new_line()
        self._write_deck_heading("Data")
        # Material cards.
        if len(self.sim.sys.materials) > 0:
            self._write_data_heading("Material")
            self._write_dictionary(self.sim.sys.materials)
        # Source cards.
        if len(self.sim.source) > 0:
            self._write_data_heading("Source")
            self._write_dictionary(self.sim.source)
        # Distribution cards.
        if len(self.sim.dists) > 0:
            self._write_data_heading("Distribution")
            self._write_dictionary(self.sim.dists)
        # Tally cards.
        if len(self.sim.tally) > 0:
            self._write_data_heading("Tally")
            self._write_dictionary(self.sim.tally)
        # Transformation cards.
        if len(self.sim.transformations) > 0:
            self._write_data_heading("Transformation")
            self._write_dictionary(self.sim.transformations)
        # Misc cards.
        if len(self.sim.misc) > 0:
            self._write_data_heading("Miscellaneous")
            self._write_dictionary(self.sim.misc)
        # User cards.
        self._write_user("data", self.user_data_comm, self.user_data_cards)
        self._write_user_literal("data", self.user_data_literal)

    def _write_dictionary(self, dictionary):
        for key, card in dictionary.iteritems():
            if self.comments:
                self._write_comment(card.comment())
            self._write_card(card)
            if self.comments:
                self._write_comment()

    def _write_user(self, block, comm, cards):
        if len(comm) > 0:
            self._write_comment("User {0} cards.".format(block))
            for idx in range(len(comm)):
                if comm[idx]: self._write_comment(comm[idx])
                self._write_cardwrap(cards[idx])
            self._new_line()

    def _write_user_literal(self, block, lit):
        if len(lit) > 0:
            self._write_comment("User literal {0} cards.".format(block))
            for string in lit:
                self.fid.write(string + '\n')

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
        # Wrap the parts between newlines, but still place in newlines where
        # they are requested. str.partition('\n')[0] is left of '\n', and
        # str.partition('\n')[2] is right of '\n'.
        cont = True
        while cont:
            currstr = string.partition('\n')
            self._write_cardwrap(currstr[0])
            self._new_line()
            if currstr[2] == '': cont = False
            else: string = currstr[2]

    def _write_cardwrap(self, string):
        strlist = self.cardwrap.wrap(string)
        counter = 0
        for entry in strlist:
            counter += 1
            self.fid.write(entry)
            if counter < len(strlist): self.fid.write(self._card_end_line)

    def _new_line(self):
        self.fid.write('\n')


class SerpentInput(IInputFile):
    """Must find the cell used for a given material, and would need to create
    more than one material if necessary.

    """
    pass

