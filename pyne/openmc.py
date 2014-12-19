"""Utilities for handling OpenMC.
"""
from __future__ import print_function
import os
import io
import sys
from warnings import warn
from collections import namedtuple

if sys.version_info[0] == 2:
    from HTMLParser import HTMLParser
else:
    from html.parser import HTMLParser

from pyne import nucname
from pyne.utils import QAWarning

warn(__name__ + " is not yet QA compliant.", QAWarning)

if sys.version_info[0] > 2:
    basestring = str


class AceTable(namedtuple('_AceTable', ['alias', 'awr', 'location', 'metastable',
                                        'name', 'path', 'temperature', 'zaid'])):
    """A simple data structure reprsenenting an <ace_table /> tag in a
    cross_sections.xml file.
    """
    def __new__(cls, alias=None, awr=None, location=None, metastable=None,
                name=None, path=None, temperature=None, zaid=None,
                cross_sections_path=None):
        return super(AceTable, cls).__new__(cls, alias=alias, awr=awr,
                                            location=location,
                                            metastable=metastable, name=name,
                                            path=path, temperature=temperature,
                                            zaid=zaid)

    def __init__(self, alias=None, awr=None, location=None, metastable=None,
                 name=None, path=None, temperature=None, zaid=None,
                 cross_sections_path=None):
        """Parameters
        ----------
        alias : str, optional
            ace_table attribute.
        awr : str, optional
            ace_table attribute.
        location : str, optional
            ace_table attribute.
        metastable : str, optional
            ace_table attribute.
        name : str, optional
            ace_table attribute.
        path : str, optional
            ace_table attribute.
        temperature : str, optional
            ace_table attribute.
        zaid : str, optional
            ace_table attribute. If set or non-zero then the nucid attribute
            will be set.
        cross_sections_path : str, optional
            If this and path are both present then the abspath attribute will be
            set.
        """
        super(AceTable, self).__init__()
        nuc = None
        if zaid is not None or zaid != '0':
            meta = "0" if metastable is None else metastable
            nuc = nucname.zzaaam_to_id(zaid + meta)
            if nuc == 0:
                pass
            elif not nucname.isnuclide(nuc):  # then it's in MCNP metastable form
                nuc = nucname.mcnp_to_id(zaid)
        self.nucid = nuc
        abspath = None
        if path is not None and cross_sections_path is not None:
            if os.path.isdir(cross_sections_path):
                d = cross_sections_path
            else:
                d = os.path.dirname(cross_sections_path)
            abspath = os.path.abspath(os.path.join(d, path))
        self.abspath = abspath

    def xml(self):
        """Creates an XML representation of the ACE Table.
        """
        s = '<ace_table '
        s += " ".join(['{0}="{1}"'.format(f, getattr(self, f)) for f in self._fields
                       if getattr(self, f) is not None])
        s += '/>'
        return s


class CrossSections(HTMLParser):
    """This class represents an OpenMC cross_sections.xml file.
    """

    def __init__(self, f=None):
        """Parameters
        ----------
        f : str, file-like, or None, optional
            This is a path to the cross_sections.xml file, a file handle, or
            None indicating an empty container.
        """
        # HTMLParser is only a new-style class in python 3
        if sys.version_info[0] > 2:
            super(CrossSections, self).__init__()
        else:
            HTMLParser.__init__(self)
        self.reset()
        self._tag = None
        self.path = None
        self.filetype = 'ascii'
        self.ace_tables = []
        if f is None:
            return
        opened_here = False
        if isinstance(f, str):
            opened_here = True
            self.path = f
            f = io.open(f, 'r')
        raw = f.read()
        self.feed(raw)
        if opened_here:
            f.close()

    def handle_starttag(self, tag, attrs):
        self._tag = tag
        if tag == 'ace_table':
            self.handle_ace_table(attrs)

    def handle_endtag(self, tag):
        self._tag = None

    def handle_startendtag(self, tag, attrs):
        if tag == 'ace_table':
            self.handle_ace_table(attrs)

    def handle_data(self, data):
        if self._tag == 'filetype':
            self.filetype = data
        elif self._tag == 'directory':
            self.path = data.strip()

    def handle_ace_table(self, attrs):
        ace_table = AceTable(cross_sections_path=self.path, **dict(attrs))
        self.ace_tables.append(ace_table)

    def xml(self):
        """Returns an XML representation of the cross sections file.
        """
        template = ('<?xml version="1.0" ?>\n'
                    '<cross_sections>\n'
                    '  <filetype>{filetype}</filetype>\n'
                    '  {ace_tables}\n'
                    '</cross_sections>\n')
        ace_tables = "\n  ".join([a.xml() for a in self.ace_tables])
        s = template.format(filetype=self.filetype, ace_tables=ace_tables)
        return s
