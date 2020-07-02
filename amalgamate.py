#!/usr/bin/env python

"""Amalgate pyne C++ library sources into a single source and a single header file.
This makes the C++ portion of pyne more portable to other projects.

Originally inspired by JsonCpp: http://svn.code.sf.net/p/jsoncpp/code/trunk/jsoncpp/amalgamate.py
"""
from __future__ import print_function, unicode_literals
import os
import sys
import io
from argparse import ArgumentParser

CODE_EXTS = {'.c', '.cpp', '.cxx', '.h', '.hpp', '.hxx'}
CODE_EXTS |= {e.upper() for e in CODE_EXTS}
SOURCE_EXTS = {'.c', '.cpp', '.cxx'}
SOURCE_EXTS |= {e.upper() for e in SOURCE_EXTS}
HEADER_EXTS = {'.h', '.hpp', '.hxx'}
HEADER_EXTS |= {e.upper() for e in HEADER_EXTS}

DEFAULT_FILES = [
    'license.txt',
    'src/utils.h',
    'src/utils.cpp',
    'src/extra_types.h',
    'src/h5wrap.h',
    'src/state_map.cpp',
    'src/nucname.h',
    'src/nucname.cpp',
    'src/rxname.h',
    'src/rxname.cpp',
    'src/_atomic_data.h',
    'src/_atomic_data.cpp',
    'src/data.h',
    'src/data.cpp',
    #'src/dagmc_bridge.cpp',
    #'src/dagmc_bridge.h',
    'src/json-forwards.h',
    'src/json.h',
    'src/jsoncpp.cpp',
    'src/jsoncustomwriter.h',
    'src/jsoncustomwriter.cpp',
    'src/material.h',
    'src/material.cpp',
    'src/material_library.h',
    'src/material_library.cpp',
    'src/enrichment_cascade.h',
    'src/enrichment_cascade.cpp',
    'src/enrichment.h',
    'src/enrichment.cpp',
    'src/enrichment_symbolic.h',
    'src/enrichment_symbolic20.cpp',
    'src/_decay.h',
    'src/_decay.cpp',
    ]

class AmalgamatedFile(object):
    def __init__(self, path):
        self.path = path
        self._blocks = []
        self._filenames = []

    def append_line(self, line):
        """Adds some text to the end of the file."""
        if not line.endswith( '\n' ):
            line += '\n'
        self._blocks.append(line)

    def append_file(self, filename, comment_out=None):
        """Adds a whole file to the end of this one."""
        if comment_out is None:
            _, ext = os.path.splitext(filename)
            comment_out = ext not in CODE_EXTS
        self._blocks.append('//\n// start of {0}\n//\n'.format(filename))
        with open(filename, 'rt') as f:
            content = f.read()
        if comment_out:
            content = '// ' + content.replace('\n', '\n// ')
        self._blocks.append(content)
        self._blocks.append('//\n// end of {0}\n//\n\n\n'.format(filename))
        self._filenames.append(filename)

    def prepend_files(self):
        """Adds a file listing to the begining of the almagamted file."""
        s = '// This file is composed of the following original files:\n\n'
        for f in self._filenames:
            s += '//   {0}\n'.format(f)
        s += '\n'
        self._blocks.insert(0, s)

    def write(self):
        self.prepend_files()
        if sys.version > '3':
            txt = ''.join(self._blocks)
        else:
            txt = ''.join([block.decode('utf-8') for block in self._blocks])
        d = os.path.dirname(self.path)
        if len(d) > 0 and not os.path.isdir(d):
            os.makedirs(d)
        with io.open(self.path, 'wb') as f:
            f.write(txt.encode('utf-8'))

def main():
    parser = ArgumentParser()
    parser.add_argument('-s', dest='source_path', action='store',
                        default='pyne.cpp', help='Output *.cpp source path.')
    parser.add_argument('-i', dest='header_path', action='store',
                        default='pyne.h', help='Output header path.')
    parser.add_argument('-f', dest='files', nargs='+', help='Files to amalgamate.',
                        default=DEFAULT_FILES)
    ns = parser.parse_args()

    # header file
    hdr = AmalgamatedFile(ns.header_path)
    hdr.append_line( '// PyNE amalgated header http://pyne.io/' )
    hdr.append_line('#ifndef PYNE_52BMSKGZ3FHG3NQI566D4I2ZLY')
    hdr.append_line('#define PYNE_52BMSKGZ3FHG3NQI566D4I2ZLY')
    hdr.append_line('')
    hdr.append_line('#define PYNE_IS_AMALGAMATED')
    hdr.append_line('')
    for f in ns.files:
        _, ext = os.path.splitext(f)
        if ext in SOURCE_EXTS:
            continue
        hdr.append_file(f)
    hdr.append_line('#endif  // PYNE_52BMSKGZ3FHG3NQI566D4I2ZLY')

    # source file
    src = AmalgamatedFile(ns.source_path)
    src.append_line('// PyNE amalgated source http://pyne.io/')
    src.append_line('#include "{0}"'.format(os.path.relpath(ns.header_path,
                                            os.path.dirname(ns.source_path))))
    src.append_line('')
    for f in ns.files:
        _, ext = os.path.splitext(f)
        if ext in HEADER_EXTS:
            continue
        src.append_file(f)

    # write both
    hdr.write()
    src.write()

if __name__ == '__main__':
    main()
