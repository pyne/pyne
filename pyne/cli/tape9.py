"""This is a command line interface for manipulating ORIGEN v2.2 TAPE9.INP files.
"""
import os
import sys
import argparse
from glob import glob

import numpy as np

from pyne import ensdf
from pyne.dbgen.api import build_dir
from pyne.dbgen.decay import grab_ensdf_decay
from pyne.origen22 import write_tape9

def parse_ensdf(files):
    """Parses a list of ensdf files for origen."""
    decays = []
    branches = []
    for f in files:
        decs, brs = pyne.ensdf.origen_data(f)
        decays.extend(decs)
        branches.extend(brs)
    return decays, branches

def getnuc(nid='241AM', meta_t=1.0):
    """This computes ORIGEN data based on ENSDF data. It is necessary to run 
    the previous cells so that dlfinal and final_list are in the global
    namespace.

    Parameters
    ----------
    nid : nuc_id
        a valid string or int that can be converted into a nuc_id
    meta_t : float
        minimum lifetime of metastable state (in seconds) to be included (default 1.0 s)

    Returns
    -------
    nuc_id : int
        nuc_id of parent
    t12 : float
        half life of parent in seconds
    fb : float
        percent of B- decays per decay of parent 
    fbx : float
        percent of B- decays that end in metastable state
    fsf : float
        percent of spontaneous fission decays per decay of parent
    fpec : float
        percent of B+ and EC decays per decay of parent
    fpecx : float
        percent of B+ and EC decays that end in metastable state
    fa : float
        percent of A decays per decay of parent
    fn : float
        percent of B- + neutron decays per decay of parent
    fit : float
        percent of internal transition decays per decay of parent
    """
    fb = 0
    fbx = 0
    fsf = 0
    fpec = 0
    fpecx = 0
    fa = 0
    fn = 0
    fit = 0
    longest = 0
    longest2 = 0
    t12 = 0
    for item in final_list:
        if item[3] is not None and item[3] > meta_t and item[4] > 0 and item[3] != np.inf:
            if pyne.nucname.id(nid) == item[0]:
                if item[2] == 0:
                    if 'B-' in item[5]:
                        fbx += item[6]
                    if 'B+' in item[5] or "EC" in item[5]:
                        fpecx += item[6]
    
    for item in final_list:
        if item[3] is not None and item[3] > meta_t and item[4] > 0 and item[3] != np.inf:
            if pyne.nucname.id(nid) - 1  == item[0]:
                if item[1] == longest2:
                    if item[2] > 0:
                        if 'B-' in item[5]:
                            fbx += item[6]*item[8]
                        if 'B+' in item[5] or "EC" in item[5]:
                            fpecx += item[6]*item[8]
                elif item[1] > longest2:
                    longest2 = item[1]
                    if item[2] > 0:
                        if 'B-' in item[5]:
                            fbx = item[6]*item[8]
                        if 'B+' in item[5] or "EC" in item[5]:
                            fpecx = item[6]*item[8]
    for item in dlfinal:
        if pyne.nucname.id(nid) == item[0]:
            if item[1] == 0:
                if item[2] > meta_t:
                    if "%SF" in item[3] and item[3]["%SF"] != '?':
                        fsf = float(item[3]["%SF"])
                    if "%EC" in item[3] and item[3]["%EC"] != '?': 
                        fpec = float(item[3]["%EC"])
                    if "%B+" in item[3] and item[3]["%B+"] != '?':
                        fpec += float(item[3]["%B+"])
                    if "%EC+%B+" in item[3] and item[3]["%EC+%B+"] != '?':
                        fpec += float(item[3]["%EC+%B+"])
                    if "%B-N" in item[3]  and item[3]["%B-N"] != '?':
                        fn = float(item[3]["%B-N"])
                    if "%A" in item[3]  and item[3]["%A"] != '?':
                        fa = float(item[3]["%A"])
                    if "%IT" in item[3]  and item[3]["%IT"] != '?':
                        fit = float(item[3]["%IT"])
                    if "%B-" in item[3]  and item[3]["%B-"] != '?':
                        fb = float(item[3]["%B-"])
                    t12 = item[2]
        if pyne.nucname.id(nid)-1 == item[0]:
            if item[2] > meta_t and item[2] > longest and item[1] != 0:
                longest = item[2]
                t12 = item[2]
                if "%SF" in item[3] and item[3]["%SF"] != '?':
                    fsf = float(item[3]["%SF"])
                if "%EC" in item[3] and item[3]["%EC"] != '?': 
                    fpec = float(item[3]["%EC"])
                if "%B+" in item[3] and item[3]["%B+"] != '?':
                    fpec += float(item[3]["%B+"])
                if "%EC+%B+" in item[3] and item[3]["%EC+%B+"] != '?':
                    fpec += float(item[3]["%EC+%B+"])
                if "%B-N" in item[3]  and item[3]["%B-N"] != '?':
                    fn = float(item[3]["%B-N"])
                if "%A" in item[3]  and item[3]["%A"] != '?':
                    fa = float(item[3]["%A"])
                if "%IT" in item[3]  and item[3]["%IT"] != '?':
                    fit = float(item[3]["%IT"])
                if "%B-" in item[3]  and item[3]["%B-"] != '?':
                    fb = float(item[3]["%B-"])
            
    return pyne.nucname.id(nid), t12, fb, fbx, fsf, fpec, fpecx, fa, fn, fit


def main_gen(ns):
    """Generates an open TAPE9.INP file. by default this only uses completely open 
    data.
    """
    files = glob(os.path.join(ns.build_dir, 'ENSDF', 'ensdf.*'))
    if len(files) == 0:
        grab_ensdf_decay(ns.build_dir)
        files = glob(os.path.join(ns.build_dir, 'ENSDF', 'ensdf.*'))
    decays, branches = parse_ensdf(files)

_cmd_mains = {
    'gen': main_gen,
    }

def main():
    parser = argparse.ArgumentParser(description='Manipulates ORIGEN v2.2 '
                                                 'TAPE9.INP files.')
    subparsers = parser.add_subparsers(title='cmd', help='available sub-commands', 
                                       description='the subcommands')

    gen = subparsers.add_parser('gen', help='Creates a TAPE9 file based '
                                            'only on open data.')
    gen.add_argument('-o', dest='filename', default='TAPE9.INP', 
                      help='output filename')
    gen.add_argument('-b', dest='build_dir', action='store', default=build_dir,
                     help='path to the build directory.')
    ns = parser.parse_args()

    if ns.cmd not in _cmd_mains:
        sys.exit('command {0!r} could not be found'.format(ns.cmd))
    _cmd_mains[ns.cmd](ns)
    sys.exit()

if __name__ == '__main__':
    main()
