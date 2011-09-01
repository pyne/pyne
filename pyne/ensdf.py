import re

import numpy as np

from pyne import nucname
from pyne.utils import time_conversion


_level_regex = re.compile('([ \d]{3}[ A-Za-z]{2})  L .{30}(.{10}).{28}([ M])([ 1-9])')

_level_cont_regex = re.compile('([ \d]{3}[ A-Za-z]{2})[0-9A-Za-z] L (.*)')

def _to_zzaaam(nuc, m, s):
    nuc_zz = nucname.zzaaam(nuc.strip())
    if m == 'M':
        state = s.strip()
        if 0 < len(state):
            state = int(state)
        else:
            state = 1
        nuc_zz += state
    return nuc_zz


def _to_float(x):
    x = x.strip()
    x = x.replace('$', '')
    x = x.replace('?', '')

    if 0 == len(x):
        x = 0.0
    else:
        x = float(x)

    return x
    


_decay_to = {
    '%EC': lambda x: (x-10000)/10*10,
    '%B+': lambda x: (x-10000)/10*10,
    '%EC+%B+': lambda x: (x-10000)/10*10,
    '%B-': lambda x: (x+10000)/10*10,
    '%IT': lambda x: x/10*10,
    '%A': lambda x: (x-20040)/10*10,
    '%P': lambda x: (x-10010)/10*10,
    '%N': lambda x: (x-10)/10*10,
    }


def half_life(ensdf):
    """Grabs the half-lives from an ENSDF file.

    Parameters
    ----------
    ensdf : str or file-like object
        ENSDF file to inspect for half-life data

    Returns
    -------
    data : list of 4-tuples
        List of tuples where the indices match

        1. from_nuclide, int (zzaaam)
        2. to_nuclide,  int (zzaaam)
        3. half_life, float (seconds)
        4. branch_ratio, float (frac)

    """
    opened_here = False
    if isinstance(ensdf, basestring):
        ensdf = open(ensdf, 'r')
        opened_here = True

    data = []
    from_nuc = 0
    half_life = 0.0
    valid_from_nuc = False

    # Run through the file
    lines = ensdf.readlines()
    for i, line in enumerate(lines):
        # See if the line matches
        m = _level_regex.match(line)
        if m is not None:
            g = m.groups()

            # grab the from nuclide
            try:
                from_nuc = _to_zzaaam(g[0], g[-2], g[-1])
                valid_from_nuc = True
            except:
                valid_from_nuc = False
                continue

            # Grab the half-lives
            time_info = g[1].replace('?', '').strip()
            if 0 == len(time_info):
                valid_from_nuc = False
            elif time_info == 'STABLE':
                half_life = np.inf
                data += [(from_nuc, from_nuc, half_life, 1.0)]
            else:
                time_unit = [s.strip() for s in time_info.split()]
                if 2 == len(time_unit):
                    hl, unit = time_unit
                    half_life = time_conversion(float(hl), unit)
            continue

        m = _level_cont_regex.match(line)
        if m is not None and valid_from_nuc:
            g = m.groups()
            dat = dict([d.split('=')[:2] for d in g[-1].replace('$', ' ').split() if '=' in d])
            dat = {_decay_to[key](from_nuc): _to_float(val)*0.01 
                        for key, val in dat.items() if key in _decay_to.keys()}
            data += [(from_nuc, to_nuc, half_life, br) for to_nuc, br in dat.items() if 0.0 < br]
            continue

    if opened_here:
        ensdf.close()

    return data
