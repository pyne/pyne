import re

import numpy as np

from pyne import nucname
from pyne.utils import to_sec


_level_regex = re.compile('([ \d]{3}[ A-Za-z]{2})  L (.{10}).{20}(.{10}).{28}([ M])([ 1-9])')

_level_cont_regex = re.compile('([ \d]{3}[ A-Za-z]{2})[0-9A-Za-z] L (.*)')

def _to_id(nuc, m, s):
    nucid = nucname.id(nuc.strip())
    if m == 'M':
        state = s.strip()
        if 0 < len(state):
            state = int(state)
        else:
            state = 1
        nucid += state
    return nucid

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
    '%EC': lambda x: (x-10000000)/10000*10000,
    '%B+': lambda x: (x-10000000)/10000*10000,
    '%EC+%B+': lambda x: (x-10000000)/10000*10000,
    '%B-': lambda x: (x+10000000)/10000*10000,
    '%IT': lambda x: x/10000*10000,
    '%A': lambda x: (x-20040000)/10000*10000,
    '%P': lambda x: (x-10010000)/10000*10000,
    '%N': lambda x: (x-10000)/10000*10000,
    }

def half_life(ensdf):
    """Grabs the half-lives from an ENSDF file.

    Parameters
    ----------
    ensdf : str or file-like object
        ENSDF file to inspect for half-life data

    Returns
    -------
    data : list of 5-tuples
        List of tuples where the indices match

        1. from_nuclide, int (id)
        2. level, float (MeV) - from_nuc's energy level
        3. to_nuclide, int (id)
        4. half_life, float (seconds)
        5. branch_ratio, float (frac)

    """
    opened_here = False
    if isinstance(ensdf, basestring):
        ensdf = open(ensdf, 'r')
        opened_here = True

    data = []
    from_nuc = 0
    half_life = 0.0
    level = 0.0
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
                from_nuc = _to_id(g[0], g[-2], g[-1])
                valid_from_nuc = True
            except:
                valid_from_nuc = False
                continue

            # parse energy level
            try:
                level = float(g[1]) * 1E-3
            except ValueError:
                pass

            # Grab the half-lives
            time_info = g[2].replace('?', '').strip()
            if 0 == len(time_info):
                valid_from_nuc = False
            elif time_info == 'STABLE':
                half_life = np.inf
                data += [(from_nuc, 0.0, from_nuc, half_life, 1.0)]
            else:
                time_unit = [s.strip() for s in time_info.split()]
                if 2 == len(time_unit):
                    hl, unit = time_unit
                    half_life = to_sec(float(hl), unit)
            continue

        m = _level_cont_regex.match(line)
        if m is not None and valid_from_nuc:
            g = m.groups()
            dat = {}
            raw_children = g[-1].replace(' AP ', '=')
            raw_children = raw_children.replace('$', ' ').split()
            for raw_child in raw_children:
                if '=' in raw_child:
                    rx, br = raw_child.split('=')[:2]
                else:
                    continue
                dat[rx] = br
            dat = dict([(_decay_to[key](from_nuc), _to_float(val)*0.01) 
                        for key, val in dat.items() if key in _decay_to.keys()])
            data += [(from_nuc, level, to_nuc, half_life, br) for to_nuc, br in \
                                                              dat.items() if 0.0 < br]
            continue

    if opened_here:
        ensdf.close()

    # Hack to calculate metastable state number, make sure it doesn't go over 10,
    # and then change the from_nuc value, and remove all other states
    # FIXME: while the renaming bases on level should still happen, the limit
    # of the 10 lowest levels should be removed when id is removed.
    nuclvl = {}
    for row in data:
        from_nuc, level, to_nuc, half_life, br = row
        if from_nuc not in nuclvl:
            nuclvl[from_nuc] = set()
        nuclvl[from_nuc].add(level)
    for nuc in nuclvl.keys():
        nuclvl[nuc] = dict([(lvl, i) for i,lvl in enumerate(sorted(nuclvl[nuc])[:10])])
    data = [(row[0] + nuclvl[row[0]][row[1]],) + row[1:] for row in data \
                if (row[0] in nuclvl) and (row[1] in nuclvl[row[0]])]

    return data
