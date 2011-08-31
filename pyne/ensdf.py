import re

from pyne import nucname
from pyne.utils import time_conversion


_level_regex = re.compile('([ \d]{3}[ A-Za-z]{2})  L .{{30}}(.{10}).{28}([ M])([ \d])')

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


_decay_to = {
    '%EC': lambda x: x-10000,
    '%B+': lambda x: x-10000,
    '%EC+%B+': lambda x: x-10000,
    '%B-': lambda x: x+10000,
    '%IT': lambda x: int((x/10)*10),
    '%A': lambda x: x-20040,
    '%P': lambda x: x-10010,
    '%N': lambda x: x-10,
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
                from_nuc = nucname.zzaaam(g[0], g[-2], g[-1])
                valid_from_nuc = True
            except:
                valid_from_nuc = False
                continue

            # Grab the half-lives
            hl, unit = [s.strip() for s in g[1].split()]
            half_life = time_conversion(hl, unit)
            continue

        m = _level_cont_regex.match(line)
        if m is not None and valid_from_nuc:
            g = m.groups()
            dat = dict([d.split('=') for d in g[-1].split() if '=' in d])
            dat = {_decay_to[key](from_nuc): float(val)*0.01 for key, val in dat.items()}
            data += [(from_nuc, to_nuc, half_life, br) for to_nuc, br in dat.items()]
            continue

    if opened_here:
        ensdf.close()

    return data
