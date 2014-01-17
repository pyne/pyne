import re

import numpy as np

from pyne import nucname
from .utils import to_sec


_valexp = re.compile('([0-9.]*)([Ee][+-]\d*)')
_val = re.compile('(\d*)[.](\d*)')
_errpm = re.compile('[+](\d*)[-](\d*)')
_err = re.compile('[ ]*(\d*)')
_base = '([ \d]{3}[ A-Za-z]{2})'
_ident = re.compile(_base + '    (.{30})(.{26})(.{7})(.{6})')
_g = re.compile(_base + '  G (.{10})(.{2})(.{8})(.{2}).{24}(.{7})(.{2})')
_p = re.compile(_base + '  P (.{10})(.{2})(.{18})(.{10})(.{6}).{9}(.{10})(.{2})(.{4})')
_norm = re.compile(_base + '  N (.{10})(.{2})(.{8})(.{2})(.{8})(.{2})(.{8})(.{6})(.{7})(.{2})')
_normp = re.compile(_base + ' PN (.{10})(.{2})(.{8})(.{2})(.{8})(.{2})(.{7})(.{2})')
_decays = [' B- ', ' B+ ', ' EC ', ' IT ', ' A ']
_level_regex = re.compile(_base + '  L (.{10}).{20}(.{10}).{28}([ M])([ 1-9])')
_level_regex2 = re.compile(_base + '  L (.{10})(.{2})(.{18})(.{10})(.{6})(.{9})(.{10})(.{2})(.{1})([ M])([ 1-9])')
_level_cont_regex = re.compile('([ \d]{3}[ A-Za-z]{2})[0-9A-Za-z] L (.*)')


def _to_id(nuc, m=None, s=None):
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


def _to_time(tstr, errstr):
    t = tstr.strip()
    # This accepts questionable levels
    t = t.replace('?', '')
    tobj = [s.strip(' ()') for s in t.split()]
    if len(tobj) == 2:
        t, t_unit = tobj
        t, terr = _get_val_err(t, errstr)
        tfinal = to_sec(t, t_unit)
        tfinalerr = None
        if type(terr) == float:
            tfinalerr = to_sec(terr, t_unit)
        elif terr is not None:
            tfinalerr = to_sec(terr[0], t_unit), to_sec(terr[1], t_unit)
    elif 'STABLE' in t:
        tfinal = np.inf
        tfinalerr = None
    else:
        tfinal = None
        tfinalerr = None
    return tfinal, tfinalerr


_decay_to = {
    '%EC': lambda x: (x - 10000000) / 10000 * 10000,
    '%B+': lambda x: (x - 10000000) / 10000 * 10000,
    '%EC+%B+': lambda x: (x - 10000000) / 10000 * 10000,
    '%B-': lambda x: (x + 10000000) / 10000 * 10000,
    '%IT': lambda x: x / 10000 * 10000,
    '%A': lambda x: (x - 20040000) / 10000 * 10000,
    '%P': lambda x: (x - 10010000) / 10000 * 10000,
    '%N': lambda x: (x - 10000) / 10000 * 10000,
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
    if isinstance(ensdf, basestring):
        with open(ensdf, 'r') as f:
            lines = f.readlines()
    else:
        lines = ensdf.readlines()

    data = []

    for line in lines:
        level_l = _level_regex2.match(line)
        if level_l is not None:
            level, half_lifev, from_nuc = _parse_level_record(level_l)
            if half_lifev == np.inf and from_nuc is not None:
                data.append((from_nuc, 0.0, from_nuc, half_lifev, 1.0))
            if level is None:
                level = 0.0
            continue
        levelc = _level_cont_regex.match(line)
        if levelc is not None and from_nuc is not None and half_lifev is not None:
            dat = _parse_level_continuation_record(levelc)
            dat = dict([(_decay_to[key](from_nuc), _to_float(val) * 0.01)
                        for key, val in dat.items() if key in _decay_to])
            data += [(from_nuc, level, to_nuc, half_lifev, br)
                     for to_nuc, br in dat.items() if 0.0 < br]
            continue


    # Hack to calculate metastable state number, make sure it doesn't go over 10,
    # and then change the from_nuc value, and remove all other states
    # FIXME: while the renaming bases on level should still happen, the limit
    # of the 10 lowest levels should be removed when id is removed.
    nuclvl = {}
    for row in data:
        from_nuc, level, to_nuc, half_lifev, br = row
        if from_nuc not in nuclvl:
            nuclvl[from_nuc] = set()
        nuclvl[from_nuc].add(level)
    for nuc in nuclvl.keys():
        nuclvl[nuc] = dict([(lvl, i) for i, lvl in enumerate(sorted(nuclvl[nuc])[:10])])
    data = [(row[0] + nuclvl[row[0]][row[1]],) + row[1:] for row in data \
            if (row[0] in nuclvl) and (row[1] in nuclvl[row[0]])]

    return data


def _getvalue(obj, fn=float):
    try:
        return fn(obj)
    except:
        return None


def _get_val_err(valstr, errstr):
    pm = _errpm.match(errstr)
    err = _err.match(errstr)
    if pm is None and err.group(1) == '':
        return _getvalue(valstr), None
    val = _valexp.match(valstr)
    if val is None:
        valexp = ''
        val = valstr
    else:
        valexp = val.group(2)
        val = val.group(1)
    punc = _val.match(val.strip())
    if pm is not None:
        if punc is None:
            errplus = _getvalue(pm.group(1) + valexp)
            errminus = _getvalue(pm.group(2) + valexp)
        else:
            errplus = _get_err(len(punc.group(2)), pm.group(1), valexp)
            errminus = _get_err(len(punc.group(2)), pm.group(2), valexp)
        return _getvalue(valstr), (errplus, errminus)
    else:
        if punc is None:
            errplus = _getvalue(errstr + valexp)
        else:
            errplus = _get_err(len(punc.group(2)), errstr, valexp)
        return _getvalue(valstr), errplus


def _get_err(plen, errstr, valexp):
    errp = list((errstr.strip()).zfill(plen))
    errp.insert(-plen, '.')
    return float(''.join(errp) + valexp)


def _parse_level_record(l_rec):
    """
    This Parses and ENSDF level record
    
    Parameters
    ----------
    g : re.MatchObject
        regular expression MatchObject
        
    Returns
    -------
    e : float
        Level energy in MeV
    tfinal : float
        Half life in seconds
    from_nuc : int
        nuc id of nuclide
    """
    e, de = _get_val_err(l_rec.group(2), l_rec.group(3))
    tfinal, tfinalerr = _to_time(l_rec.group(5), l_rec.group(6))
    try:
        from_nuc = _to_id(l_rec.group(1), l_rec.group(11), l_rec.group(12))
    except:
        from_nuc = None
    return e, tfinal, from_nuc


def _parse_level_continuation_record(lc_rec):
    """
    This Parses and ENSDF level record
    
    Parameters
    ----------
    g : re.MatchObject
        regular expression MatchObject

    Returns
    -------
    dat : dict
        dictionary of branching ratios of different reaction channels
    """
    g = lc_rec.groups()
    dat = {}
    raw_children = g[-1].replace(' AP ', '=')
    raw_children = raw_children.replace('$', ' ').split()
    for raw_child in raw_children:
        if '=' in raw_child:
            rx, br = raw_child.split('=')[:2]
        else:
            continue
        dat[rx] = br
    return dat


def _parse_gamma_record(g):
    """
    This parses an ENSDF gamma record

    Parameters
    ----------
    g : re.MatchObject
        regular expression MatchObject

    Returns
    -------
    dat : np.ndarray
        This array contains 6 floats corresponding to:
            * gamma ray energy in keV
            * uncertainty in energy
            * intensity
            * uncertainty in intensity
            * electron conversion intensity
            * uncertainty in electron conversion intensity
    """
    dat = np.zeros(6)
    en, en_err = _get_val_err(g.group(2), g.group(3))
    inten, inten_err = _get_val_err(g.group(4), g.group(5))
    conv, conv_err = _get_val_err(g.group(6), g.group(7))
    dat[:] = en, en_err, inten, inten_err, conv, conv_err
    return dat


def _parse_normalization_record(n_rec):
    """
    This parses an ENSDF normalization record

    Parameters
    ----------
    n_rec : re.MatchObject
        regular expression MatchObject

    Returns
    -------
    nr : float
        Multiplier for converting relative photon intensity to photons per 100
        decays of the parent through the decay branch or to photons per 100
        neutron captures for (n,g).
    nr_err : float
        Uncertainty in nr
    nt : float
        Multiplier for converting relative transition intensity to transitions
        per 100 decays of the parent through the decay branch or to photons per 100
        neutron captures for (n,g).
    nt_err : float
        Uncertainty in nt
    br : float
        Branching ratio multiplier for converting intensity per 100 decays
        through this decay branch to intensity per 100 decays of the parent
        nuclide.
    br_err : float
        Uncertainty in br
    nb : float
        Multiplier for converting relative B- and EC intensities to intensities
        per 100 decays through this decay branch.
    nb_err : float
        Uncertainty in nb

    """
    nr, nr_err = _get_val_err(n_rec.group(2), n_rec.group(3))
    nt, nt_err = _get_val_err(n_rec.group(4), n_rec.group(5))
    br, br_err = _get_val_err(n_rec.group(6), n_rec.group(7))
    nb, nb_err = _get_val_err(n_rec.group(8), n_rec.group(9))
    if nr is not None and br is not None:
        nrbr = nr * br
    else:
        nrbr = None
    if nr_err is not None and br_err is not None:
        nrbr_err = np.sqrt((nr_err ** 2) + (br_err ** 2))
    else:
        nrbr_err = None
    return nr, nr_err, nt, nt_err, br, br_err, nb, nb_err, nrbr, nrbr_err


def _parse_production_normalization_record(np_rec):
    """
    This parses an ENSDF production normalization record

    Parameters
    ----------
    np_rec : re.MatchObject
        regular expression MatchObject

    Returns
    -------
    nrbr : float
        Multiplier for converting relative photon intensity to photons per 100
        decays of the parent nuclide
    nrbr_err : float
        Uncertainty in nrbr
    ntbr : float
        Multiplier for converting relative transition intensity to transitions
        per 100 decays of the parent nuclide
    ntbr_err : float
        Uncertainty in ntbr
    nbbr: float
        Multiplier for converting relative B- and EC intensities to intensity
        per 100 decays of the parent nuclide
    nbbr_err : float
        Uncertainty in nbbr
    """
    nrbr, nrbr_err = _get_val_err(np_rec.group(2), np_rec.group(3))
    ntbr, ntbr_err = _get_val_err(np_rec.group(4), np_rec.group(5))
    nbbr, nbbr_err = _get_val_err(np_rec.group(6), np_rec.group(7))
    return nrbr, nrbr_err, ntbr, ntbr_err, nbbr, nbbr_err


def _parse_parent_record(p_rec):
    """
    This parses an ENSDF parent record

    Parameters
    ----------
    p_rec : re.MatchObject
        regular expression MatchObject

    Returns
    -------
    tfinal : float
        half-life in seconds
    tfinalerr : float
        Uncertainty in half-life in seconds
    """
    e, e_err = _get_val_err(p_rec.group(2), p_rec.group(3))
    j = p_rec.group(4)
    tfinal, tfinalerr = _to_time(p_rec.group(5), p_rec.group(6))
    return tfinal, tfinalerr


def _parse_decay_dataset(lines, decay_s):
    """
    This parses a gamma ray dataset It returns a tuple of the data.

    Parameters
    ----------
    lines : list of str
        list containing lines from one dataset of an ensdf file
    decay_s : str
        string of the decay type

    Returns
    -------
    int
        nuc_id of the parent
    int
        nuc_id of the daughter
    str
        decay type
    float
        half-life in seconds
    float
        half-life error in seconds
    float
        Conversion factor for gamma intensity to photons per 100 decays of the
        parent
    float
        Error in conversion factor for gamma intensity
    numpy.ndarray
        a numpy array containing information about each gamma ray:
            * energy in keV
            * uncertainty in energy
            * intensity
            * uncertainty in intensity
            * electron conversion intensity
            * uncertainty in electron conversion intensity

    """
    gammarays = []
    ident = _ident.match(lines[0])
    daughter = ident.group(1)
    parent = ident.group(2).split()[0]
    tfinal = None
    tfinalerr = None
    nrbr = None
    nrbr_err = None
    for line in lines:
        g_rec = _g.match(line)
        if g_rec is not None:
            dat = _parse_gamma_record(g_rec)
            if not np.isnan(dat[0]):
                gammarays.append(dat)
        n_rec = _norm.match(line)
        if n_rec is not None:
            nr, nr_err, nt, nt_err, br, br_err, nb, nb_err, nrbr, nrbr_err = \
                _parse_normalization_record(n_rec)
        np_rec = _normp.match(line)
        if np_rec is not None:
            nrbr2, nrbr_err2, ntbr, ntbr_err, nbbr, nbbr_err = \
                _parse_production_normalization_record(np_rec)
            if nrbr2 is not None:
                nrbr = nrbr2
                nrbr_err = nrbr_err2
        p_rec = _p.match(line)
        if p_rec is not None:
            tfinal, tfinalerr = _parse_parent_record(p_rec)
    if len(gammarays) > 0:
        gammas = np.array(gammarays)
        return _to_id(parent), _to_id(daughter), decay_s.strip(), tfinal, tfinalerr, \
               nrbr, nrbr_err, gammas
    return None


def gamma_rays(f='ensdf.001'):
    """
    This splits an ENSDF file into datasets. It then passes the dataset to the
    appropriate parser. Currently only a subset of decay datasets are
    supported. The output is a list of objects containing information
    pertaining to a particular decay. This object is described in detail in the
    _parse_decay_dataset function.

    Parameters
    ----------
    filename : str
        Name of ENSDF formatted file

    Returns
    -------
    decaylist : list of tuples
        list of objects containing information pertaining to a particular
        decay. Contents of the tuple are described in the returns of the
        _parse_decay_dataset function.

    """
    if isinstance(f, str):
        with open(f, 'r') as f:
            dat = f.read()
    else:
        dat = f.read()
    decaylist = []
    datasets = dat.split(80 * " " + "\n")[0:-1]
    for dataset in datasets:
        lines = dataset.splitlines()
        ident = re.match(_ident, lines[0])
        if ident is not None:
            for decay_s in _decays:
                if decay_s in ident.group(2):
                    decay = _parse_decay_dataset(lines, decay_s)
                    if decay is not None:
                        decaylist.append(decay)
    return decaylist
