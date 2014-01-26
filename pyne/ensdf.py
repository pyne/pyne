import re
import urllib
import os
import warnings

import numpy as np

from pyne import nucname
from .utils import to_sec


_valexp = re.compile('([0-9.]*)([Ee][+-]\d*)')
_val = re.compile('(\d*)[.](\d*)')
_errpm = re.compile('[+](\d*)[-](\d*)')
_err = re.compile('[ ]*(\d*)')
_base = '([ \d]{3}[ A-Za-z]{2})'
_ident = re.compile(_base + '    (.{30})(.{26})(.{7})(.{6})')
_g = re.compile(_base + '  G (.{10})(.{2})(.{8})(.{2}).{24}(.{7})(.{2})(.{10})'
                + '(.{2})')
_gc = re.compile(_base + '[0-9A-Za-z] G (.{70})')
_beta = re.compile(_base + '  B (.{10})(.{2})(.{8})(.{2}).{10}(.{8})(.{6})')
_ec = re.compile(_base + '  E (.{10})(.{2})(.{8})(.{2})(.{8})(.{2})(.{8})(.{6})')
_p = re.compile(_base + '  P (.{10})(.{2})(.{18})(.{10})(.{6}).{9}(.{10})(.{2})'
                + '(.{4})')
_norm = re.compile(_base + '  N (.{10})(.{2})(.{8})(.{2})(.{8})(.{2})(.{8})'
                   + '(.{6})(.{7})(.{2})')
_normp = re.compile(_base +
                    ' PN (.{10})(.{2})(.{8})(.{2})(.{8})(.{2})(.{7})(.{2})')
_decays = ['B-', 'B+A', 'EC', 'B-A', 'B+', 'B+P', 'B-N', 'ECP', 'EC2P', 'N',
           '2N', 'IT', 'B+2P', 'B-2N', 'B+3P', 'ECA', 'P', '2P', '2B-', 'SF',
           'A', '2B+', '2EC', '14C']
_level_regex = re.compile(_base + '  L (.{10}).{20}(.{10}).{28}([ M])([ 1-9])')
_level_regex2 = re.compile(_base + '  L (.{10})(.{2})(.{18})(.{10})(.{6})(.{9})'
                           + '(.{10})(.{2})(.{1})([ M])([ 1-9])')
_level_cont_regex = re.compile('([ \d]{3}[ A-Za-z]{2})[0-9A-Za-z] L (.*)')


def _getvalue(obj, fn=float, rn=None):
    x = obj.strip()
    x = x.replace('$', '')
    x = x.replace('?', '')
    try:
        return fn(x)
    except ValueError:
        return rn


def _readpoint(line, dstart, dlen):
    data = _getvalue(line[dstart:dstart + dlen])
    error = _getvalue(line[dstart + dlen:dstart + dlen + 2])
    return data, error


def _read_variablepoint(line, dstart, dlen):
    sub = line[dstart:dstart + dlen + 2].split()
    data = None
    error = None
    if len(sub) == 2:
        data = _getvalue(sub[0])
        error = _getvalue(sub[1])
    return data, error


def _build_xray_table():
    i = 0
    j = 0
    dat = np.zeros((105, 26))
    medfile = os.path.join(os.path.dirname(__file__), 'mednew.dat')
    if not os.path.isfile(medfile):
        urllib.urlretrieve('http://www.nndc.bnl.gov/nndcscr/ensdf_pgm/'
                           + 'analysis/radlst/mednew.dat', medfile)
    with open(medfile, 'r') as f:
        lines = f.readlines()
    for line in lines:
        if (-1) ** i == 1:
            Z = int(line[0:3])
            k_shell_fluor, k_shell_fluor_error = _readpoint(line, 9, 6)
            l_shell_fluor, l_shell_fluor_error = _readpoint(line, 18, 6)
            #Probability of creating L-shell vacancy by filling K-shell vacancy
            prob, prob_error = _readpoint(line, 27, 6)
            k_shell_be, k_shell_be_err = _readpoint(line, 36, 8)
            li_shell_be, li_shell_be_err = _readpoint(line, 47, 8)
            mi_shell_be, mi_shell_be_err = _readpoint(line, 58, 8)
            ni_shell_be, ni_shell_be_err = _readpoint(line, 69, 8)
        else:
            Kb_to_Ka, Kb_to_Ka_err = _read_variablepoint(line, 9, 7)
            Ka2_to_Ka1, Ka2_to_Ka1_err = _read_variablepoint(line, 19, 7)
            L_auger = _getvalue(line[29:36])
            K_auger = _getvalue(line[36:42])
            Ka1_X_ray_en, Ka1_X_ray_en_err = _readpoint(line, 43, 8)
            Ka2_X_ray_en, Ka2_X_ray_en_err = _readpoint(line, 54, 7)
            Kb_X_ray_en = _getvalue(line[65:69])
            L_X_ray_en = _getvalue(line[70:76])
            dat[j] = Z, k_shell_fluor, k_shell_fluor_error, l_shell_fluor,\
                     l_shell_fluor_error, prob, k_shell_be, k_shell_be_err, \
                     li_shell_be, li_shell_be_err, mi_shell_be, mi_shell_be_err,\
                     ni_shell_be, ni_shell_be_err, \
                     Kb_to_Ka, Kb_to_Ka_err, Ka2_to_Ka1, Ka2_to_Ka1_err,\
                     L_auger, K_auger, Ka1_X_ray_en, Ka1_X_ray_en_err, \
                     Ka2_X_ray_en, Ka2_X_ray_en_err, Kb_X_ray_en, L_X_ray_en
            j += 1
        i += 1
    return dat


def _to_id(nuc, m=None, s=None):
    if not 'NN' in nuc:
        nucid = nucname.id(nuc.strip())
    else:
        warnings.warn('Neutron data not supported!')
        return None
    if m == 'M':
        state = s.strip()
        if 0 < len(state):
            state = int(state)
        else:
            state = 1
        nucid += state
    return nucid


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
            dat = dict([(_decay_to[key](from_nuc), _getvalue(val, rn=0) * 0.01)
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
    for nuc in nuclvl:
        nuclvl[nuc] = dict([(lvl, i) for i, lvl in enumerate(sorted(nuclvl[nuc])[:10])])
    data = [(row[0] + nuclvl[row[0]][row[1]],) + row[1:] for row in data \
            if (row[0] in nuclvl) and (row[1] in nuclvl[row[0]])]

    return data


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
    return _getvalue(''.join(errp) + valexp)


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
    from_nuc = _to_id(l_rec.group(1), l_rec.group(11), l_rec.group(12))
    if e is not None:
        e *= 1.0E-3
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
        if '%' in rx:
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
    dat = np.zeros(8)
    en, en_err = _get_val_err(g.group(2), g.group(3))
    inten, inten_err = _get_val_err(g.group(4), g.group(5))
    conv, conv_err = _get_val_err(g.group(6), g.group(7))
    tti, tti_err = _get_val_err(g.group(8), g.group(9))
    dat[:] = en, en_err, inten, inten_err, conv, conv_err, tti, tti_err
    return dat


def _parse_gamma_continuation_record(g, gammaray):
    """
    This parses an ENSDF gamma continuation record

    """
    conversions = {}
    entries = g.group(2).split('$')
    for entry in entries:
        entry = entry.replace('AP', '=')
        entry = entry.replace('EL1C+EL2C', 'LC')
        if 'C+' in entry:
            continue
        tsplit = entry.split('C')
        greff = gammaray[2]
        if '/T' in entry:
            tsplit = entry.split('/T')
            greff = gammaray[6]
            if np.isnan(greff):
                greff = gammaray[2]
        if len(tsplit) == 2:
            conv = None
            err = None
            contype = tsplit[0].lstrip('E')
            eff = tsplit[1].lstrip('= ').split()
            if len(eff) == 2:
                conv, err = _get_val_err(eff[0], eff[1])
            elif len(eff) == 1:
                conv = _getvalue(eff[0])
            if conv is None and not contype in conversions:
                conversions[contype] = (None, None)
            elif not contype in conversions:
                conversions[contype] = (conv * greff, err)
    return conversions


def _parse_beta_record(b_rec):
    """
    This parses an ENSDF beta minus record

    Parameters
    ----------
    b_rec : re.MatchObject
        regular expression MatchObject

    Returns
    -------
    en : float
        b- endpoint energy in keV
    en_err : float
        error in b- endpoint energy
    ib : float
        branch intensity
    dib : float
        error in branch intensity
    logft : float
        logft of the decay
    dft : float
        error in logft
    """
    en, en_err = _get_val_err(b_rec.group(2), b_rec.group(3))
    ib, dib = _get_val_err(b_rec.group(4), b_rec.group(5))
    logft, dft = _get_val_err(b_rec.group(6), b_rec.group(7))
    return en, en_err, ib, dib, logft, dft


def _parse_ec_record(e_rec):
    """
    This parses an ENSDF electron capture + b+ record

    Parameters
    ----------
    e_rec : re.MatchObject
        regular expression MatchObject

    Returns
    -------
    en : float
        b+ endpoint energy in keV
    en_err : float
        error in b+ endpoint energy
    ib : float
        b+ branch intensity
    dib : float
        error in b+ branch intensity
    ie : float
        ec branch intensity
    die : float
        error in ec branch intensity
    logft : float
        logft of the decay
    dft : float
        error in logft
    """
    en, en_err = _get_val_err(e_rec.group(2), e_rec.group(3))
    ib, dib = _get_val_err(e_rec.group(4), e_rec.group(5))
    ie, die = _get_val_err(e_rec.group(6), e_rec.group(7))
    logft, dft = _get_val_err(e_rec.group(8), e_rec.group(9))
    return en, en_err, ib, dib, ie, die, logft, dft


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
    return tfinal, tfinalerr, e, e_err


def _update_xrays(conv, xrays, nuc_id):
    """
    Update X-ray data for a given decay
    """
    z = nucname.znum(nuc_id)
    xka1 = 0
    xka2 = 0
    xkb = 0
    xl = 0
    if 'K' in conv and conv['K'][0] is not None and not np.isnan(conv['K'][0]):
        xk = _xraydat[z - 1, 1] * conv['K'][0]
        xka = xk / (1.0 + _xraydat[z - 1, 14])
        xka1 = xka / (1.0 + _xraydat[z - 1, 16])
        xka2 = xka - xka1
        xkb = xk - xka
        if 'L' in conv and conv['L'][0] is not None and not np.isnan(conv['L'][0]):
            xl = _xraydat[z - 1, 3] * (conv['L'][0] + conv['K'][0] * _xraydat[z - 1, 5])
        else:
            xl = 0
    elif 'L' in conv and conv['L'][0] is not None and not np.isnan(conv['L'][0]):
        xl = _xraydat[z - 1, 3] * (conv['L'][0])

    xrays = np.array([_xraydat[z - 1, 20], xka1 + xrays[1], _xraydat[z - 1, 22], xka2 + xrays[3],
                      _xraydat[z - 1, 24], xkb + xrays[5], _xraydat[z - 1, 25], xl + xrays[7]])
    return xrays


def _parse_decay_dataset(lines, decay_s):
    """
    This parses a gamma ray dataset. It returns a tuple of the parsed data.

    Parameters
    ----------
    lines : list of str
        list containing lines from one dataset of an ensdf file
    decay_s : str
        string of the decay type

    Returns
    -------
    Tuple of decay parameters which is described in detail in gamma_rays docs

    """
    gammarays = []
    conv = {}
    ident = _ident.match(lines[0])
    daughter = ident.group(1)
    parent = ident.group(2).split()[0]
    tfinal = None
    tfinalerr = None
    nrbr = None
    nrbr_err = None
    xrays = np.zeros(8)
    gamma = 0
    goodgray = False

    for line in lines:
        g_rec = _g.match(line)
        if g_rec is not None:
            dat = _parse_gamma_record(g_rec)
            if gamma == 1:
                gamma = 0
                #store old atomic data
                xrays = _update_xrays(conv, xrays, _to_id(daughter))
                #clear conversion data
                conv = {}

            if not np.isnan(dat[0]):
                gammarays.append(dat)
                goodgray = True
            else:
                goodgray = False
            continue
        gc_rec = _gc.match(line)
        if gc_rec is not None and goodgray is True:
            conv_temp = conv
            conv = _parse_gamma_continuation_record(gc_rec, gammarays[-1])
            if gamma == 0:
                gamma = 1
            conv.update(conv_temp)
            continue
        n_rec = _norm.match(line)
        if n_rec is not None:
            nr, nr_err, nt, nt_err, br, br_err, nb, nb_err, nrbr, nrbr_err = \
                _parse_normalization_record(n_rec)
            continue
        np_rec = _normp.match(line)
        if np_rec is not None:
            nrbr2, nrbr_err2, ntbr, ntbr_err, nbbr, nbbr_err = \
                _parse_production_normalization_record(np_rec)
            if nrbr2 is not None:
                nrbr = nrbr2
                nrbr_err = nrbr_err2
            continue
        p_rec = _p.match(line)
        if p_rec is not None:
            tfinal, tfinalerr, e, e_err = _parse_parent_record(p_rec)
            continue
    if len(gammarays) > 0:
        gammas = np.array(gammarays)
        if gamma == 1:
            xrays = _update_xrays(conv, xrays, _to_id(daughter))
        pfinal = []
        parent = parent.split('(')[0]
        parents = parent.split(',')
        if len(parents) > 1:
            for item in parents:
                pfinal.append(_to_id(item))
        else:
            pfinal = _to_id(parents[0])
        return pfinal, _to_id(daughter), decay_s.strip(), tfinal, tfinalerr, \
               nrbr, nrbr_err, xrays, gammas
    return None


def origen_data(filename='ensdf.001'):
    """
    This function parses assorted data from an ensdf file in order to collect
    the necessary information to generate data for origen input decks.

    Parameters
    ----------
    filename : str
        Name of ENSDF formatted file

    Returns
    -------
    decaylist : list
        This is a list of tuples containing:
            * decay type
            * half life
            * parent nuc_id
            * daughter nuc_id
            * branching ratio
            * energy level of parent
            * whether the daughter has a long half life
    branchlist : list
    pid, level, half_lifev, brs
        This is a list of tuples containing:
            * parent nuc_id
            * energy level of parent
            * half life of parent
            * dictionary of branching ratios
    """
    if isinstance(filename, str):
        with open(filename, 'r') as f:
            dat = f.read()
    else:
        dat = filename.read()
    datasets = dat.split(80 * " " + "\n")[0:-1]
    decaylist = []
    branchlist = []
    for dataset in datasets:
        lines = dataset.splitlines()
        ident = re.match(_ident, lines[0])
        if ident is not None:
            if 'DECAY' in ident.group(2):
                daughter = ident.group(1)
                parent = ident.group(2).split()[0]
                dtype = ident.group(2).split()[1]
                tfinal = 0
                e = 0.0
                ie = None
                ib = None
                nb = None
                br = None
                newlevel = False
                for line in lines:
                    b_rec = _beta.match(line)
                    if b_rec is not None:
                        en, en_err, ib, dib, logft, dft = _parse_beta_record(b_rec)
                        continue
                    e_rec = _ec.match(line)
                    if e_rec is not None:
                        en, en_err, ib, dib, ie, die, logft, dft = _parse_ec_record(e_rec)
                        continue
                    n_rec = _norm.match(line)
                    if n_rec is not None:
                        nr, nr_err, nt, nt_err, br, br_err, nb, nb_err, nrbr, nrbr_err = \
                            _parse_normalization_record(n_rec)
                        continue
                    p_rec = _p.match(line)
                    if p_rec is not None:
                        tfinal, tfinalerr, e, e_err = _parse_parent_record(p_rec)
                        continue
                    level_l = _level_regex2.match(line)
                    if level_l is not None:
                        if newlevel:
                            #save old level data
                            if (ib is not None or ie is not None) and nb is not None and br is not None:
                                if ib is None:
                                    ib = 0.0
                                if ie is None:
                                    ie = 0.0
                                decaylist.append((_to_id(parent), tfinal, e, half_lifev, level, dtype, (ib+ie), nb, br))
                        level, half_lifev, from_nuc = _parse_level_record(level_l)
                        newlevel = True
                        continue
                if newlevel:
                    #save old level data
                    if (ib is not None or ie is not None) and nb is not None and br is not None:
                        if ib is None:
                            ib = 0.0
                        if ie is None:
                            ie = 0.0
                        decaylist.append((_to_id(parent), tfinal, e, half_lifev, level, dtype, (ib+ie), nb, br))
            if 'ADOPTED LEVELS' in ident.group(2):
                parent = ident.group(1)
                pid = _to_id(parent)
                brs = {}
                levelc_found = False
                for line in lines:
                    level_l = _level_regex2.match(line)
                    if level_l is not None:
                        if levelc_found and half_lifev is not None:
                            levelc_found = False
                            if len(brs) > 0:
                                branchlist.append((pid, level, half_lifev, brs))
                            brs = {}
                        level, half_lifev, from_nuc = _parse_level_record(level_l)
                        continue
                    levelc = _level_cont_regex.match(line)
                    if levelc is not None:
                        brs.update(_parse_level_continuation_record(levelc))
                        levelc_found = True
                        continue
                if levelc_found and half_lifev is not None and pid is not None:
                    if len(brs) > 0:
                        branchlist.append((pid, level, half_lifev, brs))
    return decaylist, branchlist


def _dlist_gen(f='ensdf.001'):
    """
    This compiles a list of decay types in an ensdf file

    Parameters
    ----------
    f : str
        Name of ENSDF formatted file

    Returns
    -------
    decaylist : list
        list of decay types in the ENSDF file eg. ['B+','B-','A']
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
            if 'DECAY' in ident.group(2):
                #print ident.group(2)
                fin = ident.group(2).split()[1]
                if not fin in decaylist:
                    decaylist.append(fin)

    return decaylist


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
        decay. This information is in the following format:
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
        X-ray energies and intensities in the following format:
            * K_alpha1 energy
            * K_alpha1 intensity (multiply by conversion factor for percentage)
            * K_alpha2 energy
            * K_alpha2 intensity (multiply by conversion factor for percentage)
            * K_beta energy
            * K_beta intensity (multiply by conversion factor for percentage)
            * L energy
            * L intensity (multiply by conversion factor for percentage)
    numpy.ndarray
        a numpy array containing information about each gamma ray:
            * energy in keV
            * uncertainty in energy
            * intensity (multiply by conversion factor for percentage)
            * uncertainty in intensity
            * electron conversion intensity
            * uncertainty in electron conversion intensity
            * total transition intensity
            * total transition intensity error

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
                if 'DECAY' in ident.group(2) and decay_s == ident.group(2).split()[1]:
                    decay = _parse_decay_dataset(lines, decay_s)
                    if decay is not None:
                        decaylist.append(decay)
    return decaylist


_xraydat = _build_xray_table()
