from __future__ import division
import re
import sys
import copy

try:
    from collections.abc import defaultdict
except ImportError:
    from collections import defaultdict
from warnings import warn
from pyne.utils import QA_warn, time_conv_dict

import numpy as np

from pyne import nucname, rxname, data

if sys.version_info[0] > 2:
    basestring = str

QA_warn(__name__)

_valexp = re.compile("([0-9.]*)([Ee][+-]?\d*)")
_val = re.compile("(\d*)[.](\d*)")
_specialval = re.compile("([0-9. ]*)[+]([A-Z])")
_specialval2 = re.compile("([A-Z]*)[+]([0-9.]*)")
_errpm = re.compile("[+](\d*)[-](\d*)")
_err = re.compile("[ ]*(\d*)")
_base = "([ \d]{3}[ A-Za-z]{2})"
_ident = re.compile(_base + "    (.{30})(.{26})(.{7})(.{6})")
_g = re.compile(
    _base + "  G (.{10})(.{2})(.{8})(.{2}).{24}(.{7})(.{2})(.{10})" + "(.{2})"
)
_gc = re.compile(_base + "[0-9A-Za-z] [GE] (.{70})")
_beta = re.compile(_base + "  B (.{10})(.{2})(.{8})(.{2}).{10}(.{8})(.{6})")
_betac = re.compile(_base + "[0-9A-Za-z] ([BE]) (.{70})")
_ec = re.compile(
    _base + "  E (.{10})(.{2})(.{8})(.{2})" + "(.{8})(.{2})(.{8})(.{6})(.{10})(.{2})"
)
_p = re.compile(
    _base + "  P (.{10})(.{2})(.{18})(.{10})" + "(.{6}).{9}(.{10})(.{2})(.{4})"
)
_norm = re.compile(
    _base + "  N (.{10})(.{2})(.{8})(.{2})(.{8})(.{2})(.{8})" + "(.{6})(.{7})(.{2})"
)
_normp = re.compile(_base + " PN (.{10})(.{2})(.{8})(.{2})(.{8})(.{2})(.{7})(.{2})")
_q = re.compile(_base + "  Q (.{10})(.{2})(.{8})(.{2})" + "(.{8})(.{2})(.{8})(.{6})")
_alpha = re.compile(_base + "  A (.{10})(.{2})(.{8})(.{2})(.{8})(.{2})")
_dp = re.compile(_base + "  D(.{1})(.{10})(.{2})(.{8})(.{2})(.{8})(.{10})" + "(.{6})")
_decays = [
    "B-",
    "B+A",
    "EC",
    "B-A",
    "B+",
    "B+P",
    "B-N",
    "ECP",
    "EC2P",
    "N",
    "2N",
    "IT",
    "B+2P",
    "B-2N",
    "B+3P",
    "ECA",
    "P",
    "2P",
    "2B-",
    "SF",
    "A",
    "2B+",
    "2EC",
    "14C",
]
_level_regex = re.compile(
    _base
    + "  L (.{10})(.{2})(.{18})(.{10})(.{6})"
    + "(.{9})(.{10})(.{2})(.{1})([ M])([ 1-9])"
)
_level_cont_regex = re.compile("([ \d]{3}[ A-Za-z]{2})[0-9A-Za-z] L (.*)")


def _getvalue(obj, fn=float, rn=None):
    x = obj.strip()
    x = x.replace("$", "")
    x = x.replace("?", "")
    try:
        return fn(x)
    except ValueError:
        return rn


def _to_id(nuc):
    if "NN" not in nuc:
        nucid = nucname.ensdf_to_id(nuc.strip())
    else:
        warn("Neutron data not supported!")
        return 0
    return nucid


# Energy to half-life conversion:  T1/2= ln(2) × (h/2 pi) / energy
# See http://www.nndc.bnl.gov/nudat2/help/glossary.jsp#halflife
# NIST CODATA https://physics.nist.gov/cgi-bin/cuu/Value?hbar
#    h-bar = 1.054 571 800(13) x 1e-34 J
#    1 J = 6.241 509 126(38) x 1e18 eV
HBAR_LN2 = 4.5623775832376968e-16  # h-bar ln(2) in eV s
energy_conv_dict = {
    "ev": HBAR_LN2,
    "kev": 1e-3 * HBAR_LN2,
    "mev": 1e-6 * HBAR_LN2,
}


def _halflife_to_seconds(value, err, units):
    """Converts a halflife with err and units to seconds.

    Parameters
    ----------
    value: number
        Time or energy, depending on units.
    err : number or (number, number)
        Uncertainty, or (plus, minus) uncertainty in [units].
    units : str
        Units flag, eg 'min', 'ms', 'days', or even 'MeV'.

    Returns
    -------
    sec_time : float
        Time value in [sec].
    sec_err : None or float or (float, float) in [sec].
        Time uncertainty in [sec], or (plus, minus) if asymmetric uncertainty.
    """
    if err is None:
        plus, minus = 0, 0
    elif np.isscalar(err):
        plus, minus = err, err
    else:
        plus, minus = err

    units = units.lower()
    scale = time_conv_dict.get(units, None)
    if scale is not None:
        sec_time = scale * value
        sec_err = (scale * plus, scale * minus)
    else:
        scale = energy_conv_dict[units]
        sec_time = scale / value
        sec_err = (
            scale / max(0.1 * value, value - minus) - sec_time,
            sec_time - scale / (value + plus),
        )
    if err is None:
        return sec_time, None
    elif sec_err[0] == sec_err[1]:
        return sec_time, sec_err[0]
    else:
        return sec_time, sec_err


def _to_time(tstr, errstr):
    t = tstr.strip()
    # This accepts questionable levels
    t = t.replace("?", "")
    tobj = [s.strip(" ()") for s in t.split()]
    if len(tobj) == 2:
        t, t_unit = tobj
        value, err = _get_val_err(t, errstr)
        tfinal, tfinalerr = _halflife_to_seconds(value, err, t_unit)
    elif "STABLE" in t:
        tfinal = np.inf
        tfinalerr = None
    else:
        tfinal = None
        tfinalerr = None
    return tfinal, tfinalerr


def _get_val_err(valstr, errstr):
    pm = _errpm.match(errstr)
    err = _err.match(errstr)
    if pm is None and err.group(1) == "":
        return _getvalue(valstr), None
    val = _valexp.match(valstr)
    if val is None:
        valexp = ""
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
    errp.insert(-plen, ".")
    return _getvalue("".join(errp) + valexp)


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
        Level energy in keV
    tfinal : float
        Half life in seconds
    from_nuc : int
        nuc id of nuclide
    state : int
        metastable state of level
    special : str
        A-Z character denoting a group of known levels with no reference
        to the ground state. P and N are special characters reserved for
        proton and neutron resonances given in center of mass system energy.
    """
    lm = re.match("[ ]*([A-Z]+)(?![A-Z0-9+])", l_rec.group(2))
    spv = _specialval.match(l_rec.group(2).strip())
    spv2 = _specialval2.match(l_rec.group(2).strip())
    special = " "
    if lm is not None:
        special = lm.group(1)
        if "S" in special and len(special.strip()) > 1:
            special = special.strip()[1]
        e = 0.0
        de = np.nan
    elif spv is not None:
        e, de = _get_val_err(spv.group(1), l_rec.group(3))
        special = spv.group(2)
    elif spv2 is not None:
        e, de = _get_val_err(spv2.group(2), l_rec.group(3))
        special = spv2.group(1)
        if "S" in special and len(special.strip()) > 1:
            special = special.strip()[1]
    else:
        e, de = _get_val_err(l_rec.group(2).strip("() "), l_rec.group(3))
    tfinal, tfinalerr = _to_time(l_rec.group(5), l_rec.group(6))
    from_nuc = _to_id(l_rec.group(1))
    m = l_rec.group(11)
    s = l_rec.group(12)
    state = 0
    if m == "M":
        state = s.strip()
        if 0 < len(state):
            state = int(state)
        else:
            state = 1
    return e, tfinal, from_nuc, state, special


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
    raw_children = g[-1].replace(" AP ", "=")
    raw_children = raw_children.replace("$", " ").split()
    for raw_child in raw_children:
        if "=" in raw_child:
            rx, br = raw_child.split("=")[:2]
            br = br.strip()
        else:
            continue
        if "%" in rx and "?" not in br and len(br) > 0:
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
    en, en_err = _get_val_err(g.group(2), g.group(3))
    inten, inten_err = _get_val_err(g.group(4), g.group(5))
    conv, conv_err = _get_val_err(g.group(6), g.group(7))
    tti, tti_err = _get_val_err(g.group(8), g.group(9))
    return [en, en_err, inten, inten_err, conv, conv_err, tti, tti_err]


def _parse_gamma_continuation_record(g, inten, tti):
    """
    This parses an ENSDF gamma continuation record

    """
    conversions = {}
    entries = g.group(2).split("$")
    for entry in entries:
        entry = entry.replace("AP", "=")
        entry = entry.replace("EL1C+EL2C", "LC")
        if "+=" in entry or "EAV" in entry:
            continue
        if "C=" in entry:
            tsplit = entry.split("C")
        else:
            tsplit = entry.split("=")
            tsplit[0] = tsplit[0].lstrip("C")
        greff = inten
        if "/T" in entry:
            tsplit = entry.split("/T")
            greff = tti
            if greff is None:
                greff = inten
        if greff is None:
            greff = 1.0
        if len(tsplit) == 2:
            conv = None
            err = None
            contype = tsplit[0].lstrip("E")
            eff = tsplit[1].lstrip("= ").split()
            if len(eff) == 2:
                conv, err = _get_val_err(eff[0], eff[1])
            elif len(eff) == 1:
                conv = _getvalue(eff[0])
            if conv is None and contype not in conversions:
                conversions[contype] = (None, None)
            elif contype not in conversions:
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


def _parse_beta_continuation_record(bc_rec):
    """
    This parse the beta continuation record for EAV
    """
    entries = bc_rec.group(3).split("$")
    eav = None
    eav_err = None
    for entry in entries:
        if "EAV" in entry and "=" in entry:
            dat = entry.split("=")[1]
            dat = dat.split()
            if len(dat) == 2:
                eav, eav_err = _get_val_err(dat[0], dat[1])
            elif len(dat) == 1:
                eav = _getvalue(dat[0])
    return eav, eav_err


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
    tti, dtti = _get_val_err(e_rec.group(10), e_rec.group(11))
    return en, en_err, ib, dib, ie, die, logft, dft, tti, dtti


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
        per 100 decays of the parent through the decay branch or to photons
        per 100 neutron captures for (n,g).
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
        nrbr_err = nrbr * np.sqrt((br_err / br) ** 2 * (nr_err / nr) ** 2)
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
    lm = re.match("[ ]*([A-Z]+)(?![A-Z0-9+])", p_rec.group(2))
    spv = _specialval.match(p_rec.group(2).strip())
    spv2 = _specialval2.match(p_rec.group(2).strip())
    special = " "
    if lm is not None:
        special = lm.group(1)
        if "S" in special and len(special.strip()) > 1:
            special = special.strip()[1]
        e = 0.0
        de = np.nan
    elif spv is not None:
        e, de = _get_val_err(spv.group(1), p_rec.group(3))
        special = spv.group(2)
    elif spv2 is not None:
        e, de = _get_val_err(spv2.group(2), p_rec.group(3))
        special = spv2.group(1)
        if "S" in special and len(special.strip()) > 1:
            special = special.strip()[1]
    else:
        e, de = _get_val_err(p_rec.group(2).strip("() "), p_rec.group(3))
    j = p_rec.group(4)
    tfinal, tfinalerr = _to_time(p_rec.group(5), p_rec.group(6))
    return p_rec.group(1), tfinal, tfinalerr, e, de, special


def _parse_qvalue_record(q_rec):
    """
    This parses and ENSDF q-value record

    Parameters
    ----------
    q_rec : re.MatchObject
        regular expression MatchObject

    Returns
    -------
    qminus : float
        total energy for B- decay (if qminus > 0 B- decay is possible)
    dqminus : float
        standard uncertainty in qminus
    sn : float
        neutron separation energy in keV
    dsn : float
        standard uncertainty in sn
    sp : float
        neutron separation energy in keV
    dsp : float
        standard uncertainty in sp
    qa : float
        total energy available for alpha decay of the ground state
    dqa : float
        standard uncertainty in qa
    """
    qminus, dqminus = _get_val_err(q_rec.group(2), q_rec.group(3))
    sn, dsn = _get_val_err(q_rec.group(4), q_rec.group(5))
    sp, dsp = _get_val_err(q_rec.group(5), q_rec.group(7))
    qa, dqa = _get_val_err(q_rec.group(8), q_rec.group(9))
    return qminus, dqminus, sn, dsn, sp, dsp, qa, dqa


def _parse_alpha_record(a_rec):
    """
    This parses and ENSDF alpha record

    Parameters
    ----------
    q_rec : re.MatchObject
        regular expression MatchObject

    Returns
    -------
    e : float
        energy of alpha particle
    de : float
        standard uncertainty in energy
    ia : float
        intensity of the decay branch in percent
    dia : float
        standard uncertainty in intensity
    hf : float
        hindrance factor
    dhf : float
        standard uncertainty in hindrance factor
    """
    e, de = _get_val_err(a_rec.group(2), a_rec.group(3))
    ia, dia = _get_val_err(a_rec.group(4), a_rec.group(5))
    hf, dhf = _get_val_err(a_rec.group(5), a_rec.group(7))
    return e, de, ia, dia, hf, dhf


def _parse_delayed_particle_record(dp_rec):
    """
    This parses and ENSDF delayed particle record

    Parameters
    ----------
    dp_rec : re.MatchObject
        regular expression MatchObject

    Returns
    -------
    ptype : str
        symbol for delayed particle
    e : float
        particle energy
    de : float
        standard uncertainty in energy
    ip : float
        intensity of delayed particle in percent
    dip : float
        standard uncertainty in intensity
    ei : float
        energy level of the intermediate
    t : float
        half-life of the transition (in seconds)
    dt : float
        standard uncertainty in half-life
    """
    ptype = dp_rec.group(2)
    e, de = _get_val_err(dp_rec.group(3), dp_rec.group(4))
    ip, dip = _get_val_err(dp_rec.group(5), dp_rec.group(6))
    ei = _getvalue(dp_rec.group(7))
    t, dt = _to_time(dp_rec.group(8), dp_rec.group(9))
    return ptype, e, de, ip, dip, ei, t, dt


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
    betas = []
    alphas = []
    ecbp = []
    ident = _ident.match(lines[0])
    daughter = ident.group(1)
    daughter_id = abs(_to_id(daughter))
    parent = ident.group(2).split()[0]
    parent = parent.split("(")[0]
    parents = parent.split(",")
    if len(parents) > 1:
        pfinal = abs(_to_id(parents[0]))
    else:
        pfinal = abs(_to_id(parents[0][:5]))
    tfinal = None
    tfinalerr = None
    nrbr = None
    nbbr = None
    nrbr_err = None
    nbbr_err = None
    nb_err = None
    br_err = None
    nb = None
    br = None
    level = None
    special = " "
    goodgray = False
    parent2 = None
    for line in lines:
        level_l = _level_regex.match(line)
        if level_l is not None:
            level, half_lifev, from_nuc, state, special = _parse_level_record(level_l)
            continue
        b_rec = _beta.match(line)
        if b_rec is not None:
            dat = _parse_beta_record(b_rec)
            if parent2 is None:
                bparent = pfinal
            else:
                bparent = parent2
            level = 0.0 if level is None else level
            bdaughter = abs(data.id_from_level(_to_id(daughter), level))
            betas.append([bparent, bdaughter, dat[0], 0.0, dat[2]])
        bc_rec = _betac.match(line)
        if bc_rec is not None:
            bcdat = _parse_beta_continuation_record(bc_rec)
            if bcdat[0] is not None:
                if bc_rec.group(2) == "B":
                    betas[-1][3] = bcdat[0]
                else:
                    ecbp[-1][3] = bcdat[0]
                    bggc = _gc.match(line)
                    conv = _parse_gamma_continuation_record(bggc, dat[2], dat[8])
                    if "K" in conv:
                        ecbp[-1][-3] = conv["K"][0]
                    if "L" in conv:
                        ecbp[-1][-2] = conv["L"][0]
                    if "M" in conv:
                        ecbp[-1][-1] = conv["M"][0]
        a_rec = _alpha.match(line)
        if a_rec is not None:
            dat = _parse_alpha_record(a_rec)
            if parent2 is None:
                aparent = pfinal
            else:
                aparent = parent2
            level = 0.0 if level is None else level
            adaughter = abs(data.id_from_level(_to_id(daughter), level))
            alphas.append([aparent, adaughter, dat[0], dat[2]])
        ec_rec = _ec.match(line)
        if ec_rec is not None:
            dat = _parse_ec_record(ec_rec)
            if parent2 is None:
                ecparent = pfinal
            else:
                ecparent = parent2
            level = 0.0 if level is None else level
            ecdaughter = abs(data.id_from_level(_to_id(daughter), level))
            ecbp.append([ecparent, ecdaughter, dat[0], 0.0, dat[2], dat[4], 0, 0, 0])
            continue
        g_rec = _g.match(line)
        if g_rec is not None:
            dat = _parse_gamma_record(g_rec)
            if dat[0] is not None:
                gparent = 0
                gdaughter = 0
                if level is not None:
                    gparent = abs(data.id_from_level(_to_id(daughter), level, special))
                    dlevel = level - dat[0]
                    gdaughter = abs(
                        data.id_from_level(_to_id(daughter), dlevel, special)
                    )
                if parent2 is None:
                    gp2 = pfinal
                else:
                    gp2 = parent2
                dat.insert(0, daughter_id)
                dat.insert(0, gp2)
                dat.insert(0, gdaughter)
                dat.insert(0, gparent)
                for i in range(3):
                    dat.append(0)
                gammarays.append(dat)
                goodgray = True
            else:
                goodgray = False
            continue
        gc_rec = _gc.match(line)
        if gc_rec is not None and goodgray is True:
            conv = _parse_gamma_continuation_record(
                gc_rec, gammarays[-1][6], gammarays[-1][10]
            )
            if "K" in conv:
                gammarays[-1][-3] = conv["K"][0]
            if "L" in conv:
                gammarays[-1][-2] = conv["L"][0]
            if "M" in conv:
                gammarays[-1][-1] = conv["M"][0]
            continue
        n_rec = _norm.match(line)
        if n_rec is not None:
            (
                nr,
                nr_err,
                nt,
                nt_err,
                br,
                br_err,
                nb,
                nb_err,
                nrbr,
                nrbr_err,
            ) = _parse_normalization_record(n_rec)
            if nb is not None and br is not None:
                nbbr = nb * br
            if nb_err is not None and br_err is not None and nb_err != 0:
                nbbr_err = nbbr * ((br_err / br) ** 2 * (nb_err / nb) ** 2) ** 0.5
            continue
        np_rec = _normp.match(line)
        if np_rec is not None:
            (
                nrbr2,
                nrbr_err2,
                ntbr,
                ntbr_err,
                nbbr2,
                nbbr_err2,
            ) = _parse_production_normalization_record(np_rec)
            if nrbr2 is not None and nrbr is None:
                nrbr = nrbr2
                nrbr_err = nrbr_err2
            if nbbr2 is not None and nbbr is None:
                nbbr = nbbr2
                nbbr_err = nbbr_err2
            continue
        p_rec = _p.match(line)
        if p_rec is not None:
            # only 2 parents are supported so this can be here
            multi = False
            if parent2 is not None:
                multi = True
                pfinal = [
                    parent2,
                ]
                tfinal = [
                    t,
                ]
                tfinalerr = [
                    terr,
                ]
            parent2, t, terr, e, e_err, special = _parse_parent_record(p_rec)
            parent2 = abs(data.id_from_level(_to_id(parent2), e, special))
            if terr is not None and not isinstance(terr, float):
                terr = (terr[0] + terr[1]) / 2.0
            if multi:
                tfinal.append(t)
                tfinalerr.append(terr)
                pfinal.append(parent2)
            else:
                tfinal = t
                tfinalerr = terr
                pfinal = parent2
            continue
    if len(gammarays) > 0 or len(alphas) > 0 or len(betas) > 0 or len(ecbp) > 0:
        if len(parents) > 1 and parent2 is None:
            pfinal = []
            for item in parents:
                pfinal.append(_to_id(item))
        return (
            pfinal,
            daughter_id,
            rxname.id(decay_s.strip().lower()),
            tfinal,
            tfinalerr,
            br,
            br_err,
            nrbr,
            nrbr_err,
            nbbr,
            nbbr_err,
            gammarays,
            alphas,
            betas,
            ecbp,
        )
    return None


_BAD_RX = frozenset(
    [
        # Be-6 doesn't really alpha decay (leaving He-2), rather it emits 2p
        (40060000, 1089),
        # Li-8 -> He-4 + beta- + alpha is really a shortcut for
        # Li-8 -> Be-8 + beta- -> He-4 + alpha
        (30080000, 1355894000),
    ]
)


def _adjust_ge100_branches(levellist):
    """This adjust branches that are greater than or equal to 100% to be
    100% - sum(other branches).  This helps prevent unphysical errors
    downstream.
    """
    n = len(levellist)
    brsum = defaultdict(float)
    bridx = defaultdict(lambda: (-1, -1.0))
    baddies = []
    for i, (nuc, rx, hl, lvl, br, ms, sp) in enumerate(levellist):
        if rx == 0:
            continue
        if br >= bridx[nuc][1]:
            bridx[nuc] = (i, br)
        brsum[nuc] += br
        nucrx = (nuc, rx)
        if nucrx in _BAD_RX:
            baddies.append(i)
    # adjust branch ratios
    for nuc, (i, br) in bridx.items():
        row = levellist[i]
        # this line ensures that all branches sum to 100.0 within floating point
        new_br = 100.0 - brsum[nuc] + br
        new_row = row[:4] + (new_br,) + row[5:]
        levellist[i] = new_row
    # remove bad reaction rows
    for i in baddies[::-1]:
        del levellist[i]


# State Id, Bad Metastable Number, (Replacement State ID, optional) Replacement Metastable Number
_BAD_METASTABLES = {
    # Rh-110 misreports its ground state as a first meta-stable and its first
    # metastable as its second.
    (451100000, 1): 0,
    (451100001, 2): 1,
    # Pm-154 misreports its ground state as a first metastable
    (611540000, 1): 0,
    # Ga-72M is not listed as metastable
    (310720002, 0): 1,
    # Rh-108M is not listed as metastable
    (451080004, 0): 1,
    # Pm-136 mislabels two states as both metastable or ground.
    # Replacing with what KAERI and NNDC report
    (611360001, 2): (611360000, 0),
    (611360000, 1): (611360001, 1),
}


def _adjust_metastables(levellist):
    """Adjusts misreported metastable states in place."""
    for i in range(len(levellist)):
        key = (levellist[i][0], levellist[i][5])
        if key in _BAD_METASTABLES:
            row = list(levellist[i])
            new_id = _BAD_METASTABLES[key]
            if not isinstance(new_id, int):
                row[0], new_id = new_id
            row[5] = new_id
            levellist[i] = tuple(row)


# State Id, Rx Id : New Half-lives
_BAD_HALF_LIVES = {
    # Eu-151 lists a very long half-life (5.364792e+25) even though it
    # lists no reaction, and thus no children, and no branch ratio.
    # set to infinity for consistency.
    (631510000, 0): float("inf"),
}


def _adjust_half_lives(levellist):
    """Resets misbehaving half-lives to new value."""
    for i in range(len(levellist)):
        key = levellist[i][:2]
        if key in _BAD_HALF_LIVES:
            row = list(levellist[i])
            row[2] = _BAD_HALF_LIVES[key]
            levellist[i] = tuple(row)


def levels(filename, levellist=None):
    """
    This takes an ENSDF filename or file object and parses the ADOPTED LEVELS
    records to assign level numbers by energy. It also parses the different
    reported decay types and branching ratios.

    Parameters
    ----------
    filename : str or file
        Name of ENSDF formatted file or a file-like object containing ENSDF
        formatted data
    levellist : list of tuples
        This is a list object which all newly processed levels will be added
        to. If it's None a new one will be created.

    Returns
    -------
    levellist : list of tuples
        This is a list of all the level data. Each level has base entry with a
        reaction id of 0 and additional entries for any listed decays. The
        format of each row is:
        nuc_id : int
            The state_id of the level
        rx_id : int
            The id of the decay "reaction" in PyNE reaction id form.
        half_life : float
            Half life of the state in s
        level : float
            energy of the level in keV
        branch_ratio : float
            if rx_id != 0 this is the percent of decays in that channel
        metastable : int
            metastable id number of the level (if given)
        special : string
            single character denoting levels with unknown relation to ground
            state
    """
    badlist = [
        "ecsf",
        "34si",
        "|b{+-}fission",
        "{+24}ne",
        "{+22}ne",
        "24ne",
        "b-f",
        "{+20}o",
        "2|e",
        "b++ec",
        "ecp+ec2p",
        "ecf",
        "mg",
        "ne",
        "{+20}ne",
        "{+25}ne",
        "{+28}mg",
        "sf(+ec+b+)",
    ]
    special = ""
    if levellist is None:
        levellist = []
    if isinstance(filename, str):
        with open(filename, "r") as f:
            dat = f.read()
    else:
        dat = filename.read()
    datasets = dat.split(80 * " " + "\n")[0:-1]
    for dataset in datasets:
        lines = dataset.splitlines()
        ident = re.match(_ident, lines[0])
        if ident is None:
            continue
        if "ADOPTED LEVELS" in ident.group(2):
            leveln = 0
            brs = {}
            level_found = False
            for line in lines:
                level_l = _level_regex.match(line)
                if level_l is not None:
                    if len(brs) > 0:
                        for key, val in brs.items():
                            goodkey = True
                            keystrip = key.replace("%", "").lower()
                            for item in badlist:
                                if keystrip == item:
                                    goodkey = False
                            if goodkey is True:
                                rx = rxname.id(keystrip)
                                branch_percent = float(val.split("(")[0])
                                levellist.append(
                                    (
                                        nuc_id,
                                        rx,
                                        half_lifev,
                                        level,
                                        branch_percent,
                                        state,
                                        special,
                                    )
                                )
                    if level_found is True:
                        levellist.append(
                            (nuc_id, 0, half_lifev, level, 0.0, state, special)
                        )
                    brs = {}
                    level, half_lifev, from_nuc, state, special = _parse_level_record(
                        level_l
                    )
                    if from_nuc is not None:
                        nuc_id = from_nuc + leveln
                        leveln += 1
                        level_found = True
                    else:
                        level_found = False
                    continue
                levelc = _level_cont_regex.match(line)
                if levelc is not None:
                    brs.update(_parse_level_continuation_record(levelc))
                    continue
            if len(brs) > 0:
                for key, val in brs.items():
                    goodkey = True
                    keystrip = key.replace("%", "").lower()
                    for item in badlist:
                        if keystrip == item:
                            goodkey = False
                    if goodkey is True:
                        rx = rxname.id(keystrip)
                        branch_percent = float(val.split("(")[0])
                        levellist.append(
                            (
                                nuc_id,
                                rx,
                                half_lifev,
                                level,
                                branch_percent,
                                state,
                                special,
                            )
                        )
            if level_found is True:
                levellist.append((nuc_id, 0, half_lifev, level, 0.0, state, special))
    _adjust_ge100_branches(levellist)
    _adjust_metastables(levellist)
    _adjust_half_lives(levellist)
    return levellist


def decays(filename, decaylist=None):
    """
    This splits an ENSDF file into datasets. It then passes the dataset to the
    appropriate parser. Currently only a subset of decay datasets are
    supported. The output is a list of objects containing information
    pertaining to a particular decay.

    Parameters
    ----------
    filename : str or file
        Name of ENSDF formatted file or a file-like object containing ENSDF
        formatted data
    decaylist : list of tuples
        This is a list object which all newly processed decays will be added
        to. If it's None a new one will be created.

    Returns
    -------
    decaylist : list of tuples
        list of objects containing information pertaining to a particular
        decay. This information is in the following format:

    int
        nuc_id of the parent
    int
        nuc_id of the daughter
    int
        PyNE reaction id
    float
        half-life in seconds
    float
        half-life error in seconds
    float
        branching ratio (percent)
    float
        Conversion factor for gamma intensity to photons per 100 decays of the
        parent
    float
        Error in conversion factor for gamma intensity
    float
        Conversion factor for electron capture/beta intensity to electron
        captures/betas per 100 decays of the parent
    float
        Error in conversion factor for electron capture/beta intensity
    list
        a list containing information about each gamma ray:
            * starting level of gamma transition in stats_id form
            * final level of gamma transition in state_id form
            * original parent
            * energy in keV
            * uncertainty in energy
            * intensity (multiply by conversion factor for percentage)
            * uncertainty in intensity
            * electron conversion intensity
            * uncertainty in electron conversion intensity
            * total transition intensity
            * total transition intensity error
            * k electron conversion intensity
            * l electron conversion intensity
            * m electron conversion intensity
    list
        a list containing information about each alpha:
            * parent nuclide id in state_id form
            * child nuclide id in state_id form
            * alpha energy
            * alpha intensity in percent of total alphas
    list
        a list containing information about each beta minus from the parent
        decay:
            * parent nuclide id in state_id form
            * child nuclide id in state_id form
            * beta endpoint energy
            * beta average energy
            * beta intensity (multiply by conversion factor for percentage)
    list
        a list containing information about each beta plus and electron capture
        from the parent decay:
            * parent nuclide id in state_id form
            * child nuclide id in state_id form
            * beta plus endpoint energy
            * beta plus average energy
            * beta intensity (multiply by conversion factor for percentage)
            * electron capture intensity (multiply by conversion factor for
              percentage)
            * k electron conversion intensity
            * l electron conversion intensity
            * m electron conversion intensity
    """
    if decaylist is None:
        decaylist = []
    if isinstance(filename, str):
        with open(filename, "r") as f:
            dat = f.read()
    else:
        dat = filename.read()
    datasets = dat.split(80 * " " + "\n")
    for dataset in datasets:
        lines = dataset.splitlines()
        if len(lines) == 0:
            continue
        ident = re.match(_ident, lines[0])
        if ident is None:
            continue
        if "DECAY" in ident.group(2):
            decay_s = ident.group(2).split()[1]
            decay = _parse_decay_dataset(lines, decay_s)
            if decay is not None:
                if isinstance(decay[0], list):
                    if isinstance(decay[3], list):
                        for i, parent in enumerate(decay[0]):
                            dc = copy.deepcopy(list(decay))
                            dc[0] = parent
                            dc[3] = decay[3][i]
                            dc[4] = decay[4][i]
                            for gamma in dc[11]:
                                gamma[2] = parent
                            for alpha in dc[12]:
                                alpha[0] = parent
                            for beta in dc[13]:
                                beta[0] = parent
                            for ecbp in dc[14]:
                                ecbp[0] = parent
                            decaylist.append(tuple(dc))
                    else:
                        for parent in decay[0]:
                            dc = copy.deepcopy(list(decay))
                            dc[0] = parent
                            for gamma in dc[11]:
                                gamma[2] = parent
                            for alpha in dc[12]:
                                alpha[0] = parent
                            for beta in dc[13]:
                                beta[0] = parent
                            for ecbp in dc[14]:
                                ecbp[0] = parent
                            decaylist.append(tuple(dc))
                else:
                    decaylist.append(decay)
    return decaylist


def _dlist_gen(f):
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
        with open(f, "r") as f:
            dat = f.read()
    else:
        dat = f.read()
    decaylist = []
    datasets = dat.split(80 * " " + "\n")[0:-1]
    for dataset in datasets:
        lines = dataset.splitlines()
        ident = re.match(_ident, lines[0])
        if ident is not None:
            if "DECAY" in ident.group(2):
                fin = ident.group(2).split()[1]
                if fin not in decaylist:
                    decaylist.append(fin)

    return decaylist


def _level_dlist_gen(f, keys):
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
        with open(f, "r") as f:
            dat = f.read()
    else:
        dat = f.read()
    datasets = dat.split(80 * " " + "\n")[0:-1]
    for dataset in datasets:
        lines = dataset.splitlines()
        ident = re.match(_ident, lines[0])
        if ident is not None:
            if "ADOPTED LEVELS" in ident.group(2):
                for line in lines:
                    levelc = _level_cont_regex.match(line)
                    if levelc is None:
                        continue
                    ddict = _parse_level_continuation_record(levelc)
                    for item in ddict.keys():
                        if item in keys:
                            continue
                        keys.append(item)
    return keys
