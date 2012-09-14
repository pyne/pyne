import re
from copy import deepcopy
from itertools import chain, imap, izip

import numpy as np

from pyne import nucname
from pyne.material import Material, from_atom_frac

# Table 4.2 in ORIGEN 2.2 manual
ORIGEN_TIME_UNITS = [None,              # No zero unit
                     1.0,               # seconds
                     60.0,              # minutes
                     3600.0,            # hours
                     86400.0,           # days
                     31556926.0,        # years...which are fuzzily defined.
                     np.inf,            # stable
                     31556926.0 * 1E3,  # ky
                     31556926.0 * 1E6,  # My
                     31556926.0 * 1E9,  # Gy
                    ]


def sec_to_time_unit(s):
    """Converts seconds to ORIGEN time and units.

    Parameters
    ----------
    s : float
        time in seconds

    Returns
    -------
    t : float
        time in units
    unit : int
        time unit that t is in. Represents index 
        into ORIGEN_TIME_UNITS, which matches 
        Table 4.2 in ORIGEN 2.2 manual.
    """
    for i, val in enumerate(ORIGEN_TIME_UNITS):
        if val is None:
            continue

        t = s / val
        unit = i 

        if t != 0.0 and val == np.inf:
            # Origen spec for stable nuclides
            t = 0.0
            break
        elif 0.0 < t < 1.0:
            if i == 1: 
                pass
            elif i == 7:
                unit -= 2 
            else:
                unit -= 1
            t = s / ORIGEN_TIME_UNITS[unit]
            break

    return t, unit
        

# Regex helpers
_data_format = "\d+\.\d*[EeDd]?[ +-]?\d+"



###################################
### ORIGEN Input Deck Functions ###
###################################

def write_tape4(mat, outfile="TAPE4.INP"):
    """Writes a TAPE4.INP ORIGEN input file for a material.

    Parameters
    ----------
    mat : Material 
        Material with mass weights in units of grams.
    outfile : str or file handler, optional 
        Path to tape4 file or file-like object.
    """
    lower_z = mat[:'AC']
    upper_z = mat['AC':]

    lower_lines = ["1 {0} {1:.10E}   0 0   0 0   0 0".format(nuc, mass) for nuc, mass in lower_z.mult_by_mass().items()]
    upper_lines = ["2 {0} {1:.10E}   0 0   0 0   0 0".format(nuc, mass) for nuc, mass in upper_z.mult_by_mass().items()]
    lines = lower_lines + upper_lines + ["0 0 0 0\n"]

    tape4 = "\n".join(lines)

    # Write to the file
    opened_here = False
    if isinstance(outfile, basestring):
        outfile = open(outfile, 'w')
        opened_here = True

    outfile.write(tape4)

    if opened_here:
        outfile.close()



_tape5_irradiation_template = """\
  -1
  -1
  -1
  CUT     5 {CUT_OFF} -1
  RDA     Make sure thet the library identifier numbers match those in the TAPE9.INP file
  LIB     0 {DECAY_NLB1} {DECAY_NLB2} {DECAY_NLB3} {XSFPY_NLB1} {XSFPY_NLB2} {XSFPY_NLB3} 9 3 0 4 0
  OPTL    {optl}
  OPTA    {opta}
  OPTF    {optf}
  INP     1 -1  0  -1  4  4
  RDA     All irradiation (IRF and IRP) cards must be between burnup (BUP) cards.
  BUP
  {irr_type}     {irr_time}  {irr_value}   1   2   4  2
  BUP
  OUT     2  1 1 0
  END
"""

_nes_table = np.zeros((2, 2, 2), dtype=int)
_nes_table[True, True, True]    = 1
_nes_table[True, True, False]   = 2
_nes_table[True, False, True]   = 3
_nes_table[False, True, True]   = 4
_nes_table[True, False, False]  = 5
_nes_table[False, True, False]  = 6
_nes_table[False, False, True]  = 7
_nes_table[False, False, False] = 8


def _out_table_string(out_table_nes, out_table_num):
    """Makes a string output table line from relevant information."""
    if out_table_num == None: 
        arr = np.ones(24, dtype=int)
        s = np.array2string(arr)[1:-1]
        return s

    arr = 8 * np.ones(24, dtype=int)
    idx = np.array(out_table_num, dtype=int) - 1

    arr[idx] = _nes_table[tuple(out_table_nes)]

    s = np.array2string(arr)[1:-1]

    return s


def write_tape5_irradiation(irr_type, irr_time, irr_value, 
                            outfile="TAPE5.INP",
                            decay_nlb=(1, 2, 3), 
                            xsfpy_nlb=(204, 205, 206), 
                            cut_off=1E-10, 
                            out_table_nes=(False, False, True),
                            out_table_laf=(True,  True,  True),
                            out_table_num=None):
    """Writes an irradiation TAPE5 file.

    Parameters
    ----------
    irr_type : str
        Flag that determines whether this is a constant power "IRP"
        irradiation or a constant flux "IRF" irradiation calculation.
    irr_time : float 
        Irradiation time durration in days.
    irr_value : float 
        Magnitude of the irradiation. If irr_type = "IRP", then
        this is a power.  If irr_type = "IRF", then this is a flux. 
    outfile : str or file-like object
        Path or file to write the tape5 to.
    decay_nlb : length 3 sequence
        Three tuple of library numbers from the tape9 file decay data, eg (1, 2, 3).
    xsfpy_nlb : length 3 sequence
        Three tuple of library numbers from the tape9 file for cross section and fission
        product yields, eg (204, 205, 206).
    cut_off : float, optional
        Cut-off concentration, below which reults are not recorded.
    out_table_nes :  length 3 sequence of bools, optional
        Specifies which type of output tables should be printed by ORIGEN.  The fields 
        represent (Nuclide, Element, Summary).  The default value of (False, False, True) 
        only prints the summary tables. 
    out_table_laf :  length 3 sequence of bools, optional 
        Specifies whether to print the activation products (l), actinides (a), and 
        fission products (f).  By default all three are printed.
    out_table_num : sequence of ints or None
        Specifies which tables, by number, to print according to the rules given by 
        out_table_nes and out_table_laf.  For example the list [10, 5] would print 
        tables 5 and 10.  There are 24 tables available. If None, then all tables 
        are printed.   
    """
    if irr_type not in ["IRP", "IRF"]:
        raise TypeError("Irradiation type must be either 'IRP' or 'IRF'.")
    
    # Make template fill-value dictionary
    tape5_kw = {
        'CUT_OFF': "{0:.3E}".format(cut_off),
        'DECAY_NLB1': decay_nlb[0],
        'DECAY_NLB2': decay_nlb[1],
        'DECAY_NLB3': decay_nlb[2],
        'XSFPY_NLB1': xsfpy_nlb[0],
        'XSFPY_NLB2': xsfpy_nlb[1],
        'XSFPY_NLB3': xsfpy_nlb[2],
        'irr_type': irr_type,
        'irr_time': '{0:.10E}'.format(irr_time),
        'irr_value': '{0:.10E}'.format(irr_value),
        }

    no_print_string = np.array2string(8 * np.ones(24, dtype=int))[1:-1]

    # Activation Product Print String
    if out_table_laf[0]:
        tape5_kw['optl'] = _out_table_string(out_table_nes, out_table_num)
    else:
        tape5_kw['optl'] = no_print_string

    # Actinide Print String
    if out_table_laf[1]:
        tape5_kw['opta'] = _out_table_string(out_table_nes, out_table_num)
    else:
        tape5_kw['opta'] = no_print_string

    # Fission Product Print String
    if out_table_laf[2]:
        tape5_kw['optf'] = _out_table_string(out_table_nes, out_table_num)
    else:
        tape5_kw['optf'] = no_print_string


    # Fill the template and write it to a file
    tape5 = _tape5_irradiation_template.format(**tape5_kw)

    opened_here = False
    if isinstance(outfile, basestring):
        outfile = open(outfile, 'w')
        opened_here = True

    outfile.write(tape5)

    if opened_here:
        outfile.close()




#
# Tape6 functions
#

_rx_bu_data_line = re.compile(' (TIME, SEC|NEUT. FLUX|SP POW,MW|BURNUP,MWD|K INFINITY|NEUT PRODN|NEUT DESTN|TOT BURNUP|AVG N FLUX|AVG SP POW) (.*)')

_rx_bu_key_map = {
    "TIME, SEC":  "time_sec",
    "NEUT. FLUX": "flux", 
    "SP POW,MW":  "specific_power_MW", 
    "BURNUP,MWD": "burnup_MWD", 
    "K INFINITY": "k_inf", 
    "NEUT PRODN": "neutron_production_rate",
    "NEUT DESTN": "neutron_destruction_rate",
    "TOT BURNUP": "total_burnup", 
    "AVG N FLUX": "average_flux",
    "AVG SP POW": "average_specific_power",
    }

_table_header_line = re.compile("[ 0]\s+(\d+) (NUCLIDE|ELEMENT|SUMMARY) TABLE:([ A-Z]+),([ 0-9A-Za-z*]+)")

_table_header_alpha_line = re.compile("[ 0]\s+(\d+) (NUCLIDE|ELEMENT|SUMMARY) TABLE:\s*(ALPHA RADIOACTIVITY)\s+([ 0-9A-Za-z*]+)")

_nuclide_line = re.compile(" ([ A-Z][A-Z][ \d][ \d]\d[ M])\s+(.*)")

_element_line = re.compile(" ([ A-Z][A-Z])   \s+(.*)")

_species_group_line = re.compile('[ +]\s+(ACTIVATION PRODUCTS|ACTINIDES[+]DAUGHTERS|FISSION PRODUCTS)')

_group_key_map = {
    'ACTIVATION PRODUCTS': 'activation_products',
    'ACTINIDES+DAUGHTERS': 'actinides',
    'FISSION PRODUCTS': 'fission_products',
    }

_alpha_n_header_line = re.compile('\s*(\(ALPHA,N\) NEUTRON SOURCE), (NEUTRONS/SEC)')

_spont_fiss_header_line = re.compile('\s*(SPONTANEOUS FISSION NEUTRON SOURCE), (NEUTRONS/SEC)')

_n_source_key_map ={
    '(ALPHA,N) NEUTRON SOURCE': 'alpha_neutron_source',
    'SPONTANEOUS FISSION NEUTRON SOURCE': 'spont_fiss_neutron_source',
    }

_photon_spec_header_line = re.compile("\s+PHOTON SPECTRUM FOR(.*)")

def parse_tape6(tape6="TAPE6.OUT"):
    """Parses an ORIGEN 2.2 TAPE6.OUT file. 

    Parameters
    ----------
    tape6 : str or file-like object
        Path or file to read the tape6 file from.

    Returns
    -------
    results : dict 
        Dictionary of parsed values.

    Warnings
    --------
    This method currently only functions to extract neutronic data from TAPE6
    files.  It does not yet parse out photonic data.  If you would like to see
    this feature added, please contact the developers.

    Notes
    -----
    The results dictionary that is returned is highly structured and generally
    matches the layout of the TAPE6 file.  Data is stored as 1d numpy float arrays
    which (if the TAPE6 is well-formed) will all be of the same length and match
    the time vector.  The possible layout of results is as follows::

      |- 'time_sec': time per index in [seconds]
      |- 'flux': neutron flux at this time [n/cm^2/s]
      |- 'specific_power_MW': recator specific power at this time [MW]
      |- 'burnup_MWD': reactor burnup since last time step [MWd/input mass [g] from TAPE4]
      |- 'k_inf': infinite multiplication factor [unitless]
      |- 'neutron_production_rate': Total reactor neutron production rate [n/s]
      |- 'neutron_destruction_rate: Total reactor neutron destruction rate [n/s]
      |- 'total_burnup': Cummulative burnup over all time [MWd/input mass [g] from TAPE4]
      |- 'average_flux': average neutron flux over preceeding time interval [n/cm^2/s]
      |- 'average_specific_power: recator specific power over preceeding time interval [MW]
      |- 'materials': list of Materials of same length as 'time_sec', only present if 
      |               'table_3' or 'table_5' exist and have 'nuclide' output.
      |- 'alpha_neutron_source': dict
      |                          |- 'title': str
      |                          |- 'units': str
      |                          |- nuclide or element str: (alpha, n) neutron source [n/s] 
      |- 'spont_fiss_neutron_source': dict
      |                          |- 'title': str
      |                          |- 'units': str
      |                          |- nuclide or element str: spontaneous fission neutron source [n/s] 
      |- 'table_{n}': dict
      |               |- 'nuclide': dict
      |               |             |- 'title': str
      |               |             |- 'units': str
      |               |             |- 'activation_products': dict of (nuc-zzaaam, data) pairs
      |               |             |- 'actinides': dict of (nuc-zzaaam, data) pairs
      |               |             |- 'fission_products': dict of (nuc-zzaaam, data) pairs
      |               |- 'element': dict
      |               |             |- 'title': str
      |               |             |- 'units': str
      |               |             |- 'activation_products': dict of (elem str, data) pairs
      |               |             |- 'actinides': dict of (elem str, data) pairs
      |               |             |- 'fission_products': dict of (elem str, data) pairs
      |               |- 'summary': dict
      |               |             |- 'title': str
      |               |             |- 'units': str
      |               |             |- 'activation_products': dict of (elem or nuc str, data) pairs
      |               |             |- 'actinides': dict of (elem or nuc str, data) pairs
      |               |             |- 'fission_products': dict of (elem or nuc str, data) pairs
 
    """
    # Read the TAPE6 file
    opened_here = False
    if isinstance(tape6, basestring):
        tape6 = open(tape6, 'r')
        opened_here = True

    lines = tape6.readlines()

    if opened_here:
        tape6.close()

    # Prep to parse the file
    results = {}

    # Defaults
    table_key = None
    table_type = None
    table_group = None

    # Read in the file line-by-line
    for i, line in enumerate(lines):
        # Get reactivity and burnup data
        m = _rx_bu_data_line.match(line)
        if m is not None:
            key, data = m.groups()
            new_key = _rx_bu_key_map[key]
            arr_data = np.array(data.split(), dtype=float)
            curr_data = results.get(new_key, [])
            results[new_key] = np.append(curr_data, arr_data)
            continue

        # Get table spcies group
        m = _species_group_line.match(line)
        if m is not None:
            table_group = _group_key_map[m.group(1)]
            continue

        # Get table header info
        m = _table_header_line.match(line) or _table_header_alpha_line.match(line)
        if m is not None:
            tnum, ttype, ttitle, tunits = m.groups()

            table_key = "table_{0}".format(tnum)
            if table_key not in results.keys():
                results[table_key] = {}

            table_type = ttype.lower()
            if table_type not in results[table_key]:
                results[table_key][table_type] = {}

            results[table_key][table_type]["title"] = ttitle.strip().lower()
            results[table_key][table_type]["units"] = tunits.strip().lower()
            if table_group not in results[table_key][table_type]:
                results[table_key][table_type][table_group] = {}
            continue


        # Grab nuclide data lines
        m = _nuclide_line.match(line)
        if (m is not None) and (table_key is not None):
            nuc, data = m.groups()
            nuc_name = nuc.replace(' ', '')

            # Don't know WTF element 'SF' is suppossed to be! (Spent fuel, spontaneous fission)
            if nuc_name == 'SF250':
                continue

            nuc_zz = nucname.zzaaam(nuc_name)
            nuc_key = nuc_zz if table_type == 'nuclide' else nuc_name
            nuc_data = np.array(data.split(), dtype=float)

            if table_key.startswith('table_'): 
                curr_data = results[table_key][table_type][table_group].get(nuc_key, [])
                results[table_key][table_type][table_group][nuc_key] = np.append(curr_data, nuc_data)
            else:
                curr_data = results[table_key].get(nuc_key, [])
                results[table_key][nuc_key] = np.append(curr_data, nuc_data)
            continue

        # Grab element data line
        m = _element_line.match(line)
        if (m is not None) and (table_key is not None):
            elem, data = m.groups()
            elem = elem.replace(' ', '')

            # Still don't know WTF element 'SF' is suppossed to be! (Spent fuel, spontaneous fission)
            if elem == 'SF':
                continue

            elem_data = np.array(data.split(), dtype=float)

            if table_key.startswith('table_'): 
                curr_data = results[table_key][table_type][table_group].get(elem, [])
                results[table_key][table_type][table_group][elem] = np.append(curr_data, elem_data)
            else:
                curr_data = results[table_key].get(elem, [])
                results[table_key][elem] = np.append(curr_data, elem_data)
            continue

        # Grab (alpha, n) and spontaneous fission headers
        m = _alpha_n_header_line.match(line) or _spont_fiss_header_line.match(line)
        if m is not None:
            ttitle, tunits = m.groups()

            table_key = _n_source_key_map[ttitle]
            if table_key not in results.keys():
                results[table_key] = {}

            table_type = None
            table_group = None

            results[table_key]["title"] = ttitle.strip().lower()
            results[table_key]["units"] = tunits.strip().lower()
            continue

        # Photon spectra parsing is not yet supported
        m = _photon_spec_header_line.match(line)
        if m is not None:
            table_key = None
            table_type = None
            table_group = None

    # Done with parsing, try to convert to material
    tbl = None
    if ('table_5' in results) and ('nuclide' in results['table_5']):
        tbl = 'table_5'
        mat_gen = Material
    elif ('table_3' in results) and ('nuclide' in results['table_3']):
        tbl = 'table_3'
        mat_gen = from_atom_frac

    if tbl is not None:
        T = len(results['time_sec'])
        mats = [Material() for t in range(T)]

        for grp in _group_key_map.values():
            if grp in results[tbl]['nuclide']:
                mats = [m + mat_gen(dict([(nuc, arr[i]) for nuc, arr in \
                                            results[tbl]['nuclide'][grp].items()])) \
                                            for i, m in enumerate(mats)]

        results['materials'] = mats

    return results


#
# Tape9 functions
#

title_card_re = re.compile("(\d+)\s+(\S.*)")

# Decay library regex
decay_card1_re = re.compile("(\d+)\s+(\d{{5,7}})\s+(\d)\s+({num})\s+({num})\s+({num})\s+({num})\s+({num})\s+({num})".format(num=_data_format))
decay_card2_re = re.compile("(\d+)\s+({num})\s+({num})\s+({num})\s+({num})\s+({num})\s+({num})".format(num=_data_format))

# Cross section and fission product yeild library regex
xsfpy_card1_re = re.compile("(\d+)\s+(\d{{5,7}})\s+({num})\s+({num})\s+({num})\s+({num})\s+({num})\s+({num})\s+([+-]?{num})".format(num=_data_format))
xsfpy_card2_re = re.compile("(\d+)\s+({num})\s+({num})\s+({num})\s+({num})\s+({num})\s+({num})\s+({num})\s+({num})".format(num=_data_format))


def _parse_tape9_decay(deck):
    pdeck = {'_type': 'decay'}
    pdeck['title'] = title_card_re.match(deck[0]).group(2).strip()

    # Parse the cards into a structured arrau
    cards = [m.groups()[1:] + n.groups()[1:] for m, n in 
             izip(imap(decay_card1_re.match, deck[1::2]), imap(decay_card2_re.match, deck[2::2]))]
    cards = [tuple(d.replace(' ', '') for d in card) for card in cards]
    cards = np.array(cards, dtype='i4,i4' + ',f8'*12)
    pdeck['_cards'] = cards

    # Add the first cards
    pdeck['half_life'] = dict([(nuc, ORIGEN_TIME_UNITS[unit]*(val or 1.0)) for nuc, unit, val in cards[['f0', 'f1', 'f2']] ])
    pdeck['frac_beta_minus_x'] = dict([(nuc, val) for nuc, val in cards[['f0', 'f3']] ])
    pdeck['frac_beta_plus_or_electron_capture'] = dict([(nuc, val) for nuc, val in cards[['f0', 'f4']] ])
    pdeck['frac_beta_plus_or_electron_capture_x'] = dict([(nuc, val) for nuc, val in cards[['f0', 'f5']] ])
    pdeck['frac_alpha'] = dict([(nuc, val) for nuc, val in cards[['f0', 'f6']] ])
    pdeck['frac_isomeric_transition'] = dict([(nuc, val) for nuc, val in cards[['f0', 'f7']] ])

    # Add the second cards
    pdeck['frac_spont_fiss'] = dict([(nuc, val) for nuc, val in cards[['f0', 'f8']] ])
    pdeck['frac_beta_n'] = dict([(nuc, val) for nuc, val in cards[['f0', 'f9']] ])
    pdeck['recoverable_energy'] = dict([(nuc, val) for nuc, val in cards[['f0', 'f10']] ])
    pdeck['frac_natural_abund'] = dict([(nuc, val*0.01) for nuc, val in cards[['f0', 'f11']] ])
    pdeck['inhilation_concentration'] = dict([(nuc, val) for nuc, val in cards[['f0', 'f12']] ])
    pdeck['ingestion_concentration'] = dict([(nuc, val) for nuc, val in cards[['f0', 'f13']] ])

    return pdeck


def _parse_tape9_xsfpy(deck):
    pdeck = {'_type': 'xsfpy'}
    pdeck['title'] = title_card_re.match(deck[0]).group(2).strip()

    # Pasre the deck
    cards = []
    no_fpy = ('0.0', ) * 8

    i = 1
    deck_size = len(deck)
    while (i < deck_size):
        first_card = xsfpy_card1_re.match(deck[i]).groups()[1:]
        if 0.0 < float(first_card[-1]):
            i += 1
            second_card = xsfpy_card2_re.match(deck[i]).groups()[1:]
        else:
            second_card = no_fpy
        cards.append(first_card + second_card)
        i += 1

    cards = [tuple(d.replace(' ', '') for d in card) for card in cards]
    cards = np.array(cards, dtype='i4' + ',f8'*15)
    pdeck['_cards'] = cards

    # Try to determine subtype
    if (0.0 < cards['f7']).any():
        subtype = 'fission_products'
    elif 890000 < cards['f0'].mean():
        subtype = 'actinides'
    else:
        subtype = 'activation_products'
    pdeck['_subtype'] = subtype


    # Parse first cards
    pdeck['sigma_gamma'] = dict([(nuc, val) for nuc, val in cards[['f0', 'f1']] ])
    pdeck['sigma_2n'] = dict([(nuc, val) for nuc, val in cards[['f0', 'f2']] ])

    f3_keys = {'fission_products': 'sigma_alpha', 'actinides': 'sigma_3n', 'activation_products': 'sigma_alpha'}
    pdeck[f3_keys[subtype]] = dict([(nuc, val) for nuc, val in cards[['f0', 'f3']] ])

    f4_keys = {'fission_products': 'sigma_p', 'actinides': 'sigma_f', 'activation_products': 'sigma_p'}
    pdeck[f4_keys[subtype]] = dict([(nuc, val) for nuc, val in cards[['f0', 'f4']] ])

    pdeck['sigma_gamma_x'] = dict([(nuc, val) for nuc, val in cards[['f0', 'f5']] ])
    pdeck['sigma_2n_x'] = dict([(nuc, val) for nuc, val in cards[['f0', 'f6']] ])

    pdeck['fiss_yields_present'] = dict([(nuc, 0.0 < val) for nuc, val in cards[['f0', 'f7']] ])

    # parse second cards if of correct subtype 
    if subtype == 'fission_products':
        pdeck['TH232_fiss_yield'] = dict([(nuc, val) for nuc, val in cards[['f0', 'f8']] ])
        pdeck['U233_fiss_yield'] = dict([(nuc, val) for nuc, val in cards[['f0', 'f9']] ])
        pdeck['U235_fiss_yield'] = dict([(nuc, val) for nuc, val in cards[['f0', 'f10']] ])
        pdeck['U238_fiss_yield'] = dict([(nuc, val) for nuc, val in cards[['f0', 'f11']] ])
        pdeck['PU239_fiss_yield'] = dict([(nuc, val) for nuc, val in cards[['f0', 'f12']] ])
        pdeck['PU241_fiss_yield'] = dict([(nuc, val) for nuc, val in cards[['f0', 'f13']] ])
        pdeck['CM245_fiss_yield'] = dict([(nuc, val) for nuc, val in cards[['f0', 'f14']] ])
        pdeck['CF249_fiss_yield'] = dict([(nuc, val) for nuc, val in cards[['f0', 'f15']] ])

    return pdeck


def parse_tape9(tape9="TAPE9.INP"):
    """Parses an ORIGEN 2.2 TAPE9 file and returns the data as a dictionary of nuclide dictionaries.
    
    Parameters
    ----------
    tape9 : str or file-like object, optional
        Path to the tape9 file.

    Returns
    -------
    parsed : dict
        A dictionary of the data from the TAPE9 file.

    Notes
    -----
    The TAPE9 format is highly structured. Therefore the in-memory representation contains a 
    non-trivial amount of nesting.  At the top level, the dictionary keys are the library 
    numbers::

        tape9
          |- keys : deck number (1, 2, 3, 241, ...) 
          |- values : sub-dictionaries for each deck.

    Each deck contains keys which vary by deck type (and subtype).  All dictionary-typed data
    maps zzaaam-nuclide integers to the appropriate value::

        all decks
          |- '_type' : str in ['decay', 'xsfpy'] # photon libs not yet supported
          |- '_subtype' : str for 'xsfpy' in ['activation_products', 'actinides', 'fission_products']
          |- 'title' : str, deck name
          |- '_cards' : optional, numpy structrued array of deck data

        decay decks
          |- 'half_life' : float-valued dict [seconds]
          |- 'frac_beta_minus_x' : float-valued dict [fraction of decays via beta minus 
          |                        which leave an excited nucleus]
          |- 'frac_beta_plus_or_electron_capture' : float-valued dict [fraction of decays
          |                                         via positron emission or electron capture]
          |- 'frac_beta_plus_or_electron_capture_x' : float-valued dict [fraction of decays
          |                                           via positron emission or electron capture
          |                                           which leave an excited nucleus]
          |- 'frac_alpha' : float-valued dict [fraction of decays via alpha emission]
          |- 'frac_isomeric_transition' : float-valued dict [fraction of decays from an excitied 
          |                               state to the ground state]
          |- 'frac_spont_fiss' : float-valued dict [fraction of decays via spontateous fission]
          |- 'frac_beta_n' : float-valued dict [fraction of decays via beta plus a neutron]
          |- 'recoverable_energy' : float-valued dict, Total recoverable energy [MeV / decay]
          |- 'frac_natural_abund' : float-valued dict, natrual occuring abundance [atom fraction]
          |- 'inhilation_concentration' : float-valued dict, continuous inhilation [RCG]
          |- 'ingestion_concentration' : float-valued dict, continuous ingestion [RCG]

        cross section and fission product yield decks
          |- 'sigma_gamma' : float-valued dict, (n, gamma) cross section [barns]
          |- 'sigma_2n' : float-valued dict, (n, 2n) cross section [barns]
          |- 'sigma_gamma_x' : float-valued dict, (n, gamma *) cross section [barns]
          |- 'sigma_2n_x' : float-valued dict, (n, 2n *) cross section [barns]
          |- 'fiss_yields_present' : bool-valued dict, Whether fission product yields are 
                                     included for this nuclide.

        activation product cross section decks
          |- '_subtype' : 'activation_products'
          |- 'sigma_3n' : float-valued dict, (n, 3n) cross section [barns]
          |- 'sigma_p' : float-valued dict, (n, proton) cross section [barns]

        actinide cross section decks
          |- '_subtype' : 'actinides'
          |- 'sigma_alpha' : float-valued dict, (n, alpha) cross section [barns]
          |- 'sigma_f' : float-valued dict, (n, fission) cross section [barns]

        fission product cross section and yield decks
          |- '_subtype' : 'fission_products'
          |- 'sigma_3n' : float-valued dict, (n, 3n) cross section [barns]
          |- 'sigma_p' : float-valued dict, (n, proton) cross section [barns]
          |- 'TH232_fiss_yield' : float-valued dict, yield from Th-232 fission [frac]
          |- 'U233_fiss_yield' : float-valued dict, yield from U-233 fission [frac]
          |- 'U235_fiss_yield' : float-valued dict, yield from U-235 fission [frac]
          |- 'U238_fiss_yield' : float-valued dict, yield from U-238 fission [frac]
          |- 'PU239_fiss_yield' : float-valued dict, yield from Pu-239 fission [frac]
          |- 'PU241_fiss_yield' : float-valued dict, yield from Pu-241 fission [frac]
          |- 'CM245_fiss_yield' : float-valued dict, yield from Cm-245 fission [frac]
          |- 'CF249_fiss_yield' : float-valued dict, yield from Cf-249 fission [frac]
    """
    # Read and strip lines
    opened_here = False
    if isinstance(tape9, basestring):
        tape9 = open(tape9, 'r')
        opened_here = True

    tape9_lines = [line.strip() for line in tape9]

    if opened_here:
        tape9.close()

    # Split lines into various decks.
    decks = []
    while 0 < tape9_lines.count('-1'):
        n = tape9_lines.index('-1')
        decks.append(tape9_lines[:n])
        tape9_lines = tape9_lines[n+1:]

    # parse the individual decks.
    parsed = {}
    for deck in decks:
        # Decay deck
        m = decay_card1_re.match(deck[1])
        if m is not None:
            deck_num = int(m.group(1))
            parsed[deck_num] = _parse_tape9_decay(deck)
            continue

        # Cross section deck
        m = xsfpy_card1_re.match(deck[1])
        if m is not None:
            deck_num = int(m.group(1))
            parsed[deck_num] = _parse_tape9_xsfpy(deck)
            continue

    return parsed



def merge_tape9(tape9s):
    """Merges a sequence of full or partial TAPE9s into a single tape9 dictionary.
    Data from the first tape9 has precednce over the second, the second over the
    third, etc.

    Parameters
    ----------
    tape9s : list or tuple
        A sequence of full or partial tape9 dictionaries.  See parse_tape9()
        for more information on the structure of these dictionaires.

    Returns
    -------
    tape9 : dictionary
        A tape9 file which is the merger of the all of the tape9s.  See 
        parse_tape9() for more information on the structure of this dictionary.
    """
    tape9 = {}

    for t9 in tape9s[::-1]:
        for nlb, deck in t9.items():
            if nlb in tape9:
                # Make sure the decks are of the same type
                assert tape9[nlb]['_type'] == deck['_type']
                if ('_subtype' in tape9[nlb]) and ('_subtype' in deck):
                    assert tape9[nlb]['_subtype'] == deck['_subtype']

                # _cards has been invalidated... remove it
                tape9[nlb].pop('_cards', None)

                # Update all of the keys, except _cards
                for key, value in deck.items():
                    if key in tape9[nlb] and hasattr(value, 'keys'):
                        tape9[nlb][key].update(value)
                    elif key != '_cards':
                        tape9[nlb][key] = deepcopy(value)
            else:
                # New library number, make a copy
                tape9[nlb] = deepcopy(deck)

    return tape9



def _double_get(dict, key1, key2, default=0.0):
    if key1 in dict:
        return dict[key1].get(key2, default)
    else:
        return default


_deck_title_fmt = "{nlb:>4}    {title:^72}\n"

_decay_card_fmt = ("{nlb:>4}{nuc:>8}  {unit}     {time:<9.{p}E} {fbx:<9.{p}E} {fpec:<9.{p}E} {fpecx:<9.{p}E} {fa:<9.{p}E} {fit:<9.{p}E}\n"
                   "{nlb:>4}                {fsf:<9.{p}E} {fn:<9.{p}E} {qrec:<9.{p}E} {abund:<9.{p}E} {arcg:<9.{p}E} {wrcg:<9.{p}E}\n")

_xs_card_fmt = "{nlb:>4}{nuc:>8} {sg:<9.{p}E} {s2n:<9.{p}E} {s3n_or_a:<9.{p}E} {sf_or_p:<9.{p}E} {sg_x:<9.{p}E} {s2n_x:<9.{p}E} {fpy_flag:>6.1F} \n"

_fpy_card_fmt = "{nlb:>4}     {y1:<9.{p}E} {y2:<9.{p}E} {y3:<9.{p}E} {y4:<9.{p}E} {y5:<9.{p}E} {y6:<9.{p}E} {y7:<9.{p}E} {y8:<9.{p}E}\n"

def _decay_deck_2_str(nlb, deck, precision):
    # Get unique isotopes 
    nucset = set([nuc for nuc in chain(*[v.keys() for k, v in deck.items() if hasattr(v, 'keys')]) ])

    s = ""
    for nuc in nucset:
        t, unit = sec_to_time_unit(_double_get(deck, 'half_life', nuc))
        s += _decay_card_fmt.format(nlb=nlb,
                                    nuc=nuc,
                                    unit=unit,
                                    time=t,
                                    fbx=_double_get(deck, 'frac_beta_minus_x', nuc),
                                    fpec=_double_get(deck, 'frac_beta_plus_or_electron_capture', nuc),
                                    fpecx=_double_get(deck, 'frac_beta_plus_or_electron_capture_x', nuc),
                                    fa=_double_get(deck, 'frac_alpha', nuc),
                                    fit=_double_get(deck, 'frac_isomeric_transition', nuc),
                                    fsf=_double_get(deck, 'frac_spont_fiss', nuc),
                                    fn=_double_get(deck, 'frac_beta_n', nuc),
                                    qrec=_double_get(deck, 'recoverable_energy', nuc),
                                    abund=_double_get(deck, 'frac_natural_abund', nuc),
                                    arcg=_double_get(deck, 'inhilation_concentration', nuc),
                                    wrcg=_double_get(deck, 'ingestion_concentration', nuc),
                                    p=precision,
                                    )
    return s


def _xs_deck_2_str(nlb, deck, precision):
    # Get unique isotopes 
    nucset = set([nuc for nuc in chain(*[v.keys() for k, v in deck.items() if hasattr(v, 'keys')]) ])

    is_actinides = deck['_subtype'] == 'actinides'

    s = ""
    for nuc in nucset:
        fpy_flag = -1.0
        fpy_present = _double_get(deck, 'fiss_yields_present', nuc, False)
        if fpy_present:
            fpy_flag = 1.0

        s += _xs_card_fmt.format(nlb=nlb,
                                 nuc=nuc,
                                 sg=_double_get(deck, 'sigma_gamma', nuc),
                                 s2n=_double_get(deck, 'sigma_2n', nuc),
                                 s3n_or_a=_double_get(deck, 'sigma_alpha', nuc) if is_actinides else _double_get(deck, 'sigma_3n', nuc),
                                 sf_or_p=_double_get(deck, 'sigma_f', nuc) if is_actinides else _double_get(deck, 'sigma_p', nuc),
                                 sg_x=_double_get(deck, 'sigma_gamma_x', nuc),
                                 s2n_x=_double_get(deck, 'sigma_2n_x', nuc),
                                 fpy_flag=fpy_flag,
                                 p=precision,
                                 )
    return s


def _xsfpy_deck_2_str(nlb, deck, precision):
    # Get unique isotopes 
    nucset = set([nuc for nuc in chain(*[v.keys() for k, v in deck.items() if hasattr(v, 'keys')]) ])

    is_actinides = deck['_subtype'] == 'actinides'

    s = ""
    for nuc in nucset:
        fpy_flag = -1.0
        fpy_present = _double_get(deck, 'fiss_yields_present', nuc, False)
        if fpy_present:
            fpy_flag = 1.0

        s += _xs_card_fmt.format(nlb=nlb,
                                 nuc=nuc,
                                 sg=_double_get(deck, 'sigma_gamma', nuc),
                                 s2n=_double_get(deck, 'sigma_2n', nuc),
                                 s3n_or_a=_double_get(deck, 'sigma_alpha', nuc) if is_actinides else _double_get(deck, 'sigma_3n', nuc),
                                 sf_or_p=_double_get(deck, 'sigma_f', nuc) if is_actinides else _double_get(deck, 'sigma_p', nuc),
                                 sg_x=_double_get(deck, 'sigma_gamma_x', nuc),
                                 s2n_x=_double_get(deck, 'sigma_2n_x', nuc),
                                 fpy_flag=fpy_flag,
                                 p=precision,
                                 )
        if fpy_present:
            s += _fpy_card_fmt.format(nlb=nlb,
                                 y1=_double_get(deck, 'TH232_fiss_yield', nuc),
                                 y2=_double_get(deck, 'U233_fiss_yield', nuc),
                                 y3=_double_get(deck, 'U235_fiss_yield', nuc),
                                 y4=_double_get(deck, 'U238_fiss_yield', nuc),
                                 y5=_double_get(deck, 'PU239_fiss_yield', nuc),
                                 y6=_double_get(deck, 'PU241_fiss_yield', nuc),
                                 y7=_double_get(deck, 'CM245_fiss_yield', nuc),
                                 y8=_double_get(deck, 'CF249_fiss_yield', nuc),
                                 p=precision,
                                 )
    return s


_DECK_2_STR_MAP = {
    ('decay', None): _decay_deck_2_str,
    ('xsfpy', 'activation_products'): _xs_deck_2_str,
    ('xsfpy', 'actinides'): _xs_deck_2_str,
    ('xsfpy', 'fission_products'): _xsfpy_deck_2_str,
    }


def write_tape9(tape9, outfile="TAPE9.INP", precision=3):
    """Writes an ORIGEN 2.2 TAPE9.INP file given a tape9 dictionary of values.

    Parameters
    ----------
    tape9 : dict
        A tape9 dictionary. See parse_tape9() for more information on the structure.
    outfile : str or file-like object, optional
        Path to the new tape9 file.
    precision :  int, optional 
        The number of significant figures that all output data is given to beyond 
        the decimal point.
    """
    t9 = ""

    for nlb, deck in tape9.items():
        t9 += _deck_title_fmt.format(nlb=nlb, title=deck['title'])
        t9 += _DECK_2_STR_MAP[deck['_type'], deck.get('_subtype', None)](nlb, deck, precision)
        t9 += "  -1\n"

    opened_here = False
    if isinstance(outfile, basestring):
        outfile = open(outfile, 'w')
        opened_here = True

    outfile.write(t9)

    if opened_here:
        outfile.close()
