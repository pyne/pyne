import re

import numpy as np

from pyne import nucname

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
  RDA     FIND CROSS SECTION LIBRARY IDENTIFIER NUMBERS IN YOUR LIBRARY FILE
  LIB     0 1 2 3 {NLB1} {NLB2} {NLB3} 9 3 0 3 0
  OPTL    {optl}
  OPTA    {opta}
  OPTF    {optf}
  INP     1 -1  0  -1  4  4
  HED     1     IN FUEL
  RDA     ALL IRRADIATION (IRF and IRP) CARDS MUST TAKE PLACE IN BETWEEN BURNUP (BUP) CARDS
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


def write_tape5_irradiation(irr_type, irr_time, irr_value, nlb, 
                            cut_off=1E-10, outfile="TAPE5.INP",
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
    nlb : length 3 sequence
        Three tuple of library numbers from the tape9 file, eg (204, 205, 206).
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
        'NLB1': nlb[0],
        'NLB2': nlb[1],
        'NLB3': nlb[2],
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



def write_tape9(outfile="TAPE9.INP", origfile="tape9.original", sigma_gamma=None, sigma_2n=None, 
                sigma_3n=None, sigma_f=None, sigma_alpha=None, sigma_p=None, sigma_gamma_x=None, 
                sigma_2n_x=None, nuclist=None, precision=3):
    """Writes an ORIGEN 2.2 TAPE9.INP file given cross section values, and optionally based on an 
    original TAPE9 template file.

    Parameters
    ----------
    outfile : str or file-like object, optional
        Path to the new tape9 file.
    origfile : str or file-like object, optional
        Path to the tape9 template file.
    sigma_gamma : dict, optional
        Microscopic (n, gamma) cross-section nuclide map of form {nuc: xs_(n,g) [barns]}.
    sigma_2n : dict, optional 
        Microscopic (n, 2n) cross-section nuclide map of form {iso: xs_(n,2n) [barns]}.
    sigma_3n : dict, optional
        Microscopic (n, 3n) cross-section nuclide map of form {iso: xs_(n,3n) [barns]}.
    sigma_f : dict, optional
        Microscopic (n, fission) cross-section nuclide map of form {iso: xs_(n,f) [barns]}.
    sigma_alpha : dict, optional
        Microscopic (n, alpha) cross-section nuclide map of form {iso: xs_(n,alpha) [barns]}.
    sigma_p : dict, optional 
        Microscopic (n, proton) cross-section nuclide map of form {iso: xs_(n,p) [barns]}.
    sigma_gamma_x : dict, optional
        Microscopic (n, gamma*) excited state cross-section nuclide map of form {iso: xs_(n,g*) [barns]}.
    sigma_2n_x : dict, optional
        Microscopic (n, 2n*) excited state cross-section nuclide map of form {iso: xs_(n,2n*) [barns]}.
    nuclist : sequence, optional
        List of nuclides to write over.
    precision :  int, optional 
        The number of significant figures that all output data is given to beyond the decimal point.

    Notes
    -----
    All nuclides must be given in zzaaam form.
    """

    #perform some basic set-up
    zero_space = "0.0       "

    if nuclist is None:
        nucset = set()
        nucset.update(sigma_gamma.keys())
        nucset.update(sigma_2n.keys())
        nucset.update(sigma_3n.keys())
        nucset.update(sigma_f.keys())
        nucset.update(sigma_alpha.keys())
        nucset.update(sigma_p.keys())
        nucset.update(sigma_gamma_x.keys())
        nucset.update(sigma_2n_x.keys())
    else:
        nucset = set(nuclist)

    #Open the first tape9 for reading and the new tape9 for writing.
    tape9_o = open(name_org, 'r')   #Original File
    tape9_n = open(name_new, 'w')   #New File

    #Write the new file...
    for line in tape9_o:
        ls = line.split()
        if ls == []:
            #Rewrites blank lines
            tape9_n.write(line)
            continue
        elif int(ls[0]) <= 3:
            #Rewrites Decay data and '-1' spacer lines
            tape9_n.write(line)
            continue

        #Rewrites Title Lines
        try:
            int(ls[1])
            iso = isoname.mixed_2_zzaaam(ls[1])
        except ValueError:
            tape9_n.write(line)
            continue

        #If we don't care about the isotope, just rewrite the line...
        if iso not in iso_list:
            tape9_n.write(line)
            continue

        #Fixes messed up exponential data in the original file...
        orig_data = []
        SkipNext = False
        for n in range(len(ls)):
            if SkipNext:
                SkipNext = False
            elif 'E' in ls[n]:
                try:
                    orig_data.append(float(ls[n]))
                except ValueError:
                    orig_data.append(float(ls[n] + ls[n+1]))
                    SkipNext = True
            else:
                orig_data.append(float(ls[n]))
        ls = orig_data	#This is what ls was suppossed to look like! Stupid ORIGEN...			 

        newline = line[:13]

        #(n, gamma) XS
        if iso in SNG.keys(): 
            if SNG[iso] == 0.0:
                newline = newline + zero_space
            else:
                newline = newline + "{0:.{1}E} ".format(SNG[iso], sigfig)
        else:
            newline = newline + "{0:.{1}E} ".format(ls[2], sigfig)

        #(n, 2n) XS
        if iso in SN2N.keys():
            if SN2N[iso] == 0.0:
                newline = newline + zero_space
            else:
                newline = newline + "{0:.{1}E} ".format(SN2N[iso], sigfig)
        else:
            newline = newline + "{0:.{1}E} ".format(ls[3], sigfig)

        #Check if in actinide set, because positional arguments mean different things.
        if (iso/10000) in isoname.act:
            #(n, 3n) XS
            if iso in SN3N.keys():
                if SN3N[iso] == 0.0:
                    newline = newline + zero_space
                else:
                    newline = newline + "{0:.{1}E} ".format(SN3N[iso], sigfig)
            else:
                newline = newline + "{0:.{1}E} ".format(ls[4], sigfig)

            #(n, fission) XS
            if iso in SNF.keys():
                if SNF[iso] == 0.0:
                    newline = newline + zero_space
                else:
                    newline = newline + "{0:.{1}E} ".format(SNF[iso], sigfig)
            else:
                newline = newline + "{0:.{1}E} ".format(ls[5], sigfig)

        #If not an actinide, then...
        else:
            #(n, alpha) XS
            if iso in SNA.keys():
                if SNA[iso] == 0.0:
                    newline = newline + zero_space
                else:
                    newline = newline + "{0:.{1}E} ".format(SNA[iso], sigfig)
            else:
                newline = newline + "{0:.{1}E} ".format(ls[4], sigfig)

            #(n, proton) XS
            if iso in SNP.keys():
                if SNP[iso] == 0.0:
                    newline = newline + zero_space
                else:
                    newline = newline + "{0:.{1}E} ".format(SNP[iso], sigfig)
            else:
                newline = newline + "{0:.{1}E} ".format(ls[5], sigfig)

        #(n, g*) XS
        if iso in SNGX.keys(): 
            if SNGX[iso] == 0.0:
                newline = newline + zero_space
            else:
                newline = newline + "{0:.{1}E} ".format(SNGX[iso], sigfig)
        elif ls[2] == 0.0:
            newline = newline + zero_space
        elif iso in SNG.keys():
            sngx = SNG[iso] * ls[6] / ls[2]
            if sngx == 0.0:
                newline = newline + zero_space
            else:
                newline = newline + "{0:.{1}E} ".format(sngx, sigfig)
        else:
            newline = newline + "{0:.{1}E} ".format(ls[6], sigfig)

        #(n, 2n*) XS
        if iso in SN2NX.keys(): 
            if SN2NX[iso] == 0.0:
                newline = newline + zero_space
            else:
                newline = newline + "{0:.{1}E} ".format(SN2NX[iso], sigfig)
        elif ls[3] == 0.0:
            newline = newline + zero_space
        elif iso in SN2N.keys():
            sn2nx = SN2N[iso] * ls[7] / ls[3]
            if sn2nx == 0.0:
                newline = newline + zero_space
            else:
                newline = newline + "{0:.{1}E} ".format(sn2nx, sigfig)
        else:
            newline = newline + "{0:.{1}E} ".format(ls[7], sigfig)

        #End the line
        newline = newline + line[-8:]
        tape9_n.write(newline)

    #Close out the files and return
    tape9_o.close()
    tape9_n.close()





def parse_tape6(name = "TAPE6.OUT"):
    """Parses an ORIGEN TAPE6.OUT file that is in the current directory + path p.

    Keyword Args:
        * name (str): Path to the tape6 file to parse:

    Returns:
        * results (dict): Dictionary of parsed values.
    """
    results = {}

    # Defaults
    in_table = False
    table_key = None
    table_type = None

    # Read the TAPE6 file
    with open(name, 'r') as tape6:
        for line in tape6:
            # Skip trivial lines
            if len(line) == 0:
                continue

            # Spliut the line
            ls = line.split()

            # Grab Basis lines
            if "TIME, SEC" in line:
                results["time_sec"] = float(ls[-1])

            elif "NEUT. FLUX" in line:
                results["flux"] = float(ls[-1])

            elif "SP POW,MW" in line:
                results["specific_power_MW"] = float(ls[-1])

            elif "BURNUP,MWD" in line:
                results["burnup_MWD"] = float(ls[-1])

            elif "K INFINITY" in line:
                results["k_inf"] = float(ls[-1])

            elif "NEUT PRODN" in line:
                results["neutron_production_rate"] = float(ls[-1])

            elif "NEUT DESTN" in line:
                results["neutron_destruction_rate"] = float(ls[-1])

            elif "TOT BURNUP" in line:
                results["total_burnup"] = float(ls[-1])

            elif "AVG N FLUX" in line:
                results["average_flux"] = float(ls[-1])

            elif "AVG SP POW" in line:
                results["average_specific_power"] = float(ls[-1])

            elif ("TABLE:" in line):
                in_table = True

                # Set table key
                if line[0] == "0":
                    table_key = "table_{0}".format(ls[1])
                else:
                    table_key = "table_{0}".format(ls[0])

                if table_key not in results.keys():
                    results[table_key] = {}

                # Set table type
                if line[0] == "0":
                    table_type = ls[2].lower()
                else:
                    table_type = ls[1].lower()

                if table_type not in results[table_key].keys():
                    results[table_key][table_type] = {}

                    pline = line.partition(":")[2].partition(",")
                    title = pline[0].strip()
                    units = pline[2].strip()

                    results[table_key][table_type]["title"] = title
                    results[table_key][table_type]["units"] = units
                    results[table_key][table_type]["data"]  = {}
                    

            elif in_table and ("OUTPUT UNIT = " in line):
                # restore defaults
                in_table = False
                table_key = None
                table_type = None

            elif in_table:
                ind = 0
                try:
                    iso = isoname.LLAAAM_2_zzaaam(ls[0])
                    ind = 1
                except:
                    try:
                        iso = isoname.LLAAAM_2_zzaaam(ls[0] + ls[1])
                        ind = 2
                    except:
                        continue

                #results[table_key][table_type]["data"][iso] = float(ls[-1])
                if iso not in results[table_key][table_type]["data"]:
                    results[table_key][table_type]["data"][iso] = []
                results[table_key][table_type]["data"][iso].append(np.array(ls[ind:], dtype=float))
            else:
                continue

    return results


data_format = "\d+\.\d*[EeDd]?[+-]?\d+"
title_card_re = re.compile("(\d+)\s+(\S.*)")

# Decay library regex
decay_card1_re = re.compile("(\d+)\s+(\d{{5,7}})\s+(\d)\s+({num})\s+({num})\s+({num})\s+({num})\s+({num})\s+({num})".format(num=data_format))
decay_card2_re = re.compile("(\d+)\s+({num})\s+({num})\s+({num})\s+({num})\s+({num})\s+({num})".format(num=data_format))

# Cross section and fission product yeild library regex
xsfpy_card1_re = re.compile("(\d+)\s+(\d{{5,7}})\s+({num})\s+({num})\s+({num})\s+({num})\s+({num})\s+({num})\s+([+-]?{num})".format(num=data_format))
xsfpy_card2_re = re.compile("(\d+)\s+({num})\s+({num})\s+({num})\s+({num})\s+({num})\s+({num})\s+({num})\s+({num})".format(num=data_format))


def _parse_tape9_decay(deck):
    pdeck = {'_type': 'decay'}
    pdeck['title'] = title_card_re.match(deck[0]).group(2).strip()
    return pdeck


def _parse_tape9_xsfpy(deck):
    pdeck = {'_type': 'xsfpy'}
    pdeck['title'] = title_card_re.match(deck[0]).group(2).strip()
    return pdeck


def parse_tape9(tape9="TAPE9.INP"):
    """Parses an ORIGEN 2.2 TAPE9 file and returns the data as a dictionary of nuclide dictionaries.
    
    Parameters
    ----------
    tape9 : str or file-like object, optional
        Path to the tape9 file.
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
