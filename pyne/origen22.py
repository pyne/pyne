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

    lower_lines = ["1 {0} {1:.10E}   0 0   0 0   0 0".format(nuc, mass) for nuc, mass in lower_z.mult_by_mass()]
    upper_lines = ["2 {0} {1:.10E}   0 0   0 0   0 0".format(nuc, mass) for nuc, mass in upper_z.mult_by_mass()]
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


tape5_irradiation_template = ("  -1\n"
                              "  -1\n"
                              "  -1\n"
                              "  CUT     5 {CUT_OFF} -1\n"
                              "  RDA     FIND CROSS SECTION LIBRARY IDENTIFIER NUMBERS IN YOUR LIBRARY FILE\n"
                              "  LIB     0 1 2 3 {NLB1} {NLB2} {NLB3} 9 3 0 3 0\n"
                              "  OPTL    {optl}\n"
                              "  OPTA    {opta}\n"
                              "  OPTF    {optf}\n"
                              "  INP     1 -1  0  -1  4  4\n"
                              "  HED     1     IN FUEL\n"
                              "  RDA     ALL IRRADIATION (IRF and IRP) CARDS MUST TAKE PLACE IN BETWEEN BURNUP (BUP) CARDS\n"
                              "  BUP\n"
                              "  {ir_type}     {ir_time}  {ir_value}   1   2   4  2\n"
                              "  BUP\n"
                              "  OUT     2  1 1 0\n"
                              "  END\n")

nes_table = np.zeros((2, 2, 2), dtype=int)
nes_table[True, True, True]    = 1
nes_table[True, True, False]   = 2
nes_table[True, False, True]   = 3
nes_table[False, True, True]   = 4
nes_table[True, False, False]  = 5
nes_table[False, True, False]  = 6
nes_table[False, False, True]  = 7
nes_table[False, False, False] = 8

def out_table_string(out_table_nes, out_table_num):
    """Makes a string output table line from relevant information."""
    if out_table_num == None: 
        arr = np.ones(24, dtype=int)
        s = np.array2string(arr)[1:-1]
        return s

    arr = 8 * np.ones(24, dtype=int)
    idx = np.array(out_table_num, dtype=int) - 1

    arr[idx] = nes_table[tuple(out_table_nes)]

    s = np.array2string(arr)[1:-1]

    return s

def write_tape5_irradiation(ir_type, ir_time, ir_value, nlb, 
                            cut_off=10.0**-10, name="TAPE5.INP",
                            out_table_nes=(False, False, True),
                            out_table_laf=(True,  True,  True),
                            out_table_num=None):
    """Writes a TAPE5 files to `name`.

    Args:
        * ir_type (str): Flag that determines whether this is a constant power "IRP"
          irradiation or a constant flux "IRF" irradiation calculation.
        * ir_time (float): Irradiation time durration in days.
        * ir_value (float): Magnitude of the irradiation. If ir_type = "IRP", then
          this is a power.  If ir_type = "IRF", then this is a flux. 
        * nlb: Three tuple of library numbers from the tape9 file, eg (204, 205, 206).

    Keyword Args:
        * cut_off (float): Cut-off concentration, below which reults are not recorded.
        * out_table_nes (three-tuple of bools): Specifies which type of output
          tables should be printed by ORIGEN.  The fields represent 
          (Nuclide, Element, Summary).  Thus, the default value of (False, False, True) 
          only prints the summary tables. 
        * out_table_laf (three-tuple of bools): Specifies whether to print the 
          activation products (l), actinides (a), and fission products (f).
          By default all three are printed.
        * out_table_num (sequence of ints or None):  Specifies which tables, by number,
          to print according to the rules given by out_table_nes and 
          out_table_laf.  For example the list [10, 5] would print tables 5 
          and 10.  There are 24 tables available.
          If None, then all tables are printed.   
    """

    if ir_type not in ["IRP", "IRF"]:
        raise TypeError("Irradiation type must be either 'IRP' or 'IRF'.")
    
    # Make template fill-value dictionary
    tape5_kw = {
        'CUT_OFF':  "{0:.3E}".format(cut_off),
        'NLB1':     nlb[0],
        'NLB2':     nlb[1],
        'NLB3':     nlb[2],
        'ir_type':  ir_type,
        'ir_time':  '{0:.10E}'.format(ir_time),
        'ir_value': '{0:.10E}'.format(ir_value),
        }

    # Activation Product Print String
    if out_table_laf[0]:
        tape5_kw['optl'] = out_table_string(out_table_nes, out_table_num)
    else:
        tape5_kw['optl'] = np.array2string(8 * np.ones(24, dtype=int))[1:-1]

    # Actinide Print String
    if out_table_laf[1]:
        tape5_kw['opta'] = out_table_string(out_table_nes, out_table_num)
    else:
        tape5_kw['opta'] = np.array2string(8 * np.ones(24, dtype=int))[1:-1]

    # Fission Product Print String
    if out_table_laf[2]:
        tape5_kw['optf'] = out_table_string(out_table_nes, out_table_num)
    else:
        tape5_kw['optf'] = np.array2string(8 * np.ones(24, dtype=int))[1:-1]

    # Fill the template and write it to a file
    with open(name, 'w') as f:
        f.write(tape5_irradiation_template.format(**tape5_kw))

    return

def run_origen(self):
    """Runs the ORIGEN Burnup Calculations."""
    os.chdir('libs/ORIGEN/')

    # Grab General Data from the HDF5 File
    libfile = tb.openFile("../{0}.h5".format(reactor), 'r')
    CoreLoadIsos = list(libfile.root.CoreLoad_zzaaam)
    libfile.close()

    if 0 < verbosity:
        print(message("Preping the ORIGEN Directories..."))
    t1 = time.time()
    for t in FineTime[1:]:
        self.make_input_tape5(t)
        self.make_input_tape9(t)
    for iso in CoreLoadIsos: 
        os.mkdir("{0}".format(iso))
    t2 = time.time()
    if 0 < verbosity:
        print(message("...Done!  That only took {0:time} min.\n", "{0:.3G}".format((t2-t1)/60.0) ))

    if 0 < verbosity:
        print(message("  ~~~~~  Starting ORIGEN Runs  ~~~~~  "))
    orit1 = time.time()

    # Initialize the data structures
    self.BU  = {}
    self.k   = {}
    self.Pro = {}
    self.Des = {}
    self.Tij = {}

    for iso in CoreLoadIsos:
        isoLL = isoname.zzaaam_2_LLAAAM(iso)
        if 0 < verbosity:
            print(message("  ~~~~~  Now on {0:iso}  ~~~~~  \n", "Isotope {0}".format(isoLL)))
        isot1 = time.time()

        # Initilize iso data, for t = 0
        self.BU[iso]  = [0.0]
        self.k[iso]   = [0.0]
        self.Pro[iso] = [0.0]
        self.Des[iso] = [0.0]
        self.Tij[iso] = [{iso: 1000.0}]

        for t in FineTime[1:]:
            if 0 < verbosity:
                print(message("Starting ORIGEN run for {0:iso} at {1:time}...", isoLL, "Time {0}".format(t)))
            t1 = time.time()

            os.chdir("{0}".format(iso))

            # Make/Get Input Decks
            self.make_input_tape4(Tij[iso][-1])
            shutil.copy("../{0}_T{1}.tape5".format(reactor, t), "TAPE5.INP")
            shutil.copy("../{0}_T{1}.tape9".format(reactor, t), "TAPE9.INP")

            # Run ORIGEN
            subprocess.call(self.run_str, shell=True)

            # Parse Output
            self.parse_tape6()
            self.BU[iso].append(  BU[iso][-1] + self.tape6_BU )
            self.k[iso].append(   self.tape6_k )
            self.Pro[iso].append( self.tape6_Pro )
            self.Des[iso].append( self.tape6_Des )
            self.Tij[iso].append( self.tape6_outvec )

            # Clean up the directory
            for f in os.listdir('.'):
                if f[-4:] in ['.INP', '.OUT']:
                    metasci.SafeRemove(f)
            os.chdir('../') # Back to ORIGEN Directory

            t2 = time.time()
            if 0 < verbosity:
                print(message("ORIGEN run completed in {0:time} min!", "{0:.3G} min".format((t2-t1)/60.0) ))
    
        isot2 = time.time()
        if 0 < verbosity:
            print(message("  ~~~~~  Isotope {0:iso} took {1:time} min!  ~~~~~  \n", isoLL, "{0:.3G} min".format((isot2-isot1)/60.0) ))


    # Kludge to put Tij in the right units and form
    allORIGENisoList = []
    for iso in CoreLoadIsos:
        for t in Tij[iso]:
            for j in t.keys():
                if (j not in allORIGENisoList):
                    allORIGENisoList.append(j)
    for iso in CoreLoadIsos:
        for n_t in range(len(Tij[iso])):
            for j in allORIGENisoList:
                if j in Tij[iso][n_t].keys():
                    Tij[iso][n_t][j] = Tij[iso][n_t][j] / (10.0**3)
                else:
                    Tij[iso][n_t][j] = 0.0
    
    orit2 = time.time()
    if 0 < verbosity:
        print(message("  ~~~~~  ORIGEN took {0:time} to run!  ~~~~~  ", "{0:.3G} min".format((orit2-orit1)/60.0) ))

    os.chdir('../../') #Back to 'reactor' root
    return 

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


def write_tape9(name_new="TAPE9.INP", name_org="original.tape9", SNG={}, SN2N={}, 
    SN3N={}, SNF={}, SNA={}, SNP={}, SNGX={}, SN2NX={}, iso_list=[], sigfig=3):
    """Writes a new ORIGEN TAPE9.INP file based on an original TAPE9 template 
    file and given cross section values.

    Keyword Args:
    	* `name_new` (str): Path to the new tape9 file.
    	* `name_org` (str): Path to the tape9 template file.
        * `SNG` (dict): (n, gamma) cross-section isotopic vector, of form {iso: xs_(n,g) [barns]}.
        * `SN2N` (dict): (n, 2n) cross-section isotopic vector, of form {iso: xs_(n,2n) [barns]}.
        * `SN3N` (dict): (n, 3n) cross-section isotopic vector, of form {iso: xs_(n,3n) [barns]}.
        * `SNF` (dict): (n, fission) cross-section isotopic vector, of form {iso: xs_(n,f) [barns]}.
        * `SNA` (dict): (n, alpha) cross-section isotopic vector, of form {iso: xs_(n,alpha) [barns]}.
        * `SNP` (dict): (n, proton) cross-section isotopic vector, of form {iso: xs_(n,p) [barns]}.
        * `SNGX` (dict): (n, gamma*) excited state cross-section isotopic vector, of form {iso: xs_(n,g*) [barns]}.
        * `SN2NX` (dict): (n, 2n*) excited state cross-section isotopic vector, of form {iso: xs_(n,2n*) [barns]}.
        * `iso_list` (list): List of isotopes to write over.
        * `sigfig` (int): Ensures that all overwritten data is given to this many digits beyond the 
          decimal point via "{xs:.{sigfig}E}".format().
    """

    #perform some basic set-up
    zero_space = "0.0       "

    SNG   = isoname.isovec_keys_2_zzaaam(SNG)
    SN2N  = isoname.isovec_keys_2_zzaaam(SN2N)
    SN3N  = isoname.isovec_keys_2_zzaaam(SN3N)
    SNF   = isoname.isovec_keys_2_zzaaam(SNF)
    SNA   = isoname.isovec_keys_2_zzaaam(SNA)
    SNP   = isoname.isovec_keys_2_zzaaam(SNP)
    SNGX  = isoname.isovec_keys_2_zzaaam(SNGX)
    SN2NX = isoname.isovec_keys_2_zzaaam(SN2NX)

    if iso_list == []:
        iso_set = set()
        iso_set.update(SNG.keys())
        iso_set.update(SN2N.keys())
        iso_set.update(SN3N.keys())
        iso_set.update(SNF.keys())
        iso_set.update(SNA.keys())
        iso_set.update(SNP.keys())
        iso_set.update(SNGX.keys())
        iso_set.update(SN2NX.keys())
        iso_list = list(iso_set)
    else:
        iso_list = isoname.mixed_2_zzaaam_List(iso_list)

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
    return
