"""
Module for parsing FISPACT output data
FISPACT and FISPACT-II are bateman equation solvers for transmutation 
and fission product yield calculations.  FISPACT-II is developed and maintained
 by the UKAEA. 
it supports not only neutron irradiaiton but gamma, proton, deuteron and triton
irradiation. it has support for self shielding and can read endf format nuclear 
data

this module has methods for parsing the fispact output file, 
extracting data, processing the data
"""

import numpy as np

class fispact_output():
    """ fispact output data"""
    def __init__(self):
        """ define data"""
        self.file_name = ""
        self.sumdat = []
        self.timestep_data = []
        self.pathway_data = []
        self.pathway_targets = []
        self.num_cool_step = 0   # number of steps after zero keyword
        self.num_irrad_step = 0  # number of steps with flux > 0
        self.version = ""
        self.isFisII = False
        self.cpu_time = 0.0
        self.tot_irrad_time = 0.0
        self.tot_fluence = 0.0
        self.ave_flux = 0.0
        self.time_days = []


class fispact_timestep():
    """ data for an individual time step can be heating or cooling """
    def __init__(self):
        self.step_num = 1
        self.step_length = 0
        self.flux_amp = 0
        self.is_cooling = False
        self.num_nuclides = 0
        self.alpha_act = 0 # bq
        self.beta_act = 0  # bq
        self.gamma_act = 0 # bq
        self.total_act = 0 # bq
        self.total_act_no_trit = 0 # bq
        self.alpha_heat = 0 # kw
        self.beta_heat = 0  # kw
        self.gamma_heat = 0 # kw
        self.total_heat = 0 # kw
        self.total_heat_no_trit = 0 # kw
        self.num_fissions = 0
        self.neutron_flux = 0 # n/cm**2/s
        self.initial_mass = 0 # kg
        self.total_mass = 0   # kg
        self.density = 0      # g/cc
        self.actinide_burn = 0 # %
        self.appm_h1 = 0 
        self.appm_h2 = 0
        self.appm_h3 = 0
        self.appm_he3 = 0
        self.appm_he4 = 0
        
        self.dom_data = []
        self.inventory = []
        self.gspec = []
        self.composition = []


def read_fis_out(path):
    """ parse a fispact output file
        returns fo, a fispact output object
    """

    fo = fispact_output()
    fo.file_name = path
    with open(path) as f:
        lines = f.read().splitlines()
    f.close()

    fo.version = check_fisp_version(lines)
    fo.isFisII = isFisII(lines)
    fo.sumdat = read_summary_data(lines)
    #fo.pathway_data = read_pathways(lines)
    #fo.pathway_targets = read_pathway_targets(lines)
    
    fo.ave_flux = read_parameter(lines, "Mean flux")
    fo.tot_irrad_time = read_parameter(lines, "Total irradiation time")
    fo.tot_fluence = read_parameter(lines, "Total fluence")
    fo.num_irrad_step = read_parameter(lines, "Number of on-times")
    
    if isFisII:
        search_string = "fispact run time"
    else:
        search_string = "CPU Time used for case"
    fo.cpu_time = read_parameter(lines, search_string)
    
    # find where each time step starts
    time_step_inds = []
    for line in lines:
        if len(line) > 0:
            if line[0:7] == "1 * * *":
                time_step_inds.append(lines.index(line))
    
    # parse all time steps
    i = 0
    while i < len(time_step_inds)-1:
        data = lines[time_step_inds[i]:time_step_inds[i+1]]
        fo.timestep_data.append(read_time_step(data, i))
        i = i + 1
    # final timestep
    data = lines[time_step_inds[-1]:]
    fo.timestep_data.append(read_time_step(data, i))
       
    return fo


def read_time_step(lines, i):
    """ reads a particular time step """
    ts = fispact_timestep()
    ts.step_num = i + 1
    
    ts.step_length = float(lines[0][50:60])

    
    ind = find_ind(lines, "TOTAL NUMBER OF NUCLIDES PRINTED IN INVENTORY")
    ts.num_nuclides = int(lines[ind][50:])

    ind = find_ind(lines, "ALPHA BECQUERELS")
    ts.alpha_act = float(lines[ind][22:34])
    ts.beta_act = float(lines[ind][54:66])
    ts.gamma_act = float(lines[ind][87:99])
    
    ind = find_ind(lines, "TOTAL ACTIVITY FOR ALL MATERIALS ")
    ts.total_act = float(lines[ind][40:51])    

    ind = find_ind(lines, "TOTAL ACTIVITY EXCLUDING TRITIUM ")
    ts.total_act_no_trit = float(lines[ind][40:51])

    ind = find_ind(lines, "TOTAL ALPHA HEAT")
    ts.alpha_heat = float(lines[ind][40:51])
    ts.beta_heat = float(lines[ind + 1][40:51])
    ts.gamma_heat = float(lines[ind + 2][40:51])
    ts.total_heat = float(lines[ind + 2][90:101])
    ts.initial_mass = float(lines[ind + 3][40:51])
    ts.total_heat_no_trit = float(lines[ind + 3][90:101])
    ts.total_mass = float(lines[ind + 4][40:51])
    ts.neutron_flux = float(lines[ind + 5][40:51])
    ts.num_fissions = float(lines[ind + 6][40:51])
    ts.actinide_burn = float(lines[ind + 6][90:101])

    ind = find_ind(lines, "DENSITY")
    ts.density = float(lines[ind][78:86])

    if ts.total_act > 0.0:
        ind = find_ind(lines, "APPM OF He  4 ")
        ts.appm_he4 = float(lines[ind][23:33])
        ts.appm_he3 = float(lines[ind+1][23:33])
        ts.appm_h3 = float(lines[ind+2][23:33])
        ts.appm_h2 = float(lines[ind+3][23:33])
        ts.appm_h1 = float(lines[ind+4][23:33])
        ind=1

        ts.dom_data = parse_dominant(lines)
        ts.composition = parse_composition(lines)
        ts.gspec = parse_spectra(lines)
    ts.inventory = parse_inventory(lines)

    return ts


def check_fisp_version(data):
    """ Checks which version of fispact was used to produced data

        Parameters
        ---------
        data: list with each element being a line from the fispact output file

        Returns
        ------
        v: a string of the version name"""

    sub = "FISPACT VERSION 07.0/0"
    data = data[:50]
    if next((s for s in data if sub in s), None):
        v = "FISP07"
    else:
        v = "FISPACT-II"
    return v


def isFisII(data):
    """boolean check if file is fispact-ii output """
    v = check_fisp_version(data)
    if v == "FISPACT-II":
        return True
    else:
        return False


def read_summary_data(data):
    """ Processes the summary block at the end of the file"""

    if isFisII(data):
        cool_str = " -----Irradiation Phase-----"
    else:
        cool_str = "  COOLING STEPS"

    start_ind = data.index(cool_str)
    end_ind = [i for i, line in enumerate(data) if "0 Mass" in line]
    sum_lines = data[start_ind+1:end_ind[0]]
    sum_data = []
    time_yrs = []
    act = []
    act_un = []
    dr = []
    dr_un = []
    heat = []
    heat_un = []
    ing = []
    ing_un = []
    inhal = []
    inhal_un = []
    trit = []
    to = 0

    for l in sum_lines:
        if isFisII(data):
            if l[1] == "-":
                to = time_yrs[-1] 
            else:
                
                time_yrs.append(float(l[24:32]) + to)
                act.append(l[35:43])
                dr.append(l[58:66])
                heat.append(l[81:89])
                ing.append(l[104:112])
                inhal.append(l[127:135])
                trit.append(l[150:158])
                # act_un.append(l[])
                # dr_un.append(l[])
                # heat_un.append(l[])
                # ing_un.append(l[])
                # inhal_un.append(l[])
        else:
            time_yrs.append(l[20:28])
            act.append(l[31:39])
            dr.append(l[54:62])
            heat.append(l[77:85])
            ing.append(l[100:108])
            inhal.append(l[123:131])
            trit.append(l[146:154])
            # act_un.append(l[])
            # dr_un.append(l[])
            # heat_un.append(l[])
            # ing_un.append(l[])
            # inhal_un.append(l[])
    sum_data.append(time_yrs)
    sum_data.append(act)
    sum_data.append(dr)
    sum_data.append(heat)
    sum_data.append(ing)
    sum_data.append(inhal)
    sum_data.append(trit)
    sum_data.append(act_un)
    sum_data.append(dr_un)
    sum_data.append(heat_un)
    sum_data.append(ing_un)
    sum_data.append(inhal_un)

    return sum_data


def parse_dominant(data):
    """parse dominant nuclides section and return a list of lists """
    p1_ind = find_ind(data, "DOMINANT NUCLIDES")
    data = data[p1_ind:]
    d1_ind = find_ind(data, "(Bq) ")
    d2_ind = find_ind(data, "GAMMA HEAT")
    topset = data[d1_ind+2:d2_ind-1]
    topset = np.array(topset)
    lowerset = data[d2_ind+3:]

    act_nuc = []
    act = []
    act_percent = []
    heat_nuc = []
    heat = []
    heat_percent = []
    dr_nuc = []
    dr = []
    dr_percent = []
    gheat_nuc = []
    gheat = []
    gheat_percent = []
    bheat_nuc = []
    bheat = []
    bheat_percent = []    

    for l in topset:
        act_nuc.append(l[7:13])
        act.append(l[15:25])
        act_percent.append(l[27:36])
        heat_nuc.append(l[38:44])
        heat.append(l[46:56])
        heat_percent.append(l[58:67])
        dr_nuc.append(l[69:75])
        dr.append(l[77:87])
        dr_percent.append(l[89:98])

    for l in lowerset:
        gheat_nuc.append(l[7:13])
        gheat.append(l[15:25])
        gheat_percent.append(l[27:36])
        bheat_nuc.append(l[38:44])
        bheat.append(l[46:56])
        bheat_percent.append(l[58:67])

    dom_data = []
    dom_data.append(act_nuc)
    dom_data.append(act)
    dom_data.append(act_percent)
    dom_data.append(heat_nuc)
    dom_data.append(heat)
    dom_data.append(heat_percent)
    dom_data.append(dr_nuc)
    dom_data.append(dr)
    dom_data.append(dr_percent)
    dom_data.append(gheat_nuc)
    dom_data.append(gheat)
    dom_data.append(gheat_percent)
    dom_data.append(bheat_nuc)
    dom_data.append(bheat)
    dom_data.append(bheat_percent)

    return dom_data


def parse_composition(data):
    """ parse compostions section """
    p1 = find_ind(data, "COMPOSITION  OF  MATERIAL  BY  ELEMENT")
    p2 = find_ind(data, "GAMMA SPECTRUM AND ENERGIES/SECOND")
    data = data[p1+5:p2-3]
    ele_list = []
    atoms = []

    for l in data:
        ele_list.append(l[12:14])
        atoms.append(float(l[20:30]))

    composition = []
    composition.append(ele_list)
    composition.append(atoms)

    return composition 


def parse_spectra(data):
    """ reads gamma spectra data for each timestep
        returns list of length 24 corresponding to 24 gamma energy groups
        data is in gamma/s/cc
    """
    p1 = find_ind(data, "GAMMA SPECTRUM AND ENERGIES/SECOND")
    data = data[p1+7:p1+31]
    spectra = []
    for l in data:
        spectra.append(float(l[130:141]))
    return spectra


def parse_inventory(data):
    """ """
    inv = []
    return inv


def find_ind(data, sub):
    """ """
    for i, s in enumerate(data):
        if sub in s:
            ind = i
    return ind


def read_parameter(data, sub):
    """ """
    ind=find_ind(data, sub)
    line = data[ind]
    line = line.split("=")
    line = line[1].strip()
    line=line.split(" ")
    param=float(line[0])
    return param


