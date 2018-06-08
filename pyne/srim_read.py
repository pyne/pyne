"""This module is for SRIM files. Currently defines a
SRIM class, reads a srim file.
"""

from warnings import warn
from pyne.utils import QAWarning
import numpy as np

warn(__name__ + " is not yet QA compliant.", QAWarning)


class srim(object):
    """srim class includes srim specific variables"""

    def __init__(self, ion='', ion_mass=0, target_density=0,
                 target_comp=None, date='', bragg_corr=0,
                 filename='',ions=None,stop_pow_ele=None,stop_pow_nuc=None,
                 ion_units=None,proj_ranges=None,proj_range_units=None,
                 long_straggs=None,long_stragg_units=None,lat_straggs=None,
                 lat_stragg_units=None):
        """Define srim specific variables """
        super(srim, self).__init__()
        self.ion = ion
        self.ion_mass = ion_mass
        self.target_density = target_density
        self.target_comp = [] if target_comp is None else target_comp
        self.date = date
        self.bragg_corr = bragg_corr
        self.filename = filename
        self.ions = ions
        self.stop_pow_ele = stop_pow_ele
        self.stop_pow_nuc = stop_pow_nuc
        self.ion_units = ion_units
        self.proj_ranges = proj_ranges
        self.proj_range_units = proj_range_units
        self.long_straggs = long_straggs
        self.long_stragg_units = long_stragg_units
        self.lat_straggs = lat_straggs
        self.lat_stragg_units = lat_stragg_units

    def __str__(self):
        """Print debug information"""
        print_string = ('Debug print of all header variables\n'
                        'The ion is: {x.ion}\n'
                        'The ion mass is: {x.ion_mass}\n'
                        'The target density is: {x.target_density}\n'
                        'Target compostion is {x.target_comp}\n'
                        'Start date: {x.date}\n'
                        'Bragg Correction: {x.bragg_corr}\n'
                        'Stopping power electric: {x.stop_pow_ele}\n'
                        'Stopping power nuclear: {x.stop_pow_nuc}\n'
                        'Projected Range: {x.proj_ranges}\n'
                        'File name: {x.file_name}'
                        'Ions Energy: {x.ions}\n'
                        'Projected ranges: {x.proj_ranges}\n'
                        'Projected ranges units: {x.proj_range_units}\n'
                        'Longitudinal Straggling: {x.long_straggs}\n'
                        'Longitudinal Straggling units: {x.long_stragg_units}\n'
                        'Lateral Straggling: {x.lat_straggs}\n'
                        'Lateral Straggling units: {x.lat_stragg_units}\n').format(x=self)
        return print_string

def read_srim_output(filepath):
    #def read_srim_output(file_path):
    """Reads a srim file 
    
    Parameters
    ----------
    filepath : str
        Path to srim file

    Returns
    -------
    info: a srim object
        Contains all information from srim file

    """

    with open(filepath, "r") as srim_file:
        full_file_text = srim_file.read()
    file_split = full_file_text.splitlines()
    srim_file.close()
    
    info = srim()
    
    for item in file_split:
        if "Disk" in item:
            info.filename = item.split('= ')[1]
        elif "Bragg Correction" in item:
            info.bragg_corr = float((item.split('= ')[1]).strip('%'))
        elif "Ion =" in item:
            info.ion = (item.split('= ')[1]).split( )[0]
            info.ion_mass = float((item.split('= ')[2]).strip(' amu'))
        elif "Stopping Units =" in item:
            info.stopping_units = item.split('= ')[1]
        elif "Calc. date" in item:
            info.date = item.split('---> ')[1].strip( )
        elif "Target Density" in item:
            info.target_density_grams = item.split('=')[1].strip( )
            info.target_density_atoms = item.split('=')[2].strip( )
    
    data_start = file_split.index("   Ion        dE/dx      dE/dx     Projected  Longitudinal   Lateral") +3
    data_end = file_split.index("-----------------------------------------------------------") -1
    ions = []
    stop_pow_ele = []
    stop_pow_nuc = []
    proj_ranges = []
    ion_units = []
    proj_range_units = []
    long_straggs = []
    long_stragg_units = []
    lat_straggs = []
    lat_stragg_units = []
    for i in range(data_start,data_end+1):
        data = file_split[i].split( )
        ion = data[0]
        ion_unit = data[1]
        ele_dedx = data[2] 
        nuc_dedx = data[3]
        proj_range = data[4]
        proj_range_unit = data[5]
        long_stragg = data[6]
        long_stragg_unit = data[7]
        lat_stragg = data[8]
        lat_stragg_unit = data[9]
        ions.append(float(ion))
        ion_units.append(ion_unit)
        stop_pow_ele.append(float(ele_dedx))
        stop_pow_nuc.append(float(nuc_dedx))
        proj_ranges.append(float(proj_range))
        proj_range_units.append(proj_range_unit)
        long_straggs.append(float(long_stragg))
        long_stragg_units.append(long_stragg_unit)
        lat_straggs.append(float(lat_stragg))
        lat_stragg_units.append(lat_stragg_unit)
        
    for j in range(len(ions)):
        if "MeV" in ion_units[j]:
            ions[j] = ions[j] *1000 #If MeV, convert to keV
            ion_units[j] = 'keV'
    for j in range(len(proj_range)):
        if "um" in proj_range_units[j]:
            proj_range[j] = proj_range[j] *10000 #If um, convert to A
            proj_range_units[j] = 'A'
    for j in range(len(long_straggs)):
        if "um" in long_stragg_units[j]:
            long_straggs[j] = long_straggs[j] *10000 #If um, convert to A
            long_stragg_units[j] = 'A'
    for j in range(len(lat_straggs)):
        if "um" in lat_stragg_units[j]:
            lat_straggs[j] = lat_straggs[j] *10000 #If um, convert to A
            lat_stragg_units[j] = 'A'
    info.ions = ions
    info.stop_pow_ele = stop_pow_ele
    info.stop_pow_nuc = stop_pow_nuc
    info.ion_units = ion_units
    info.proj_ranges = proj_ranges
    info.proj_range_units = proj_range_units
    info.long_stragg_units = long_stragg_units
    info.long_straggs = long_straggs
    info.long_stragg_units = long_stragg_units
    info.lat_straggs = lat_straggs
    info.lat_stragg_units = lat_stragg_units
    return info

