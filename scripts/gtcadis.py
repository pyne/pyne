#!/usr/bin/env python
import argparse
import ConfigParser
from os.path import isfile

import numpy as np
from pyne.mesh import Mesh
from pyne.partisn import write_partisn_input, isotropic_vol_source
from pyne.dagmc import discretize_geom, load
from pyne import nucname
from pyne.bins import pointwise_collapse


config_filename = 'config.ini'

config = \
"""
# Optional step to assess all materials in geometry for compatibility with 
# SNILB criteria
[step0]

# Prepare PARTISN input file for adjoint photon transport
[step1]
# Path to hdf5 geometry file for photon transport
geom_file:
# Volume ID of adjoint photon source cell on 
# DAGMC input [Trelis/Cubit .sat file]
src_cell:
# Volume [cm^3] of source cell (detector)
src_vol:
# Define mesh that covers entire geometry:
# xmin, xmax, number of intervals
xmesh:
# ymin, ymax, number of intervals
ymesh:
# zmin, zmax, number of intervals
zmesh:

# Calculate T matrix for each material
[step2]

# Calculate adjoint neutron source
[step3]

# Prepare PARTISN input for adjoint neutron transport
[step4]

# Generate Monte Carlo variance reduction parameters 
# (biased source and weight windows)
[step5]


"""



def setup():
    with open(config_filename, 'w') as f:
        f.write(config)
    print('File "{}" has been written'.format(config_filename))
    print('Fill out the fields in these files then run ">> gtcadis.py step1"')

def _names_dict():
    names = ['h1', 'd', 'h3', 'he3', 'he4', 'li6', 'li7', 'be9', 'b10', 'b11',
    'c12', 'n14', 'n15', 'o16', 'f19', 'na23', 'mgnat', 'al27', 'si28', 'si29',
    'si30', 'p31', 'snat', 'cl35', 'cl37', 'knat', 'canat', 'ti46', 'ti47', 'ti48',
    'ti49', 'ti50', 'vnat', 'cr50', 'cr52', 'cr53', 'cr54', 'mn55', 'fe54', 'fe56',
    'fe57', 'fe58', 'co59', 'ni58', 'ni60', 'ni61', 'ni62', 'ni64', 'cu63', 'cu65',
    'ganat', 'zrnat', 'nb93', 'mo92', 'mo94', 'mo95', 'mo96', 'mo97', 'mo98',
    'mo100', 'snnat', 'ta181', 'w182', 'w183', 'w184', 'w186', 'au197', 'pb206',
    'pb207', 'pb208', 'bi209']

    names_formatted = ['h1', 'h2', 'h3', 'he3', 'he4', 'li6', 'li7', 'be9', 'b10', 'b11',
    'c12', 'n14', 'n15', 'o16', 'f19', 'na23', 'mg', 'al27', 'si28', 'si29',
    'si30', 'p31', 's', 'cl35', 'cl37', 'k', 'ca', 'ti46', 'ti47', 'ti48',
    'ti49', 'ti50', 'v', 'cr50', 'cr52', 'cr53', 'cr54', 'mn55', 'fe54', 'fe56',
    'fe57', 'fe58', 'co59', 'ni58', 'ni60', 'ni61', 'ni62', 'ni64', 'cu63', 'cu65',
    'ga', 'zr', 'nb93', 'mo92', 'mo94', 'mo95', 'mo96', 'mo97', 'mo98',
    'mo100', 'sn', 'ta181', 'w182', 'w183', 'w184', 'w186', 'au197', 'pb206',
    'pb207', 'pb208', 'bi209']

    names_dict = {nucname.id(x):y for x, y in zip(names_formatted, names)}

    return names_dict
 
def _cards():
    cards = {"block1": {"isn": 16,
                        "maxscm": '3E8',
                        "maxlcm": '6E8',
                       },
             "block3": {"lib": "xsf21-71", # name of cross section library
                       "lng":175,
                       "maxord": 5,
                       "ihm": 227,
                       "iht": 10,
                       "ihs": 11,
                       "ifido": 1,
                       "ititl": 1,
                       "i2lp1": 0,
                       "savbxs": 0,
                       "kwikrd": 0
                       },
            "block5": {"source": source,
                       "ith":1,
                       "isct":5}
            }
    return cards
   
def step1():
    # Parse config file
    config = ConfigParser.ConfigParser()
    config.read(config_filename)

    # Get user-input from config file
    geom = config.get('step1', 'geom_file') 
    cells = [config.getint('step1', 'src_cell')]
    src_vol = [config.getfloat('step1', 'src_vol')]
    xmesh = config.get('step1', 'xmesh').split(',')
    ymesh = config.get('step1', 'ymesh').split(',')
    zmesh = config.get('step1', 'zmesh').split(',')
    
    # Create structured mesh
    sc = [np.linspace(float(xmesh[0]), float(xmesh[1]), float(xmesh[2])+1),
          np.linspace(float(ymesh[0]), float(ymesh[1]), float(ymesh[2])+1),
          np.linspace(float(zmesh[0]), float(zmesh[1]), float(zmesh[2])+1)]
    mesh = Mesh(structured=True, structured_coords=sc)
   
    # Generate 42 photon energy bins [eV] 
    #  First bin has been replaced with 1 for log interpolation
    photon_bins =  [1, 1e4, 2e4, 3e4, 4.5e4, 6e4, 7e4, 7.5e4, 1e5, 1.5e5, 2e5, 3e5, 4e5,
                    4.5e5, 5.1e5, 5.12e5, 6e5, 7e5, 8e5, 1e6, 1.33e6, 1.34e6, 1.5e6, 1.66e6, 2e6,
                    2.5e6, 3e6, 3.5e6, 4e6, 4.5e6, 5e6, 5.5e6, 6e6, 6.5e6, 7e6, 7.5e6, 8e6, 1e7,
                    1.2e7, 1.4e7, 2e7, 3e7, 5e7]
    # Convert to MEV
    photon_bins = np.array([x/1E6 for x in photon_bins])
    # ICRP 74 flux-to-dose conversion factors
    de = np.array([0.01, 0.015, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.1, 0.15, 0.2, 0.3,
                   0.4, 0.5, 0.6, 0.8, 1, 2, 4, 6, 8, 10])
    df = np.array([0.0485, 0.1254, 0.205, 0.2999, 0.3381, 0.3572, 0.378, 0.4066, 0.4399, 0.5172,
                   0.7523, 1.0041, 1.5083, 1.9958, 2.4657, 2.9082, 3.7269, 4.4834, 7.4896,
                   12.0153, 15.9873, 19.9191, 23.76])
    # Convert to Sv/s per photon FLUX (not fluence)
    df = np.array([x*1E-12 for x in df])
    # Convert pointwise data to group data for log interpolation
    photon_spectrum = pointwise_collapse(photon_bins, de, df, logx=True, logy=True)
    #  Anything below 0.01 MeV should be assigned the DF value of 0.01 MeV
    photon_spectrum[0] = df[0]
    # Total number of groups is 217 (42 photon + 175 neutron)
    spectra = [np.append(photon_spectrum, np.zeros(175))]
    # The spectrum is normalized by PyNE, so we need to mutliply by the sum of intensities in the spectrum.
    # Additionally, we divide by the volume of the source cell in order to get source density.
    intensities = [np.sum(spectra)/src_vol]
    
    # Load geometry into DAGMC 
    load(geom)
    # Generate isotropic photon volume source
    source, dg = isotropic_vol_source(geom, mesh, cells, spectra, intensities)
   
    # PARTISN input 
    ngroup = 217
    cards = _cards()
    names_dict = _names_dict()
    
    write_partisn_input(mesh, geom, ngroup, cards=cards, dg=dg, names_dict=names_dict, data_hdf5path="/materials", nuc_hdf5path="/nucid", fine_per_coarse=1)

def main():

    gtcadis_help = ('This script automates the GT-CADIS process of \n'
                    'producing variance reduction parameters to optimize the\n'
                    'neutron transport step of the Rigorous 2-Step (R2S) method.\n')
    setup_help = ('Prints the file "config.ini" to be\n'
                  'filled in by the user.\n')
    step1_help = 'Creates the PARTISN input file for adjoint photon transport.'
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(help=gtcadis_help, dest='command')

    setup_parser = subparsers.add_parser('setup', help=setup_help)
    step1_parser = subparsers.add_parser('step1', help=step1_help)

    args, other = parser.parse_known_args()
    if args.command == 'setup':
        setup()
    elif args.command == 'step1':
        step1()

if __name__ == '__main__':
    main()
