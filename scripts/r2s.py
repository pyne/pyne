#!/usr/bin/env python
import argparse
try:
    import ConfigParser
except ImportError:
    import configparser as ConfigParser

from os.path import isfile

from pyne.mesh import Mesh, NativeMeshTag
from pyne.dagmc import cell_materials, load, discretize_geom
from pyne.r2s import resolve_mesh, irradiation_setup, photon_sampling_setup,\
    total_photon_source_intensity
from pyne.alara import photon_source_to_hdf5, photon_source_hdf5_to_mesh,\
    phtn_src_energy_bounds
from pyne.mcnp import Meshtal

config_filename = 'config.ini'
alara_params_filename = 'alara_params.txt'

config = \
    """[general]
# Specify whether this problem uses structured or unstructured mesh
structured: True
# Specify whether this problem uses sub-voxel r2s
sub_voxel: False

[step1]
# Path to MCNP MESHTAL file containing neutron fluxes or a DAG-MCNP5
# unstructured mesh tally .h5m file.
meshtal: meshtal
# Tally number within the meshtal file containing the fluxes for activation.
tally_num: 4
# The name of the tag used to store flux data on the mesh. For unstructured
# mesh this tag must already exist within the file specified in <meshtal>.
flux_tag: n_flux
# Path to the DAGMC material-laden geometry file (.h5m).
geom: geom.h5m
# If True the fluxes in the fluxin file will be printed in the reverse
# order of how they appear within the flux vector tag. Since MCNP and
# the Meshtal class order fluxes from low energy to high energy, this
# option should be true if the transmutation data being used is
# ordered from high-energy to low-energy.
reverse: True
# Number of rays to fire down each mesh row in each direction to calculate
# cell volume fractions.
num_rays: 10
# If true, rays will be fired down mesh rows in evenly spaced intervals.
# In this case <num_rays> must be a perfect square. If false, rays are fired
# down mesh rows in random intervals.
grid: False

[step2]
# List of decays times, seperated by commas. These strings much match exactly
# with their counterparts in the phtn_src file produced in step1. No spaces
# should appear in this line except the space between the time and the time unit
# for each entry.
decay_times:1E3 s,12 h,3.0 d
# The prefix of the .h5m files containing the source density distributations for
# each decay time.
output: source
# The name of the output files containing the total photon source intensities for
# each decay time
tot_phtn_src_intensities : total_photon_source_intensities.txt
"""

alara_params =\
    """material_lib alara_matlib
element_lib nuclib
data_library alaralib fendl3bin

output zone
       integrate_energy
       # Energy group upper bounds. The lower bound is always zero.
       photon_source  fendl3bin  phtn_src 24 1.00E4 2.00E4 5.00E4 1.00E5
       2.00E5 3.00E5 4.00E5 6.00E5 8.00E5 1.00E6 1.22E6 1.44E6 1.66E6
       2.00E6 2.50E6 3.00E6 4.00E6 5.00E6 6.50E6 8.00E6 1.00E7 1.20E7
       1.40E7 2.00E7
end

#     flux name    fluxin file   norm   shift   unused
flux  my_flux     alara_fluxin  1e10     0      default

# Specify the irradiation schedule below.
# Syntax is found in the ALARA user manual
# This example is for a single 3.5 h pulse
schedule    my_schedule
    3.5 d my_flux my_pulse_history 0  s
end
pulsehistory  my_pulse_history
    1    0.0    s
end

#other parameters
truncation 1e-12
impurity 5e-6 1e-3
dump_file dump.file
"""


def setup():
    with open(config_filename, 'w') as f:
        f.write(config)
    with open(alara_params_filename, 'w') as f:
        f.write(alara_params)
    print('File "{}" has been written'.format(config_filename))
    print('File "{}" has been written'.format(alara_params_filename))
    print('Fill out the fields in these filse then run ">> r2s.py step1"')


def step1():
    config = ConfigParser.ConfigParser()
    config.read(config_filename)

    structured = config.getboolean('general', 'structured')
    sub_voxel = config.getboolean('general', 'sub_voxel')
    meshtal = config.get('step1', 'meshtal')
    tally_num = config.getint('step1', 'tally_num')
    flux_tag = config.get('step1', 'flux_tag')
    decay_times = config.get('step2', 'decay_times').split(',')
    geom = config.get('step1', 'geom')
    reverse = config.getboolean('step1', 'reverse')
    num_rays = config.getint('step1', 'num_rays')
    grid = config.getboolean('step1', 'grid')

    load(geom)

    # get meshtal info from meshtal file
    flux_mesh = resolve_mesh(meshtal, tally_num, flux_tag)

    # create the cell_fracs array before irradiation_steup
    if flux_mesh.structured:
        cell_fracs = discretize_geom(flux_mesh, num_rays=num_rays, grid=grid)
        # tag cell fracs for both default and subvoxel r2s modes
        flux_mesh.tag_cell_fracs(cell_fracs)
    else:
        cell_fracs = discretize_geom(flux_mesh)

    cell_mats = cell_materials(geom)
    irradiation_setup(flux_mesh, cell_mats, cell_fracs, alara_params_filename, tally_num,
                      num_rays=num_rays, grid=grid, reverse=reverse,
                      flux_tag=flux_tag, decay_times=decay_times,
                      sub_voxel=sub_voxel)

    # create a blank mesh for step 2:
    ves = list(flux_mesh.iter_ve())
    tags_keep = ("cell_number", "cell_fracs",
                 "cell_largest_frac_number", "cell_largest_frac")
    for tag in flux_mesh.get_all_tags():
        if tag.name not in tags_keep and isinstance(tag, NativeMeshTag):
            tag.delete()
    flux_mesh.write_hdf5('blank_mesh.h5m')
    print('The file blank_mesh.h5m has been saved to disk.')
    print('Do not delete this file; it is needed by r2s.py step2.\n')

    print('R2S step1 complete, run ALARA with the command:')
    print('>> alara alara_inp > output.txt')


def step2():
    config = ConfigParser.ConfigParser()
    config.read(config_filename)
    structured = config.getboolean('general', 'structured')
    sub_voxel = config.getboolean('general', 'sub_voxel')
    decay_times = config.get('step2', 'decay_times').split(',')
    output = config.get('step2', 'output')
    tot_phtn_src_intensities = config.get('step2', 'tot_phtn_src_intensities')
    tag_name = "source_density"

    if sub_voxel:
        geom = config.get('step1', 'geom')
        load(geom)
        cell_mats = cell_materials(geom)
    else:
        cell_mats = None
    h5_file = 'phtn_src.h5'
    if not isfile(h5_file):
        photon_source_to_hdf5('phtn_src')
    intensities = "Total photon source intensities (p/s)\n"
    for i, dc in enumerate(decay_times):
        print('Writing source for decay time: {0}'.format(dc))
        mesh = Mesh(structured=structured, mesh='blank_mesh.h5m')
        tags = {('TOTAL', dc): tag_name}
        photon_source_hdf5_to_mesh(mesh, h5_file, tags, sub_voxel=sub_voxel,
                                   cell_mats=cell_mats)
        mesh.write_hdf5('{0}_{1}.h5m'.format(output, i+1))
        intensity = total_photon_source_intensity(mesh, tag_name,
                                                  sub_voxel=sub_voxel)
        intensities += "{0}: {1}\n".format(dc, intensity)

    with open(tot_phtn_src_intensities, 'w') as f:
        f.write(intensities)

    e_bounds = phtn_src_energy_bounds("alara_inp")
    e_bounds_str = ""
    for e in e_bounds:
        e = e/1e6  # convert unit to MeV
        e_bounds_str += "{0}\n".format(e)
    with open("e_bounds", 'w') as f:
        f.write(e_bounds_str)

    print('R2S step2 complete.')


def main():

    r2s_help = ('This script automates the process of preforming Rigorous Two-\n'
                'Step (R2S) analysis using DAG-MCNP5 and the ALARA activation code.\n'
                'Infomation on how to use this script can be found at:\n'
                'http://pyne.io/usersguide/r2s.html\n')
    setup_help = ('Prints the files "config.ini" and "alara_params.txt, to be\n'
                  'filled in by the user.\n')
    step1_help = 'Creates the necessary files for running ALARA.'
    step2_help = 'Creates mesh-based photon sources for photon transport.'
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(help=r2s_help, dest='command')

    setup_parser = subparsers.add_parser('setup', help=setup_help)
    step1_parser = subparsers.add_parser('step1', help=step1_help)
    step2_parser = subparsers.add_parser('step2', help=step2_help)

    args, other = parser.parse_known_args()
    if args.command == 'setup':
        setup()
    elif args.command == 'step1':
        step1()
    elif args.command == 'step2':
        step2()


if __name__ == '__main__':
    main()
