#!/usr/bin/env python
import argparse

try:
    import ConfigParser
except ImportError:
    import configparser as ConfigParser

from os.path import isfile

from pyne.mesh import Mesh, NativeMeshTag
from pyne.dagmc import cell_materials, load, discretize_geom
from pyne.r2s import resolve_mesh, irradiation_setup
from pyne.alara import (
    photon_source_to_hdf5,
    photon_source_hdf5_to_mesh,
    phtn_src_energy_bounds,
    response_to_hdf5,
    response_hdf5_to_mesh,
)
from pyne.mcnp import Meshtal

config_filename = "config.ini"
alara_params_filename = "alara_params.txt"

config = """[general]
# Specify whether this problem uses structured or unstructured mesh
structured: True
# List of requested response, seperated by commas.
# Available options: 'decay_heat', 'specific_activity',
# 'alpha_heat', 'beta_heat', 'gamma_heat', 'wdr', 'photon_source'.
responses: decay_heat
# wdr file for response wdr.
wdr_file: IAEA.clearance

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
# ALARA output filename
alara_out: output.txt
"""

alara_params = """material_lib alara_matlib
element_lib nuclib
data_library alaralib fendl3bin

# Specify the flux condition below.
#     flux name    fluxin file   norm   shift   unused
flux  my_flux     alara_fluxin  1e10     0      default

# Specify the irradiation schedule below.
# Syntax is found in the ALARA user manual
# This example is for a single 3.5 d pulse
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
    with open(config_filename, "w") as f:
        f.write(config)
    with open(alara_params_filename, "w") as f:
        f.write(alara_params)
    print('File "{}" has been written'.format(config_filename))
    print('File "{}" has been written'.format(alara_params_filename))
    print(
        'Fill out the fields in these files then run ">> activation_response.py step1"'
    )


def step1():
    config = ConfigParser.ConfigParser()
    config.read(config_filename)

    structured = config.getboolean("general", "structured")
    meshtal = config.get("step1", "meshtal")
    tally_num = config.getint("step1", "tally_num")
    flux_tag = config.get("step1", "flux_tag")
    decay_times = config.get("step2", "decay_times").split(",")
    geom = config.get("step1", "geom")
    reverse = config.getboolean("step1", "reverse")
    num_rays = config.getint("step1", "num_rays")
    grid = config.getboolean("step1", "grid")
    responses = config.get("general", "responses").split(",")
    wdr_file = config.get("general", "wdr_file")

    load(geom)

    # get meshtal info from meshtal file
    flux_mesh = resolve_mesh(meshtal, tally_num, flux_tag)

    # create the cell_fracs array before irradiation_steup
    if flux_mesh.structured:
        cell_fracs = discretize_geom(flux_mesh, num_rays=num_rays, grid=grid)
    else:
        cell_fracs = discretize_geom(flux_mesh)

    cell_mats = cell_materials(geom)
    irradiation_setup(
        flux_mesh,
        cell_mats,
        cell_fracs,
        alara_params_filename,
        tally_num,
        num_rays=num_rays,
        grid=grid,
        reverse=reverse,
        output_mesh="activation_responses_step1.h5m",
        flux_tag=flux_tag,
        decay_times=decay_times,
        sub_voxel=False,
        responses=responses,
        wdr_file=wdr_file,
    )

    # create a blank mesh for step 2:
    ves = list(flux_mesh.iter_ve())
    for tag in flux_mesh.get_all_tags():
        tag.delete()
    flux_mesh.write_hdf5("blank_mesh.h5m")
    print("The file blank_mesh.h5m has been saved to disk.")
    print("Do not delete this file; it is needed by activation_responses.py step2.\n")
    print("Decay_heat step1 complete, run ALARA with the command:")
    print(">> alara alara_inp > output.txt")


def step2():
    config = ConfigParser.ConfigParser()
    config.read(config_filename)
    structured = config.getboolean("general", "structured")
    decay_times = config.get("step2", "decay_times").split(",")
    try:
        alara_out = config.get("step2", "alara_out")
    except:
        alara_out = "output.txt"
    responses = config.get("general", "responses").split(",")

    for response in responses:
        response_to_hdf5(alara_out, response)
        tag_name = response
        for i, dt in enumerate(decay_times):
            print("Writing {0} for decay time: {1}".format(response, dt))
            mesh = Mesh(structured=structured, mesh="blank_mesh.h5m")
            tags = {("TOTAL", dt): tag_name}
            response_hdf5_to_mesh(mesh, "".join([response, ".h5"]), tags, response)
            mesh.write_hdf5("{0}_{1}.h5m".format(response, i + 1))

    print("Activation_response step2 complete.")


def main():
    """
    Setup the activation responses calculation.
    """
    activation_response_help = (
        "This script automates the process of preforming\n"
        "activation analysis using DAG-MCNP5 and the ALARA activation code.\n"
    )
    setup_help = (
        'Prints the files "config.ini" and "alara_params.txt, to be\n'
        "filled in by the user.\n"
    )
    step1_help = "Creates the necessary files for running ALARA."
    step2_help = "Creates mesh-based activation response data for visulization."
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(help=activation_response_help, dest="command")

    setup_parser = subparsers.add_parser("setup", help=setup_help)
    step1_parser = subparsers.add_parser("step1", help=step1_help)
    step2_parser = subparsers.add_parser("step2", help=step2_help)

    args, other = parser.parse_known_args()
    if args.command == "setup":
        setup()
    elif args.command == "step1":
        step1()
    elif args.command == "step2":
        step2()


if __name__ == "__main__":
    main()
