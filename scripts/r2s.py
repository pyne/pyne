import argparse
import ConfigParser

from pyne.dagmc import cell_materials
from pyne.r2s import irradiation_setup, photon_sampling_setup

config_filename = 'config.ini'
alara_params_filename = 'alara_params.txt'

config = \
"""[step1]
# Path to MCNP MESHTAL file containing neutron fluxes
meshtal: meshtal
# Tally number within the meshtal file containing the fluxes for activation.
tally_num: 4
# Path to the DAGMC material-laden geometry.
geom: geom.h5m
reverse: True
# Number of rays to fire down each mesh row in each direction to calculate
# cell volume fractions
num_rays: 10
grid: False

[step2]
# List of decays times, seperated by commas. These strings much match exactly with
# their counterparts in the "cooling" block within alara_params.txt
decay_times = "1E3 s", "12 h", "3.0 d" 
"""

alara_params =\
"""material_lib alara_matlib
element_lib nuclib
data_library alaralib fendl3bin

cooling
    1E3 s
    12 h
    3.0 d
end

output zone
       integrate_energy
       photon_source  fendl3bin  phtn_src1 24 1.00E4 2.00E4 5.00E4 1.00E5
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

    meshtal = config.get('step1', 'meshtal')
    tally_num = config.get('step1', 'tally_num')
    geom = config.get('step1', 'geom')
    reverse = config.get('step1', 'reverse')
    num_rays = config.get('step1', 'num_rays')
    grid = config.get('step1', 'grid')

    cell_mats = cell_materials(geom)
    irradiation_setup(meshtal, cell_mats, alara_params_filename, tally_num,
                      num_rays=num_rays, grid=grid, reverse=reverse)
  
    print("R2S step1 complete, run ALARA with the command:")
    print(">> alara alara_inp > output.txt")

def step2():
    config = ConfigParser.ConfigParser()
    config.read(config_filename)

    decay_times = config.get('step2', 'decay_times')
    print(decay_times)

def main():

    r2s_help = ""
    setup_help = ""
    step1_help = ""
    step2_help = ""
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(help=r2s_help, dest='command')

    setup_parser = subparsers.add_parser("setup", help=setup_help)
    step1_parser = subparsers.add_parser("step1", help=step1_help)
    step2_parser = subparsers.add_parser("step2", help=step2_help)

    args, other = parser.parse_known_args()
    if args.command == 'setup':
        setup()
    if args.command == 'step1':
        step1()
    elif args.command == 'step2':
        step2()

if __name__ == '__main__':
    main()

