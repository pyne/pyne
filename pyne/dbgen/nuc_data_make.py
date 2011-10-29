import os
import argparse
from distutils.file_util import copy_file, move_file
from distutils.dir_util import mkpath, remove_tree

from pyne.api import nuc_data
from pyne.dbgen.api import build_dir
from pyne.dbgen.decay import make_decay
from pyne.dbgen.atomic_weight import make_atomic_weight
from pyne.dbgen.scattering_lengths import make_scattering_lengths
from pyne.dbgen.simple_xs import make_simple_xs
from pyne.dbgen.cinder import make_cinder

# Thanks to http://patorjk.com/software/taag/
# and http://www.chris.com/ascii/index.php?art=creatures/dragons (Jeff Ferris)
# for ASCII art inspiriation

pyne_logo = """\

                                  /   \       
 _                        )      ((   ))     (                          
(@)                      /|\      ))_((     /|\                          
|-|                     / | \    (/\|/\)   / | \                      (@) 
| | -------------------/--|-voV---\`|'/--Vov-|--\---------------------|-|
|-|                         '^`   (o o)  '^`                          | |
| |                               `\Y/'                               |-|
|-|                                                                   | |
| |        /\             ___           __  __             /\         |-|
|-|       /^~\           / _ \_   _  /\ \ \/__\           /^~\        | |  
| |       /^~\          / /_)/ | | |/  \/ /_\             /^~\        |-|
|-|       /^~\         / ___/| |_| / /\  //__             /^~\        | | 
| |       ^||`         \/     \__, \_\ \/\__/             ^||`        |-|  
|-|        ||                |____/                        ||         | | 
| |       ====                                            ====        |-|
|-|                                                                   | |
| |                                                                   |-|
|-|___________________________________________________________________| |
(@)              l   /\ /         ( (       \ /\   l                `\|-|
                 l /   V           \ \       V   \ l                  (@)
                 l/                _) )_          \I                   
                                   `\ /'
                                     `  
"""


def main():
    """Entry point for nuc_data_make utility."""
    print pyne_logo

    make_funcs = [('atomic_weight', make_atomic_weight),
                  ('scattering_lengths', make_scattering_lengths),
                  ('decay', make_decay), 
                  ('simple_xs', make_simple_xs), 
                  ('cinder', make_cinder), 
                 ]
    make_map = dict(make_funcs)

    # Parse the command line arguments
    parser = argparse.ArgumentParser(description='Make a nuclear data library.')
    parser.add_argument('-o', dest='nuc_data', action='store', default=nuc_data,
                        help='path to the output database file.')
    parser.add_argument('-b', dest='build_dir', action='store', default=build_dir,
                        help='path to the build directory.')
    parser.add_argument('-m', dest='make', action='store', default='all',
                        help='comma-separated parts of nuc_data to make: ' + \
                        ", ".join([mf[0] for mf in make_funcs]) + ', all, and none.')
    parser.add_argument('--clean', dest='clean', type=int, default=0,
                        help="""level to clean up existing files.
                                0: no cleaning (default).
                                1: clean nuc_data.
                                2: clean nuc_data and build_dir.""")
    args = parser.parse_args()

    # clean nuc data
    if args.clean in [1, 2]:
        print "removing nuc_data from {0}".format(args.nuc_data)
        try:
            os.remove(args.nuc_data)
        except OSError:
            pass

    # Make the build dir
    if args.clean == 2:
        print "removing build_dir from {0}".format(args.build_dir)
        remove_tree(args.build_dir)
    mkpath(args.build_dir)

    # Determine what to make
    if args.make == 'none':
        make_order = []
    elif args.make == 'all':
        make_order = [mf[0] for mf in make_funcs]
    else:
        make_order = args.make.replace(' ', "").split(',')

    print "Making nuc_data at {0}".format(args.nuc_data)

    # Make the various tables
    for mo in make_order:
        make_map[mo](args.nuc_data, args.build_dir)


if __name__ == '__main__':
    main()
