import os
import argparse
from distutils.file_util import copy_file, move_file
from distutils.dir_util import mkpath, remove_tree

from pyne.api import nuc_data
from pyne.dbgen.api import build_dir
from pyne.dbgen.atomic_weight import make_atomic_weight

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
    print pyne_logo

    # Parse the command line arguments
    parser = argparse.ArgumentParser(description='Make a nuclear data library.')
    parser.add_argument('-o', dest='nuc_data', action='store', default=nuc_data,
                        help='path to the output database file.')
    parser.add_argument('-b', dest='build_dir', action='store', default=build_dir,
                        help='path to the build directory.')
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

    # Stop after cleaning
    if 0 < args.clean:
        return

    print "Making nuc_data at {0}".format(args.nuc_data)

    # make atomic weight table
    make_atomic_weight(args.nuc_data, args.build_dir)



if __name__ == '__main__':
    main()
