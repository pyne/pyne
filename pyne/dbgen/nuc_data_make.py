import os
import argparse
import urllib2
import shutil
from distutils.file_util import copy_file, move_file
from distutils.dir_util import mkpath, remove_tree

from pyne.api import nuc_data
from pyne.utils import message
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

def _fetch_prebuilt(args):
    nuc_data, build_dir = args.nuc_data, args.build_dir
    prebuilt_nuc_data = os.path.join(build_dir, 'prebuilt_nuc_data.h5')
    prebuilt_nuc_data_url = "http://s3.amazonaws.com/pyne/prebuilt_nuc_data.h5"

    if not os.path.exists(prebuilt_nuc_data):
        print "Fetching pre-built nuc_data.h5 from " + prebuilt_nuc_data_url
        pnd = urllib2.urlopen(prebuilt_nuc_data_url)
        with open(prebuilt_nuc_data, 'wb') as f:
            f.write(pnd.read())

    if not os.path.exists(nuc_data):
        shutil.copyfile(prebuilt_nuc_data, nuc_data)


def main():
    """Entry point for nuc_data_make utility."""
    print message(pyne_logo)

    make_funcs = [('atomic_weight', make_atomic_weight),
                  ('scattering_lengths', make_scattering_lengths),
                  ('decay', make_decay), 
                  ('simple_xs', make_simple_xs), 
                  ('cinder', make_cinder), 
                  ]
    make_map = dict(make_funcs)
    make_open = set(['atomic_weight', 'scattering_lengths', 'simple_xs'])

    # Parse the command line arguments
    parser = argparse.ArgumentParser(description='Make a nuclear data library.')
    parser.add_argument('-o', dest='nuc_data', action='store', default=nuc_data,
                        help='path to the output database file.')
    parser.add_argument('-b', dest='build_dir', action='store', default=build_dir,
                        help='path to the build directory.')
    parser.add_argument('--datapath', dest='datapath', action='store', default="",
                        help='MCNP DATAPATH.')
    parser.add_argument('--fetch-prebuilt', dest='fetch_prebuilt', action='store', 
                        type=lambda s: 't' in s.lower() or 'y' in s.lower(), 
                        default=True, help='grab partially assembled file [y/n].')
    parser.add_argument('--make-open-only', dest='make_open_only', action='store', 
                        type=lambda s: 't' in s.lower() or 'y' in s.lower(), 
                        default=False, help='only add open data to file [y/n].')
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
        print "Removing nuc_data from {0}".format(args.nuc_data)
        try:
            os.remove(args.nuc_data)
        except OSError:
            pass

    # Make the build dir
    if args.clean == 2:
        print "Removing build_dir from {0}".format(args.build_dir)
        remove_tree(args.build_dir)
    mkpath(args.build_dir)

    # Determine what to make
    if args.make == 'none':
        make_order = []
    elif args.make == 'all':
        make_order = [mf[0] for mf in make_funcs]
    else:   
        make_order = args.make.replace(' ', "").split(',')

    if args.make_open_only:
        make_order = [mo for mo in make_order if mo in make_open]

    # fetch prebuilt data library if possible
    if args.fetch_prebuilt:
        _fetch_prebuilt(args)

    # Make the various tables
    print "Making nuc_data at {0}".format(args.nuc_data)
    for mo in make_order:
        make_map[mo](args)


if __name__ == '__main__':
    main()
