import os

from pyne.api import nuc_data
from pyne.dbgen.api import build_dir
from pyne.dbgen.atomic_weight import make_atomic_weight

def main():
    # Make the build dir
    try:
        os.makedirs(build_dir)
    except OSError:
        pass

    print "Making nuc_data.h5 at: {0}".format(nuc_data)
    make_atomic_weight(nuc_data, build_dir)



if __name__ == '__main__':
    main()
