#!/usr/bin/env python
from pyne.mcnp import Meshtal
import argparse

def main():
    parser = argparse.ArgumentParser(description=(
             'Reads an MCNP meshtal file and creates a an h5m mesh file '
             'for each meshtally within the file. The output mesh file are'
             'named <filename>_tally_<tally_num>.h5m'))
    parser.add_argument('filename', help='Name of the MCNP meshtal file')
    args = parser.parse_args()
    
    m = Meshtal(args.filename)
    for num, tal in m.tally.items():
        tal.mesh.save("{0}_tally_{1}.h5m".format(args.filename, num))

if __name__ == '__main__':
    main()
