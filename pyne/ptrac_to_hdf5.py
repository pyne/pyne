#!/usr/bin/env python
"""Read a MCNP Ptrac file and save it in HDF5 format."""
from pyne.utils import QA_warn

import tables

from . import mcnp

try:
    import argparse
except ImportError:
    from . import _argparse as argparse

QA_warn(__name__)


def main():
    argparser = argparse.ArgumentParser(
        description="write the contents of a MCNP PTRAC file to a HDF5 table"
    )
    argparser.add_argument("ptrac_file", help="MCNP PTRAC file to read from")
    argparser.add_argument(
        "hdf5_file", help="HDF5 file to write to (will be created if it does not exist)"
    )
    argparser.add_argument(
        "-n",
        "--table-name",
        default="ptrac",
        help='name of the HDF5 table (default is "ptrac")',
    )
    argparser.add_argument(
        "-t",
        "--table-title",
        default="Ptrac data",
        help='title of the HDF5 table (default is "Ptrac data")',
    )
    argparser.add_argument(
        "-s", "--show-progress", action="store_true", help="show progress indicator"
    )
    args = argparser.parse_args()

    ptrac_filename = args.ptrac_file
    hdf5_filename = args.hdf5_file
    table_name = args.table_name
    table_title = args.table_title
    print_progress = 1000000 if args.show_progress else 0

    ptrac = mcnp.PtracReader(ptrac_filename)

    # open HDF5 file and create table if it doesn't exist yet
    h5file = tables.open_file(hdf5_filename, mode="a", title=ptrac.problem_title)
    table_path = "/" + table_name
    if table_path in h5file:
        table = h5file.get_node(table_path)
    else:
        table = h5file.create_table("/", table_name, mcnp.PtracEvent, table_title)

    ptrac.write_to_hdf5_table(table, print_progress=print_progress)

    table.flush()
    h5file.close()


if __name__ == "__main__":
    main()
