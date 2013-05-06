#!/usr/bin/env python
"""Read a MCNP Ptrac file and save it in HDF5 format."""
import sys
import tables
from pyne import mcnp

if len(sys.argv) <= 3:
    sys.stderr.write("Not enough command line arguments.\n")
    sys.stderr.write("\n")
    sys.stderr.write("Usage:\n")
    sys.stderr.write("{0} ptrac_file hdf5_file tablename\n".format(sys.argv[0]))
    exit(-1)

ptrac_filename = sys.argv[1]
hdf5_filename = sys.argv[2]
tablename = sys.argv[3]

ptrac = mcnp.PtracReader(ptrac_filename)

# open HDF5 file and create table if it doesn't exist yet
h5file = tables.openFile(hdf5_filename, mode="a", title=problem_title)
tablepath = "/" + tablename
if tablepath in h5file:
    table = h5file.getNode(tablepath)
else:
    table = h5file.createTable(group, tablename, mcnp.PtracEvent, "Ptrac data")

ptrac.write_to_hdf5_table(table, print_progress=1000000)

table.flush()
h5file.close()