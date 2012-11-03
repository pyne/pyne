#!/usr/bin/env python
 
import os
import sys
import glob
import json
from distutils.file_util import copy_file, move_file
from distutils.dir_util import mkpath, remove_tree
from copy import deepcopy

from Cython.Compiler.Version import version as CYTHON_VERSION



INFO = {
    'version': '0.1',
    }


def main():
    "Run functions specified on the command line"
    if len(sys.argv) <= 1:
        raise SystemExit("no command(s) specified")
    cmds = sys.argv[1:]
    if '-h' in cmds or '--help' in cmds:
        raise SystemExit("usage: " + sys.argv[0] + " <func-name> [<func-name>]")
    glbs = globals()
    for cmd in cmds:
        if cmd not in glbs:
            raise SystemExit(cmd + " not found")
    for cmd in cmds:
        if callable(glbs[cmd]):
            glbs[cmd]()
        else:
            raise SystemExit(cmd + " not callable")


def metadata(path):
    """Build a metadata file."""
    md = {"HDF5_DIR": HDF5_DIR,
          }
    md.update(INFO)

    # write the metadata file
    with open(path, 'w') as f:
        json.dump(md, f, indent=2)

    return md

def cleanup(mdpath="<None>"):
    # Clean includes after setup has run
    if os.path.exists('pyne/includes'):
        remove_tree('pyne/includes')

    # clean up metadata file
    if os.path.exists(mdpath):
        os.remove(mdpath)


def final_message(setup_success=True, metadata=None):
    if setup_success:
        return
        
    if metadata is not None:
        msg = "\n\nCURRENT METADATA:\n"
        for k, v in sorted(metadata.items()):
            msg += "  {0} = {1}\n".format(k, repr(v))
        print msg[:-1]

    if os.name != 'nt':
        return

    try: 
        import tables as tb
        h5ver = tb.getHDF5Version()
    except ImportError:
        h5ver = '1.8.5-patch1'

    msg = ("\n\nIf compilation is failing with HDF5 issues please try the "
           "following steps:\n\n"
           "    1. Install EPD [1].\n"
           "    2. Download the HDF5 Windows binarys from [2].\n"
           "    3. Unzip them to the C-drive (C:\\hdf5-{h5ver}).\n"
           "    4. Re-run setup with the '--hdf5' option:\n\n"
           "        python setup.py install --user --hdf5=C:\\hdf5-{h5ver}\n\n"
           "Should this still fail, please report your problem to pyne-dev@googlegroups.com\n\n"
           "[1] http://www.enthought.com/products/epd.php\n"
           "[2] http://www.hdfgroup.org/ftp/HDF5/releases/hdf5-{h5ver}/bin/windows/\n"
           ).format(h5ver=h5ver)
    print msg


def cython_version():
    pxi = ("# Cython compile-time version information\n"
           "DEF CYTHON_VERSION_MAJOR = {major}\n"
           "DEF CYTHON_VERSION_MINOR = {minor}\n"
           "DEF CYTHON_VERSION_MICRO = {micro}")
    cyver = CYTHON_VERSION.split('-')[0].split('.')
    while len(cyver) < 3:
        cyver = cyver + [0]
    cyver = dict([(k, int(cv)) for k, cv in zip(['major', 'minor', 'micro'], cyver)])
    pxi = pxi.format(**cyver)
    if not os.path.exists('pyne/includes'):
        os.mkdir('pyne/includes')
    with open('pyne/includes/cython_version.pxi', 'w') as f:
        f.write(pxi)


def includes():
    ds = ['pyne/includes', 'pyne/includes/pyne']
    cpfs = []

    for root, dirs, files in os.walk('cpp'):
        incroot = root.replace('cpp', 'pyne/includes')
        ds += [os.path.join(incroot, d) for d in dirs]
        cpfs += [(os.path.join(root, f), incroot) for f in files if f.endswith('.h')]

    for root, dirs, files in os.walk('pyne'):
        incroot = root.replace('pyne', 'pyne/includes/pyne')
        ds += [os.path.join(incroot, d) for d in dirs]
        cpfs += [(os.path.join(root, f), incroot) for f in files if f.endswith('.pxd')]

    for d in ds:
        mkpath(d)

    for src, dst in cpfs:
        copy_file(src, dst)



if __name__ == "__main__":
    main()
