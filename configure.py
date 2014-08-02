#!/usr/bin/env python

from __future__ import print_function
import os
import sys
import json

INFO = {
    'version': '0.4-dev',
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


def metadata(path="pyne/metadata.json"):
    """Build a metadata file."""
    md = {}
    md.update(INFO)

    # FIXME: Add the contents of CMakeCache.txt to the metadata dictionary

    # write the metadata file
    with open(path, 'w') as f:
        json.dump(md, f, indent=2)

    return md


def final_message(success=True):
    if success:
        return

    metadata = None
    mdpath = os.path.join('pyne', 'metadata.json')
    if os.path.exists(mdpath):
        with open(mdpath) as f:
            metadata = json.load(f)
    if metadata is not None:
        msg = "\n\nCURRENT METADATA:\n"
        for k, v in sorted(metadata.items()):
            msg += "  {0} = {1}\n".format(k, repr(v))
        print(msg[:-1])

    if os.name != 'nt':
        return

    try: 
        import tables as tb
        h5ver = tb.getHDF5Version()
    except ImportError:
        h5ver = '1.8.5-patch1'

    msg = ("\n\nUSAGE: "
           "python setup.py <distutils-args> [-- <cmake-arg>] [-- <make-args>]\n"
           "CMake and make command line arguments are optional, but must be preceeded "
           "by '--'.\n"
           "\n\nIf compilation is failing with HDF5 issues please try the "
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
    print(msg)


def setup():
    from distutils import core
    scripts = [os.path.join('scripts', f) for f in os.listdir('scripts')]
    scripts = [s for s in scripts if (os.name == 'nt' and s.endswith('.bat')) or 
                                     (os.name != 'nt' and not s.endswith('.bat'))]
    packages = ['pyne', 'pyne.lib', 'pyne.dbgen', 'pyne.apigen', 'pyne.xs', 
                'pyne.transmute', 'pyne.gui', 'pyne.cli']
    pack_dir = {
        'pyne': 'pyne',
        'pyne.xs': 'pyne/xs',
        'pyne.lib': 'pyne/lib',
        'pyne.gui': 'pyne/gui',
        'pyne.cli': 'pyne/cli',
        'pyne.dbgen': 'pyne/dbgen',
        'pyne.apigen': 'pyne/apigen',
        'pyne.transmute': 'pyne/transmute',
        }
    extpttn = ['*.dll', '*.so', '*.dylib', '*.pyd', '*.pyo']
    pack_data = {
        'pyne': ['*.pxd', 'include/*.h', 'include/*.pxi', 'include/*/*.h', '*.inp',
                 'include/*/*/*.h', 'include/*/*/*/*.h', '*.json', '_includes/*.txt',
                 '_includes/*.pxd', '_includes/*/*', '_includes/*/*/*'] + extpttn,
        'pyne.xs': ['*.pxd'] + extpttn,
        'pyne.lib': extpttn,
        'pyne.gui': ['*.pyw'],
        'pyne.dbgen': ['*.html', '*.csv', 'abundances.txt', 'mass.mas12'],
        }
    setup_kwargs = {
        "name": "pyne",
        "version": INFO['version'],
        "description": 'The Nuclear Engineering Toolkit',
        "author": 'PyNE Development Team',
        "author_email": 'scopatz@gmail.com',
        "url": 'http://pyne.github.com/',
        "packages": packages,
        "package_dir": pack_dir,
        "package_data": pack_data,
        "scripts": scripts,
        }
    rtn = core.setup(**setup_kwargs)


if __name__ == "__main__":
    main()
