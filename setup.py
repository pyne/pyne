#!/usr/bin/env python
from __future__ import print_function

import os
import sys
import imp
import argparse
import subprocess
from glob import glob

import numpy as np

# Thanks to http://patorjk.com/software/taag/
# and http://www.chris.com/ascii/index.php?art=creatures/dragons
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

def assert_np_version():
    low = (1, 8, 0)
    v = np.version.short_version
    cur = tuple(map(int, v.split('.')))
    if cur < low:
        msg = "numpy version too low! {0} (have) < 1.8.0 (min)".format(v)
        raise ValueError(msg)

def assert_ipython_version():
    try:
        import IPython
        low = (1, 2, 1)
        v = IPython.__version__.split('-')[0]
        cur = tuple(map(int, v.split('.')))
        if cur < low:
            msg = "ipython version is too low! {0} (have) < 2.0.0 (min)".format(v)
            raise ValueError(msg)
    except ImportError:
        pass;

def assert_ubuntu_version():
    import platform,warnings
    v = platform.uname()
    for itm in v:
        if 'precise' in itm:
            msg = ("ubuntu 12/precise packages may be outdated, it is highly "
                   "recommended to update to ubuntu 14 LTS.")
            warnings.warn(msg, Warning)

def assert_dep_versions():
    assert_np_version()
    assert_ubuntu_version()
    assert_ipython_version()
        
def parse_args():
    distutils_args = []
    cmake = []
    make = []
    argsets = [distutils_args, cmake, make]
    i = 0
    for arg in sys.argv:
        if arg == '--':
            i += 1
        else:
            argsets[i].append(arg)
    # handle HDF5
    hdf5opt = [o.split('=')[1] for o in distutils_args \
               if o.startswith('--hdf5=')]
    if 0 < len(hdf5opt):
        os.environ['HDF5_ROOT'] = hdf5opt[0]  # Expose to CMake
        distutils_args = [o for o in distutils_args \
                          if not o.startswith('--hdf5=')]

    # handle build type
    btopt = [o.split('=')[1] for o in distutils_args \
             if o.startswith('--build-type=')]
    if 0 < len(btopt):
        cmake.append('-DCMAKE_BUILD_TYPE=' + btopt[0])
        distutils_args = [o for o in distutils_args \
                          if not o.startswith('--build-type=')]

    # Change egg-base entry to absolute path so it behaves as expected
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--egg-base')
    res, distutils_args = parser.parse_known_args(distutils_args)
    if res.egg_base is not None:
        local_path = os.path.dirname(os.path.abspath(sys.argv[0]))
        distutils_args.append('--egg-base='+os.path.join(local_path,
                                                         res.egg_base))
    return distutils_args, cmake, make

INFO = {'version': '0.5-dev'}

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
           "python setup.py <distutils-args> [-- <cmake-arg>] "
           "[-- <make-args>]\n CMake and make command line "
           "arguments are optional, but must be preceeded by '--'.\n"
           "Should this still fail, please report your problem to "
           "pyne-dev@googlegroups.com\n\n"
           ).format(h5ver=h5ver)
    print(msg)


def setup():
    from distutils import core
    scripts = [os.path.join('scripts', f) for f in os.listdir('scripts')]
    scripts = [s for s in scripts if (os.name == 'nt' and s.endswith('.bat'))
                                     or (os.name != 'nt' and
                                         not s.endswith('.bat'))]
    packages = ['pyne', 'pyne.dbgen', 'pyne.apigen', 'pyne.xs',
                'pyne.transmute', 'pyne.gui', 'pyne.cli']
    pack_dir = {
        'pyne': 'pyne',
        'pyne.xs': 'pyne/xs',
        'pyne.gui': 'pyne/gui',
        'pyne.cli': 'pyne/cli',
        'pyne.dbgen': 'pyne/dbgen',
        'pyne.apigen': 'pyne/apigen',
        'pyne.transmute': 'pyne/transmute',
        }
    extpttn = ['*.dll', '*.so', '*.dylib', '*.pyd', '*.pyo']
    pack_data = {
        'lib': extpttn,
        'pyne': ['*.pxd', 'include/*.h', 'include/*.pxi', 'include/*/*.h',
                 '*.inp', 'include/*/*/*.h', 'include/*/*/*/*.h', '*.json',
                 '_includes/*.txt', '_includes/*.pxd', '_includes/*/*',
                 '_includes/*/*/*'] + extpttn,
        'pyne.xs': ['*.pxd'] + extpttn,
        'pyne.gui': ['*.pyw'],
        'pyne.dbgen': ['*.html', '*.csv', 'abundances.txt', 'mass.mas12'],
        }
    libpynes = set()
    for ext in extpttn:
        libpynes |= set(glob('src/' + ext))
    data_files = [
        ('lib', libpynes),
        ('include/pyne', glob('../src/*.h')),
        ]
    setup_kwargs = {
        "name": "pyne",
        "version": INFO['version'],
        "description": 'The Nuclear Engineering Toolkit',
        "author": 'PyNE Development Team',
        "author_email": 'pyne-dev@googlegroups.com',
        "url": 'http://pyne.github.com/',
        "packages": packages,
        "package_dir": pack_dir,
        "package_data": pack_data,
        "data_files": data_files,
        "scripts": scripts,
        }
    rtn = core.setup(**setup_kwargs)


def main_body():
    if not os.path.exists('build'):
        os.mkdir('build')
    sys.argv, cmake_args, make_args = parse_args()
    makefile = os.path.join('build', 'Makefile')
    if not os.path.exists(makefile):
        if os.name != 'nt':
            rtn = subprocess.call(['which', 'cmake'])
            if rtn != 0:
                sys.exit('CMake is not installed, aborting PyNE build.')
        cmake_cmd = ['cmake', '..'] + cmake_args
        cmake_cmd += ['-DPYTHON_EXECUTABLE=' + sys.executable]
        if os.name == 'nt':
            files_on_path = set()
            for p in os.environ['PATH'].split(';')[::-1]:
                if os.path.exists(p):
                    files_on_path.update(os.listdir(p))
            if 'cl.exe' in files_on_path:
                pass
            elif 'sh.exe' in files_on_path:
                cmake_cmd += ['-G "MSYS Makefiles"']
            elif 'gcc.exe' in files_on_path:
                cmake_cmd += ['-G "MinGW Makefiles"']
            cmake_cmd = ' '.join(cmake_cmd)
        is_nt = os.name == 'nt'
        cmake_cmdstr = cmake_cmd if isinstance(cmake_cmd, str) else ' '.join(cmake_cmd)
        print("CMake command is\n", cmake_cmdstr, sep="")
        rtn = subprocess.check_call(cmake_cmd, cwd='build', shell=is_nt)
    rtn = subprocess.check_call(['make'] + make_args, cwd='build')
    cwd = os.getcwd()
    os.chdir('build')
    configure.setup()
    os.chdir(cwd)

def old_main():
    success = False
    try:
        main_body()
        success = True
    finally:
        configure.final_message(success)
    # trick to get install path
    abspath = os.path.abspath
    joinpath = os.path.join
    cwd = abspath(os.getcwd())
    pypath = [p for p in sys.path if len(p) > 0 and cwd != abspath(p)]
    try:
        _, pynepath, _ = imp.find_module('pyne', pypath)
    except ImportError:
        pynepath = "${HOME}/.local/python2.7/site-packages"
    libpath = abspath(joinpath(pynepath, '..', '..', '..'))
    binpath = abspath(joinpath(libpath, '..', 'bin'))
    msg = ("\nNOTE: If you have not done so already, please be sure that your PATH and "
           "LD_LIBRARY_PATH (or DYLD_FALLBACK_LIBRARY_PATH on Mac OSX) has been "
           "appropriately set to the install prefix of pyne. For this install of "
           "pyne you may add the following lines to your '~/.bashrc' file or "
           "equivalent:\n\n"
           "# PyNE Environment Settings\n"
           'export PATH="{binpath}:${{PATH}}"\n'
           'export LD_LIBRARY_PATH="{libpath}:${{LD_LIBRARY_PATH}}"'
           ).format(binpath=binpath, libpath=libpath)
    print(msg, file=sys.stderr)

def main():
    assert_dep_versions()
    parse_args()


if __name__ == "__main__":
    main()
