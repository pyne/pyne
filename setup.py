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

VERSION = '0.5-dev'
IS_NT = os.name == 'nt'

CMAKE_BUILD_TYPES = {
    'none': 'None', 
    'debug': 'Debug',
    'release': 'Release',
    'relwithdebinfo': 'RelWithDebInfo',
    'minsizerel': 'MinSizeRel',
    }

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

def parse_setup(ns):
    a = [sys.argv[0], ns.cmd]
    if ns.user:
        a.append('--user')
    if ns.egg_base is not None:
        local_path = os.path.dirname(os.path.abspath(sys.argv[0]))
        a.append('--egg-base='+os.path.join(local_path, ns.egg_base))
    return a

def parse_cmake(ns):
    a = []
    if ns.D is not None:
        a += ['-D' + x for x in ns.D]
    if ns.build_type is not None:
        a.append('-DCMAKE_BUILD_TYPE=' + CMAKE_BUILD_TYPES[ns.build_type.lower()])
    return a

def parse_make(ns):
    a = []
    if ns.j is not None:
        a.append('-j' + ns.j)
    return a

def parse_others(ns):
    if ns.hdf5 is not None:
        os.environ['HDF5_ROOT'] = ns.hdf5

def parse_args():
    parser = argparse.ArgumentParser()

    setup = parser.add_argument_group('setup', 'Group for normal setup.py arguments')
    setup.add_argument('cmd', help="command to send to normal setup, e.g. "
                       "install or build.")
    parser.add_argument('--user', nargs='?', const=True, default=False)
    parser.add_argument('--egg-base')

    cmake = parser.add_argument_group('cmake', 'Group for CMake arguments.')
    cmake.add_argument('-D', metavar='VAR', action='append', 
                       help='Set enviornment variable.')
    cmake.add_argument('--build-type', metavar='BT', 
                       help='Set build type via CMAKE_BUILD_TYPE, '
                            'e.g. Release or Debug.')

    make = parser.add_argument_group('make', 'Group for make arguments.')
    make.add_argument('-j', help='Degree of parallelism for build.')

    other = parser.add_argument_group('other', 'Group for miscellaneous arguments.')
    other.add_argument('--hdf5', help='Path to HDF5 root directory.')

    ns = parser.parse_args()
    sys.argv = parse_setup(ns)
    cmake_args = parse_cmake(ns)
    make_args = parse_make(ns)
    parse_others(ns)

    return cmake_args, make_args


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
        "version": VERSION,
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

def cmake_cli(cmake_args):
    if not IS_NT:
        rtn = subprocess.call(['which', 'cmake'])
        if rtn != 0:
            sys.exit('CMake is not installed, aborting PyNE build.')
    cmake_cmd = ['cmake', '..'] + cmake_args
    cmake_cmd += ['-DPYTHON_EXECUTABLE=' + sys.executable]
    if IS_NT:
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
    cmake_cmdstr = cmake_cmd if isinstance(cmake_cmd, str) else ' '.join(cmake_cmd)
    print("CMake command is\n", cmake_cmdstr, sep="")
    return cmake_cmd

def main_body():
    assert_dep_versions()
    cmake_args, make_args = parse_args()
    if not os.path.exists('build'):
        os.mkdir('build')
    cmake_cmd = cmake_cli(cmake_args)
    rtn = subprocess.check_call(cmake_cmd, cwd='build', shell=IS_NT)
    
    rtn = subprocess.check_call(['make'] + make_args, cwd='build')

    cwd = os.getcwd()
    os.chdir('build')
    setup()
    os.chdir(cwd)

def final_message(success=True):
    if success:
        return
    msg = ("\n\nYou seem to be having issues building pyne. please report your problem "
           "to pyne-dev@googlegroups.com or look for help at http://pyne.io\n\n"
           )
    print("-"*20 + msg + '-'* 20)


def main():
    success = False
    try:
        main_body()
        success = True
    finally:
        final_message(success)
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

if __name__ == "__main__":
    main()
