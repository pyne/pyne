#!/usr/bin/env python
"""Welcome to PyNE's setup.py script. This is a little non-standard because pyne
is a multilanguage projects.  Still this script follows a predicatable ordering:

1. Parse command line arguments,
2. Call cmake from the 'build' directory
3. Call make from the 'build' directory
4. Use distuitls/setuptools from the 'build' directory

This gives us the best of both worlds. Compiled code is installed with cmake/make
and Cython/Python code is installed with normal Python tools. The only trick here is
how the various command line arguments are handed off to the three sub-processes.

To acomplish this we use argparser groups to group command line arguments based on
whether they go to:

1. the setup() function,
2. cmake,
3. make, or
4. other - typically used for args that apply to multiple other groups or
   modify the environment in some way.

To add a new command line argument, first add it to the appropriate group in the
``parse_args()`` function.  Then, modify the logic in the cooresponding
``parse_setup()``, ``parse_cmake()``, ``parse_make()``, or ``parse_others()``
functions to consume your new command line argument.  It is OK for more than
one of the parser functions to comsume the argument. Where appropriate,
ensure the that argument is appended to the argument list that is returned by these
functions.
"""
from __future__ import print_function

import io
import os
import sys
import imp
import shutil
import tarfile
import argparse
import platform
import warnings
import subprocess
from glob import glob
from distutils import core, dir_util
from contextlib import contextmanager
if sys.version_info[0] < 3:
    from urllib import urlopen
else:
    from urllib.request import urlopen
try:
    from setuptools import setup as _setup
    have_setuptools = True
except ImportError:
    from distutils.core import setup as _setup
    have_setuptools = False

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

VERSION = '0.5.0-rc1'
IS_NT = os.name == 'nt'

CMAKE_BUILD_TYPES = {
    'none': 'None',
    'debug': 'Debug',
    'release': 'Release',
    'relwithdebinfo': 'RelWithDebInfo',
    'minsizerel': 'MinSizeRel',
    }


@contextmanager
def indir(path):
    orig = os.getcwd()
    os.chdir(path)
    yield
    os.chdir(orig)


@contextmanager
def cleanpypath(path):
    orig = sys.path
    sys.path = [p for p in sys.path if p != path]
    yield
    sys.path = orig


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
    except ImportError:
        return
    low = (1, 2, 1)
    v = IPython.__version__.split('-')[0]
    cur = tuple(map(int, v.split('.')))
    if cur < low:
        msg = "ipython version is too low! {0} (have) < 2.0.0 (min)".format(v)
        raise ValueError(msg)


def assert_ubuntu_version():
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


DECAY_H = os.path.join('src', 'decay.h')
DECAY_CPP = os.path.join('src', 'decay.cpp')
DECAY_H_REP = os.path.join('src', '_decay.h')
DECAY_CPP_REP = os.path.join('src', '_decay.cpp')
DECAY_URL = 'http://data.pyne.io/decay.tar.gz'


def download_decay():
    print('Downloading ' + DECAY_URL)
    try:
        durl = urlopen(DECAY_URL)
        d = durl.read()
        durl.close()
    except IOError:
        print('...failed!')
        return False
    f = io.BytesIO(d)
    tar = tarfile.open(fileobj=f, mode='r:gz')
    tar.extractall('src')
    tar.close()
    durl.close()
    return True


def generate_decay():
    with indir('src'):
        try:
            import decaygen
        except ImportError:
            return False
        try:
            decaygen.build()
        except Exception:
            return False
    return True


def ensure_decay():
    mb = 1024**2
    if os.path.isfile(DECAY_H) and os.path.isfile(DECAY_CPP) and \
       os.stat(DECAY_CPP).st_size > mb:
        return
    downloaded = download_decay()
    if downloaded:
        return
    generated = generate_decay()
    if generated:
        return
    print('!'*42)
    print('Decay files could not be downloaded or generated, using surrogates instead.')
    print('Please consider using the --bootstrap command line argument.')
    print('!'*42 + '\n')
    shutil.copy(DECAY_H_REP, DECAY_H)
    shutil.copy(DECAY_CPP_REP, DECAY_CPP)


def ensure_nuc_data():
    import tempfile
    tdir = tempfile.gettempdir()
    with cleanpypath('.'), cleanpypath(os.getcwd()), indir(tdir):
        from pyne.dbgen import nuc_data_make
        from pyne.dbgen.api import build_dir
        bdir = os.path.join(os.getcwd(), 'build', build_dir)
        nuc_data_make.main(args=['-b', bdir])

def parse_setup(ns):
    a = [sys.argv[0], ns.cmd]
    if ns.user:
        a.append('--user')
    if ns.prefix is not None:
        a.append('--prefix=' + ns.prefix)
    if ns.egg_base is not None:
        local_path = os.path.dirname(os.path.abspath(sys.argv[0]))
        a.append('--egg-base=' + os.path.join(local_path, ns.egg_base))
    if ns.cmd == 'clean':
        if os.path.exists('build'):
            dir_util.remove_tree('build')
        print('build directory cleaned ... exiting')
        sys.exit()
    if ns.clean:
        if os.path.exists('build'):
            dir_util.remove_tree('build')
    return a


def parse_cmake(ns):
    a = []
    if ns.D is not None:
        a += ['-D' + x for x in ns.D]
    if ns.build_type is not None:
        a.append('-DCMAKE_BUILD_TYPE=' + CMAKE_BUILD_TYPES[ns.build_type.lower()])
    if ns.prefix is not None:
        a.append('-DCMAKE_INSTALL_PREFIX=' + ns.prefix)
    if have_setuptools:
        a.append('-DHAVE_SETUPTOOLS=TRUE')
    return a


def parse_make(ns):
    a = []
    if ns.j is not None:
        a.append('-j' + ns.j)
    return a


def parse_others(ns):
    if ns.hdf5 is not None:
        os.environ['HDF5_ROOT'] = ns.hdf5
    if ns.moab is not None:
        os.environ['MOAB_ROOT'] = ns.moab
    if ns.clean:
        if os.path.isfile(DECAY_H):
            os.remove(DECAY_H)
        if os.path.isfile(DECAY_CPP):
            os.remove(DECAY_CPP)


def parse_args():
    argv = [a for a in sys.argv[1:] if a != '--']  # needed for backwards compat.
    parser = argparse.ArgumentParser()

    setup = parser.add_argument_group('setup', 'Group for normal setup.py arguments')
    setup.add_argument('cmd', help="command to send to normal setup, e.g. "
                       "install or build.")
    parser.add_argument('--clean', nargs='?', const=True, default=False)
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
    other.add_argument('--moab', help='Path to MOAB root directory.')
    other.add_argument('--prefix', help='Prefix for install location.')
    other.add_argument('--bootstrap', default=False, action='store_true', 
                       help='Bootstraps the PyNE installation, including '
                            'nuc_data_make and possibly decaygen.')

    ns = parser.parse_args(argv)
    sys.argv = parse_setup(ns)
    cmake_args = parse_cmake(ns)
    make_args = parse_make(ns)
    parse_others(ns)

    return cmake_args, make_args, ns.bootstrap


def setup():
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
        #"zip_safe": False,
        }
    rtn = _setup(**setup_kwargs)


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


def main_body(cmake_args, make_args):
    assert_dep_versions()
    ensure_decay()
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
    msg = ("\n\nIf you are having issues building pyne, please report your problem "
           "to pyne-dev@googlegroups.com or look for help at http://pyne.io\n\n"
           )
    print('\n' + '-'*20 + msg + '-'*20)


def main_safe(cmake_args, make_args):
    success = False
    try:
        main_body(cmake_args, make_args)
        success = True
    finally:
        final_message(success)


def main():
    cmake_args, make_args, bootstrap = parse_args()
    main_safe(cmake_args, make_args)
    if bootstrap:
        ensure_nuc_data()
        main_safe(cmake_args, make_args)
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
