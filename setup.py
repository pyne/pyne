#!/usr/bin/env python
"""Welcome to PyNE's setup.py script. This is a little non-standard because pyne
is a multilanguage project.  Still this script follows a predictable ordering:

1. Parse command line arguments,
2. Call cmake from the 'build' directory
3. Call make from the 'build' directory
4. Use distutils/setuptools from the 'build' directory

This gives us the best of both worlds. Compiled code is installed with cmake/make
and Cython/Python code is installed with normal Python tools. The only trick here is
how the various command line arguments are handed off to the three sub-processes.

To accomplish this we use argparser groups to group command line arguments based on
whether they go to:

1. the setup() function,
2. cmake,
3. make, or
4. other - typically used for args that apply to multiple other groups or
   modify the environment in some way.

To add a new command line argument, first add it to the appropriate group in the
``parse_args()`` function.  Then, modify the logic in the corresponding
``parse_setup()``, ``parse_cmake()``, ``parse_make()``, or ``parse_others()``
functions to consume your new command line argument.  It is OK for more than
one of the parser functions to consume the argument. Where appropriate,
ensure that the argument is appended to the argument list that is returned by these
functions.
"""

from __future__ import print_function
import numpy as np
import io
import os
import re
import sys
import imp
import ssl
import shutil
import tarfile
import argparse
import platform
import warnings
import subprocess
from glob import glob
from distutils import core, dir_util, sysconfig
from contextlib import contextmanager
if sys.version_info[0] < 3:
    from urllib import urlopen
else:
    from urllib.request import urlopen


# import src into pythonpath - needed to actually run decaygen/atomicgen
if '.' not in sys.path:
    sys.path.append(os.getcwd() + '/src')


def absexpanduser(x): return os.path.abspath(os.path.expanduser(x))


VERSION = '0.5.11'
IS_NT = os.name == 'nt'
LOCALDIR = absexpanduser('~/.local')
CMAKE_BUILD_TYPES = {
    'none': 'None',
    'debug': 'Debug',
    'release': 'Release',
    'relwithdebinfo': 'RelWithDebInfo',
    'minsizerel': 'MinSizeRel',
}
ON_DARWIN = platform.system() == 'Darwin'
LIBEXT = '.dylib' if ON_DARWIN else '.so'

SKIP_OPTION = "SKIP"


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
    low = (1, 8)
    v = np.version.short_version
    cur = tuple(map(int, v.split('.')[:2]))
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


def ssl_context():
    # this is compitble for both Python 2 & 3
    # on Python 3, you can do just ssl.SSLContext()
    ctx = ssl.create_default_context()
    ctx.check_hostname = False
    ctx.verify_mode = ssl.CERT_NONE
    return ctx


DECAY_H = os.path.join('src', 'decay.h')
DECAY_CPP = os.path.join('src', 'decay.cpp')
DECAY_S = glob('src/decay*.s')
DECAY_H_REP = os.path.join('src', '_decay.h')
DECAY_CPP_REP = os.path.join('src', '_decay.cpp')
DECAY_URL = 'http://raw.githubusercontent.com/pyne/data/master/decay.tar.gz'


CRAM_H = os.path.join('src', 'cram.h')
CRAM_C = os.path.join('src', 'cram.c')
CRAM_S = glob('src/cram*.s')
CRAM_URL = 'http://raw.githubusercontent.com/pyne/data/master/cram.tar.gz'


local_ensdf_evaluators = ['alphad', 'delta', 'gtol', 'bldhst', 'hsicc', 'hsmrg',
                          'seqhst', 'logft', 'radd', 'ruler']
local_ensdf_tools = [['ensdf_processing/RADD/98AK04.in', '98AK04.in'],
                     ['ensdf_processing/RADD/ELE.in', 'ELE.in']]


def copy_ensdf_executables(exe_dest):
    print('Copying ENSDF Executables to install directory')
    # Hack for copying the executables the first time PyNE is installed, before
    # pyne has been added to the python path.
    if exe_dest[-4:] != 'pyne':
        exe_dest = sysconfig.get_python_lib()
        for f in os.listdir(sysconfig.get_python_lib()):
            if re.match('pyne', f):
                exe_dest = exe_dest + '/' + f
        exe_dest = exe_dest + '/pyne'
    for tool in local_ensdf_evaluators:
        try:
            local_path = os.path.join('build', os.path.join('src', tool))
            dest_path = os.path.join(exe_dest, tool)
            shutil.copy(local_path, dest_path)
        except Exception:
            print('Some ENSDF processing executables were unable to be copied to the \
                   install directory.')
    for tool in local_ensdf_tools:
        try:
            local_path = os.path.join('src', tool[0])
            dest_path = os.path.join(exe_dest, tool[1])
            shutil.copy(local_path, dest_path)
        except Exception:
            print('Some ENSDF processing executables were unable to be copied to the \
                   install directory.')


ATOMIC_H = os.path.join('src', 'atomic_data.h')
ATOMIC_CPP = os.path.join('src', 'atomic_data.cpp')
ATOMIC_H_UNDER = os.path.join('src', '_atomic_data.h')
ATOMIC_CPP_UNDER = os.path.join('src', '_atomic_data.cpp')


def generate_atomic():
    with indir('src'):
        try:
            import atomicgen
        except ImportError:
            return False
        try:
            atomicgen.main()
        except Exception:
            return False
    return True


def ensure_atomic():
    # generate the data
    generated = generate_atomic()
    if generated:
        return
    # last resort - if generate atomic failed, use the backup
    if not os.path.isfile(ATOMIC_H) and not os.path.isfile(ATOMIC_CPP):
        print("!!! Could not generate atomic data, using backup.")
        shutil.copy(ATOMIC_H_UNDER, ATOMIC_H)
        shutil.copy(ATOMIC_CPP_UNDER, ATOMIC_CPP)


def ensure_nuc_data():
    import tempfile
    tdir = tempfile.gettempdir()
    with cleanpypath('.'), cleanpypath(os.getcwd()), indir(tdir):
        from pyne.dbgen import nuc_data_make
        from pyne.dbgen.api import build_dir
        bdir = os.path.join(os.getcwd(), 'build', build_dir)
        nuc_data_make.main(args=['-b', bdir])


def update_setup_args(ns):
    if ns.user and ns.prefix is None:
        ns.prefix = LOCALDIR
    elif ns.prefix is not None:
        pass
    else:
        ns.prefix = sys.prefix

    files = [DECAY_H, DECAY_CPP, CRAM_H, CRAM_C] + DECAY_S + CRAM_S
    if ns.cmd == 'clean':
        if os.path.exists(ns.build_dir):
            dir_util.remove_tree(ns.build_dir)
        for f in files:
            if os.path.isfile(f):
                print('deleting ' + f)
                os.remove(f)
        print('build directory cleaned ... exiting')
        sys.exit()
    if ns.clean:
        if os.path.exists(ns.build_dir):
            dir_util.remove_tree(ns.build_dir)
        for f in files:
            if os.path.isfile(f):
                print('deleting ' + f)
                os.remove(f)


def update_cmake_args(ns):
    ns.cmake_args = ['-DCMAKE_INSTALL_PREFIX=' + ns.prefix]
    if ns.D is not None:
        ns.cmake_args += ['-D' + x for x in ns.D]
    if ns.build_type is not None:
        bt = CMAKE_BUILD_TYPES[ns.build_type.lower()]
        ns.cmake_args.append('-DCMAKE_BUILD_TYPE=' + bt)
    if ns.hdf5 is not None:
        h5root = absexpanduser(ns.hdf5)
        ns.cmake_args += [
            '-DHDF5_ROOT=' + h5root,
            '-DHDF5_LIBRARIES={0}/lib/libhdf5{1};{0}/lib/libhdf5_hl{1}'.format(
                h5root, LIBEXT),
            '-DHDF5_LIBRARY_DIRS=' + h5root + '/lib',
            '-DHDF5_INCLUDE_DIRS=' + h5root + '/include',
        ]
    if ns.moab is not SKIP_OPTION:
        ns.cmake_args.append('-DWITH_MOAB=ON')
        if ns.moab is not None:
            ns.cmake_args.append('-DMOAB_ROOT=' + absexpanduser(ns.moab))

    if ns.dagmc is not SKIP_OPTION:
        assert ns.moab is not SKIP_OPTION, "If the --dagmc option is present," \
                                           " --moab must be as well"
        ns.cmake_args.append('-DWITH_DAGMC=ON')
        if ns.dagmc is not None:
            ns.cmake_args.append('-DDAGMC_ROOT=' + absexpanduser(ns.dagmc))

    if ns.deps_root:
        ns.cmake_args.append('-DDEPS_ROOT_DIR=' + absexpanduser(ns.deps_root))
    if ns.fast is not None:
        fast = 'TRUE' if ns.fast else 'FALSE'
        ns.cmake_args.append('-DPYNE_FAST_COMPILE=' + fast)


def update_make_args(ns):
    ns.make_args = []
    if ns.j is not None:
        ns.make_args.append('-j' + ns.j)


def update_other_args(ns):
    if ns.hdf5 is not None:
        os.environ['HDF5_ROOT'] = ns.hdf5
    if ns.moab is not None:
        os.environ['MOAB_ROOT'] = ns.moab
    if ns.dagmc is not None:
        os.environ['DAGMC_ROOT'] = ns.dagmc


def parse_args():
    # needed for backwards compat.
    argv = [a for a in sys.argv[1:] if a != '--']
    parser = argparse.ArgumentParser()
    parser.add_argument('--clean', nargs='?', const=True, default=False,
                        help='removes the build directory before continuing.')
    parser.add_argument('--user', nargs='?', const=True, default=False,
                        help='Installs into ~/.local')

    setup = parser.add_argument_group('setup', 'Normal setup.py arguments')
    setup.add_argument('cmd', help="command to send to normal setup, e.g. "
                       "install or build.")
    setup.add_argument('--egg-base', help="does nothing, for compatability")

    cmake = parser.add_argument_group('cmake', 'CMake arguments.')
    cmake.add_argument('-D', metavar='VAR', action='append',
                       help='Set environment variable.')
    cmake.add_argument('--build-type', metavar='BT',
                       help='Set build type via CMAKE_BUILD_TYPE, '
                            'e.g. Release or Debug.')
    cmake.add_argument('--deps-root', default=None, dest='deps_root',
                       help="the path to the directory containing "
                            "all dependencies")
    cmake.add_argument('--fast', default=None, dest='fast',
                       action='store_true', help="Will try to compile "
                       "from assembly, if possible. This is faster than "
                       "compiling from source (default).")
    cmake.add_argument('--slow', dest='fast',
                       action='store_false', help="Will NOT try to compile "
                       "from assembly, if possible. This is slower as it "
                       "must compile from source.")

    make = parser.add_argument_group('make', 'Make arguments.')
    make.add_argument('-j', help='Degree of parallelism for build.')

    other = parser.add_argument_group('other', 'Miscellaneous arguments.')
    other.add_argument('--hdf5', help='Path to HDF5 root directory.')
    other.add_argument('--moab', help='Path to MOAB root directory.',
                       nargs='?', default=SKIP_OPTION)
    other.add_argument('--dagmc', help='Path to DAGMC root directory.',
                       nargs='?', default=SKIP_OPTION)
    other.add_argument('--prefix', help='Prefix for install location.',
                       default=None)
    other.add_argument('--build-dir', default='build', dest="build_dir",
                       help='where to place the build directory')
    other.add_argument('--bootstrap', default=False, action='store_true',
                       help='Bootstraps the PyNE installation, including '
                            'nuc_data_make and possibly decaygen.')

    ns = parser.parse_args(argv)
    update_setup_args(ns)
    update_cmake_args(ns)
    update_make_args(ns)
    update_other_args(ns)
    return ns


def cmake_cli(cmake_args):
    if not IS_NT:
        rtn = subprocess.call(['which', 'cmake'])
        if rtn != 0:
            sys.exit('CMake is not installed, aborting PyNE build.')
    cmake_cmd = ['cmake', '..'] + cmake_args
    cmake_cmd += ['-DPYTHON_EXECUTABLE=' + sys.executable]
    cmake_cmdstr = cmake_cmd if isinstance(
        cmake_cmd, str) else ' '.join(cmake_cmd)
    print("CMake command is\n", cmake_cmdstr, sep="")
    return cmake_cmd


def main_body(ns):
    assert_dep_versions()
    ensure_atomic()
    if not os.path.exists(ns.build_dir):
        os.mkdir(ns.build_dir)
    cmake_cmd = cmake_cli(ns.cmake_args)
    rtn = subprocess.check_call(cmake_cmd, cwd=ns.build_dir, shell=IS_NT)
    rtn = subprocess.check_call(['make'] + ns.make_args, cwd=ns.build_dir)
    if ns.cmd == 'install':
        rtn = subprocess.check_call(['make', 'install'], cwd=ns.build_dir)


def final_message(success=True):
    if success:
        return
    msg = ("\n\nIf you are having issues building pyne, please report your problem "
           "to pyne-dev@googlegroups.com or look for help at http://pyne.io\n\n"
           )
    print('\n' + '-' * 20 + msg + '-' * 20)


def main_safe(ns):
    success = False
    try:
        main_body(ns)
        success = True
    finally:
        final_message(success)


def main():
    ns = parse_args()
    main_safe(ns)
    if ns.bootstrap:
        ensure_nuc_data()
        main_safe(ns)
    binpath = ns.prefix + '/bin'
    msg = ("\nNOTE: If you have not done so already, please be sure that your "
           "PATH has been appropriately set to the install prefix of pyne. "
           "For this install of pyne you may add the following lines to your "
           "'~/.bashrc' file or equivalent:\n\n"
           "  # PyNE Environment Settings\n\n"
           '  export PATH="{binpath}:${{PATH}}"'
           ).format(binpath=binpath)
    print(msg, file=sys.stderr)


if __name__ == "__main__":
    main()
