#!/usr/bin/env python
from __future__ import print_function

import os
import sys
import imp
import subprocess

import numpy as np

import configure

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

def assert_dep_versions():
    assert_np_version()

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
    if btopt:
        cmake += ['-D','CMAKE_BUILD_TYPE:STRING=' + btopt[0]]
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


def main_body():
    assert_dep_versions()
    if not os.path.exists('build'):
        os.mkdir('build')
    sys.argv, cmake_args, make_args = parse_args()
    makefile = os.path.join('build', 'Makefile')
    source_dir = os.path.dirname(__file__)
    if not os.path.exists(makefile):
        if os.name != 'nt':
            rtn = subprocess.call(['which', 'cmake'])
            if rtn != 0:
                sys.exit('CMake is not installed, aborting PyNE build.')
        cmake_cmd = ['cmake'] + cmake_args
        cmake_cmd += ['-D', 'PYTHON_EXECUTABLE:FILEPATH=' + sys.executable]
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
        cmake_cmd.append(source_dir)
        cmake_cmdstr = cmake_cmd if isinstance(cmake_cmd, str) else ' '.join(cmake_cmd)
        print("CMake command is\n", cmake_cmdstr, sep="")
        rtn = subprocess.check_call(cmake_cmd, cwd='build', shell=is_nt)
    rtn = subprocess.check_call(['make'] + make_args, cwd='build')
    cwd = os.getcwd()
    os.chdir('build')
    configure.setup()
    os.chdir(cwd)

def main():
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

if __name__ == "__main__":
    main()
