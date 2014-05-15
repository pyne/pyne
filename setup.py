#!/usr/bin/env python
 
import os
import sys
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
        raise ValueError("numpy version too low! {0} (have) < 1.8.0 (min)".format(v))

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
    hdf5opt = [o.split('=')[1] for o in distutils_args if o.startswith('--hdf5=')]
    if 0 < len(hdf5opt):
        os.environ['HDF5_ROOT'] = hdf5opt[0]  # Expose to CMake
        distutils_args = [o for o in distutils_args if not o.startswith('--hdf5=')]
    # Change egg-base entry to absolute path so it behaves as expected
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--egg-base')
    res, distutils_args = parser.parse_known_args(distutils_args)
    if res.egg_base is not None:
        local_path = os.path.dirname(os.path.abspath(sys.argv[0]))
        distutils_args.append('--egg-base='+os.path.join(local_path, res.egg_base))
    return distutils_args, cmake, make


def main_body():
    assert_dep_versions()
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
        cmake_cmd += ['-DPYTHON_EXECUTABLE=' + sys.executable, ]
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
        rtn = subprocess.check_call(cmake_cmd, cwd='build', shell=(os.name=='nt'))
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

if __name__ == "__main__":
    main()
