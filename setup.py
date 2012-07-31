#!/usr/bin/env python
 
import os
import sys
import glob
import json
from copy import deepcopy

from distutils.core import setup, run_setup
from distutils import sysconfig
from distutils.ccompiler import CCompiler
from distutils.extension import Extension
from distutils.util import get_platform
from distutils.file_util import copy_file, move_file
from distutils.dir_util import mkpath, remove_tree
from distutils.sysconfig import get_python_version, get_config_vars, get_python_lib
from Cython.Distutils import build_ext


import numpy as np

INFO = {
    'version': '0.1',
    }

SITE_PACKAGES = get_python_lib()

_local_subsititues = {'darwin': 'Library'}
HOME = os.environ['HOME'] if os.name != 'nt' else os.environ['UserProfile']
PYNE_DIR = os.path.join(HOME, 
                        _local_subsititues.get(sys.platform, '.local'),
                        'pyne')
                        
HDF5_DIR = os.environ.get('HDF5_DIR', '')
args = sys.argv[:]
for arg in args:
    if arg.find('--hdf5=') == 0:
        HDF5_DIR = os.path.expanduser(arg.split('=')[1])
        sys.argv.remove(arg)


def make_metadata(path):
    """Build a metadata file."""
    md = {"HDF5_DIR": HDF5_DIR,
          }
    md.update(INFO)

    # write the metadata file
    with open(path, 'w') as f:
        json.dump(md, f, indent=2)

    return md


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

###########################################
### Set compiler options for extensions ###
###########################################
pyt_dir = os.path.join('pyne')
cpp_dir = os.path.join('cpp')
dat_dir = os.path.join('data')


# Get numpy include dir
numpy_include = np.get_include()

# HDF5 stuff
#posix_hdf5_libs = ["z", "m", "hdf5", "hdf5_hl",]
posix_hdf5_libs = ["hdf5", "hdf5_hl",]
#nt_hdf5_libs = ["szip", "zlib1", "hdf5dll", "hdf5_hldll",]
nt_hdf5_libs = ["hdf5dll", "hdf5_hldll",]
nt_hdf5_extra_compile_args = ["/EHsc"]
nt_hdf5_macros = [("_WIN32_MSVC", None), ("_HDF5USEDLL_", None),]


###############################
### Platform specific setup ###
###############################
def darwin_linker_paths():
    paths = [os.path.join(PYNE_DIR, 'lib')]
    vars = ['LD_LIBRARY_PATH', 'DYLD_FALLBACK_LIBRARY_PATH', 'DYLD_LIBRARY_PATH', 'LIBRARY_PATH']
    for v in vars:
        curvar = os.getenv(v, '')
        varpaths = paths + ([] if 0 == len(curvar) else [curvar])
        os.environ[v] = ":".join(varpaths)


def darwin_build_ext_decorator(f):
    def new_build_ext(self, ext):
        rtn = f(self, ext)
        if not ext.name.split('.')[-1].startswith('lib'):
            return rtn
        
        libpath = os.path.join(PYNE_DIR, 'lib')
        if not os.path.exists(libpath):
            mkpath(libpath)
        
        copy_file(self.get_ext_fullpath(ext.name), libpath)
        return rtn
    return new_build_ext


def darwin_setup():
    darwin_linker_paths()
    build_ext.build_extension = darwin_build_ext_decorator(build_ext.build_extension)


def win32_finalize_opts_decorator(f):
    def replace(lst, val, rep):
        if val not in lst:
            return lst
        ind = lst.index(val)
        del lst[ind]
        for r in rep[::-1]:
            lst.insert(ind, r)
        return lst

    def posix_like_ext(ext):
        replace(ext.extra_compile_args, "__COMPILER__", [])
        replace(ext.define_macros, "__COMPILER__", [])
        replace(ext.libraries, "__USE_HDF5__", nt_hdf5_libs)
        replace(ext.extra_compile_args, "__USE_HDF5__", [])
        replace(ext.define_macros, "__USE_HDF5__", [])

    def nt_like_ext(ext):
        replace(ext.extra_compile_args, "__COMPILER__", ["/EHsc"])
        replace(ext.define_macros, "__COMPILER__", [("_WIN32_MSVC", None)])
        replace(ext.libraries, "__USE_HDF5__", 
                    ["/DEFAULTLIB:" + lib + ".lib" for lib in nt_hdf5_libs])
        replace(ext.extra_compile_args, "__USE_HDF5__", nt_hdf5_extra_compile_args)
        replace(ext.define_macros, "__USE_HDF5__", nt_hdf5_macros)
    
    update_ext = {'mingw32': posix_like_ext, 
                  'cygwin': posix_like_ext, 
                  'msvc': nt_like_ext, 
                 }
    
    def new_finalize_opts(self):
        rtn = f(self)
        comp = self.compiler
        for ext in self.extensions:
            update_ext[comp](ext)
        if sysconfig.get_config_var('CC') is None:
            if comp in ['mingw32', 'cygwin']:
                # Hack to get compiler to be recognized
                sysconfig._config_vars['CC'] = 'gcc'
        return rtn
    return new_finalize_opts


def win32_build_ext_decorator(f):
    def new_build_ext(self, ext):
        dll_name = ext.name.split('.')[-1]
        islib = dll_name.startswith('lib')
        config_vars = get_config_vars()
        if islib:
            config_vars['SO'] = '.dll'
        rtn = f(self, ext)
        if islib:
            config_vars['SO'] = '.pyd'
        return rtn
    return new_build_ext


def win32_get_exp_sym_decorator(f):
    def new_get_exp_sym(self, ext):
        rtn = f(self, ext)
        dll_name = ext.name.split('.')[-1]
        islib = dll_name.startswith('lib')
        initfunc_name = 'init' + dll_name 
        if islib and initfunc_name in rtn:
            rtn.remove(initfunc_name)
        return rtn
    return new_get_exp_sym


def win32_exec_decorator(f):
    def new_exec(self, func, args, msg=None, level=1):
        if 2 == len(args) and args[0].endswith('.def') and not args[0].startswith('lib'):
            filename, contents = args
            contents = [c for c in contents if not c.startswith('LIBRARY ')]
            args = (filename, contents)
        print "Args = ", args
        rtn = f(self, func, args, msg, level)
        return rtn
    return new_exec


def win32_setup():
    build_ext.finalize_options = win32_finalize_opts_decorator(build_ext.finalize_options)
    build_ext.build_extension = win32_build_ext_decorator(build_ext.build_extension)
    build_ext.get_export_symbols = win32_get_exp_sym_decorator(build_ext.get_export_symbols)
    #build_ext.execute = win32_exec_decorator(build_ext.execute)
    CCompiler.execute = win32_exec_decorator(CCompiler.execute)
    
platform_setup = {'darwin': darwin_setup, 'win32': win32_setup}



##########################
### Exetension Creator ###
##########################

def cpp_ext(name, sources, libs=None, use_hdf5=False):
    """Helper function for setting up extension dictionary.

    Parameters
    ----------
    name : str
        Module name
    sources : list of str
        Files to compile
    libs : list of str
        Additional files to link against
    use_hdf5 : bool
        Link against hdf5?
    """
    ext = {'name': name}

    ext['sources'] = [os.path.join(cpp_dir, s) for s in sources if s.endswith('cpp')] + \
                     [os.path.join(pyt_dir, s) for s in sources if s.endswith('pyx')] + \
                     [s for s in sources if not any([s.endswith(suf) for suf in ['cpp', 'pyx']])]

    ext["libraries"] = []
    ext['include_dirs'] = [pyt_dir, cpp_dir, numpy_include]
    if 0 < len(HDF5_DIR):
        ext['include_dirs'].append(os.path.join(HDF5_DIR, 'include'))
    ext['language'] = "c++"

    # may need to be more general
    ext['library_dirs'] = ['build/lib/pyne/lib',
                           'build/lib.{0}-{1}/pyne/lib'.format(get_platform(), get_python_version()),
                           ]
    if os.name == 'nt':
        ext['library_dirs'] += [SITE_PACKAGES, os.path.join(HDF5_DIR, 'dll'),]
    ext['library_dirs'].append(os.path.join(HDF5_DIR, 'lib'))
                           

    # perfectly general, thanks to dynamic runtime linking of $ORIGIN
    #ext['runtime_library_dirs'] = ['${ORIGIN}/lib', '${ORIGIN}']
    ext['runtime_library_dirs'] = ['${ORIGIN}/lib', '${ORIGIN}', '${ORIGIN}/.', 
                                   '${ORIGIN}/../lib', '${ORIGIN}/..',]
    #ext['runtime_library_dirs'] = ['${ORIGIN}/lib', '${ORIGIN}', '${ORIGIN}/.'] + \
    #                              [os.path.abspath(p) for p in ext['library_dirs']] + \
    #                              [os.path.abspath(p + '/pyne/lib') for p in sys.path] + \
    #                              [os.path.abspath(p + '/pyne') for p in sys.path] + \
    #                              [os.path.abspath(p) for p in sys.path]

    if sys.platform == 'linux2':
        #ext["extra_compile_args"] = ["-Wno-strict-prototypes"]
        ext["undef_macros"] = ["NDEBUG"]
        if use_hdf5:
            ext["libraries"] += posix_hdf5_libs
        if libs is not None:
            ext["libraries"] += libs
    elif sys.platform == 'darwin':
        ext["undef_macros"] = ["NDEBUG"]
        if use_hdf5:
            ext["libraries"] += posix_hdf5_libs
        if libs is not None:
            ext["libraries"] += libs
        config_vars = get_config_vars()
        #config_vars['SO'] = '.dylib'
        config_vars['LDSHARED'] = config_vars['LDSHARED'].replace('-bundle', '-Wl,-x') 

        ext['library_dirs'] = []
        ext['runtime_library_dirs'] = []
        ext["extra_compile_args"] = ["-dynamiclib",
                                     "-undefined", "dynamic_lookup", 
                                     '-shared',
                                     ]
        ext["extra_link_args"] = ["-dynamiclib", 
                                  "-undefined", "dynamic_lookup", 
                                  '-shared',
                                  "-install_name" , os.path.join(PYNE_DIR, 'lib', name.split('.')[-1] + config_vars['SO']),
                                  ]
    elif sys.platform == 'win32':
        ext["extra_compile_args"] = ["__COMPILER__"]
        ext["define_macros"] = ["__COMPILER__"]

        if use_hdf5:
            ext["libraries"].append("__USE_HDF5__")
            ext["extra_compile_args"].append("__USE_HDF5__")
            ext["define_macros"].append("__USE_HDF5__")

        if libs is not None:
            ext["libraries"] += libs
    elif sys.platform == 'cygwin':
        pass

    return ext


#
# For extensions
# 
exts = []

# Pure C/C++ share libraries
# pyne lib
exts.append(cpp_ext("pyne.lib.libpyne", ['pyne.cpp']))

# nucname
exts.append(cpp_ext("pyne.lib.libpyne_nucname", ['nucname.cpp'], ['pyne']))

# data
exts.append(cpp_ext("pyne.lib.libpyne_data", ['data.cpp'], ['pyne', 'pyne_nucname'], True))

# material
exts.append(cpp_ext("pyne.lib.libpyne_material", ['material.cpp'], ['pyne', 'pyne_nucname', 'pyne_data'], True))

# enrichment
exts.append(cpp_ext("pyne.lib.libpyne_enrichment", ['enrichment.cpp'], ['pyne_material',]))


# Python extension modules
# STL converters
#exts.append(cpp_ext("pyne.stlconverters", ['stlconverters.pyx']))

# pyne_config
exts.append(cpp_ext("pyne.pyne_config", ['pyne_config.pyx'], ['pyne']))


# _utils
exts.append(cpp_ext("pyne._utils", ['_utils.pyx']))

# nucname
exts.append(cpp_ext("pyne.nucname", ['nucname.pyx'], ['pyne', 'pyne_nucname']))

# data
exts.append(cpp_ext("pyne.data", ['data.pyx'], ['pyne', 'pyne_nucname', 'pyne_data'], True))

# material
exts.append(cpp_ext("pyne.material", ['material.pyx'], ['pyne', 'pyne_nucname', 'pyne_data', 'pyne_material'], True))

# ace
exts.append(cpp_ext("pyne.ace", ['ace.pyx'], ['pyne', 'pyne_nucname']))

# xs.models
exts.append(cpp_ext("pyne.xs.models", ['xs/models.pyx'], ['pyne', 'pyne_nucname']))

# material
exts.append(cpp_ext("pyne.enrichment", ['enrichment.pyx'], ['pyne_enrichment', 'pyne_material']))


##########################
### Setup Package Data ###
##########################
packages = ['pyne', 'pyne.lib', 'pyne.dbgen', 'pyne.xs']

pack_dir = {'pyne': 'pyne', 'pyne.dbgen': 'pyne/dbgen', 'pyne.xs': 'pyne/xs'}

pack_data = {'pyne': ['includes/*.h', 'includes/pyne/*.pxd', '*.json'],
             'pyne.dbgen': ['*.html'],
            }

ext_modules=[Extension(**ext) for ext in exts]

# Compiler directives
compiler_directives = {'embedsignature': False}
for e in ext_modules:
    e.pyrex_directives = compiler_directives


# Utility scripts
scripts=['scripts/nuc_data_make']
if os.name == 'nt':
    scripts = [s + '.bat' for s in scripts]


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


###################
### Call setup! ###
###################
def pyne_setup():
    print pyne_logo

    # clean includes dir and recopy files over
    if os.path.exists('pyne/includes'):
        remove_tree('pyne/includes')

    mkpath('pyne/includes')
    for header in glob.glob('cpp/*.h'):
        copy_file(header, 'pyne/includes')

    mkpath('pyne/includes/pyne')
    for header in glob.glob('pyne/*.pxd'):
        copy_file(header, 'pyne/includes/pyne')

    # Create metadata file
    mdpath = "pyne/metadata.json"
    md = make_metadata(mdpath)

    # Platform specific setup
    platform_setup.get(sys.platform, lambda: None)()
    
    setup_kwargs = {
        "name": "pyne",
        "version": INFO['version'],
        "description": 'Python for Nuclear Engineering',
        "author": 'PyNE Development Team',
        "author_email": 'scopatz@gmail.com',
        "url": 'http://pyne.github.com/',
        "packages": packages,
        "package_dir": pack_dir,
        "package_data": pack_data,
        "cmdclass": {'build_ext': build_ext}, 
        "ext_modules": ext_modules,
        "scripts": scripts, 
        }
        
    # call setup
    setup_success= False
    try:
        rtn = setup(**setup_kwargs)
        setup_success= True
    finally:
        cleanup(mdpath)
        final_message(setup_success, metadata=md)


if __name__ == "__main__":
    pyne_setup()
