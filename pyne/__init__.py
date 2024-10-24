import os
import sys
import glob
from warnings import warn


if os.name == "nt":
    p = os.environ["PATH"].split(";")
    lib = os.path.join(os.path.split(__file__)[0], "lib")
    os.environ["PATH"] = ";".join([lib] + p)

try:
    from .pyne_version import PYNE_VERSION

    __version__ = PYNE_VERSION
    from .pyne_config import *

except ImportError:
    msg = (
        "Error importing PyNE: you should not try to import PyNE from "
        "its source directory; please exit the PyNE source tree, and relaunch "
        "your python interpreter from there."
    )
    warn(msg, Warning)
    raise

def get_path(subdir, pattern=""):
    """Helper function to return paths that match a given pattern within a subdirectory."""
    path = os.path.normpath(os.path.join(__path__[0], '..', '..', '..', '..', subdir))
    return glob.glob(os.path.join(path, pattern)) if os.path.exists(path) else []

def get_include_path():
    """Return the include directory path for PyNE headers."""
    return get_path("include")

def get_core_libraries():
    """Return library paths and library directory paths."""
    #lib_paths = [os.path.join(sys.prefix, subdir) for subdir in ["lib", "lib64"]]
    # TODO: Temporary fix for old pyne
    lib_paths = get_path("lib")
    lib = get_path("lib", "*pyne*")
    return lib, lib_paths

def get_extra_libraries():
    """List all the extra libraries of PyNE."""
    libs_path = os.path.join(__path__[0], ".dylibs") if sys.platform == "darwin" else os.path.join(__path__[0], "..", "pyne.libs")
    return (glob.glob(os.path.join(libs_path, "*")), libs_path) if os.path.exists(libs_path) else ([], [])

# Setup variables
include_path = get_include_path()
lib, lib_path = get_core_libraries()
extra_lib, extra_lib_path = get_extra_libraries()

# Export variables for easy access
__all__ = ["include_path", "lib", "lib_path", "extra_lib", "extra_lib_path"]