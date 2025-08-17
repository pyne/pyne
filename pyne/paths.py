import os
import sys
import glob
from . import __path__

PYNE_CORE_BASE_PATH = os.path.join(__path__[0], "core")


def get_paths(subdir, pattern="*", recursive=False):
    """
    Helper function to return paths that match a given pattern within a subdirectory.

    Args:
        subdir (str): The subdirectory within the 'core' directory.
        pattern (str): The pattern to match files or directories.
        recursive (bool): Whether to search recursively in subdirectories.

    Returns:
        list: A list of matched paths.
    """
    search_pattern = (
        os.path.join(PYNE_CORE_BASE_PATH, subdir, "**", pattern)
        if recursive
        else os.path.join(PYNE_CORE_BASE_PATH, subdir, pattern)
    )
    return glob.glob(search_pattern, recursive=recursive)


def get_include_path():
    """Return includes and include path for PyNE headers."""
    include = get_paths("include", "*", recursive=True)
    include_path = get_paths("include", "", recursive=False)
    return include, include_path


def get_core_libraries():
    """Return libraries and library paths for PyNE."""
    lib = [
        lib_file
        for lib in ["lib", "lib64"]
        for lib_file in get_paths(lib, "libpyne*", recursive=True)
    ]
    lib_path = [
        lib_file
        for lib in ["lib", "lib64"]
        for lib_file in get_paths(lib, "", recursive=False)
    ]
    return lib, lib_path


def get_extra_libraries():
    """Return the extra libraries installed by auditwheel or delocate."""
    libs_path = (
        os.path.join(__path__[0], ".dylibs")
        if sys.platform == "darwin"
        else os.path.normpath(os.path.join(__path__[0], "..", "pyne.libs"))
    )
    return (
        (glob.glob(os.path.join(libs_path, "*")), libs_path)
        if os.path.exists(libs_path)
        else ([], [])
    )


# Setup variables
include, include_path = get_include_path()
lib, lib_path = get_core_libraries()
extra_lib, extra_lib_path = get_extra_libraries()

# Export variables for easy access
__all__ = ["include", "include_path", "lib", "lib_path", "extra_lib", "extra_lib_path"]
