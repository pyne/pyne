#!/usr/bin/env python
"""
Amalgamate PyNE's C++ library sources into a single header and source file.

This script consolidates selected C++ components of the PyNE library into:
    - pyne.h   : Self-contained header
    - pyne.cpp : Self-contained source

Inspired by the JsonCpp amalgamation tool:
    http://svn.code.sf.net/p/jsoncpp/code/trunk/jsoncpp/amalgamate.py

Usage:
    python amalgamate.py [-s OUTPUT.cpp] [-i OUTPUT.h] [-f file1.h file2.cpp ...] [-o OUTPUT_DIR]
"""

from __future__ import print_function, unicode_literals
import os
import subprocess
from pathlib import Path
from argparse import ArgumentParser
import sys
import io

# Configuration
BASE_DIR = Path(__file__).resolve().parent

SOURCE_EXTS = {".c", ".cpp", ".cxx"} | {ext.upper() for ext in [".c", ".cpp", ".cxx"]}
HEADER_EXTS = {".h", ".hpp", ".hxx"} | {ext.upper() for ext in [".h", ".hpp", ".hxx"]}
CODE_EXTS = SOURCE_EXTS | HEADER_EXTS

DEFAULT_FILES = [
    "license.txt",
    "version.h",
    "src/utils.h",
    "src/utils.cpp",
    "src/extra_types.h",
    "src/h5wrap.h",
    "src/state_map.cpp",
    "src/nucname.h",
    "src/nucname.cpp",
    "src/rxname.h",
    "src/rxname.cpp",
    "src/_atomic_data.h",
    "src/_atomic_data.cpp",
    "src/data.h",
    "src/data.cpp",
    #'src/dagmc_bridge.cpp',
    #'src/dagmc_bridge.h',
    "src/json-forwards.h",
    "src/json.h",
    "src/jsoncpp.cpp",
    "src/jsoncustomwriter.h",
    "src/jsoncustomwriter.cpp",
    "src/material.h",
    "src/material.cpp",
    "src/material_library.h",
    "src/material_library.cpp",
    "src/enrichment_cascade.h",
    "src/enrichment_cascade.cpp",
    "src/enrichment.h",
    "src/enrichment.cpp",
    "src/enrichment_symbolic.h",
    "src/enrichment_symbolic20.cpp",
    "src/_decay.h",
    "src/_decay.cpp",
]

DEFAULT_FILES = [os.path.join(BASE_DIR, f) for f in DEFAULT_FILES]


# Version Handling
def get_version():
    def in_git_repo():
        try:
            subprocess.check_output(
                ["git", "rev-parse", "--is-inside-work-tree"], stderr=subprocess.STDOUT
            )
            return True
        except (subprocess.CalledProcessError, FileNotFoundError) as e:
            if isinstance(e, subprocess.CalledProcessError) and b"fatal" in e.output:
                raise RuntimeError(
                    "Git command failed. Ensure you are in a valid Git repository.\n"
                    f"Error: {e.output.decode('utf-8').strip()}"
                )
            return False

    if in_git_repo():
        try:
            version = (
                subprocess.check_output(
                    ["git", "describe", "--tags"], stderr=subprocess.STDOUT
                )
                .strip()
                .decode("utf-8")
            )
            if not version:
                raise RuntimeError("Empty version string from git.")
            print(f"[✓] Version from git: {version}")
            return version
        except subprocess.CalledProcessError as e:
            raise RuntimeError(
                f"Git describe failed. Output: {e.output.decode('utf-8').strip()}\n"
                "Hint: Ensure your repo has tags.\n"
                "   git fetch --tags\n"
                "If you are using forked repositories, ensure you have the correct upstream set.\n"
                "   git remote add upstream https://github.com/pyne/pyne.git\n"
                "   git fetch upstream --tags"
            )
    else:
        archival = os.path.join(BASE_DIR, ".git_archival.txt")
        if os.path.exists(archival):
            with open(archival, "r", encoding="utf-8") as f:
                for line in f:
                    if line.startswith("describe-name:"):
                        version = line.split(":", 1)[1].strip()
                        print(f"[✓] Version from .git_archival.txt: {version}")
                        return version
            raise RuntimeError("describe-name not found in .git_archival.txt")
        raise RuntimeError("Not in a Git repo and .git_archival.txt is missing.")


def create_version_header(version, file_name="version.h"):
    content = f"""\
#ifndef PYNE_VERSION_HEADER
#define PYNE_VERSION_HEADER

#include <string>

namespace pyne {{
inline std::string pyne_version() {{
    return "{version} (amalgamated)";
}}
}}  // namespace pyne

#endif  // PYNE_VERSION_HEADER
"""
    output_path = os.path.join(BASE_DIR, file_name)
    with open(output_path, "w", encoding="utf-8") as f:
        f.write(content)

class AmalgamatedFile(object):
    def __init__(self, path):
        self.path = path
        self._blocks = []
        self._filenames = []

    def append_line(self, line):
        """Adds some text to the end of the file."""
        if not line.endswith("\n"):
            line += "\n"
        self._blocks.append(line)

    def append_file(self, filename, comment_out=None):
        """Adds a whole file to the end of this one."""
        if comment_out is None:
            _, ext = os.path.splitext(filename)
            comment_out = ext not in CODE_EXTS
        self._blocks.append("//\n// start of {0}\n//\n".format(filename))
        with open(filename, "rt", encoding="utf-8") as f:
            content = f.read()
        if comment_out:
            content = "// " + content.replace("\n", "\n// ")
        self._blocks.append(content)
        self._blocks.append("//\n// end of {0}\n//\n\n\n".format(filename))
        self._filenames.append(filename)

    def prepend_files(self):
        """Adds a file listing to the begining of the almagamted file."""
        s = "// This file is composed of the following original files:\n\n"
        for f in self._filenames:
            s += "//   {0}\n".format(f)
        s += "\n"
        self._blocks.insert(0, s)

    def write(self):
        self.prepend_files()
        if sys.version > "3":
            txt = "".join(self._blocks)
        else:
            txt = "".join([block.decode("utf-8") for block in self._blocks])
        d = os.path.dirname(self.path)
        if len(d) > 0 and not os.path.isdir(d):
            os.makedirs(d)
        with io.open(self.path, "wb") as f:
            f.write(txt.encode("utf-8"))


def main():
    parser = ArgumentParser()
    parser.add_argument(
        "-s",
        dest="source_path",
        action="store",
        default="pyne.cpp",
        help="Output *.cpp source path.",
    )
    parser.add_argument(
        "-i",
        dest="header_path",
        action="store",
        default="pyne.h",
        help="Output header path.",
    )
    parser.add_argument(
        "-f",
        dest="files",
        nargs="+",
        help="Files to amalgamate.",
        default=DEFAULT_FILES,
    )
    ns = parser.parse_args()
    version = get_version()
    create_version_header(version)
    # header file
    hdr = AmalgamatedFile(ns.header_path)
    hdr.append_line("// PyNE amalgated header http://pyne.io/")
    hdr.append_line("#ifndef PYNE_52BMSKGZ3FHG3NQI566D4I2ZLY")
    hdr.append_line("#define PYNE_52BMSKGZ3FHG3NQI566D4I2ZLY")
    hdr.append_line("")
    hdr.append_line("#define PYNE_IS_AMALGAMATED")
    hdr.append_line("")
    for f in ns.files:
        _, ext = os.path.splitext(f)
        if ext in SOURCE_EXTS:
            continue
        hdr.append_file(f)
    hdr.append_line("#endif  // PYNE_52BMSKGZ3FHG3NQI566D4I2ZLY")

    # source file
    src = AmalgamatedFile(ns.source_path)
    src.append_line("// PyNE amalgated source http://pyne.io/")
    src.append_line(
        '#include "{0}"'.format(
            os.path.relpath(ns.header_path, os.path.dirname(ns.source_path))
        )
    )
    src.append_line("")
    for f in ns.files:
        _, ext = os.path.splitext(f)
        if ext in HEADER_EXTS:
            continue
        src.append_file(f)

    # write both
    hdr.write()
    src.write()


if __name__ == "__main__":
    main()
