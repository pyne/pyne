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

# Configuration
BASE_DIR = Path(__file__).resolve().parent.parent

CODE_EXTS = {".c", ".cpp", ".cxx", ".h", ".hpp", ".hxx"}
CODE_EXTS |= {ext.upper() for ext in CODE_EXTS}
SOURCE_EXTS = {".c", ".cpp", ".cxx"} | {ext.upper() for ext in [".c", ".cpp", ".cxx"]}
HEADER_EXTS = {".h", ".hpp", ".hxx"} | {ext.upper() for ext in [".h", ".hpp", ".hxx"]}

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
        except (subprocess.CalledProcessError, FileNotFoundError):
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
                "Hint: Ensure your repo has tags or fallback to .git_archival.txt."
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



# Amalgamation Logic
class AmalgamatedFile:
    def __init__(self, output_path):
        self.path = output_path
        self._blocks = []
        self._filenames = []

    def append_line(self, line):
        if not line.endswith("\n"):
            line += "\n"
        self._blocks.append(line)

    def append_file(self, filename, comment_out=None):
        _, ext = os.path.splitext(filename)
        if comment_out is None:
            comment_out = ext not in CODE_EXTS
        try:
            with open(filename, "r", encoding="utf-8") as f:
                content = f.read()
        except Exception as e:
            print(f"[!] Warning: Skipping unreadable file: {filename}\n    Reason: {e}")
            return

        header = f"//\n// Begin: {filename}\n//\n"
        footer = f"//\n// End: {filename}\n//\n\n"
        if comment_out:
            content = "// " + content.replace("\n", "\n// ")

        self._blocks.append(header + content + "\n" + footer)
        self._filenames.append(filename)

    def prepend_file_listing(self):
        listing = "// Amalgamated from the following files:\n"
        for f in self._filenames:
            listing += f"//   {f}\n"
        self._blocks.insert(0, listing + "\n")

    def write(self):
        self.prepend_file_listing()
        final = "".join(self._blocks)
        output_dir = os.path.dirname(self.path)
        if output_dir:
            os.makedirs(output_dir, exist_ok=True)
        with open(self.path, "w", encoding="utf-8") as f:
            f.write(final)
        print(f"[✓] Written: {self.path}")


# Main Entry
def main():
    parser = ArgumentParser(description="Amalgamate PyNE C++ code.")
    parser.add_argument(
        "-s",
        dest="source_name",
        default="pyne.cpp",
        help="Output C++ source file name.",
    )
    parser.add_argument(
        "-i", dest="header_name", default="pyne.h", help="Output header file name."
    )
    parser.add_argument(
        "-f",
        dest="files",
        nargs="+",
        default=DEFAULT_FILES,
        help="Input files to amalgamate.",
    )
    parser.add_argument(
        "-o",
        dest="output_dir",
        default=".",
        help="Output directory for generated files.",
    )

    args = parser.parse_args()
    output_dir = os.path.abspath(args.output_dir)

    version = get_version()
    create_version_header(version)

    header_path = os.path.join(output_dir, args.header_name)
    source_path = os.path.join(output_dir, args.source_name)

    # Header and Source Generation
    header = AmalgamatedFile(header_path)
    header.append_line("// Amalgamated PyNE header - http://pyne.io/")
    header.append_line("#ifndef PYNE_AMALGAMATED_HEADER")
    header.append_line("#define PYNE_AMALGAMATED_HEADER\n")
    header.append_line("#define PYNE_IS_AMALGAMATED\n")

    source = AmalgamatedFile(source_path)
    source.append_line("// Amalgamated PyNE source - http://pyne.io/")
    rel_header = os.path.relpath(header_path, os.path.dirname(source_path))
    source.append_line(f'#include "{rel_header}"\n')

    # Process all files in a single loop
    for file in args.files:
        _, ext = os.path.splitext(file)
        if ext in HEADER_EXTS:
            header.append_file(file)
        elif ext in SOURCE_EXTS:
            source.append_file(file)
        else:
            # Assume non-code files (like license.txt) go in both
            header.append_file(file)
            source.append_file(file)

    header.append_line("#endif  // PYNE_AMALGAMATED_HEADER")
    header.write()
    source.write()
    os.remove(os.path.join(BASE_DIR, "version.h"))


if __name__ == "__main__":
    main()
