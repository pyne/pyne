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
import subprocess
from pathlib import Path
from argparse import ArgumentParser

# Configuration
BASE_DIR = Path(__file__).resolve().parent

SOURCE_EXTS = {".c", ".cpp", ".cxx"} | {ext.upper() for ext in [".c", ".cpp", ".cxx"]}
HEADER_EXTS = {".h", ".hpp", ".hxx"} | {ext.upper() for ext in [".h", ".hpp", ".hxx"]}
CODE_EXTS = SOURCE_EXTS | HEADER_EXTS
LICENSE_FILE = "license.txt"

DEFAULT_FILES = [
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

DEFAULT_FILES = [BASE_DIR / f for f in DEFAULT_FILES]


# Version Handling
def get_version():
    """Retrieves the version string from Git or a fallback file."""

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
        archival = BASE_DIR / ".git_archival.txt"
        if archival.exists():
            with archival.open("r", encoding="utf-8") as f:
                for line in f:
                    if line.startswith("describe-name:"):
                        version = line.split(":", 1)[1].strip()
                        print(f"[✓] Version from .git_archival.txt: {version}")
                        return version
            raise RuntimeError("describe-name not found in .git_archival.txt")
        raise RuntimeError("Not in a Git repo and .git_archival.txt is missing.")


def generate_version_header_content(version):
    """Generates the C++ version header content as a string."""
    return f"""\
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


# Amalgamation Logic
class AmalgamatedFile:
    """Class to handle the amalgamation of files into a single output file."""

    def __init__(self, output_path, amalgamated_headers=None):
        self.path = output_path
        self._blocks = []
        self._filenames = []

    def append_line(self, line):
        """Appends a single line to the amalgamated file."""
        if not line.endswith("\n"):
            line += "\n"
        self._blocks.append(line)

    def append_commented_block(self, content, title):
        """Appends a block of text, formatting it as a C++ comment."""
        header = f"//\n// Begin: {title}\n//\n"
        footer = f"//\n// End: {title}\n//\n\n"
        commented_content = "// " + content.replace("\n", "\n// ")
        self._blocks.append(header + commented_content + "\n" + footer)

    def append_file(self, filename):
        """Appends the content of a file to the amalgamated output."""
        try:
            lines = filename.read_text(encoding="utf-8").splitlines(keepends=True)
        except Exception as e:
            print(f"[!] Warning: Skipping unreadable file: {filename}\n    Reason: {e}")
            return

        content = "".join(lines)
        header = f"//\n// Begin: {filename.name}\n//\n"
        footer = f"//\n// End: {filename.name}\n//\n\n"

        self._blocks.append(header + content + "\n" + footer)
        self._filenames.append(str(filename.relative_to(BASE_DIR)))

    def prepend_file_listing(self):
        """Adds a list of all files included in the amalgamation to the top."""
        listing = "// Amalgamated from the following files:\n"
        listing += "//   " + "\n//   ".join(self._filenames) + "\n"
        self._blocks.insert(0, listing + "\n")

    def write(self):
        """Writes the collected blocks to the output file."""
        self.prepend_file_listing()
        final = "".join(self._blocks)
        output_dir = self.path.parent
        output_dir.mkdir(parents=True, exist_ok=True)
        self.path.write_text(final, encoding="utf-8")
        print(f"[✓] Written: {self.path}")


# Main Entry
def main():
    parser = ArgumentParser(description="Amalgamate PyNE C++ code.")
    parser.add_argument(
        "-s",
        dest="source_name",
        default="pyne.cpp",
        help="Output C++ (.cpp) source file name.",
    )
    parser.add_argument(
        "-i", dest="header_name", default="pyne.h", help="Output header (.h) file name."
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
    output_dir = Path(args.output_dir).resolve()
    input_files = [Path(f) for f in args.files]

    # Get version and generate header
    version = get_version()
    version_header_content = generate_version_header_content(version)

    # Read license content
    license_path = BASE_DIR / LICENSE_FILE
    try:
        license_content = license_path.read_text(encoding="utf-8")
    except FileNotFoundError:
        license_content = f"{LICENSE_FILE} not found."
        print(f"[!] Warning: {license_content}")

    header_path = output_dir / args.header_name
    source_path = output_dir / args.source_name

    # Header File Generation
    header = AmalgamatedFile(header_path)
    header.append_line("// Amalgamated PyNE header - http://pyne.io/")
    header.append_commented_block(license_content, "License")
    header.append_line("#ifndef PYNE_AMALGAMATED_HEADER")
    header.append_line("#define PYNE_AMALGAMATED_HEADER")
    header.append_line("#define PYNE_IS_AMALGAMATED\n")

    # Add version header content
    header.append_line("//\n// Begin: version.h\n//")
    header.append_line(version_header_content)
    header.append_line("//\n// End: version.h\n//\n")

    # Source File Generation
    source = AmalgamatedFile(source_path)
    source.append_line("// Amalgamated PyNE source - http://pyne.io/")
    source.append_commented_block(license_content, "License")
    rel_header_path = header_path.relative_to(source_path.parent).as_posix()
    source.append_line(f'#include "{rel_header_path}"\n\n')

    # Process all user-specified files
    for file_path in input_files:
        if file_path.suffix in HEADER_EXTS:
            header.append_file(file_path)
        elif file_path.suffix in SOURCE_EXTS:
            source.append_file(file_path)
        else:
            print(f"[!] Warning: Skipping file with unknown extension: {file_path}")

    header.append_line("#endif  // PYNE_AMALGAMATED_HEADER")
    header.write()
    source.write()


if __name__ == "__main__":
    main()
