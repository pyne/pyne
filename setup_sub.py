#!/usr/bin/env python
"""Welcome to PyNE's setup.py sub script."""
from __future__ import print_function

import io
import os
import re
import sys
import imp
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
from distutils.core import setup

from pyne.pyne_version import PYNE_VERSION


IS_NT = os.name == "nt"


def main():
    scripts = [os.path.join("scripts", f) for f in os.listdir("scripts")]
    scripts = [
        s
        for s in scripts
        if (os.name == "nt" and s.endswith(".bat"))
        or (os.name != "nt" and not s.endswith(".bat"))
    ]
    packages = [
        "pyne",
        "pyne.dbgen",
        "pyne.apigen",
        "pyne.xs",
        "pyne.transmute",
        "pyne.gui",
        "pyne.cli",
        "pyne.fortranformat",
    ]
    pack_dir = {
        "pyne": "pyne",
        "pyne.xs": "pyne/xs",
        "pyne.gui": "pyne/gui",
        "pyne.cli": "pyne/cli",
        "pyne.dbgen": "pyne/dbgen",
        "pyne.apigen": "pyne/apigen",
        "pyne.transmute": "pyne/transmute",
        "pyne.fortranformat": "pyne/fortranformat",
    }
    extpttn = ["*.dll", "*.so", "*.dylib", "*.pyd", "*.pyo"]
    pack_data = {
        "lib": extpttn,
        "pyne": [
            "*.pxd",
            #'include/*.h', 'include/*.pxi', 'include/*/*.h',
            #'include/*/*/*.h', 'include/*/*/*/*.h',
            "*.json",
            "*.inp",
            #'_includes/*.txt', '_includes/*.pxd', '_includes/*/*',
            #'_includes/*/*/*'
        ]
        + extpttn,
        "pyne.xs": ["*.pxd"] + extpttn,
        "pyne.gui": ["*.pyw"],
        "pyne.dbgen": ["*.html", "*.csv", "abundances.txt", "mass.mas16", "*.dat"],
    }
    setup_kwargs = {
        "name": "pyne",
        "version": PYNE_VERSION,
        "description": "The Nuclear Engineering Toolkit",
        "author": "PyNE Development Team",
        "author_email": "pyne-dev@googlegroups.com",
        "url": "http://pyne.github.com/",
        "packages": packages,
        "package_dir": pack_dir,
        "package_data": pack_data,
        "scripts": scripts,
    }
    rtn = setup(**setup_kwargs)


if __name__ == "__main__":
    main()
