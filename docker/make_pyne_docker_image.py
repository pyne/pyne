#! /usr/bin/env python
from __future__ import print_function, unicode_literals
import os
import sys
import tarfile
import platform
import subprocess
import shutil
import io
import argparse as ap


def absexpanduser(x): return os.path.abspath(os.path.expanduser(x))


def consistancy_check(args):
    args.moab = args.moab or args.dagmc or args.pymoab

    return args


def build_name(args):
    if args.name != '':
        return args.name

    name = 'ubuntu_18.04'

    if args.moab and not args.dagmc and not args.pymoab:
        name += "_moab"
    else:
        if args.dagmc:
            name += '_dagmc'
        if args.pymoab:
            name += '_pymoab'

    if args.deps:
        name += '_pyne_deps'
    return name


def build_docker(args):

    dockerfile = ['-f', 'ubuntu_18.04-dev.dockerfile']

    tag = ['-t', build_name(args) ]
    docker_args = []
    if args.moab:
        docker_args += ["--build-arg", "build_moab=YES"]
    if args.dagmc:
        docker_args += ["--build-arg", "build_dagmc=YES"]
    if args.pymoab:
        docker_args += ["--build-arg", "enable_pymoab=YES"]
    if args.deps:
        docker_args += ["--build-arg", "build_pyne=NO"]

    rtn = subprocess.check_call(["docker",  "build" ] + tag + dockerfile + docker_args + ["."], shell=(os.name == 'nt'))

def main():
    localdir = absexpanduser('~/.local')
    description = "Build a docker image for PyNE"
    parser = ap.ArgumentParser(description=description)

    moab = 'Build and install MOAB'
    parser.add_argument('--moab', help=moab,
                        action='store_true', default=False)

    dagmc = 'Build and install DAGMC'
    parser.add_argument('--dagmc', help=dagmc,
                        action='store_true', default=False)

    pymoab = 'Enable pymoab'
    parser.add_argument('--pymoab', help=pymoab,
                        action='store_true', default=False)

    deps = 'Depdendencies only'
    parser.add_argument('--deps', help=deps,
                        action='store_true', default=False)

    name = "Set docker imgae name"
    parser.add_argument('--name', help=name, default='')

    args = parser.parse_args()

    args = consistancy_check(args)

    build_docker(args)


if __name__ == "__main__":
    main()
