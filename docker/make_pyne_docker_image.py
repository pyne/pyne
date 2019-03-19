#! /usr/bin/env python
from __future__ import print_function, unicode_literals
import os
import subprocess
import argparse as ap


def consistancy_check(args):
    """Checks arguments consistencies: for example, if DAGMC is request, MOAB
    needs to be installed

    Parameters
    ----------
    args : list of arguments 

    """
    args.moab = args.moab or args.dagmc or args.pymoab or args.all
    args.pymoab = args.pymoab or args.all
    args.dagmc = args.dagmc or args.all

    return args


def build_name(args):
    """Build the container name, if not specified by the user. Container name
    will be "pyne/ubuntu_18.04{_moab/_dagmc{_pymoab}}{_pyne-deps/}"

    Parameters
    ----------
    args : list of arguments 

    """
    if args.name != '':
        return args.name

    name = 'pyne/ubuntu_18.04'
    name += '_py'+str(args.py_version)
    if args.moab and not args.dagmc and not args.pymoab:
        name += "_moab"
    else:
        if args.dagmc:
            name += '_dagmc'
        if args.pymoab:
            name += '_pymoab'

    if args.deps:
        name += '_pyne-deps'
    return name


def build_docker(args):
    """Forms the docker command line and executes it. IF the push flag is
    provided also pushes it to DockerHub which requires PyNE admin right on
    DockerHub(Docker login is assumed...)

    Parameters
    ----------
    args : list of arguments 

    """
    dockerfile = ['-f', 'ubuntu_18.04-dev.dockerfile']
    tag = build_name(args)
    tag_flag = ['-t', tag]
    docker_args = []
    if args.moab:
        docker_args += ["--build-arg", "build_moab=YES"]
    if args.dagmc:
        docker_args += ["--build-arg", "build_dagmc=YES"]
    if args.pymoab:
        docker_args += ["--build-arg", "enable_pymoab=YES"]
    if args.deps:
        docker_args += ["--build-arg", "build_pyne=NO"]
    if args.py_version:
        if args.py_version == 2:
            docker_args += ["--build-arg", "py_version=2.7"]
        elif args.py_version == 3:
            docker_args += ["--build-arg", "py_version=3.6"]
        else:
            print("Can only deal with python 2 or 3")
            return

    rtn = subprocess.check_call(
        ["docker",  "build"] + tag_flag + dockerfile + docker_args + ["."], shell=(os.name == 'nt'))

    if args.push:
        rtn = subprocess.check_call(
            ["docker", "push",  tag], shell=(os.name == 'nt'))


def main():
    """Parse the different arguments and call the proper methods.
    """
    description = 'Build a docker image for PyNE'
    parser = ap.ArgumentParser(description=description)

    py_version = 'Require a specific python version'
    parser.add_argument('--py_version', type=int, help=py_version)

    moab = 'Build and install MOAB'
    parser.add_argument('--moab', help=moab,
                        action='store_true', default=False)

    dagmc = 'Build and install DAGMC'
    parser.add_argument('--dagmc', help=dagmc,
                        action='store_true', default=False)

    pymoab = 'Enable pymoab'
    parser.add_argument('--pymoab', help=pymoab,
                        action='store_true', default=False)

    all_deps = 'Add all dependencies'
    parser.add_argument('--all', '-a', '-all', help=all_deps,
                        action='store_true', default=False)

    deps = 'Depdendencies only'
    parser.add_argument('--deps', help=deps,
                        action='store_true', default=False)

    name = 'Set docker image name'
    parser.add_argument('--name', help=name, default='')

    push = 'Push docker image on dockerhub'
    parser.add_argument('--push', '-p', '-push', help=push,
                        action='store_true', default=False)

    args = parser.parse_args()

    args = consistancy_check(args)

    build_docker(args)


if __name__ == "__main__":
    main()
