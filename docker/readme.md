Docker Script Documentation
===========================

The files in this folder are intended to be mainly used by PyNE maintainers,
they provide an "easy" mechanism to build and push the various containers
required for the multiple configurations tested by the PyNE CI.

This folder contains:

- A custom Dockerfile: "ubuntu_18.04-dev.dockerfile", allowing the user 
to build docker containers with a custom set of
  dependencies required by the CI. The configuration of the Docker containers
  can be specified using the `--build-args <flag>` argument. The different flags are:
    
    - `build_moab=YES` add MOAB in the docker container (default: `NO`)
    - `build_pymoab=YES` install MOAB with pymoab (default: `NO`)
    - `build_dagmc=YES` add DAGMC (default: `NO`)
    - `install_openmc=YES` install OpenMC Python3 API (default: `NO`)
    - `build_pyne=NO` do not build PyNE, only add the PyNE dependencies to the container (default: `YES`)
    - `py_version=X.Y` specify the python version to install (options are 2.7 or
      3.6) (default: `2.7`)


Note, if using multiple build arguments `--build-args` has to be repeated, i.e. :
to request DAGMC dependencies only with Python 2.7, the docker build arguments will be: 
`--build-args build_dagmc=YES --build-args build_pyne=NO --build-args py_version=2.7` 

- A python script, `make_pyne_docker_image.py`, is also present to simplify the usage of the docker file. It
 allows maintainers to build, name, and push the docker container into
  DockerHub. The available options are:
    
    - `--py_version=X` with "X" = major python version (e.g. 2 or 3)
    - `--moab` install MOAB in the container
    - `--dagmc` install DAGMC in the container
    - `--pymoab` install pyMOAB in the container (along with MOAB)
    - `--openmc` install OpenMC Python3 API in the container
    - `--all/-a/-all` install all the optional dependencies (MOAB/pyMOAB/DAGMC) 
    - `--deps` don't install PyNE, only install the dependencies (both required and selected optional dependencies)
    - `--push` push the docker container to DockerHub after building it
      (requires right to push on the PyNE DockerHub account)
    - `--name=` manually set the docker container name (NOT recommended to be used
      with the push option)
The default value will be the one defined in the `Dockerfile` (Python 2.7, no
optionnal dependencies, PyNE built).

The Python script will check for consistency among the required dependencies to ensure
working build. For example, pyMOAB requires MOAB to be install (MOAB will be installed if
the pyMOAB flag is provided, even if the MOAB flag is not), DAGMC requires pyMOAB (and
MOAB).

The Python script will also name the container according to the
convention, allowing the automatic CI to use the pushed docker image. We
**STRONGLY** advise against manually setting the container name when updating the CI
docker containers in the PyNE DockerHuB. 


We recommend maintainers use the python script to generate and push the docker
container in order to form container name correctly and to include the proper
dependency tree.
