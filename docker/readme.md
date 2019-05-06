Docker Script Documentation
===========================

The files in this folder are intended to be mainly used by PyNE maintainers,
they provide an "easy" mechanism to build and push the various containers
required for the multiple configurations tested by the PyNE CI.

This folder contains:

- A custom Dockerfile: "ubuntu_18.04-dev.dockerfile", allowing the user 
to build docker containers with a custom set of
  dependencies required by the CI, the configuration of the Docker containers
  can be specified using the "--build-args" flags. The different flags are:
    
    - "build_moab=YES" add MOAB in the docker container
    - "build_pymoab=YES" install MOAB with pymoab
    - "build_dagmc=YES" add DAGMC 
    - "build_pyne=NO" only add the PyNE dependencies to the container
    - "py_version=X.Y" specify the python version to install (option are 2.7 or
      3.6)

- A python script, "make_pyne_docker_image.py", is also present to simplify the usage of the docker file, it
  will allows maintainers to build, name and push the docker container into
  DockerHub. The different option are:
    
    - "--py_version=X" with "X" = 2 or 3, specify the python version
    - "--moab" install MOAB in the container
    - "--dagmc" install DAGMC in the container
    - "--pymoab" install pyMOAB in the container (along with MOAB)
    - "--all/-a/-all" install all the optional dependencies (MOAB/pyMOAB/DAGMC)
    - "--deps" don't install PyNE only install the selected optional PyNE dependencies (along with the required one)
    - "--push" push the docker container to DockerHub after building it
      (requires right to push on the PyNE DockerHub account)
    - "--name=" manually set the docker container name (NOT recommended with used
      with the push option)


The Python script will check consistency in the required dependency to ensure
working build: pyMOAB requires MOAB to be install (MOAB will be install if
pyMOAB flag is provided, even if MOAB flag is not), DAGMC requires pyMOAB (and
MOAB).

The Python script will also name container according to some
convention allowing the automatic CI to use the pushed docker image, One
**STRONGLY** advices against manually set the container name when updating the CI
docker containers in the PyNE DockerHuB. 


We recommend maintainers use the python script to generate and push the docker
container in order to form container name correctly, and include the proper
dependency tree.
