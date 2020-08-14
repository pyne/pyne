.. _ubuntu:

=================================
Ubuntu build script
=================================

Scripts for installing PyNE and all its dependencies from scratch on Ubuntu
14.04 - 15.04 are found `here
<https://github.com/pyne/install_scripts/>`_.

In this repository, you will find both bash scripts and dockerfiles for Ubuntu.
These scripts will install PyNE and its dependencies to either your device
or an `image <https://docs.docker.com/get-started/#images-and-containers>`_. Docker
allows users to build these images and to operate them separate from the rest of
your device.

You can download Docker `here <https://docs.docker.com/get-docker/>`_.


-------------------
Bash Scripts (*.sh)
-------------------

The script used to install PyNE should correspond
to the user's operating system and version.

The intention of these
scripts is that PyNE can be ready to use on a clean install of any of
the supported operating systems. Furthermore, the user should choose either
to build a stable version of PyNE or the current develop
branch of PyNE by supplying a secondary argument. Choosing between dev and stable is
a matter of preference, if you intend to contribute to PyNE it is recommended that
you use dev.

Example for installing the most recent stable branch on Ubuntu 16.04::

    ./ubuntu_16.04.sh stable
    
Example for installing the develop branch on Mint 18.01::
    
    ./mint_18.01.sh dev
    

----------------------------
Docker Builds (*.dockerfile)
----------------------------

The user may choose either
to build a stable version of PyNE ("-stable") or the current develop
branch of PyNE ("-dev") for each.

Example for building a docker image of the latest stable branch of PyNE based on
Ubuntu 16.04::

    docker build -f ubuntu_16.04-stable.dockerfile -t pyne-ubuntu-16.04-stable .
