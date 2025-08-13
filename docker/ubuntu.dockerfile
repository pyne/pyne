# ------------------------------
# Stage 1: Base Python Environment
# ------------------------------
ARG PYNE_TEST_BASE=openmc
ARG UBUNTU_VERSION=24.04
ARG MAKE_CORES=4
ARG INSTALL_PATH=/opt

FROM ubuntu:${UBUNTU_VERSION} AS base

ARG INSTALL_PATH

# Set environment variables
ENV TZ=America/Chicago \
    HOME=/root \
    VENV_PATH=${INSTALL_PATH}/venv

# Set working directory
WORKDIR ${HOME}

# Use bash shell
SHELL ["/bin/bash", "-c"]

# System setup
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone && \
    apt-get update && \
    apt-get install -y --no-install-recommends \
        software-properties-common \
        python3-dev \
        python3-venv \
        wget \
        build-essential \
        git \
        cmake \
        gfortran \
        libeigen3-dev \
        libblas-dev \
        liblapack-dev && \
    apt-get clean && rm -rf /var/lib/apt/lists/*

# Create and activate virtual environment
RUN python3 -m venv ${VENV_PATH} && \
    source ${VENV_PATH}/bin/activate && \
    pip install --upgrade pip && \
    pip install \
        numpy \
        scipy \
        cython \
        pytest \
        tables \
        matplotlib \
        jinja2 \
        setuptools \
        future \
        progress

ENV PATH=${VENV_PATH}/bin:$PATH \
    LD_LIBRARY_PATH=${VENV_PATH}/lib

# HDF5 setup
ARG HDF5_VERSION=1.14.6
ARG MAKE_CORES

RUN git clone --depth 1 --branch hdf5-${HDF5_VERSION} https://github.com/HDFGroup/hdf5.git && \
    cd hdf5 && \
    mkdir build && cd build && \
    cmake .. \
        -DCMAKE_INSTALL_PREFIX=${VENV_PATH} \
        -DBUILD_SHARED_LIBS=ON && \
    make -j ${MAKE_CORES} && \
    make install && \
    rm -rf $HOME/hdf5

# ------------------------------
# Stage 2: MOAB Setup
# ------------------------------
FROM base AS moab

ARG MOAB_VERSION=5.5.1
ARG MAKE_CORES

RUN git clone --depth 1 --branch ${MOAB_VERSION} https://bitbucket.org/fathomteam/moab && \
    cd moab && \
    mkdir build && cd build && \
    cmake .. \
        -DENABLE_PYMOAB=ON \
        -DCMAKE_INSTALL_PREFIX=${VENV_PATH} \
        -DENABLE_HDF5=ON \
        -DBUILD_SHARED_LIBS=ON \
        -DENABLE_BLASLAPACK=OFF \
        -DENABLE_FORTRAN=OFF && \
    make -j ${MAKE_CORES} && \
    make install && \
    rm -rf $HOME/moab

ENV PYNE_MOAB_ARGS="--moab ${VENV_PATH}"

# ------------------------------
# Stage 3: DAGMC Setup
# ------------------------------
FROM moab AS dagmc

ARG DAGMC_VERSION=v3.2.4
ARG MAKE_CORES

RUN git clone --depth 1 --branch ${DAGMC_VERSION} https://github.com/svalinn/DAGMC.git && \
    cd DAGMC && \
    mkdir build && cd build && \
    cmake .. \
        -DMOAB_DIR=${VENV_PATH} \
        -DCMAKE_INSTALL_PREFIX=${VENV_PATH} \
        -DBUILD_STATIC_LIBS=OFF \
        -DBUILD_UWUW=OFF \
        -DBUILD_TALLY=OFF \
        -DBUILD_MAKE_WATERTIGHT=OFF \
        -DBUILD_OVERLAP_CHECK=OFF \
        -DBUILD_TESTS=OFF && \
    make -j ${MAKE_CORES} && \
    make install && \
    rm -rf $HOME/DAGMC

ENV PYNE_DAGMC_ARGS="--dagmc ${VENV_PATH}"

# ------------------------------
# Stage 4: OpenMC Setup
# ------------------------------
FROM dagmc AS openmc

ARG OPENMC_VERSION=v0.15.2
ARG MAKE_CORES

RUN git clone --depth 1 --branch ${OPENMC_VERSION} https://github.com/openmc-dev/openmc.git && \
    cd openmc && \
    pip install . && \
    rm -rf $HOME/openmc

# ------------------------------
# Stage 5: PyNE Build & Test
# ------------------------------
FROM ${PYNE_TEST_BASE} AS pyne
ARG MAKE_CORES

COPY . ${HOME}/pyne
WORKDIR ${HOME}/pyne

RUN python setup.py install $PYNE_MOAB_ARGS $PYNE_DAGMC_ARGS \
        --clean -j ${MAKE_CORES} && \
    cd tests && \
    nuc_data_make && \
    pytest -ra
