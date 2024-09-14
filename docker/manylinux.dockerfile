# Setup base image
ARG PyNE_STAGE=openmc
ARG MANYLINUX_IMAGE=manylinux_2_28_x86_64
ARG Python_ABI="cp312-cp312"
ARG HDF5_VERSION="1.14.3"
ARG EIGEN3_VERSION="3.4.0"
ARG LAPACK_VERSION="3.12.0"
ARG MOAB_VERSION="5.5.1"
ARG DAGMC_VERSION="2.6.0"
ARG OpenMC_VERSION="0.14.0"

# Build base stage
FROM quay.io/pypa/${MANYLINUX_IMAGE} AS base

ARG Python_ABI
ARG HDF5_VERSION
ARG EIGEN3_VERSION
ARG LAPACK_VERSION

# Set timezone
ENV TZ=America/Chicago
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone

# Install basic dependencies
RUN yum install -y \ 
        wget \
        git \
        gcc \
        gcc-gfortran \
        make && \
    yum clean all

# Use Python from manylinux as the default Python
ENV PATH="/opt/python/${Python_ABI}/bin:${PATH}"
RUN ln -sf /opt/python/${Python_ABI}/bin/python3 /usr/bin/python

# Set environment variables for installation paths
ENV HDF5_ROOT=/opt/hdf5
ENV EIGEN3_ROOT=/opt/eigen3
ENV MOAB_ROOT=/opt/moab
ENV DAGMC_ROOT=/opt/dagmc
ENV LAPACK_ROOT=/opt/lapack

# Build and install HDF5
RUN HDF5_MAJOR_VERSION=$(echo ${HDF5_VERSION} | cut -d'.' -f1) && \
    HDF5_MINOR_VERSION=$(echo ${HDF5_VERSION} | cut -d'.' -f2) && \
    wget https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-${HDF5_MAJOR_VERSION}.${HDF5_MINOR_VERSION}/hdf5-${HDF5_VERSION}/src/hdf5-${HDF5_VERSION}.tar.gz && \
    tar -xzf hdf5-${HDF5_VERSION}.tar.gz && \
    cd hdf5-${HDF5_VERSION} && \
    mkdir -p build && cd build && \
    cmake .. \
        -DCMAKE_INSTALL_PREFIX=${HDF5_ROOT} \
        -DBUILD_SHARED_LIBS=ON \
        -DBUILD_STATIC_LIBS=OFF \
        -DBUILD_TESTING=OFF \
        -DHDF5_BUILD_TOOLS=OFF \
        -DHDF5_BUILD_EXAMPLES=OFF && \
    make -j$(nproc) && \
    make install && \
    cd ../.. && \
    rm -rf hdf5-${HDF5_VERSION} hdf5-${HDF5_VERSION}.tar.gz

# Add HDF5 to the system path
ENV PATH="${HDF5_ROOT}/bin:${PATH}"
ENV LD_LIBRARY_PATH="${HDF5_ROOT}/lib:${LD_LIBRARY_PATH}"

# Build and install Eigen3
RUN wget https://gitlab.com/libeigen/eigen/-/archive/${EIGEN3_VERSION}/eigen-${EIGEN3_VERSION}.tar.bz2 && \
    tar -xjf eigen-${EIGEN3_VERSION}.tar.bz2 && \
    cd eigen-${EIGEN3_VERSION} && \
    mkdir -p build && cd build && \
    cmake .. \
        -DCMAKE_INSTALL_PREFIX=${EIGEN3_ROOT} && \
    make install && \
    cd ../.. && \
    rm -rf eigen-${EIGEN3_VERSION} eigen-${EIGEN3_VERSION}.tar.bz2

# Add Eigen3 to the system path
ENV CPLUS_INCLUDE_PATH="${EIGEN3_ROOT}/include/eigen3:${CPLUS_INCLUDE_PATH}"

# Build and install LAPACK
RUN wget https://github.com/Reference-LAPACK/lapack/archive/refs/tags/v${LAPACK_VERSION}.tar.gz && \
    tar -xzf v${LAPACK_VERSION}.tar.gz && \
    cd lapack-${LAPACK_VERSION} && \
    mkdir -p build && cd build && \
    cmake .. \
        -DCMAKE_INSTALL_PREFIX=${LAPACK_ROOT} \
        -DBUILD_SHARED_LIBS=ON \
        -DUSE_OPTIMIZED_BLAS=ON && \
    make -j$(nproc) && make install && \
    cd ../.. && \
    rm -rf lapack-${LAPACK_VERSION} v${LAPACK_VERSION}.tar.gz

# Add LAPACK to the system path
ENV LD_LIBRARY_PATH="${LAPACK_ROOT}/lib:${LD_LIBRARY_PATH}"

# Build MOAB stage
FROM base AS moab

ARG MOAB_VERSION

# Build and install MOAB
RUN git clone --depth 1 -b ${MOAB_VERSION} https://bitbucket.org/fathomteam/moab.git moab && \
    cd moab && \
    mkdir build && cd build && \
    cmake .. \
        -DCMAKE_INSTALL_PREFIX=${MOAB_ROOT} \
        -DBUILD_SHARED_LIBS=ON \
        -DENABLE_HDF5=ON \
        -DHDF5_ROOT=${HDF5_ROOT} \
        -DEIGEN3_DIR=${EIGEN3_ROOT}/include/eigen3 \
        -DENABLE_BLASLAPACK=OFF \
        -DENABLE_FORTRAN=OFF \
        -DENABLE_PYMOAB=ON && \
    make -j$(nproc) && make install && \
    cd ../.. && \
    rm -rf moab

# Add MOAB to the system path
ENV PATH="${MOAB_ROOT}/bin:${PATH}"
ENV LD_LIBRARY_PATH="${MOAB_ROOT}/lib:${LD_LIBRARY_PATH}"
ENV PYNE_MOAB_ARGS "-DENABLE_MOAB=ON;-DMOAB_ROOT=${MOAB_ROOT}"


# Build DAGMC stage
FROM moab AS dagmc

ARG DAGMC_VERSION

# Build and install DAGMC
RUN git clone --depth 1 -b v${DAGMC_VERSION} https://github.com/svalinn/DAGMC.git dagmc && \
    cd dagmc && \
    mkdir -p build && cd build && \
    cmake .. \
        -DCMAKE_INSTALL_PREFIX=${DAGMC_ROOT} \
        -DMOAB_DIR=${MOAB_ROOT} \
        -DBUILD_STATIC_LIBS=OFF \
        -DBUILD_UWUW=OFF \
        -DBUILD_TALLY=OFF \
        -DBUILD_MAKE_WATERTIGHT=OFF \
        -DBUILD_OVERLAP_CHECK=OFF \
        -DBUILD_TESTS=OFF \
        -DBUILD_RPATH=OFF && \
    make -j$(nproc) && make install && \
    cd ../.. && \
    rm -rf dagmc


# Add DAGMC to the system path
ENV PATH="${DAGMC_ROOT}/bin:${PATH}"
ENV LD_LIBRARY_PATH="${DAGMC_ROOT}/lib:${LD_LIBRARY_PATH}"
ENV PYNE_DAGMC_ARGS "-DENABLE_DAGMC=ON;-DDAGMC_ROOT=${DAGMC_ROOT}"


# Build OpenMC stage
FROM dagmc AS openmc

ARG OpenMC_VERSION

# Build and install OpenMC
RUN git clone --depth 1 -b v${OpenMC_VERSION} https://github.com/openmc-dev/openmc.git openmc && \
    cd openmc && \
    mkdir -p build && cd build && \
    cmake .. \
        -DCMAKE_INSTALL_PREFIX=${OPENMC_ROOT} \
        -DHDF5_ROOT=${HDF5_ROOT} && \
    make -j$(nproc) && make install && \
    cd .. && \
    python -m pip install . && \
    cd .. && \
    rm -rf openmc

# Add OpenMC to the system path
ENV PATH="${OPENMC_ROOT}/bin:${PATH}"
ENV LD_LIBRARY_PATH="${OPENMC_ROOT}/lib:${LD_LIBRARY_PATH}"


# Build PyNE stage
FROM ${PyNE_STAGE} AS pyne

# Copy PyNE sources
COPY . /opt/pyne

# Configure SKBUILD CMake arguments
ENV SKBUILD_CMAKE_ARGS "-DDOWNLOAD_HDF5=OFF;\
                        -DHDF5_ROOT=$HDF5_INSTALL_PATH;\
                        -DDOWNLOAD_EIGEN3=OFF;\
                        -DDOWNLOAD_LAPACK=OFF;\
                        -DENABLE_MOAB=OFF;\
                        -DDOWNLOAD_MOAB=OFF;\
                        -DENABLE_DAGMC=OFF;\
                        -DDOWNLOAD_DAGMC=OFF;\
                        $PYNE_MOAB_ARGS;\
                        $PYNE_DAGMC_ARGS"

# Build and install PyNE
RUN cd /opt/pyne && python -m pip install .

# Test PyNE
RUN cd /opt/pyne/tests && \
    nuc_data_make && \
    pytest -ra
