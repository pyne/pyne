ARG pyne_test_base=openmc
ARG ubuntu_version=22.04

FROM ubuntu:${ubuntu_version} AS base_python

# Ubuntu Setup
ENV TZ=America/Chicago
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone

ENV HOME /root
RUN apt-get update \
    && apt-get install -y --fix-missing \
            software-properties-common \
            python3-pip \
            wget \
            build-essential \
            git \
            cmake \
            gfortran \
            libblas-dev \
            liblapack-dev \
            libeigen3-dev \
    && apt-get clean -y; \
    update-alternatives --install /usr/bin/python python /usr/bin/python3 10; \
    update-alternatives --install /usr/bin/pip pip /usr/bin/pip3 10; \
    pip install --upgrade pip; \
    pip install numpy==1.23 \
            scipy \
            cython \
            pytest \
            tables \
            matplotlib \
            jinja2 \
            setuptools \
            future \
            progress

# make starting directory
RUN mkdir -p $HOME/opt
RUN echo "export PATH=$HOME/.local/bin:\$PATH" >> ~/.bashrc

# build HDF5
ARG build_hdf5="hdf5-1_14_3"
ENV HDF5_INSTALL_PATH=$HOME/opt/hdf5/$build_hdf5
RUN cd $HOME/opt \
    && mkdir hdf5 \
    && cd hdf5 \
    && git clone --single-branch --branch $build_hdf5 https://github.com/HDFGroup/hdf5.git \
    && cd hdf5 \
    && ./configure --prefix=$HDF5_INSTALL_PATH --enable-shared \
    && make -j 3 \
    && make install \
    && cd .. \
    && rm -rf hdf5;
# put HDF5 on the path
ENV LD_LIBRARY_PATH $HDF5_INSTALL_PATH/lib:$LD_LIBRARY_PATH
ENV LIBRARY_PATH $HDF5_INSTALL_PATH/lib:$LIBRARY_PATH
RUN echo "export PATH=$PATH:$HDF5_INSTALL_PATH" >> ~/.bashrc

FROM base_python AS moab
ARG moab_version="5.5.1"
ENV INSTALL_PATH=$HOME/opt/moab

# build MOAB
RUN export MOAB_HDF5_ARGS="-DHDF5_ROOT=$HDF5_INSTALL_PATH"; \
    cd $HOME/opt \
    && mkdir moab \
    && cd moab \
    && git clone --depth 1 --single-branch -b $moab_version https://bitbucket.org/fathomteam/moab \
    && cd moab \
    && mkdir build \
    && cd build \
    && ls ..\
    # build/install shared lib
    && cmake .. \
            -DENABLE_PYMOAB=ON \
            -DCMAKE_INSTALL_PREFIX=$INSTALL_PATH \
            -DENABLE_HDF5=ON $MOAB_HDF5_ARGS \
            -DBUILD_SHARED_LIBS=ON \
            -DENABLE_BLASLAPACK=OFF \
            -DENABLE_FORTRAN=OFF \
    && make -j 3 \
    && cd pymoab \
    && pip install . \
    && cd .. \
    && make install \
    && cd .. \
    && rm -rf moab ;

# put MOAB on the path
ENV LD_LIBRARY_PATH $HOME/opt/moab/lib:$LD_LIBRARY_PATH
ENV LIBRARY_PATH $HOME/opt/moab/lib:$LIBRARY_PATH
ENV PYNE_MOAB_ARGS "-DENABLE_MOAB=ON;-DMOAB_ROOT=$HOME/opt/moab"

FROM moab AS dagmc
# build/install DAGMC
ENV INSTALL_PATH=$HOME/opt/dagmc
RUN cd /root \
    && git clone --depth 1 --branch stable https://github.com/svalinn/DAGMC.git \
    && cd DAGMC \
    && mkdir bld \
    && cd bld \
    && cmake .. -DMOAB_DIR=$HOME/opt/moab \
                -DCMAKE_INSTALL_PREFIX=$INSTALL_PATH \
                -DBUILD_STATIC_LIBS=OFF \
                -DBUILD_UWUW=OFF \
                -DBUILD_TALLY=OFF \
                -DBUILD_MAKE_WATERTIGHT=OFF \
                -DBUILD_OVERLAP_CHECK=OFF \
                -DBUILD_TESTS=OFF \
    && make -j 3\
    && make install \
    && cd ../.. \
    && rm -rf DAGMC
ENV PYNE_DAGMC_ARGS "-DENABLE_DAGMC=ON;-DDAGMC_ROOT=$HOME/opt/dagmc"

FROM dagmc AS openmc
ARG openmc_version="v0.14.0"
# build/install OpenMC Python API
RUN export HDF5_ROOT="$HDF5_INSTALL_PATH" ; \
    git clone --depth 1 --branch $openmc_version https://github.com/openmc-dev/openmc.git $HOME/opt/openmc \
    && cd  $HOME/opt/openmc \
    && pip install .

# Build/Install PyNE from release branch
FROM ${pyne_test_base} AS pyne

ENV SKBUILD_CMAKE_ARGS "-DDOWNLOAD_HDF5=OFF;-DDOWNLOAD_EIGEN3=OFF;-DDOWNLOAD_LAPACK=OFF;-DENABLE_MOAB=OFF;-DDOWNLOAD_MOAB=OFF;-DENABLE_DAGMC=OFF;-DDOWNLOAD_DAGMC=OFF;$PYNE_MOAB_ARGS;$PYNE_DAGMC_ARGS"

COPY . $HOME/opt/pyne
RUN cd $HOME/pyne \
    && python -m pip -v install .
ENV PATH $HOME/.local/bin:$PATH
RUN cd $HOME \
    && nuc_data_make \
    && cd $HOME/pyne/tests \
    && pytest -ra