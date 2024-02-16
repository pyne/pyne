ARG pkg_mgr=apt
ARG build_hdf5="NO"
ARG pyne_test_base=openmc
ARG ubuntu_version=22.04

FROM ubuntu:${ubuntu_version} AS common_base

# Ubuntu Setup
ENV TZ=America/Chicago
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone

FROM common_base AS apt_deps
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
            libhdf5-dev \
            hdf5-tools \
    && apt-get clean -y; \
    update-alternatives --install /usr/bin/python python /usr/bin/python3 10; \
    update-alternatives --install /usr/bin/pip pip /usr/bin/pip3 10; \
    pip install --upgrade pip; \
    pip install numpy==1.23 \
            scipy \
            'cython<3' \
            nose \
            pytest \
            tables \
            matplotlib \
            jinja2 \
            setuptools \
            future \
            progress

FROM common_base AS conda_deps
RUN apt-get update \
    && apt-get install -y --fix-missing \
        wget \
        bzip2 \
        ca-certificates \
    && apt-get clean -y

RUN echo 'export PATH=/opt/conda/bin:$PATH' > /etc/profile.d/conda.sh && \
    wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh && \
    /bin/bash ~/miniconda.sh -b -p /opt/conda && \
    rm ~/miniconda.sh

ENV PATH /opt/conda/bin:$PATH

RUN conda config --add channels conda-forge
RUN conda update -n base -c defaults conda
RUN conda install -y conda-libmamba-solver
RUN conda config --set solver libmamba
RUN conda install -y mamba
RUN conda uninstall -y conda-libmamba-solver
RUN conda config --set solver classic
RUN conda update -y --all && \
    mamba install -y \
                cmake \
                git \
                libblas \
                liblapack \
                hdf5 \
                setuptools \
                pytest \
                pytables \
                jinja2 \
                "cython<3" \
                && \
    mamba install -y --force-reinstall libsqlite && \
    conda clean -y --all
RUN mkdir -p `python -m site --user-site`

FROM ${pkg_mgr}_deps AS base_python
# make starting directory
RUN mkdir -p $HOME/opt
RUN echo "export PATH=$HOME/.local/bin:\$PATH" >> ~/.bashrc

# build HDF5
ARG build_hdf5="NO"
ENV HDF5_INSTALL_PATH=$HOME/opt/hdf5/$build_hdf5
RUN if [ "$build_hdf5" != "NO" ]; then \
        cd $HOME/opt \
        && mkdir hdf5 \
        && cd hdf5 \
        && git clone --single-branch --branch $build_hdf5 https://github.com/HDFGroup/hdf5.git \
        && cd hdf5 \
        && ./configure --prefix=$HDF5_INSTALL_PATH --enable-shared \
        && make -j 3 \
        && make install \
        && cd .. \
        && rm -rf hdf5; \
    fi
# put HDF5 on the path
ENV LD_LIBRARY_PATH $HDF5_INSTALL_PATH/lib:$LD_LIBRARY_PATH
ENV LIBRARY_PATH $HDF5_INSTALL_PATH/lib:$LIBRARY_PATH


FROM base_python AS moab
ARG build_hdf5
ENV INSTALL_PATH=$HOME/opt/moab

# build MOAB
RUN export MOAB_HDF5_ARGS=""; \
    if [ "$build_hdf5" != "NO" ] ; \
    then \
            export MOAB_HDF5_ARGS="-DHDF5_ROOT=$HDF5_INSTALL_PATH"; \
    fi \
    && cd $HOME/opt \
    && mkdir moab \
    && cd moab \
    && git clone --depth 1 --single-branch -b 5.3.0 https://bitbucket.org/fathomteam/moab \
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
ENV PYNE_MOAB_ARGS "--moab $HOME/opt/moab"

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
ENV PYNE_DAGMC_ARGS "--dagmc $HOME/opt/dagmc"

FROM dagmc AS openmc
ARG build_hdf5
# build/install OpenMC Python API
RUN if [ "$build_hdf5" != "NO" ]; then \
            export HDF5_ROOT="$HDF5_INSTALL_PATH" ; \
    fi ;\
    git clone https://github.com/openmc-dev/openmc.git $HOME/opt/openmc \
    && cd  $HOME/opt/openmc \
    && git checkout tags/v0.13.0 \
    && pip install .

# Build/Install PyNE from develop branch
FROM ${pyne_test_base} AS pyne-dev
ARG build_hdf5

RUN export PYNE_HDF5_ARGS="" ;\
    if [ "$build_hdf5" != "NO" ]; then \
            export PYNE_HDF5_ARGS="--hdf5 $HDF5_INSTALL_PATH" ; \
    fi \
    && cd $HOME/opt \
    && git clone -b develop --single-branch https://github.com/pyne/pyne.git \
    && cd pyne \
    && python setup.py install --user \
                                $PYNE_MOAB_ARGS $PYNE_DAGMC_ARGS \
                                $PYNE_HDF5_ARGS \
                                --clean -j 3;
ENV PATH $HOME/.local/bin:$PATH
RUN cd $HOME \
    && nuc_data_make \
    && cd $HOME/opt/pyne/tests \
    && ./ci-run-tests.sh python3

# Build/Install PyNE from release branch
FROM ${pyne_test_base} AS pyne
ARG build_hdf5

RUN export PYNE_HDF5_ARGS="" ;\
    if [ "$build_hdf5" != "NO" ]; then \
            export PYNE_HDF5_ARGS="--hdf5 $HDF5_INSTALL_PATH" ; \
    fi;
COPY . $HOME/opt/pyne
RUN cd $HOME/opt/pyne \
    && python setup.py install --user \
                                $PYNE_MOAB_ARGS $PYNE_DAGMC_ARGS \
                                $PYNE_HDF5_ARGS \
                                --clean -j 3;
ENV PATH $HOME/.local/bin:$PATH
RUN cd $HOME \
    && nuc_data_make \
    && cd $HOME/opt/pyne/tests \
    && ./ci-run-tests.sh python3
