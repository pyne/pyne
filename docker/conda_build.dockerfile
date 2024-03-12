ARG build_hdf5="NO"
ARG pyne_test_base=openmc
ARG ubuntu_version=22.04

FROM ubuntu:${ubuntu_version} AS base_conda

# Ubuntu Setup
ENV TZ=America/Chicago
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone

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
                cython \
                && \
    mamba install -y --force-reinstall libsqlite && \
    conda clean -y --all
RUN mkdir -p `python -m site --user-site`

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

FROM base_conda AS moab
RUN conda install conda-forge::moab

FROM moab AS dagmc
RUN conda install conda-forge::dagmc

FROM dagmc AS openmc
RUN conda install conda-forge::openmc

# Build/Install PyNE from release branch
FROM ${pyne_test_base} AS pyne
RUN conda install conda-forge::pyne

ENV PATH $HOME/.local/bin:$PATH
RUN cd $HOME \
    && nuc_data_make \
    && cd $HOME/opt/pyne/tests \
    && ./ci-run-tests.sh python3

 