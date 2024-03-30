ARG build_hdf5="NO"
ARG pyne_test_base=openmc
ARG ubuntu_version=22.04

FROM ubuntu:${ubuntu_version} AS base_conda

# Ubuntu Setup
ENV TZ=America/Chicago
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone
ENV HOME /root

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

# install python 3.10 to match apt dockerfile
RUN conda update conda
RUN conda install "python<3.10"
RUN conda config --add channels conda-forge
RUN conda update -n base -c defaults conda
RUN conda install -y conda-libmamba-solver
RUN conda config --set solver libmamba
RUN conda install -y mamba
RUN conda uninstall -y conda-libmamba-solver
RUN conda config --set solver classic
RUN conda update -y --all && \
    mamba install -y \
                gxx_linux-64 \
                gcc_linux-64 \
                cmake \
                make \
                gfortran \
                libblas \
                liblapack \
                eigen \
                numpy \
                scipy \
                nose \
                matplotlib \
                git \
                setuptools \
                pytest \
                pytables \
                jinja2 \
                "cython<3" \
                future \
                progress \
                && \
    mamba install -y --force-reinstall libsqlite && \
    conda clean -y --all
RUN mkdir -p `python -m site --user-site`

ENV CC /opt/conda/bin/x86_64-conda_cos6-linux-gnu-gcc
ENV CXX /opt/conda/bin/x86_64-conda_cos6-linux-gnu-g++
ENV CPP /opt/conda/bin/x86_64-conda_cos6-linux-gnu-cpp

# put conda on the path
ENV LD_LIBRARY_PATH /opt/conda/lib:$LD_LIBRARY_PATH

# make starting directory
RUN mkdir -p $HOME/opt
RUN echo "export PATH=$HOME/.local/bin:\$PATH" >> ~/.bashrc

FROM base_conda AS moab
RUN conda install "conda-forge::moab=5.3.0"

FROM moab AS dagmc
RUN conda install conda-forge::dagmc

FROM dagmc AS openmc
RUN conda install conda-forge::openmc

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
    