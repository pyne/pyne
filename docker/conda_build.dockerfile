ARG ubuntu_version=22.04

FROM ubuntu:${ubuntu_version} AS pyne-deps

# Ubuntu Setup
ENV TZ=America/Chicago
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone

ENV HOME=/root
RUN apt-get update \
    && apt-get install -y --fix-missing \
        wget \
        bzip2 \
        ca-certificates \
    && apt-get clean -y

RUN echo 'export PATH=/opt/conda/bin:$PATH' > /etc/profile.d/conda.sh && \
    wget --quiet "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh" -O ~/miniforge.sh && \
    /bin/bash ~/miniforge.sh -b -p /opt/conda && \
    rm ~/miniforge.sh
    
ENV PATH=/opt/conda/bin:$PATH

FROM pyne-deps AS pyne-conda

RUN mamba update -n base conda mamba
RUN mamba install -y conda-forge::openmc 
# required by PyNE but should already be installed by OpenMC
RUN mamba install -y \
                libblas \
                liblapack \
                eigen \
                numpy \
                scipy \
                future \
                "conda-forge::moab=5.5.1" \
                conda-forge::dagmc
# required by PyNE but not installed by OpenMC
RUN mamba install -y \
                setuptools \
                expat \
                gxx_linux-64 gcc_linux-64 \
                cmake \
                make \
                gfortran \
                matplotlib \
                git \
                pytest \
                pytables \
                jinja2 \
                cython \
                progress \
                meson \
                && \
    mamba clean -y --all
RUN mkdir -p $(python3 -m site --user-site)
ENV CC=/opt/conda/bin/x86_64-conda-linux-gnu-gcc
ENV CXX=/opt/conda/bin/x86_64-conda-linux-gnu-g++
ENV CPP=/opt/conda/bin/x86_64-conda-linux-gnu-cpp

# Build/Install PyNE from release branch
FROM pyne-conda AS pyne

# put conda on the path
ENV LD_LIBRARY_PATH=/opt/conda/lib

# make starting directory
RUN mkdir -p $HOME/opt
RUN echo "export PATH=$HOME/.local/bin:\$PATH" >> ~/.bashrc

ENV PYNE_MOAB_ARGS="--moab"
ENV PYNE_DAGMC_ARGS="--dagmc"

COPY . $HOME/opt/pyne
RUN cd $HOME/opt/pyne \
    && python3 setup.py install --user \
                            $PYNE_MOAB_ARGS $PYNE_DAGMC_ARGS \
                            --clean -j 8;

FROM pyne AS pyne-test

ENV PATH=$HOME/.local/bin:$PATH
RUN cd $HOME \
    && python $HOME/.local/.bin/nuc_data_make \
    && cd $HOME/opt/pyne/tests \
    && ./ci-run-tests.sh python3
    