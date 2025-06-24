ARG ubuntu_version=24.04
FROM ubuntu:${ubuntu_version} AS pyne-deps

# Set environment variables
ENV TZ=America/Chicago \
    HOME=/root \
    PATH=/opt/conda/bin:$PATH \
    LD_LIBRARY_PATH=/opt/conda/lib

# Use bash for RUN shell
SHELL ["/bin/bash", "-c"]

# Base system setup and timezone
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ >/etc/timezone && \
    apt-get update && \
    apt-get install -y --no-install-recommends \
        wget \
        bzip2 \
        ca-certificates && \
    apt-get clean -y  && \
    rm -rf /var/lib/apt/lists/*

# Install Miniforge (Conda) and setup environment
RUN echo 'export PATH=/opt/conda/bin:$PATH' >/etc/profile.d/conda.sh && \
    wget --quiet "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh" -O ~/miniforge.sh && \
    bash ~/miniforge.sh -b -p /opt/conda && \
    rm ~/miniforge.sh && \
    conda install -y --freeze-installed "python=3.12" && \
    mamba update -y --all && \
    mamba install -y -c conda-forge \
        expat \
        gxx \
        gcc \
        cmake \
        make \
        gfortran \
        libblas \
        liblapack \
        eigen \
        numpy \
        scipy \
        matplotlib \
        git \
        setuptools \
        pytest \
        pytables \
        jinja2 \
        cython \
        future \
        progress \
        meson \
        moab \
        dagmc \
        openmc && \
    mamba clean --all -f -y && \
    rm -rf ~/.cache ~/.conda

# ------------------------------
# Stage 2: Build PyNE
# ------------------------------
FROM pyne-deps AS pyne

# Arguments for PyNE build
ENV PYNE_MOAB_ARGS="--moab" \
    PYNE_DAGMC_ARGS="--dagmc"

# Copy PyNE source and build
COPY . $HOME/pyne
WORKDIR $HOME/pyne

RUN python setup.py install --prefix /opt/conda $PYNE_MOAB_ARGS $PYNE_DAGMC_ARGS --clean -j 4 && \
    cd tests && \
    nuc_data_make && \
    pytest -ra
