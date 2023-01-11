FROM ubuntu:20.04 AS base_python

# Ubuntu Setup
ENV TZ=America/Chicago
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone
ARG py_version=3.8

ENV HOME /root
RUN if [ "${py_version%.?}" -eq 3 ] ; \
    then \ 
            export PY_SUFIX=${py_version%.?}; \
    fi;\
    apt-get update \
    && apt-get install -y --fix-missing \
            software-properties-common \
            python${PY_SUFIX}-pip \
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
    if [ "${py_version%.?}" -eq 3 ] ; \
       then \ 
            update-alternatives --install /usr/bin/python python /usr/bin/python3 10; \
            update-alternatives --install /usr/bin/pip pip /usr/bin/pip3 10; \
    fi;\
    pip install --upgrade pip; \
    pip install numpy \
            scipy \
            cython \
            nose \
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
ENV INSTALL_PATH=$HOME/opt/moab

# build MOAB
RUN apt-get update \
    && apt-get install -y --fix-missing \
            libhdf5-dev \
            hdf5-tools \
    && apt-get clean -y;
    
    export PYMOAB_FLAG="-DENABLE_PYMOAB=ON"; \
    echo $PYMOAB_FLAG ;\
    export MOAB_HDF5_ARGS=""; \
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
            ${PYMOAB_FLAG} \
            -DCMAKE_INSTALL_PREFIX=$INSTALL_PATH \
            -DENABLE_HDF5=ON $MOAB_HDF5_ARGS \
            -DBUILD_SHARED_LIBS=ON \
            -DENABLE_BLASLAPACK=OFF \
            -DENABLE_FORTRAN=OFF \
    && make -j 3 \
    && make install \
    && cd .. \
    && rm -rf moab ;

# put MOAB on the path
ENV LD_LIBRARY_PATH $HOME/opt/moab/lib:$LD_LIBRARY_PATH
ENV LIBRARY_PATH $HOME/opt/moab/lib:$LIBRARY_PATH
ENV PYTHONPATH=$HOME/opt/moab/lib/python${py_version}/site-packages/

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

# Build/Install PyNE
RUN export PYNE_HDF5_ARGS="" ;\
    if [ "$build_hdf5" != "NO" ]; then \
            export PYNE_HDF5_ARGS="--hdf5 $HDF5_INSTALL_PATH" ; \
    fi \
    && cd $HOME/opt \
    && git clone -b develop --single-branch https://github.com/pyne/pyne.git \
    && cd pyne \
    && python setup.py install --user \
                                --moab $HOME/opt/moab --dagmc $HOME/opt/dagmc \
                                $PYNE_HDF5_ARGS \
                                --clean -j 3; \
ENV PATH $HOME/.local/bin:$PATH
RUN if [ "$build_pyne" = "YES" ]; then \
        cd $HOME \
        && nuc_data_make ; \
    fi

FROM dagmc AS openmc
# build/install OpenMC Python API
RUN if [ "$build_hdf5" != "NO" ]; then \
            export HDF5_ROOT="$HDF5_INSTALL_PATH" ; \
    fi ;\
    git clone https://github.com/openmc-dev/openmc.git $HOME/opt/openmc \
    && cd  $HOME/opt/openmc \
    && pip install . 

