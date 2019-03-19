FROM ubuntu:18.04

# Ubuntu Setup
ENV TZ=America/Chicago
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone
ARG py_version=2.7

ENV HOME /root
RUN if [ "${py_version%.?}" -eq 3 ] ; \
    then \ 
            export PY_SUFIX=${py_version%.?}; \
    fi;\
    apt-get update \
    && apt-get install -y --fix-missing \
            software-properties-common \
            wget \
            g++ \
            build-essential \
            python${PY_SUFIX}-setuptools \
            python${PY_SUFIX}-pip \
            python${PY_SUFIX}-setuptools \
            python${PY_SUFIX}-dev \
            libpython${PY_SUFIX}-dev \
            python${PY_SUFIX}-nose \
            python${PY_SUFIX}-matplotlib \
            python${PY_SUFIX}-tables \
            python${PY_SUFIX}-scipy \
            python${PY_SUFIX}-jinja2 \
            gfortran \
            git \
            cmake \
            gfortran \
            vim emacs \
            libblas-dev \
            liblapack-dev \
            libhdf5-dev \
            libhdf5-serial-dev \
            autoconf \
            libtool \
            doxygen \
    && apt-get clean -y; \
    if [ "${py_version%.?}" -eq 3 ] ; \
       then \ 
            update-alternatives --install /usr/bin/python python /usr/bin/python3 10; \
            update-alternatives --install /usr/bin/pip pip /usr/bin/pip3 10; \
    fi;\
    pip install --force-reinstall \
            sphinx \
            cloud_sptheme \
            prettytable \
            sphinxcontrib_bibtex \
            numpydoc \
            nbconvert \
            numpy \
            nose \
            cython


# make starting directory
ARG build_moab=NO
ARG enable_pymoab=NO
RUN mkdir -p $HOME/opt
RUN echo "export PATH=$HOME/.local/bin:\$PATH" >> ~/.bashrc

# build MOAB
RUN if [ "$build_moab" = "YES" ] || [ "$enable_pymoab" = "YES" ] ; then \
        if [ "$enable_pymoab" = "YES" ] ; \
        then \ 
            export PYMOAB_FLAG="-DENABLE_PYMOAB=ON"; \
        fi;\
        echo $PYMOAB_FLAG \
        && cd $HOME/opt \
        && mkdir moab \
        && cd moab \
        && git clone https://bitbucket.org/fathomteam/moab \
        && cd moab \
        && git checkout -b Version5.1.0 origin/Version5.1.0 \
        && cd .. \
        && mkdir build \
        && cd build \
        && ls ../moab/ \
      # build/install static lib
        && cmake ../moab/ \
              -DCMAKE_INSTALL_PREFIX=$HOME/opt/moab \
              -DENABLE_HDF5=ON \
              -DBUILD_SHARED_LIBS=OFF \
        && make -j 3 \
        && make install \
      # build/install shared lib
        && cmake ../moab/ \
              ${PYMOAB_FLAG} \
              -DCMAKE_INSTALL_PREFIX=$HOME/opt/moab \
              -DENABLE_HDF5=ON \
              -DBUILD_SHARED_LIBS=ON \
        && make -j 3 \
        && make install \
        && cd .. \
        && rm -rf build moab ; \
    fi

# put MOAB on the path
ENV LD_LIBRARY_PATH $HOME/opt/moab/lib:$LD_LIBRARY_PATH
ENV LIBRARY_PATH $HOME/opt/moab/lib:$LIBRARY_PATH
ENV PYTHONPATH=$HOME/opt/moab/lib/python${py_version}/site-packages/

# build/install DAGMC
ARG build_dagmc=NO
ENV INSTALL_PATH=$HOME/opt/dagmc
RUN if [ "$build_dagmc" = "YES" ]; then \
        cd /root \
        && git clone https://github.com/svalinn/DAGMC.git \
        && cd DAGMC \
        && git checkout develop \
        && mkdir bld \
        && cd bld \
        && cmake .. -DMOAB_DIR=$HOME/opt/moab \
                 -DCMAKE_INSTALL_PREFIX=$INSTALL_PATH \
        && make \
        && make install ; \
    fi

ARG build_pyne=YES
# Build/Install PyNE
RUN if [ "$build_pyne" = "YES" ]; then \
        cd $HOME/opt \
        && git clone https://github.com/cnerg/pyne.git \
        && cd pyne \
        && git checkout pymoab_cleanup \
        && python setup.py install --user \
                                    --moab $HOME/opt/moab --dagmc $HOME/opt/dagmc \
                                    --clean ; \
    fi
ENV PATH $HOME/.local/bin:$PATH
RUN if [ "$build_pyne" = "YES" ]; then \
        cd $HOME \
        && nuc_data_make ; \
    fi

