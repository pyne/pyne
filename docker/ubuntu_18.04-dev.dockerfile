FROM ubuntu:18.04

# Ubuntu Setup
ENV TZ=America/Chicago
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone

ENV HOME /root
RUN apt-get update \
    && apt-get install -y --fix-missing \
            software-properties-common wget g++ \
            build-essential python-numpy python-scipy cython python-setuptools \
            python-nose git cmake vim emacs gfortran libblas-dev \
            liblapack-dev libhdf5-dev libhdf5-serial-dev gfortran python-tables \
            python-matplotlib python-jinja2 python-dev libpython-dev \
            autoconf libtool python-setuptools python-pip doxygen \
    && apt-get clean -y \
    && pip install --force-reinstall \
            sphinx \
            cloud_sptheme \
            prettytable \
            sphinxcontrib_bibtex \
            numpydoc \
            nbconvert \
            numpy \
            cython

# Script conditional setup: Default PyNE alone
ARG build_moab=NO
ARG enable_pymoab=NO
ARG build_dagmc=NO
ARG build_pyne=YES
RUN echo "Configuration: \n" \
         "MOAB: $build_moab\n" \
         "pyMOAB: $enable_pymoab\n" \
         "DAGMC: $build_dagmc\n" \
         "PyNE: $build_pyne\n"


# make starting directory
RUN mkdir -p $HOME/opt
RUN echo "export PATH=$HOME/.local/bin:\$PATH" >> ~/.bashrc

# build MOAB
RUN if [ "$build_moab" = "YES" ] || [ "$build_dagmc" = "YES" ] || [ "$enable_pymoab" = "YES" ] ; then \
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
        && cmake ../moab/ \
              ${PYMOAB_FLAG} \
              -DCMAKE_INSTALL_PREFIX=$HOME/opt/moab \
              -DENABLE_HDF5=ON \
        && make -j 3 \
        && make install \
      # build/install static lib
        && cmake ../moab/ \
              -DCMAKE_INSTALL_PREFIX=$HOME/opt/moab \
              -DENABLE_HDF5=ON \
              -DBUILD_SHARED_LIBS=OFF \
        && make -j 3 \
        && make install \
        && cd .. \
        && rm -rf build moab ; \
    fi

# put MOAB on the path
ENV LD_LIBRARY_PATH $HOME/opt/moab/lib:$LD_LIBRARY_PATH
ENV LIBRARY_PATH $HOME/opt/moab/lib:$LIBRARY_PATH
ENV PYTHONPATH=$HOME/opt/moab/lib/python2.7/site-packages/

# build/install DAGMC
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

