#
# Ubuntu Dockerfile
#
# https://github.com/dockerfile/ubuntu
#

# pull the base image
FROM ubuntu:14.04

# install
RUN \
    apt-get -y update && \
    apt-get -y install build-essential python-numpy python-scipy cython  && \
    apt-get -y install python-nose git cmake vim emacs gfortran libblas-dev && \
    apt-get -y install liblapack-dev libhdf5-dev gfortran python-tables && \
    apt-get -y install python-matplotlib python-jinja2 autoconf libtool && \
    apt-get -y install python-pip wget 
RUN \ 
    cd /opt && \
    mkdir moab && \
    git clone https://bitbucket.org/fathomteam/moab && \
    cd moab && \
    git checkout Version5.0 && \
    autoreconf -fi && \
    mkdir bld && \
    cd bld && \
    ../configure --enable-optimize --enable-shared --disable-ahf --with-hdf5=/usr/lib/x86_64-linux-gnu/ --prefix=/opt/moab/ && \
    make -j4 && \ 
    make check && \ 
    make install && \
    cd .. && \ 
    rm -rf bld moab && \
# this line patches the current moab version for pytaps comptability
    sed -i '/IMESH_LIBS/i IMSQ_LIBS = \n MSQ_LIBS = ' lib/iMesh-Defs.inc 
RUN \
    echo "export LD_LIBRARY_PATH=/opt/moab/lib" >> /opt/env.sh && \
    echo 'export PATH=$PATH:/opt/moab/bin' >> /opt/env.sh 
RUN \
    cd /opt/ && \
    pip install numpy && \
    wget https://pypi.python.org/packages/7a/51/c5915ca128da8ab178c3f2a52eae9abcb70ac9b9f5e0d9dd917007730dc8/PyTAPS-1.4.tar.gz && \
    tar -zxf PyTAPS-1.4.tar.gz && \
    mv PyTAPS-1.4 pytaps && \
    cd pytaps && \
    python setup.py --iMesh-path=/opt/moab --without-iRel --without-iGeom install --user
RUN touch /root/.bashrc
RUN echo 'source /opt/env.sh' >> /root/.bashrc
# Set environment variables.
ENV HOME /root
# Define working directory.
WORKDIR /root
# Define default command.
CMD ["bash"]
    