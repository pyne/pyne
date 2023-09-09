#!/bin/bash
TZ=America/Chicago
ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone

# list of package installed through apt-get 
# (required to run this scripts and/or as dependancies for PyNE and its depedancies)
apt_package_list="software-properties-common \
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
                  hdf5-tools"

# list of python package required for PyNE and its depedencies (installed using pip3 python package manager)
pip_package_list="numpy==1.23 \
                  scipy \
                  \'cython<3\' \
                  nose \
                  pytest \
                  tables \
                  matplotlib \
                  jinja2 \
                  setuptools \
                  future \
                  progress"

# function to check if build folder already exists
function check_repo() {

    repo_name=$1

    if [ -d ${repo_name} ] ; then
        read -p "Delete the existing ${repo_name} directory and all contents? (y/n) " -n 1 -r
        if [[ $REPLY =~ ^[Yy]$ ]] ; then
            rm -rf ${repo_name}
        fi
    fi

}

# system update
sudo apt-get -y update
sudo apt-get install -y --fix-missing ${apt_package_list}
apt-get clean -y;
update-alternatives --install /usr/bin/python python /usr/bin/python3 10;
update-alternatives --install /usr/bin/pip pip /usr/bin/pip3 10;
pip install --user --upgrade pip;
pip install --user ${pip_package_list}

install_dir=${HOME}/opt
mkdir -p ${install_dir}
echo "export PATH=$HOME/.local/bin:\$PATH" >> ~/.bashrc

# build HDF5 if necessary
build_hdf5=$2
hdf5_libdir=$HOME/opt/hdf5/$build_hdf5
if [ $build_hdf5 != "NO" ]; then \
  cd $install_dir \
  && mkdir hdf5 \
  && cd hdf5 \
  && git clone --single-branch --branch $build_hdf5 https://github.com/HDFGroup/hdf5.git \
  && cd hdf5 \
  && ./configure --prefix=$hdf5_libdir --enable-shared \
  && make -j 3 \
  && make install \
  && cd .. \
  && rm -rf hdf5; \
fi

# put HDF5 on the path
export LD_LIBRARY_PATH=$hdf5_libdir/lib:$LD_LIBRARY_PATH
export LIBRARY_PATH=$hdf5_libdir/lib:$LIBRARY_PATH
# # need to put libhdf5.so on LD_LIBRARY_PATH (Making sure that LD_LIBRARY_PATH is defined first)
# if [ -z $LD_LIBRARY_PATH ]; then
#   export LD_LIBRARY_PATH="${hdf5_libdir}"
# else
#   export LD_LIBRARY_PATH="${hdf5_libdir}:$LD_LIBRARY_PATH"
# fi


############
### MOAB ###
############
# pre-setup
export MOAB_HDF5_ARGS=""
if [ $build_hdf5 != "NO" ]; then \
  export MOAB_HDF5_ARGS="-DHDF5_ROOT=$hdf5_libdir"; \
fi
cd ${install_dir}
check_repo moab
mkdir -p moab
cd moab

# clone and version
git clone --depth 1 --single-branch -b 5.3.0 https://bitbucket.org/fathomteam/moab moab-repo
cd moab-repo
mkdir -p build
cd build

# cmake, build and install
cmake ../ -DENABLE_PYMOAB=ON \
          -DCMAKE_INSTALL_PREFIX=${install_dir}/moab \
          -DENABLE_HDF5=ON $MOAB_HDF5_ARGS \
          -DBUILD_SHARED_LIBS=ON \
          -DENABLE_BLASLAPACK=OFF \
          -DENABLE_FORTRAN=OFF    
make -j 3 
cd pymoab
pip install --user .
cd ..
make install
cd ..
rm -rf moab-repo

# Adding MOAB/lib to $LD_LIBRARY_PATH and $LIBRARY_PATH
export LD_LIBRARY_PATH="${install_dir}/moab/lib:$LD_LIBRARY_PATH"
export LIBRARY_PATH="${install_dir}/moab/lib:$LIBRARY_PATH"

#############
### DAGMC ###
#############
# pre-setup check that the directory we need are in place
cd ${install_dir}
check_repo dagmc
mkdir -p dagmc
cd dagmc

# clone and version
git clone --depth 1 --branch stable https://github.com/svalinn/DAGMC.git dagmc-repo
cd dagmc-repo
git checkout develop
mkdir build
cd build

# cmake, build and install
cmake .. -DMOAB_DIR=${install_dir}/moab \
         -DCMAKE_INSTALL_PREFIX=${install_dir}/dagmc \
         -DBUILD_STATIC_LIBS=OFF \
         -DBUILD_UWUW=OFF \
         -DBUILD_TALLY=OFF \
         -DBUILD_MAKE_WATERTIGHT=OFF \
         -DBUILD_OVERLAP_CHECK=OFF \
         -DBUILD_TESTS=OFF 
make -j 3
make install
cd ../..
rm -rf dagmc-repo

# Adding DAGMC/lib to $LD_LIBRARY_PATH and $LIBRARY_PATH
export LD_LIBRARY_PATH="${install_dir}/dagmc/lib:$LD_LIBRARY_PATH"

# Adding dagmc/bin to $PATH
export PATH="${install_dir}/dagmc/bin:$PATH"

####################
### OpenMC API #####
####################
if [ $build_hdf5 != "NO" ]; then \
  export HDF5_ROOT="$hdf5_libdir"; \
fi
cd ${install_dir}
git clone https://github.com/openmc-dev/openmc.git
cd openmc
git checkout tags/v0.13.0
pip install --user .

############
### PyNE ###
############

# pre-setup
export PYNE_HDF5_ARGS=""
if [ $build_hdf5 != "NO" ]; then \
  export PYNE_HDF5_ARGS="--hdf5 ${hdf5_libdir}"; \
fi
cd ${install_dir}
check_repo pyne
mkdir -p pyne
cd pyne
# clone and version
git clone -b develop --single-branch https://github.com/pyne/pyne.git pyne-repo
cd pyne-repo
if [ $1 == 'stable' ] ; then
  TAG=$(git describe --abbrev=0 --tags)
  git checkout tags/`echo ${TAG}` -b `echo ${TAG}`
fi

python setup.py install --user \
                        --moab ${install_dir}/moab \
                        --dagmc ${install_dir}/dagmc \
                        ${PYNE_HDF5_ARGS} \
                        --clean -j 3

# Adding .local/lib to $LD_LIBRARY_PATH and $LIBRARY_PATH
export LD_LIBRARY_PATH="${HOME}/.local/lib:$LD_LIBRARY_PATH"
# Adding .local//bin to $PATH
export PATH="${HOME}/.local/bin:$PATH"

# Make Pyne Nuclear Data
cd  # cd without argument will take you back to your $HOME directory
nuc_data_make

# Run tests
cd ${install_dir}/pyne/pyne-repo/tests
if [ $1 == 'stable' ] ; then
  ./travis-run-tests.sh
else
  ./ci-run-tests.sh
fi



echo " \
# Add HDF5 
if [ -z \$LD_LIBRARY_PATH ]; then 
  export LD_LIBRARY_PATH=\"${hdf5_libdir}\"
else 
  export LD_LIBRARY_PATH=\"${hdf5_libdir}:\$LD_LIBRARY_PATH\" 
fi 
# Adding MOAB/lib to \$LD_LIBRARY_PATH and \$LIBRARY_PATH
export LD_LIBRARY_PATH=\"${install_dir}/moab/lib:\$LD_LIBRARY_PATH\"

# Adding pymoab to \$PYTHONPATH
PYTHON_VERSION=\$(python -c 'import sys; print(sys.version.split('')[0][0:3])')
if [ -z \$PYTHONPATH ]; then
  export PYTHONPATH=\"${install_dir}/moab/lib/python\${PYTHON_VERSION}/site-packages\"
else
  export PYTHONPATH=\"${install_dir}/moab/lib/python\${PYTHON_VERSION}/site-packages:\$PYTHONPATH\"
fi

export LD_LIBRARY_PATH=\"${install_dir}/dagmc/lib:\$LD_LIBRARY_PATH\" 
# Adding dagmc/bin to \$PATH 
export PATH=\"${install_dir}/dagmc/bin:\$PATH\" 
" >> .bashrc
echo "Run 'source ~/.bashrc' to update environment variables. PyNE may not function correctly without doing so."