#! /bin/bash


function test_img_name()
{
    if [[    "$1" != *"pyne/ubuntu_18.04_py2_pyne-deps"* \
          && "$1" != *"pyne/ubuntu_18.04_py2_dagmc_pyne-deps"* \
          && "$1" != *"pyne/ubuntu_18.04_py2_pymoab_pyne-deps"* \
          && "$1" != *"pyne/ubuntu_18.04_py2_dagmc_pymoab_pyne-deps"* \
          && "$1" != *"pyne/ubuntu_18.04_py3_pyne-deps"* \
          && "$1" != *"pyne/ubuntu_18.04_py3_dagmc_pyne-deps"* \
          && "$1" != *"pyne/ubuntu_18.04_py3_pymoab_pyne-deps"* \
          && "$1" != *"pyne/ubuntu_18.04_py3_dagmc_pymoab_pyne-deps"* \
          && "$1" != *"pyne/ubuntu_18.04_py3_dagmc_pymoab_openmc_pyne-deps"* ]];then
    echo Unknow docker image !
    echo PyNE supported docker images are:
    echo "pyne/ubuntu_18.04_py2_pyne-deps"
    echo "pyne/ubuntu_18.04_py2_dagmc_pyne-deps"
    echo "pyne/ubuntu_18.04_py2_pymoab_pyne-deps"
    echo "pyne/ubuntu_18.04_py2_dagmc_pymoab_pyne-deps"
    echo "pyne/ubuntu_18.04_py3_pyne-deps"
    echo "pyne/ubuntu_18.04_py3_dagmc_pyne-deps"
    echo "pyne/ubuntu_18.04_py3_pymoab_pyne-deps"
    echo "pyne/ubuntu_18.04_py3_dagmc_pymoab_pyne-deps"
    echo "pyne/ubuntu_18.04_py3_dagmc_pymoab_openmc_pyne-deps"
  fi
}

test_img_name $1

docker pull $1
docker run -v $PWD:/pyne -i -w "/pyne" -t $1 /bin/bash  -c "python setup.py install --prefix=~/.local -j18 --clean"
commit_id=$(docker ps -lq | tail -n1)
docker commit ${commit_id} $1
docker run -v $PWD:/pyne -i -w "/pyne" -t $1 /bin/bash