#!/bin/bash

# based on solution in https://gist.github.com/gunnarx/527fbc385a2e76f89609d837b6447f85

docker_img=$1
docker_ver=$2
vm_name=$3

docker_info=`docker images | grep $docker_img | grep $docker_ver`
docker_size=`echo ${docker_info##* }`
docker_img_size=`echo 1000*${docker_size%??} + 1000 | bc`
docker_img_size=`echo $docker_img_size/1 | bc`

dd if=/dev/zero of=$vm_name.img bs=1 count=0 seek=${docker_img_size}M
mkfs.ext2 -F $vm_name.img
sudo mkdir -p /mnt/$vm_name
sudo mount -o loop $vm_name.img /mnt/$vm_name

container=`docker run -d $docker_img:$docker_ver /bin/bash`
docker export $container | sudo tar x -C /mnt/$vm_name
docker stop $container

sudo umount /mnt/$vm_name
qemu-img convert -p -O vmdk $vm_name.img $vm_name.vmdk
rm $vm_name.img



