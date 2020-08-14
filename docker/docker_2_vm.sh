#!/bin/bash

# based on solution in https://gist.github.com/gunnarx/527fbc385a2e76f89609d837b6447f85
# with additions from https://github.com/iximiuz/docker-to-linux

set -e

LOOPDEVICE=`losetup -f`

docker_img=$1
docker_ver=$2
vm_name=$3

docker_info=`docker images | grep $docker_img | grep $docker_ver`
docker_size=`echo ${docker_info##* }`
docker_img_size=`echo 1000*${docker_size%??} + 1000 | bc`
docker_img_size=`echo $docker_img_size/1 | bc`

# create emtpy disk image
dd if=/dev/zero of=$vm_name.img bs=1 count=0 seek=${docker_img_size}M

# make partition onf disk image
sfdisk $vm_name.img < partition.txt

# format partition with ext3
losetup -D
losetup -o $(expr 512 \* 2048) $LOOPDEVICE $vm_name.img
mkfs.ext3 $LOOPDEVICE

# copy docker filesystem
sudo mkdir -p /mnt/$vm_name
sudo mount -t auto $LOOPDEVICE /mnt/$vm_name

container=`docker run -d $docker_img:$docker_ver /bin/bash`
docker export $container | sudo tar x -C /mnt/$vm_name
docker stop $container

#setup extlinux
extlinux --install /mnt/$vm_name/boot
cp syslinux.cfg /mnt/$vm_name/boot/syslinux.cfg

# unmount
sudo umount /mnt/$vm_name

# write master boot record
dd if=/usr/lib/syslinux/mbr/mbr.bin of=$vm_name.img bs=440 count=1 conv=notrunc

## convert to VM virtual disk
qemu-img convert -p -O vmdk $vm_name.img $vm_name.vmdk
rm $vm_name.img



