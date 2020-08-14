#!/bin/bash

# based on solution in https://gist.github.com/gunnarx/527fbc385a2e76f89609d837b6447f85
# and from https://iximiuz.com/en/posts/from-docker-container-to-bootable-linux-disk-image/
# with additions from https://github.com/iximiuz/docker-to-linux

set -e

LOOPDEVICE=`losetup -f`

docker_img=$1
docker_ver=$2
vm_name=$3

banner=" ############# "

echo -e "$banner Creating VM image $vm_name.vmdk from $docker_img:$docker_ver"

docker_info=`docker images | grep $docker_img | grep $docker_ver`
docker_size=`echo ${docker_info##* }`
docker_img_size=`echo 1000*${docker_size%??} + 1000 | bc`
docker_img_size=`echo $docker_img_size/1 | bc`

echo -e "$banner 1/9 - create disk image '$vm_name.img' of size $docker_img_size MB"
dd if=/dev/zero of=$vm_name.img bs=1 count=0 seek=${docker_img_size}M

echo -e "$banner 2/9 - create partition on disk image"
sfdisk $vm_name.img < partition.txt

echo -e "$banner 3/9 - format partition as ext3 filesystem"
sudo losetup -D
sudo losetup -o $(expr 512 \* 2048) $LOOPDEVICE $vm_name.img
sudo mkfs.ext3 $LOOPDEVICE

echo -e "$banner 4/9 - mount disk image in '/mnt/$vm_name'"
sudo mkdir -p /mnt/$vm_name
sudo mount -t auto $LOOPDEVICE /mnt/$vm_name

echo -e "$banner 5/9 - run docker image '$docker_img:$docker_ver' and export filesystem"
container=`docker run -d $docker_img:$docker_ver /bin/bash`
docker export $container | sudo tar x -C /mnt/$vm_name
docker stop $container

echo -e "$banner 6/9 - install booloader"
sudo extlinux --install /mnt/$vm_name/boot
sudo cp syslinux.cfg /mnt/$vm_name/boot/syslinux.cfg

echo -e "$banner 7/9 - unmount filesystem"
sudo umount /mnt/$vm_name

echo -e "$banner 8/9 - copy Master Boot Record"
dd if=/usr/lib/syslinux/mbr/mbr.bin of=$vm_name.img bs=440 count=1 conv=notrunc

echo -e "$banner 9/9 - convert disk image to VMDK virtual machine drive"
qemu-img convert -p -O vmdk $vm_name.img $vm_name.vmdk
rm $vm_name.img



