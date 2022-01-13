#!/bin/bash

#fname=$1

#folder="t01_H2/"
folder="t02_H4/"
#folder="t02_CH4/"

path1="./tests/"
path2="../Fireball-progs/TESTS/"

#pre="H_2c_"
pre="H_3c_"
post=".log"
#fnames="H_2c_t H_2c_vna H_2c_vxc H_2c_vca H_2c_vxcca H_2c_ewaldLR H_2c_ewaldSR"
fnames="t vna vxc vca vxcca ewaldLR ewaldSR"

for fname in $fnames; do
    echo ""
    echo "=========== ", $fname
    f1=$path1$folder$pre$fname.log
    f2=$path2$folder$pre$fname.log
    #echo $fname
    #echo $f1
    #echo $f2
    diff $f1 $f2
done