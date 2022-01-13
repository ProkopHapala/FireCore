#!/bin/bash

fname=$1

#folder="t01_H2/"
folder="t02_CH4/"

path1="./tests/"
path2="../Fireball-progs/TESTS/"

f1=$path1$folder$fname
f2=$path2$folder$fname

#echo $fname
echo $f1
echo $f2

diff $f1 $f2