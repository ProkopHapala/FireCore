#!/bin/bash

dfile="assemble_F.log"
fname=${1-$dfile}

#dfolder="t01_H2/"
#dfolder="t02_H4/"
dfolder="t02_CH4/"

folder=${2-$dfolder}


path1="./tests/"
path2="../Fireball-progs/TESTS/"

f1=$path1$folder$fname
f2=$path2$folder$fname



echo "#### DIFF:"
#echo $fname
echo $f1
echo $f2
echo ">>>>"
diff $f1 $f2
echo "<<<<"