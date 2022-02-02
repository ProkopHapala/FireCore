#!/bin/bash

./make.sh

echo ###############################

bkdir=`pwd`
#cd tests/t01_H2
cd tests/t02_CH4
./run.sh
cd $bkdir

#./diff.sh assemble_F.log t02_CH4/