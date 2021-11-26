#!/bin/bash

./make.sh

echo ###############################

bkdir=`pwd`
cd tests/t01_H2
./run.sh
cd $bkdir