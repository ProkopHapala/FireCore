#!/bin/bash

set -e   # stop on error  ; see  https://stackoverflow.com/questions/3474526/stop-on-first-error

dr=`pwd`


cd cpp/Build
make
cd $dr

python pyBall/gen_makefile.py VERY_OPT
#python pyBall/gen_makefile.py DEBUG
#ln -s Makefile ../build/Makefile

#mkdir build  || true
#cd    build

mkdir build_opt || true
cd    build_opt


#ln -s ../Makefile      || true
cp ../Makefile . || true
#make | tee make.log
rm make.log                        || true
echo "==== START make "
pwd
make 2>&1 | tee -a make.log
cd $dr
#ln -s build/make.log .             || true