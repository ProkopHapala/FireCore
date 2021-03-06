#!/bin/bash

set -e   # stop on error  ; see  https://stackoverflow.com/questions/3474526/stop-on-first-error

dr=`pwd`
python pyBall/gen_makefile.py
#ln -s Makefile ../build/Makefile


mkdir build  || true
cd    build
ln -s ../Makefile      || true
#make | tee make.log
rm make.log                        || true
echo "==== START make "
pwd
make 2>&1 | tee -a make.log
cd $dr
ln -s build/make.log .             || true