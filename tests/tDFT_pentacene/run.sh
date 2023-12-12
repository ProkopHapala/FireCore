# /bin/bash

wd=`pwd`
cd ../cpp/
./compile.sh
cd $wd

#LPATH=/usr/local/lib64/
LPATH2=usr/lib/x86_64-linux-gnu/
#LPATH3=`pwd`
#LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$LPATH:$LPATH2:$LPATH3
LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$LPATH2
echo  $LD_LIBRARY_PATH
export $LD_LIBRARY_PATH

python3 run.py
