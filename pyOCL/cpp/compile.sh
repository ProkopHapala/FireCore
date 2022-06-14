#! /bin/bash

#LPATH2=`pwd`
#LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$LPATH:$LPATH2

LPATH=/usr/local/lib64/
LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$LPATH
#echo  $LD_LIBRARY_PATH
#export $LD_LIBRARY_PATH

wb="../cpp_build"
wd=`pwd`

FFLAGS="-Werror -Wno-write-strings -Wno-reorder -Wno-deprecated-declarations"
g++ $FFLAGS -c -Wall -fpic -o $wb/libOCLfft.o libOCLfft.cpp

cd $wb
g++ -shared -o libOCLfft.so libOCLfft.o -lOpenCL -lclFFT
#gcc -shared -o libfft3d.so lib_fft3d.o

#echo "========= C  use_fft3d.x "
#g++ use_fft3d.cpp -o use_fft3d.x -L./  -lOpenCL -lclFFT -Bdynamic -lfft3d
#g++ use_OCLfft.c -o use_OCLfft.x -L./ -lfft3d
#./use_fft3d.x