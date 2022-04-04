# /bin/bash

LPATH=/usr/local/lib64/
LPATH2=`pwd`
LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$LPATH:$LPATH2
echo  $LD_LIBRARY_PATH
export $LD_LIBRARY_PATH

FFLAGS="-Werror -Wno-write-strings -Wno-reorder -Wno-deprecated-declarations"

g++ $FFLAGS -c -Wall -fpic -o libOCLfft.o libOCLfft.cpp
g++ -shared -o libOCLfft.so libOCLfft.o -lOpenCL -lclFFT
#gcc -shared -o libfft3d.so lib_fft3d.o

#echo "========= Python py_clfft.py "
python3 py_oclfft.py

#echo "========= C  use_fft3d.x "
#g++ use_fft3d.cpp -o use_fft3d.x -L./  -lOpenCL -lclFFT -Bdynamic -lfft3d
#g++ use_OCLfft.c -o use_OCLfft.x -L./ -lfft3d
#./use_fft3d.x