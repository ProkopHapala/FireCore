# /bin/bash

LPATH=/usr/local/lib64/
LPATH2=`pwd`
LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$LPATH:$LPATH2
echo  $LD_LIBRARY_PATH
export $LD_LIBRARY_PATH

#gcc fft3d.c -o fft3d.x -lOpenCL -lclFFT
#./fft3d.x

gcc -c -Wall -fpic -o lib_fft3d.o lib_fft3d.c
gcc -shared -o libfft3d.so lib_fft3d.o -lOpenCL -lclFFT
#gcc -shared -o libfft3d.so lib_fft3d.o

echo "========= C  use_fft3d.x "
#gcc use_fft3d.cpp -o use_fft3d.x -L./  -lOpenCL -lclFFT -Bdynamic -lfft3d
gcc use_fft3d.c -o use_fft3d.x -L./ -lfft3d
#./use_fft3d.x

echo "========= Python py_clfft.py "

python3 py_clfft.py