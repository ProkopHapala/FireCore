name=FitREQ_lib
dir=../../cpp/Build/libs/Molecular
ln -s ../../cpp/common_resources common_resources
ln -s ../../cpp/common_resources data


wd=`pwd`
cd $dir
pwd
rm lib$name
make -j4 $name
cd $wd

# ------- asan (Memory Sanitizer)
LD_PRELOAD=$(g++ -print-file-name=libasan.so)
echo   $LD_PRELOAD
export LD_PRELOAD

> FitREQ_debug.xyz
python3 run.py


