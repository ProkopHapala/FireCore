name=FitREQ_lib
dir=../../cpp/Build/libs/Molecular
ln -s ../../cpp/common_resources common_resources 

wd=`pwd`
cd $dir
pwd
rm lib$name
make -j4 $name
cd $wd

python3 run.py


