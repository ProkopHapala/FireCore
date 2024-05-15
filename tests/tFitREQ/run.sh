name=FitREQ_lib
dir=../../cpp/Build/libs/Molecular
ln -sf ../../cpp/common_resources common_resources
ln -sf ../../cpp/common_resources data


wd=`pwd`
cd $dir
#pwd
#rm -f lib$name.so
make -j4 $name
cd $wd

# ------- asan (Memory Sanitizer)
#LD_PRELOAD=$(g++ -print-file-name=libasan.so)
#echo   $LD_PRELOAD
#export LD_PRELOAD

#> FitREQ_debug.xyz
#python3 run.py
#python3 fit_manual.py
#python3 fit_manual_2.py
#python3 fit_manual_2d.py
#python3 fit_manual_OH.py

#python3 fit_manual_samp.py
python3 fit_manual_test.py


