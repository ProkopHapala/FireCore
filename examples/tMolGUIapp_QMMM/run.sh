name=MolGUIapp_QMMM
dir=../../cpp/Build/apps/MolecularEditor
ln -s ../../cpp/common_resources data
ln -s ../../cpp/common_resources common_resources 

## ---- Multiprocesing
#ncpu=`nproc`
#ncpu=$(($ncpu - 1))     # let one CPU free for user interaction
#echo "compile using ncpu="$ncpu
#OMP_NUM_THREADS=$ncpu
#export OMP_NUM_THREADS


export OMP_NUM_THREADS=1


# ---- Compilation
wd=`pwd`
cd $dir
pwd
rm $name
make -j$ncpu $name   # 2>$wd/compile_err.log
cd $wd
ln -s $dir/$name .


#export MKL_VERBOSE=1

# ------- asan (Memory Sanitizer)
#LD_PRELOAD=$(g++ -print-file-name=libasan.so)
#echo   $LD_PRELOAD
#export LD_PRELOAD


# ---- Run

#rm *.bin *.xsf

# ====== Small Molecules On Substrate

#./$name -x common_resources/xyz/H2O               -g common_resources/xyz/NaCl_1x1_L2
#./$name -x common_resources/xyz/pyridine         -g common_resources/xyz/NaCl_1x1_L2
#./$name -x common_resources/xyz/CH2O -g common_resources/xyz/NaCl_1x1_L2
./$name -x common_resources/xyz/PTCDA -g common_resources/xyz/NaCl_1x1_L2

