name=EFFapp.x
dir=../../Build/apps/EFF
ln -s ../../common_resources data
ln -s ../../common_resources common_resources 

# ---- Multiprocesing
ncpu=`nproc`
ncpu=$(($ncpu - 1))     # let one CPU free for user interaction
echo "compile using ncpu="$ncpu
OMP_NUM_THREADS=$ncpu
export OMP_NUM_THREADS

# ---- Compilation
wd=`pwd`
cd $dir
pwd
rm $name
make -j$ncpu $name   # 2>$wd/compile_err.log
cd $wd

rm $name
ln -s $dir/$name .

ln -s ../../cpp/sketches_SDL/Molecular/data data
ln -s ../../cpp/sketches_SDL/Molecular/data inputs

# ------- asan (Memory Sanitizer)
#LD_PRELOAD=$(g++ -print-file-name=libasan.so)
#echo   $LD_PRELOAD
#export LD_PRELOAD


# ---- Run

#rm *.bin *.xsf

# ====== Small Molecules

./$name # -x data/C2H4




