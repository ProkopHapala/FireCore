name=EFFapp.x
dir=../../cpp/Build/apps/EFF
#ln -s ../../cpp/common_resources data
ln -s ../../cpp/common_resources common_resources 
ln -s ../../cpp/sketches_SDL/Molecular/data data
ln -s ../../cpp/sketches_SDL/Molecular/data inputs



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


# ------- asan (Memory Sanitizer)
#LD_PRELOAD=$(g++ -print-file-name=libasan.so)
#echo   $LD_PRELOAD
#export LD_PRELOAD


# ---- Run

#rm *.bin *.xsf

# ====== Small Molecules


#./$name -f data/e.fgo
#./$name -f data/e2_singlet_far.fgo
#./$name -f data/e2_1g_2o_singlet.fgo
#./$name -f data/e2_1g_2o_triplet.fgo
./$name -f data/H2_eFF.fgo
#./$name -f data/H2O.fgo
#./$name -f data/C2H4.fgo




