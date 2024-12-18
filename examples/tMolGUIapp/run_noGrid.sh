name=MolGUIapp
dir=../../cpp/Build/apps/MolecularEditor
ln -s ../../cpp/common_resources data
ln -s ../../cpp/common_resources common_resources 

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

ln -s $dir/$name .



# ------- asan (Memory Sanitizer)
#LD_PRELOAD=$(g++ -print-file-name=libasan.so)
#echo   $LD_PRELOAD
#export LD_PRELOAD


# ---- Run

#rm *.bin *.xsf

# ====== Small Molecules

#./$name -x common_resources/H2O
#./$name -x common_resources/HCOOH
#./$name -x common_resources/formic_dimer
#./$name -x common_resources/pyridine
#./$name -x common_resources/propandiol

# ====== Polymers

#./$name -x common_resources/polydiacetylene
#./$name -x common_resources/polydiacetylene_OH
#./$name -x common_resources/polymer-2
./$name -x common_resources/polymer-2_new






#valgrind --log-file="valgrind.log" --leak-check=yes ./$name -x common_resources/H2O
