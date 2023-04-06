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
#./$name -x common_resources/polydiacetylene     -subs 4,common_resources/-COOH.xyz
#./$name -x common_resources/polydiacetylene_OH
#./$name -x common_resources/polymer-2_new

# ====== Small Molecules On Substrate

./$name -x common_resources/pyridine         -g common_resources/NaCl_1x1_L2
#./$name -x common_resources/pyridine -n 110 -g common_resources/NaCl_1x1_L2

# ====== Polymers On Substrate

#./$name -x common_resources/polydiacetylene           -g common_resources/NaCl_1x1_L2
#./$name -x common_resources/polydiacetylene           -g common_resources/NaCl_1x1_L2 -subs 4,common_resources/-COOH.xyz
#./$name -x common_resources/polydiacetylene    -n 221 -g common_resources/NaCl_1x1_L2
#./$name -x common_resources/polydiacetylene    -n 241 -g common_resources/NaCl_1x1_L2
#./$name -x common_resources/polydiacetylene    -n 241 -g common_resources/NaCl_1x1_L2 -ng 4,-2,2,3
#./$name -x common_resources/polydiacetylene    -n 141 -g common_resources/NaCl_1x1_L2 -ng 2,-1,2,3
#./$name -x common_resources/polydiacetylene    -n 141 -g common_resources/NaCl_1x1_L2 -ng 2,-1,3,4
#./$name -x common_resources/polydiacetylene    -n 141 -g common_resources/NaCl_1x1_L2 -subs 4,common_resources/-COOH.xyz  -ng 2,-1,3,4 -q 0.05

#./$name -x common_resources/polydiacetylene_OH -n 141 -g common_resources/NaCl_1x1_L2 -ng 2,-1,3,4
#./$name -x common_resources/polydiacetylene_OH        -g common_resources/NaCl_1x1_L2
#./$name -x common_resources/polydiacetylene_OH        -g common_resources/NaCl_1x1_L2

#./$name -x common_resources/polymer-2          -n 221 -g common_resources/NaCl_sym-center


#valgrind --log-file="valgrind.log" --leak-check=yes ./$name -x common_resources/H2O


