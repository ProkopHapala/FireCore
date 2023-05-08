name=MolGUIapp_multi
dir=../../cpp/Build/apps_OCL/MolecularEditorOCL
ln -s ../../cpp/common_resources data
ln -s ../../cpp/common_resources common_resources 

wd=`pwd`
cd $dir
pwd
rm $name
make -j4 $name
cd $wd

ln -s $dir/$name .

# ---- Multiprocesing
ncpu=`nproc`
ncpu=$(($ncpu - 1))     # let one CPU free for user interaction
echo "compile using ncpu="$ncpu
OMP_NUM_THREADS=$ncpu
export OMP_NUM_THREADS

rm DEBUG*

#nsys=10
nsys=50
#nsys=200


# ====== Small Molecules

#./$name -m 10 -x common_resources/H2O
#./$name -t 2 -m 10 -x common_resources/pyridine
#./$name -t 1 -m 10 -x common_resources/pyridine
#./$name -t 1 -m $nsys -x common_resources/pyridine
#./$name -t 1 -m $nsys -x common_resources/pyridine
#./$name -t 1 -m $nsys -x common_resources/PTCDA
#./$name -t 2 -m $nsys -x common_resources/PTCDA
#./$name -t 3 -m $nsys -x common_resources/PTCDA


#./$name -m 10 -x common_resources/H2O
#./$name -m 10 -x common_resources/HCOOH
#./$name -m 10 -x common_resources/formic_dimer
#./$name -m 10 -x common_resources/pyridine
#./$name -m 10 -x common_resources/propandiol

# ====== Polymers

#./$name -m $nsys -x common_resources/polydiacetylene
#./$name -m $nsys -x common_resources/polydiacetylene     -subs 4,common_resources/-COOH.xyz
#./$name -m $nsys -x common_resources/polydiacetylene_OH
#./$name -m $nsys -x common_resources/polymer-2_new

#./$name -m 1    -x common_resources/polymer-2_new
#./$name -m 4    -x common_resources/polymer-2_new
#./$name -m 1    -x common_resources/polymer-2_new_crash
#./$name -m 1    -x common_resources/polymer-2_new_crash_33


#./$name -m 1  -x common_resources/polymer-2_new
#./$name -m 5  -x common_resources/polymer-2_new
./$name -m 40  -x common_resources/polymer-2_new
#./$name -m 40 -c 10 -x common_resources/polymer-2_new
#./$name -m 2     -x common_resources/polymer-2_new
#./$name -m 1     -x common_resources/polymer-2_new

#./$name -m 10 -t 1  -x common_resources/polymer-2_new
#./$name -m $nsys -t 1  -x common_resources/polymer-2_new

# ====== Small Molecules On Substrate

#valgrind --log-file="valgrind.log" --leak-check=yes ./$name -m 10 -x common_resources/H2O -n     -g common_resources/NaCl_1x1_L2
#valgrind --leak-check=yes ./$name -m 10 -x common_resources/H2O -n     -g common_resources/NaCl_1x1_L2

#./$name -m 10    -x common_resources/H2O           -g common_resources/NaCl_1x1_L2
#./$name -m $nsys -x common_resources/pyridine      -g common_resources/NaCl_1x1_L2
#./$name -m $nsys -x common_resources/pyridine  110 -g common_resources/NaCl_1x1_L2
#./$name -m $nsys -x common_resources/PTCDA   -g common_resources/NaCl_1x1_L2
#./$name -m $nsys -x common_resources/PTCDA_SAM   -g common_resources/NaCl_1x1_L2
#./$name -m $nsys -x common_resources/BPBA   -g common_resources/NaCl_1x1_L2


# ====== Polymers On Substrate

#./$name -m $nsys -x common_resources/polydiacetylene           -g common_resources/NaCl_1x1_L2
#./$name -m $nsys -x common_resources/polydiacetylene           -g common_resources/NaCl_1x1_L2 -subs 4,common_resources/-COOH.xyz
#./$name -m $nsys -x common_resources/polydiacetylene    -n 221 -g common_resources/NaCl_1x1_L2
#./$name -m $nsys -x common_resources/polydiacetylene    -n 241 -g common_resources/NaCl_1x1_L2
#./$name -m $nsys -x common_resources/polydiacetylene    -n 241 -g common_resources/NaCl_1x1_L2 -ng 4,-2,2,3
#./$name -m $nsys -x common_resources/polydiacetylene    -n 141 -g common_resources/NaCl_1x1_L2 -ng 2,-1,2,3
#./$name -m $nsys -x common_resources/polydiacetylene    -n 141 -g common_resources/NaCl_1x1_L2 -ng 2,-1,3,4
#./$name -m $nsys -x common_resources/polydiacetylene    -n 141 -g common_resources/NaCl_1x1_L2 -subs 4,common_resources/-COOH.xyz  -ng 2,-1,3,4 -q 0.05

#./$name -m $nsys -x common_resources/polydiacetylene_OH -n 141 -g common_resources/NaCl_1x1_L2 -ng 2,-1,3,4
#./$name -m $nsys -x common_resources/polydiacetylene_OH        -g common_resources/NaCl_1x1_L2
#./$name -m $nsys -x common_resources/polydiacetylene_OH        -g common_resources/NaCl_1x1_L2

#./$name -m $nsys -x common_resources/polymer-2_new              -g common_resources/NaCl_1x1_L2
#./$name -m 10 -x common_resources/polymer-2_new                 -g common_resources/NaCl_1x1_L2   | tee out.log
#./$name -m 1 -x common_resources/polymer-2_new                  -g common_resources/NaCl_1x1_L2
#./$name -m 2 -x common_resources/polymer-2_new                  -g common_resources/NaCl_1x1_L2

#./$name -m 40 -x common_resources/polymer-2_new                  -g common_resources/NaCl_1x1_L2


#valgrind --log-file="valgrind.log" --leak-check=yes ./$name -x common_resources/H2O
#valgrind --leak-check=yes --log-file="valgrind.log" ./$name -x common_resources/HCOOH -m 10
#valgrind --leak-check=yes --track-origins=yes  --log-file="valgrind.log" ./$name -x common_resources/HCOOH -m 10


grep GPU_GFF_z out.log | cut -c 11-  > GPU_makeGridFF.log
