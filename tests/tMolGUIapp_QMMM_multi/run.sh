name=MolGUIapp_QMMM_multi
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

lscpu

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

#./$name -m 10         -x common_resources/xyz/H2O
#./$name -t 2 -m 10    -x common_resources/xyz/pyridine
#./$name -t 1 -m 10    -x common_resources/xyz/pyridine
#./$name -t 1 -m $nsys -x common_resources/xyz/pyridine
#./$name -t 1 -m $nsys -x common_resources/xyz/pyridine
#./$name -t 1 -m $nsys -x common_resources/xyz/PTCDA
#./$name -t 2 -m $nsys -x common_resources/xyz/PTCDA
#./$name -t 3 -m $nsys -x common_resources/xyz/PTCDA

#./$name -m 10 -x common_resources/xyz/H2O
#./$name -m 10 -x common_resources/xyz/HCOOH
#./$name -m 10 -x common_resources/xyz/formic_dimer
#./$name -m 10 -x common_resources/xyz/pyridine
#./$name -m 10 -x common_resources/xyz/propandiol

# ====== Polymers

#./$name -m $nsys -x common_resources/xyz/polydiacetylene
#./$name -m $nsys -x common_resources/xyz/polydiacetylene     -subs 4,common_resources/xyz/-COOH.xyz
#./$name -m $nsys -x common_resources/xyz/polydiacetylene_OH
#./$name -m $nsys -x common_resources/xyz/polymer-2_new

#./$name -m 1     -x common_resources/xyz/polymer-2_new
#./$name -m 4     -x common_resources/xyz/polymer-2_new
#./$name -m 1     -x common_resources/xyz/polymer-2_new_crash
#./$name -m 1     -x common_resources/xyz/polymer-2_new_crash_33


#./$name -m 1   -x common_resources/xyz/polymer-2_new
#./$name -m 5   -x common_resources/xyz/polymer-2_new
#./$name -m 40  -x common_resources/xyz/polymer-2_new
#./$name -m 40  -x common_resources/xyz/polymer-2_new  -c 10
#./$name -m 2   -x common_resources/xyz/polymer-2_new
#./$name -m 1   -x common_resources/xyz/polymer-2_new

#./$name -m 10 -t 1     -x common_resources/xyz/polymer-2_new
#./$name -m $nsys -t 1  -x common_resources/xyz/polymer-2_new

# ------- lattice scan

#./$name -x BB.HNH-h.NHO-hh                          -iParalel 1 -perframe 10
#./$name -x BB.HNH-h.NHO-hh -b BB.HNH-h.NHO-hh.hbonds -iParalel 1 -perframe 50 
#./$name -m 40 -x BB.HNH-h.NHO-hh -b BB.HNH-h.NHO-hh.hbonds -iParalel 1 -perframe 50 -pop lattice_scan_1d_all.xyz
#./$name  -m 30 -x BB.HNH-h.NHO-hh -b BB.HNH-h.NHO-hh.hbonds -iParalel 1 -perframe 50 -pop lattice_scan_1d_all.xyz -latscan 40,0 0.2,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0
#./$name -x BB.HNH-h.NHO-hh -b BB.HNH-h.NHO-hh.hbonds -iParalel 1
#./$name -x BB.HNH-h.NHO-hh -b BB.HNH-h.NHO-hh.hbonds -dlvec -2.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0
#./$name -x BB.HNH-h.NHO-hh -b BB.HNH-h.NHO-hh.hbonds -dlvec -0.5,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0
#./$name  -x BB.HNH-h.NHO-hh -b BB.HNH-h.NHO-hh.hbonds -dlvec -0.5,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0 -latscan 10,10 0.2,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0

# ====== Small Molecules On Substrate

#valgrind --log-file="valgrind.log" --leak-check=yes ./$name -m 10 -x common_resources/xyz/H2O -n     -g common_resources/xyz/NaCl_1x1_L2
#valgrind --leak-check=yes ./$name -m 10 -x common_resources/xyz/H2O -n     -g common_resources/xyz/NaCl_1x1_L2

#./$name -m 10    -x common_resources/xyz/H2O           -g common_resources/xyz/NaCl_1x1_L2
#./$name -m $nsys -x common_resources/xyz/pyridine      -g common_resources/xyz/NaCl_1x1_L2
#./$name -m $nsys -x common_resources/xyz/pyridine  110 -g common_resources/xyz/NaCl_1x1_L2
./$name -m $nsys -x common_resources/xyz/PTCDA         -g common_resources/xyz/NaCl_1x1_L2
#./$name -m $nsys -x common_resources/xyz/PTCDA_SAM     -g common_resources/xyz/NaCl_1x1_L2
#./$name -m $nsys -x common_resources/xyz/BPBA          -g common_resources/xyz/NaCl_1x1_L2
#./$name -m $nsys -x common_resources/xyz/BPBA          -g common_resources/xyz/NaCl_1x1_L2 -e

#./$name -m 20 -x common_resources/xyz/PTCDA_SAM  -c 59      -g common_resources/xyz/NaCl_1x1_L2   -ManipulAnim
#./$name -m 20 -x common_resources/xyz/PTCDA_SAM  -c 59      -g common_resources/xyz/NaCl_1x1_L2

# ====== Polymers On Substrate

#./$name -m $nsys -x common_resources/xyz/polydiacetylene           -g common_resources/xyz/NaCl_1x1_L2
#./$name -m $nsys -x common_resources/xyz/polydiacetylene           -g common_resources/xyz/NaCl_1x1_L2 -subs 4,common_resources/-COOH.xyz
#./$name -m $nsys -x common_resources/xyz/polydiacetylene    -n 221 -g common_resources/xyz/NaCl_1x1_L2
#./$name -m $nsys -x common_resources/xyz/polydiacetylene    -n 241 -g common_resources/xyz/NaCl_1x1_L2
#./$name -m $nsys -x common_resources/xyz/polydiacetylene    -n 241 -g common_resources/xyz/NaCl_1x1_L2 -ng 4,-2,2,3
#./$name -m $nsys -x common_resources/xyz/polydiacetylene    -n 141 -g common_resources/xyz/NaCl_1x1_L2 -ng 2,-1,2,3
#./$name -m $nsys -x common_resources/xyz/polydiacetylene    -n 141 -g common_resources/xyz/NaCl_1x1_L2 -ng 2,-1,3,4
#./$name -m $nsys -x common_resources/xyz/polydiacetylene    -n 141 -g common_resources/xyz/NaCl_1x1_L2 -subs 4,common_resources/xyz/-COOH.xyz  -ng 2,-1,3,4 -q 0.05

#./$name -m $nsys -x common_resources/xyz/polydiacetylene_OH -n 141 -g common_resources/xyz/NaCl_1x1_L2 -ng 2,-1,3,4
#./$name -m $nsys -x common_resources/xyz/polydiacetylene_OH        -g common_resources/xyz/NaCl_1x1_L2
#./$name -m $nsys -x common_resources/xyz/polydiacetylene_OH        -g common_resources/xyz/NaCl_1x1_L2

#./$name -m $nsys -x common_resources/xyz/polymer-2_new              -g common_resources/xyz/NaCl_1x1_L2
#./$name -m 10 -x common_resources/xyz/polymer-2_new                 -g common_resources/xyz/NaCl_1x1_L2   | tee out.log
#./$name -m 1 -x common_resources/xyz/polymer-2_new                  -g common_resources/xyz/NaCl_1x1_L2
#./$name -m 2 -x common_resources/xyz/polymer-2_new                  -g common_resources/xyz/NaCl_1x1_L2

#./$name -m 40 -x common_resources/xyz/polymer-2_new                  -g common_resources/NaCl_1x1_L2


#valgrind --log-file="valgrind.log" --leak-check=yes ./$name -x common_resources/xyz/H2O
#valgrind --leak-check=yes --log-file="valgrind.log" ./$name -x common_resources/xyz/HCOOH -m 10
#valgrind --leak-check=yes --track-origins=yes  --log-file="valgrind.log" ./$name -x common_resources/xyz/HCOOH -m 10


#grep GPU_GFF_z out.log | cut -c 11-  > GPU_makeGridFF.log
