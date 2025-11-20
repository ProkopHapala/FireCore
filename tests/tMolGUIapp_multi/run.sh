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
rm gopt.xyz
rm minima.dat
touch minima.dat

# ====== Small Molecules

#./$name      -m 10    -x common_resources/xyz/H2O
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

#./$name -m 1 -x common_resources/xyz/CH4  
#./$name -m 1 -x common_resources/xyz/CH4 -ex2 
# ./$name -m 10 -x common_resources/xyz/CH4 -uff -verb 4 -dt 0.005

#./$name       -x common_resources/xyz/nHexadecan_dicarboxylic  -T 100 0.01
#./$name       -x common_resources/xyz/nHexadecan_dicarboxylic -b nHexadecan_dicarboxylic.cons -T 100 0.01
#./$name       -x common_resources/xyz/nHexadecan_dicarboxylic -b nHexadecan_dicarboxylic.cons -perframe 1000 -gopt 1000,1000 0.25,1.0
#./$name -m 300 -x common_resources/xyz/butandiol-2 -iParalel 3 -perframe 10 -gopt 100,1000 0.25,1.0 -verb 0  -T 1000 0.02 -Ftol 1e-4
#./$name -m 30 -x common_resources/xyz/nHexadecan_dicarboxylic -drive nHexadecan_dicarboxylic.cons -iParalel 1 -perframe 100 -gopt 1000,1000 0.25,1.0 -verb 0  -T 1000 0.02

#./$name -m 60  -x common_resources/xyz/nHexadecan_dicarboxylic -drive nHexadecan_dicarboxylic.cons -iParalel 1 -perframe 100 -gopt 1000,1000 0.25,1.0 -verb 0  -T 1000 0.02
#./$name -m 120 -x common_resources/xyz/nHexadecan_dicarboxylic -drive nHexadecan_dicarboxylic.cons -iParalel 2 -perframe 100 -gopt 1000,1000 0.25,1.0 -verb 0  -T 1000 0.02
#./$name -m 60  -x common_resources/xyz/nHexadecan_dicarboxylic                                     -iParalel 1 -perframe 100 -gopt 1000,1000 0.0,0.0 -T 10000 0.01

#./$name -m 60  -x common_resources/xyz/nHexadecan_dicarboxylic -iParalel 1 -perframe 100 -T 1000 0.02
#./$name -m 120 -x common_resources/xyz/nHexadecan_dicarboxylic -iParalel 2 -perframe 100 -T 1000 0.02


# ====== Polymers

#./$name -m $nsys -x common_resources/xyz/polydiacetylene
#./$name -m $nsys -x common_resources/xyz/polydiacetylene     -subs 4,common_resources/xyz/-COOH.xyz
#./$name -m $nsys -x common_resources/xyz/polydiacetylene_OH
#./$name -m $nsys -x common_resources/xyz/polymer-2_new

#./$name -m 1    -x common_resources/xyz/polymer-2_new
#./$name -m 4    -x common_resources/xyz/polymer-2_new
#./$name -m 1    -x common_resources/xyz/polymer-2_new_crash
#./$name -m 1    -x common_resources/xyz/polymer-2_new_crash_33


#./$name -m 30    -x common_resources/xyz/polymer-2_new    -iParalel 1   -perframe 100  -gopt 1000,1000 0.25,1.0    -T 1000 0.02

#./$name -m 1   -x common_resources/xyz/polymer-2_new
#./$name -m 5   -x common_resources/xyz/polymer-2_new
#./$name -m 40  -x common_resources/xyz/polymer-2_new
#./$name -m 40  -x common_resources/xyz/polymer-2_new  -c 10
#./$name -m 2   -x common_resources/xyz/polymer-2_new
#./$name -m 1   -x common_resources/xyz/polymer-2_new

#./$name -m 10    -x common_resources/xyz/polymer-2_new -t 1 
#./$name -m $nsys -x common_resources/xyz/polymer-2_new -t 1 

# ------- lattice scan

#./$name -x BB.HNH-h.NHO-hh                          -iParalel 1 -perframe 10
#./$name -x BB.HNH-h.NHO-hh -b BB.HNH-h.NHO-hh.hbonds -iParalel 1 -perframe 50 
#./$name -m 40 -x BB.HNH-h.NHO-hh -b BB.HNH-h.NHO-hh.hbonds -iParalel 1 -perframe 50 -pop lattice_scan_1d_all.xyz
#./$name  -m 30 -x BB.HNH-h.NHO-hh -b BB.HNH-h.NHO-hh.hbonds -iParalel 1 -perframe 50 -pop lattice_scan_1d_all.xyz -latscan 40,0 0.2,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0
#./$name -x BB.HNH-h.NHO-hh -b BB.HNH-h.NHO-hh.hbonds -iParalel 1
#./$name -x BB.HNH-h.NHO-hh -b BB.HNH-h.NHO-hh.hbonds -iParalel 1 -dlvec -2.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0
#./$name -x BB.HNH-h.NHO-hh -b BB.HNH-h.NHO-hh.hbonds -dlvec -0.5,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0
#./$name  -x BB.HNH-h.NHO-hh -b BB.HNH-h.NHO-hh.hbonds -dlvec -0.5,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0 -latscan 10,10 0.2,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0

# ====== Small Molecules On Substrate

#valgrind --log-file="valgrind.log" --leak-check=yes ./$name -m 10 -x common_resources/xyz/H2O -n     -g common_resources/xyz/NaCl_1x1_L2
#valgrind --leak-check=yes ./$name -m 10 -x common_resources/xyz/H2O -n     -g common_resources/xyz/NaCl_1x1_L2

#./$name -m 10    -x common_resources/xyz/H2O           -g common_resources/xyz/NaCl_1x1_L2
#./$name -m 10    -x common_resources/xyz/H2O           -g common_resources/xyz/NaCl_1x1_L2 -tricubic

#./$name -m $nsys -x common_resources/xyz/pyridine      -g common_resources/xyz/NaCl_1x1_L2
#./$name -m $nsys -x common_resources/xyz/pyridine  110 -g common_resources/xyz/NaCl_1x1_L2
#./$name -m $nsys -x common_resources/xyz/PTCDA         -g common_resources/xyz/NaCl_1x1_L2

#./$name -m 1 -x common_resources/xyz/PTCDA         -g common_resources/xyz/NaCl_1x1_L2 -iParalel 3  -verb 4 -perframe 1
#./$name -m 1 -x common_resources/xyz/PTCDA         -g common_resources/xyz/NaCl_1x1_L2 -iParalel 3  -verb 4
#./$name -m 1 -x common_resources/xyz/PTCDA         -g common_resources/xyz/NaCl_1x1_L2 -iParalel 3
#./$name -m 1 -x common_resources/xyz/PTCDA         -g common_resources/xyz/NaCl_1x1_L2 -iParalel 3 -tex 1


#./$name -m 1 -x common_resources/xyz/PTCDA         -g common_resources/xyz/NaCl_1x1_L2 -iParalel 3
#./$name -m 1 -x common_resources/xyz/PTCDA         -g common_resources/xyz/NaCl_1x1_L2 -iParalel 3 -tex 1
#./$name -m 1 -x common_resources/xyz/PTCDA         -g common_resources/xyz/NaCl_1x1_L2 -iParalel -1
#./$name -m 1 -x common_resources/xyz/PTCDA         -g common_resources/xyz/NaCl_1x1_L2 -iParalel 0 -verb 3

#./$name -m $nsys -x common_resources/xyz/PTCDA_SAM     -g common_resources/xyz/NaCl_1x1_L2
#./$name -m $nsys -x common_resources/xyz/BPBA          -g common_resources/xyz/NaCl_1x1_L2
#./$name -m $nsys -x common_resources/xyz/BPBA          -g common_resources/xyz/NaCl_1x1_L2 -e

#./$name -m 20 -x common_resources/xyz/PTCDA_SAM  -c 59  -g common_resources/xyz/NaCl_1x1_L2   -ManipulAnim
#./$name -m 20 -x common_resources/xyz/PTCDA_SAM  -c 59  -g common_resources/xyz/NaCl_1x1_L2




#./$name -m 100 -x common_resources/xyz/PTCDA  -c 29,1.0      -g common_resources/xyz/NaCl_1x1_L2     -iParalel 3   -perframe 100    -Ftol 1e-4    -lua0 scan_replicas.lua

# ====== Polymers On Substrate

#./$name -m $nsys -x common_resources/xyz/polydiacetylene           -g common_resources/xyz/NaCl_1x1_L2
#./$name -m $nsys -x common_resources/xyz/polydiacetylene           -g common_resources/xyz/NaCl_1x1_L2 -subs 4,common_resources/xyz/-COOH.xyz
#./$name -m $nsys -x common_resources/xyz/polydiacetylene    -n 221 -g common_resources/xyz/NaCl_1x1_L2
#./$name -m $nsys -x common_resources/xyz/polydiacetylene    -n 241 -g common_resources/xyz/NaCl_1x1_L2
#./$name -m $nsys -x common_resources/xyz/polydiacetylene    -n 241 -g common_resources/xyz/NaCl_1x1_L2 -ng 4,-2,2,3
#./$name -m $nsys -x common_resources/xyz/polydiacetylene    -n 141 -g common_resources/xyz/NaCl_1x1_L2 -ng 2,-1,2,3
#./$name -m $nsys -x common_resources/xyz/polydiacetylene    -n 141 -g common_resources/xyz/NaCl_1x1_L2 -ng 2,-1,3,4
#./$name -m $nsys -x common_resources/xyz/polydiacetylene    -n 141 -g common_resources/xyz/NaCl_1x1_L2 -subs 4,common_resources/xyz/-COOH.xyz  -ng 2,-1,3,4 -q 0.05

#./$name -m $nsys -x common_resources/xyz/polydiacetylene_OH -n 141 -g common_resources/xyz/NaCl_1x1_L2 -ng 2,-1,3,4
#./$name -m $nsys -x common_resources/xyz/polydiacetylene_OH        -g common_resources/xyz/NaCl_1x1_L2;
#./$name -m $nsys -x common_resources/xyz/polydiacetylene_OH        -g common_resources/xyz/NaCl_1x1_L2

#./$name -m $nsys -x common_resources/xyz/polymer-2_new             -g common_resources/xyz/NaCl_1x1_L2
#./$name -m 10    -x common_resources/xyz/polymer-2_new             -g common_resources/xyz/NaCl_1x1_L2   | tee out.log
#./$name -m 1     -x common_resources/xyz/polymer-2_new             -g common_resources/xyz/NaCl_1x1_L2
#./$name -m 2     -x common_resources/xyz/polymer-2_new             -g common_resources/xyz/NaCl_1x1_L2
#./$name -m 40    -x common_resources/xyz/polymer-2_new             -g common_resources/xyz/NaCl_1x1_L2  -nogridff
#./$name -m 40    -x common_resources/xyz/polymer-2_new             -g common_resources/xyz/NaCl_1x1_L2    -Ftol 1e-12   -perframe 500



#./$name -m 3    -x common_resources/xyz/polymer-2_new              -g common_resources/xyz/NaCl_1x1_L2    -Ftol 1e-4   -perframe 1  -iParalel 3   -group  9,20,17,1,9  0,1,7,9,15,17  -bbox -30.0,-30.0,-10.0,30.0,30.0,10.0,1.0,1.0,1.0
#./$name -m 3    -x common_resources/xyz/polymer-2_new              -g common_resources/xyz/NaCl_1x1_L2    -Ftol 1e-4   -perframe 1  -iParalel 3   -group  9,20,17,1,9  0,1,7,9,15,17  -c 10,1.,1.,0.0
#./$name -m 3    -x common_resources/xyz/polymer-2_new              -g common_resources/xyz/NaCl_1x1_L2    -Ftol 1e-4   -perframe 1  -iParalel 3   -group  9,20,17,1,9  0,1,7,9,15,17  -c 10,0.,0.,1.0
#./$name -m 100  -x common_resources/xyz/polymer-2_new              -g common_resources/xyz/NaCl_1x1_L2    -Ftol 1e-4   -perframe 100  -iParalel 3   -group  9,20,17,1,9  0,1,7,9,15,17   -c 10,1.,1.,0.0
#./$name -m 100  -x common_resources/xyz/polymer-2_new              -g common_resources/xyz/NaCl_1x1_L2    -Ftol 1e-4   -perframe 100  -iParalel 3   -group  9,20,17,1,9  0,1,7,9,15,17   -bbox -30.0,-30.0,-10.0,30.0,30.0,10.0,1.0,1.0,1.0
#./$name -m 1000 -x common_resources/xyz/polymer-2_new              -g common_resources/xyz/NaCl_1x1_L2    -Ftol 1e-4   -perframe 200  -iParalel 3   -group  9,20,17,1,9  0,1,7,9,15,17   -bbox -30.0,-30.0,-10.0,30.0,30.0,10.0,1.0,1.0,1.0
#./$name -m 10   -x common_resources/xyz/polymer-2_new              -g common_resources/xyz/NaCl_1x1_L2    -Ftol 1e-4   -perframe 200  -iParalel 3   -group  9,20,17,1,9  0,1,7,9,15,17  -bbox -30.0,-30.0,-10.0,30.0,30.0,10.0,1.0,1.0,1.0  -nogridff -dt 0.05 

#./$name -m 1000 -x common_resources/xyz/polymer-2_new-OH         -g common_resources/xyz/NaCl_1x1_L2   -verb 0  -dt 0.1  -Ftol 1e-4   -perframe 200  -iParalel 3   -group  9,20,17,1,9  0,1,7,9,15,17  -group  12,19,16,3,12  12,3,2,6,14,16 -bbox -4.02,-4.02,-10.0,20.10,20.20,20.0,1e-3,1e-3,1.0

#./$name -m 10   -x common_resources/xyz/polymer-2_new          -dt 0.05   -g common_resources/xyz/NaCl_1x1_L2    -Ftol 1e-4   -perframe 200  -iParalel 3 -group  9,20,17,1,9  0,1,7,9,15,17   -c 10
#./$name -m 1000 -x common_resources/xyz/polymer-2_new                     -g common_resources/xyz/NaCl_1x1_L2    -Ftol 1e-4   -perframe 200  -iParalel 3   -group  9,20,17,1,9  0,1,7,9,15,17  -c 10
#./$name -m 500  -x common_resources/xyz/polymer-2_new                     -g common_resources/xyz/NaCl_1x1_L2    -Ftol 1e-4   -perframe 200  -iParalel 3   -group  9,20,17,1,9  0,1,7,9,15,17  -c 10

#./$name -m 10 -x common_resources/polymer-2_new                           -g common_resources/xyz/NaCl_1x1_L2    -lua0 scan_replicas.lua     -Ftol 1e-4   -perframe 200  -iParalel 3   -group  9,20,17,1,9  0,1,7,9,15,17  0.02,0.3


#./$name -m 5    -x common_resources/xyz/polymer-2_new                  -g common_resources/xyz/NaCl_1x1_L2    -Ftol 1e-12   -perframe 100  -iParalel 3    -group 0 0,1,9,17,15,7,21,27,28,25,26,36 -group 1 2,3,6,12,14,16,18,24,30,29,35
#./$name -m 100 -x common_resources/xyz/polymer-2_new                  -g common_resources/xyz/NaCl_1x1_L2    -Ftol 1e-12   -perframe 200  -iParalel 3 -T 300 0.2 -gopt 1000,100000 0.25,1.0 -verb 0 -perframe 100
#./$name -m 200  -x common_resources/xyz/polymer-2_new                  -g common_resources/xyz/NaCl_1x1_L2    -Ftol 1e-12   -perframe 100  -iParalel 3
#./$name -m 200  -x common_resources/xyz/polymer-2_new                  -g common_resources/xyz/NaCl_1x1_L2    -Ftol 1e-12   -perframe 100  -iParalel 3
#./$name -m 200  -x common_resources/xyz/polymer-2_new                  -g common_resources/xyz/NaCl_1x1_L2    -Ftol 1e-12   -perframe 100  -iParalel 2


#./$name -m 30  -x common_resources/xyz/polymer-2_new   -g common_resources/xyz/NaCl_1x1_L2    -iParalel 1   -perframe 100  -gopt 1000,1000 0.25,1.0    -T 10000 0.02  -Ftol 1e-2   -verb 0



#valgrind --log-file="valgrind.log" --leak-check=yes ./$name -x common_resources/H2O
#valgrind --leak-check=yes --log-file="valgrind.log" ./$name -x common_resources/HCOOH -m 10
#valgrind --leak-check=yes --track-origins=yes  --log-file="valgrind.log" ./$name -x common_resources/HCOOH -m 10


#grep GPU_GFF_z out.log | cut -c 11-  > GPU_makeGridFF.log

# GPU global optimization
#./$name -m 400  -x common_resources/xyz/deoxyglucose          -g common_resources/xyz/NaCl_8x8_L3_step -iParalel 3 -T 300 0.2   -gopt 1000,100000 0.25,1.0 -verb 0 -perframe 100
#./$name -m 400  -x common_resources/xyz/1,20-eicosanediol     -g common_resources/xyz/NaCl_8x8_L3_step -iParalel 3 -T 1000 0.2  -gopt 1000,100000 0.25,1.0 -verb 0 -perframe 100
#./$name -m 600  -x common_resources/xyz/deoxyglucose          -g common_resources/xyz/NaCl_8x8_L3      -iParalel 3 -T 300 0.2   -gopt 1000,100000 0.25,1.0 -verb 0 -perframe 100
#./$name -m 400  -x common_resources/xyz/1,20-eicosanediol     -g common_resources/xyz/NaCl_8x8_L3      -iParalel 3 -T 300 0.2   -gopt 1000,100000 0.25,1.0 -verb 0 -perframe 100
#./$name -m 1000    -x common_resources/xyz/xylitol_for_gridFF               -g common_resources/xyz/NaCl_8x8_L3      -iParalel 3 -T 300 0.2   -gopt 1000,100000 0.25,1.0 -verb 0 -perframe 100
#./$name -m 1000    -x common_resources/xyz/xylitol_for_gridFF                                                        -iParalel 3 -T 300 0.2   -gopt 1000,100000 0.25,1.0 -verb 0 -perframe 100
#./$name -m 1000    -x common_resources/xyz/xylitol_WO_gridFF              -g common_resources/xyz/NaCl_8x8_L3      -iParalel 3 -T 300 0.2   -gopt 1000,100000 0.25,1.0 -verb 0 -perframe 100 -nogridff
#./$name -m 1000    -x common_resources/xyz/xylitol               -g common_resources/xyz/NaCl_8x8_L3      -iParalel 3 -T 300 0.2   -gopt 1000,100000 0.25,1.0 -verb 0 -perframe 100
# ./$name -m 100 -x common_resources/xyz/xylitol               -g common_resources/xyz/NaCl_8x8_L3      -iParalel 3 -T 300 0.2   -gopt 1000,100000 0.25,1.0 -verb 0 -perframe 100
#./$name -m 20 -x common_resources/xyz/xylitol               -g common_resources/xyz/NaCl_8x8_L3_step      -iParalel 3 -T 300 0.2   -gopt 1000,100000 0.25,1.0 -verb 0 -perframe 100
#./$name -m 200 -x common_resources/xyz/xylitol               -g common_resources/xyz/NaCl_8x8_L3_NaClHole      -iParalel 3 -T 300 0.2   -gopt 1000,100000 0.25,1.0 -verb 0 -perframe 100
#./$name -m 25 -x common_resources/xyz/xylitol               -g common_resources/xyz/NaCl_8x8_L3_ClHole      -iParalel 3 -T 150 0.2   -gopt 1000,100000 0.25,1.0 -verb 0 -perframe 100
#./$name -m 150  -x common_resources/xyz/xylitol_THREE         -g common_resources/xyz/NaCl_8x8_L3      -iParalel 3 -T 300 0.2   -gopt 1000,100000 0.25,1.0 -verb 0 -perframe 100
#./$name -m 1500 -x common_resources/xyz/deoxyribose           -g common_resources/xyz/NaCl_8x8_L3      -iParalel 3 -T 300 0.2   -gopt 1000,100000 0.25,1.0 -verb 0 -perframe 100
#./$name -m 1000 -x common_resources/xyz/PTCDA_SAM             -g common_resources/xyz/NaCl_8x8_L3      -iParalel 3 -T 1000 0.2  -gopt 500,10000 0.25,1.0   -verb 0 -perframe 100
#./$name -m 1500 -x common_resources/xyz/pyridine              -g common_resources/xyz/NaCl_8x8_L3      -iParalel 3 -T 1000 0.2  -gopt 1000,100000 0.25,1.0 -verb 0 -perframe 100
#./$name -m 300   -x common_resources/xyz/nHexadecan_dicarboxylic -g common_resources/xyz/NaCl_8x8_L3    -iParalel 3 -T 1000 0.02 -gopt 1000,1000 0.25,1.0   -verb 0 -perframe 100 


./$name -m 100    -x common_resources/xyz/xylitol_WO_gridFF     -g common_resources/xyz/surfaces_for_throughput/NaCl_3x3_Cl_hole         -iParalel 3 -T 300 0.2   -gopt 1000,100000 0.25,1.0 -verb 0 -perframe 100 
#./$name -m 2000    -x common_resources/xyz/xylitol_WO_gridFF                 -iParalel 3 -T 300 0.2   -gopt 1000,100000 0.25,1.0 -verb 0 -perframe 100 -grid_nPBC 2,2,0 # -nogridff
