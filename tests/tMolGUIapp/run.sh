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

no_compile=false
if [[ "$*" == *"no"* ]]; then
    no_compile=true
fi

if [ "$no_compile" = false ]; then
# ---- Compilation
wd=`pwd`
cd $dir
pwd
rm $name
make -j$ncpu $name   # 2>$wd/compile_err.log
cd $wd
rm $name
ln -s $dir/$name .
fi


# ------- asan (Memory Sanitizer)
#LD_PRELOAD=$(g++ -print-file-name=libasan.so)
#echo   $LD_PRELOAD
#export LD_PRELOAD


# ---- Run

#rm *.bin *.xsf

# ====== Small Molecules

#./$name -x common_resources/C2H4          -iParalel 0

#./$name -x common_resources/H2O      -uff  -iParalel 0  -verb 2   -perframe 1 -dt 0.0001
#./$name -x common_resources/H2O      -uff  -iParalel 0  -verb 2   -perframe 1 -dt 0.01
#./$name -x common_resources/butandiol -uff  -iParalel 0  -verb 2   -perframe 1 -dt 0.01
#./$name -x common_resources/butandiol   -iParalel 0  -verb 2   -perframe 1 -dt 0.01


#./$name -x common_resources/H2O            -iParalel 0  -verb 0   -perframe 1 -dt 0.001
#./$name -x common_resources/H2O            -iParalel 0  -verb 0   -perframe 1 -dt 0.05
#./$name -x common_resources/H2O        -e  -iParalel 0  -verb 0

#./$name -x common_resources/HCOOH         -iParalel 0
#./$name -x common_resources/HCOOH     -e  -iParalel 0

#./$name -x common_resources/HCOOH     -e  -iParalel 0


#./$name -x common_resources/HCOOH_xy  -e  -iParalel 0
#./$name -x common_resources/CH2O      -e  -iParalel 0

#./$name -x common_resources/formic_dimer
#./$name -x common_resources/formic_dimer -e -iParalel 0

#./$name -x common_resources/pyridine
#./$name -x common_resources/pyridine -lua script.lua
#./$name -x common_resources/propandiol
#./$name -x common_resources/butandiol

#./$name -x common_resources/nHexadecan_dicarboxylic -uff -iParalel 0 -T 100 0.01 -verb 2 -perframe 1
#./$name -x common_resources/nHexadecan_dicarboxylic -uff -iParalel 0 -T 100 0.01 -verb 2 -perframe 500
#./$name -x common_resources/nHexadecan_dicarboxylic -uff -iParalel 0 -T 100 0.01 -verb 2 -perframe 2000
#./$name -x common_resources/nHexadecan_dicarboxylic -uff -iParalel 0 -T 100 0.01 -verb 2 -perframe 2000 -NBneigh


#./$name -x common_resources/nHexadecan_dicarboxylic -iParalel 0 -T 100 0.01 -verb 2 -perframe 2000
#./$name -x common_resources/nHexadecan_dicarboxylic -iParalel 0 -T 100 0.01 -verb 2 -perframe 500 -NBneigh
#./$name -x common_resources/nHexadecan_dicarboxylic -iParalel 0 -T 100 0.01 -NBneigh
#./$name -x common_resources/nHexadecan_dicarboxylic -iParalel 0
#./$name -x common_resources/nHexadecan_dicarboxylic -iParalel 0 -NBneigh


#./$name -x common_resources/nHexadecan_dicarboxylic -perframe 1 -T 10 0.1
#./$name -x common_resources/nHexadecan_dicarboxylic -perframe 1 -T 100 0.1
#./$name -x common_resources/nHexadecan_dicarboxylic -perframe 1 -T 1000 0.1
#./$name -x common_resources/nHexadecan_dicarboxylic -b nHexadecan_dicarboxylic.cons -perframe 1 -T 10 0.1

#./$name -x common_resources/nHexadecan_dicarboxylic  -T 100 0.01
#./$name -x common_resources/nHexadecan_dicarboxylic -b nHexadecan_dicarboxylic.cons -T 100 0.01
#./$name -x common_resources/nHexadecan_dicarboxylic -b nHexadecan_dicarboxylic.cons -perframe 1000 -gopt 1000,1000 0.25,1.0
#./$name -x common_resources/nHexadecan_dicarboxylic -b nHexadecan_dicarboxylic.cons -perframe 10000 -gopt 1000,1000 0.25,1.0
#./$name -x common_resources/nHexadecan_dicarboxylic -drive nHexadecan_dicarboxylic.cons -perframe 10000 -gopt 1000,1000 0.25,1.0   T 1000 0.01
#./$name -x data/hexan-dicarboxylic -b nHexadecan_dicarboxylic.cons -perframe 10000 -gopt 1000,1000 0.25,1.0 -verb 0

#./$name -x common_resources/nHexadecan_dicarboxylic -b nHexadecan_dicarboxylic.cons
#./$name -x common_resources/nHexadecan_dicarboxylic -T 500.0
#./$name -x common_resources/nHexadecan_dicarboxylic -T 500.0


#./$name -uff -x common_resources/C2H4


#./$name -x common_resources/Si10_H
#.$name -x common_resources/Si10_H -iParalel 0

#./$name -x Si255_H_relaxed
#./$name -x Si405_H_relaxed
#./$name -x si_111_surface_4x6-
#./$name -x Si705_relaxed
#./$name -x Si2505_111
#./$name -x Si2647_100
#./$name -x Si4930_110
#./$name -x Si2505_111-
#./$name -x Si2647_100-
#./$name -x Si4930_110-
#./$name -x Si2505_111-H       -perframe 1
#./$name -x Si2647_100-H      -perframe 1
#./$name -x Si4930_110-H      -perframe 1

#./$name -x Si2505_111-H-relaxed        -perframe 1
#./$name -x Si2505_111-H        -perframe 1  # -noNB
#./$name -x Si2505_111-noH-SiH3  -perframe 1   
#./$name -x Si2505_111-noH-SiH3-relaxed  -perframe 1   

#./$name -x Si2505_111-H-brak-110-relaxed  -perframe 1
#./$name -x Si2505_111-H-SiH3-relaxed      -perframe 1
#./$name -x Si2505_111-H-relaxed           -perframe 1




# ====== Small Molecules On Substrate

#./$name -x common_resources/H2O       -g common_resources/NaCl_1x1_L2
#./$name -x common_resources/H2O       -g common_resources/NaCl_1x1_L2 -tricubic

#./$name -x common_resources/H2O       -g common_resources/NaCl_1x1_L2  -lua makeGUI.lua
#./$name -x  common_resources/H2O2     -g common_resources/NaCl_1x1_L2 -e
#./$name -x common_resources/Molekuly  -g common_resources/NaCl_1x1_L2 -e


#./$name -x common_resources/NH3      -g common_resources/NaCl_1x1_L2
#./$name -x common_resources/HCOOH    -g common_resources/NaCl_1x1_L2
#./$name -x common_resources/pyridine -g common_resources/NaCl_1x1_L2
#./$name -x common_resources/pyridine -g common_resources/NaCl_1x1_L2 -lua test_add_mols.lua
#./$name -x common_resources/pyridine -n 110 -g common_resources/NaCl_1x1_L2

#./$name -x common_resources/N_2edge_assay_1  -g common_resources/NaCl_1x1_L2
#./$name -x common_resources/N_2edge_assay_1  -g common_resources/NaCl_1x1_L2 -e


#./$name -x common_resources/PTCDA -g common_resources/NaCl_1x1_L2
#./$name -x common_resources/PTCDA -g common_resources/NaCl_1x1_L2 -lua test_add_mols.lua

# ====== Polymers

#./$name -x common_resources/polydiacetylene
#./$name -x common_resources/polydiacetylene     -subs 4,common_resources/-COOH.xyz
#./$name -x common_resources/polydiacetylene_OH
#./$name -x common_resources/polymer-2_new
#./$name -x common_resources/polymer-2_new -verb 0  -lua script.lua

#./$name -x common_resources/polymer-2_new -EachAngle
#./$name -x common_resources/polymer-2_new -EachAngle -torsions


#./$name -x common_resources/polymer-2_new   -iParalel 0 -T 100 0.01 -verb 2 -perframe 500
#./$name -x common_resources/polymer-2_new   -iParalel 0 -T 100 0.01 -verb 2 -perframe 500 -NBneigh
#./$name -x common_resources/polymer-2_new   -iParalel 0 -T 100 0.01 -verb 2 -perframe 500 -noNB



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

#./$name -x common_resources/polymer-2_new            -g common_resources/NaCl_1x1_L2
#./$name -x common_resources/polymer-2_new            -g common_resources/NaCl_1x1_L2 -tricubic


#./$name -x common_resources/polymer-2_new            -g common_resources/NaCl_1x1_L2  -e
#./$name -x common_resources/polymer-2          -n 221 -g common_resources/NaCl_1x1_L2
#./$name -x common_resources/polymer-2_new   -c 10     -g common_resources/NaCl_1x1_L2

#./$name -x common_resources/polymer-2_new   -g common_resources/NaCl_1x1_L2  -perframe 100  -Ftol 1e-2 -seed 654654   -verb 0
#./$name -x common_resources/polymer-2_new   -g common_resources/NaCl_1x1_L2  -perframe 100  -T 10000 0.02  -Ftol 1e-2 -seed 654654   -verb 0
#./$name -x common_resources/polymer-2_new   -g common_resources/NaCl_1x1_L2  -perframe 100  -gopt 1000,1000 0.0,0.0  -T 1000 0.05  -Ftol 1e-4 -seed 654654   -verb 0
#./$name -x common_resources/polymer-2_new   -g common_resources/NaCl_1x1_L2  -perframe 100  -gopt 1000,1000 0.0,0.0    -T 1000 0.1  -Ftol 1e-2 -seed 654654   -verb 0
#./$name -x common_resources/polymer-2_new   -g common_resources/NaCl_1x1_L2  -drive polymer-2_new.cons -perframe 100  -gopt 1000,1000 0.0,0.0    -T 100 0.001  -Ftol 1e-2 -seed 654654   -verb 0
#./$name -x common_resources/polymer-2_new   -g common_resources/NaCl_1x1_L2  -drive polymer-2_new.cons -perframe 100  -gopt 1000,1000 0.0,0.0    -T 1000 0.01  -Ftol 1e-6 -seed 654654   -verb 0 
#./$name -x common_resources/polymer-2_new   -g common_resources/NaCl_1x1_L2  -drive polymer-2_new.cons -perframe 100  -gopt 1000,1000 0.0,0.0    -T 1000 0.1  -Ftol 1e-6 -seed 654654   -verb 0  -zspring -10.0,5.0,0.2
#./$name -x common_resources/polymer-2_new   -g common_resources/NaCl_1x1_L2   -drive polymer-2_new.cons -perframe 100  -gopt 1000,1000 0.0,0.0    -T 1000 0.1  -Ftol 1e-6 -seed 654654   -verb 0  -stuck 100,0.2


#./$name -x common_resources/polymer-2_new   -g common_resources/NaCl_1x1_L2   -Ftol 1e-6 -iParalel 0  -dt 0.05 -nogridff -perframe 1
#./$name -x common_resources/polymer-2_new   -g common_resources/NaCl_1x1_L2   -Ftol 1e-12 -iParalel 0  -dt 0.05 -nogridff -perframe 100
#./$name -x common_resources/polymer-2_new   -g common_resources/NaCl_1x1_L2   -Ftol 1e-6 -iParalel 0  -dt 0.05 -perframe 1 -gridffmode 2
./$name -x common_resources/polymer-2_new   -g common_resources/NaCl_1x1_L2   -Ftol 1e-6 -iParalel 0  -dt 0.05 -perframe 100 -gridffmode 6


#./$name  -x common_resources/polymer-2_new-OH                                                -Ftol 1e-4 -iParalel 0 -perframe 1
#./$name -verb 2 -e -x common_resources/polymer-2_new-OH   -g common_resources/NaCl_1x1_L2   -Ftol 1e-4 -iParalel 0 -perframe 200
#./$name  -x common_resources/polymer-2_new-COOH           -g common_resources/NaCl_1x1_L2   -Ftol 1e-4 -iParalel 0 -perframe 10


#./$name  -x common_resources/polymer-2_new   -g common_resources/NaCl_1x1_L2   -Ftol 1e-12 -iParalel 0 -perframe 10  -group  9,20,17,1,9  0,1,9,17,15,7,21,27,28,25,26,36 -group 12,19,16,3,12 2,3,6,12,14,16,18,24,30,29,35
#./$name  -x common_resources/polymer-2_new   -g common_resources/NaCl_1x1_L2   -Ftol 1e-12 -iParalel 0 -perframe 10  -group 12,19,16,3,12 2,3,6,12,14,16,18,24,30,29,35
#./$name  -x common_resources/polymer-2_new   -g common_resources/NaCl_1x1_L2   -Ftol 1e-12 -iParalel 0 -perframe 10  -group  9,20,17,1,9  0,1,9,17,15,7,21,27,28,25,26,36

#./$name -x common_resources/polymer-2_new   -g common_resources/NaCl_1x1_L2   -drive polymer-2_new.cons -perframe 100  -gopt 1000,1000 0.0,0.0    -T 1000 0.1  -Ftol 1e-6 -seed 654654  -stuck 300,0.2 -iParalel 0   -dt 0.05   -nogridff

#./$name -x common_resources/polymer-2_new   -g common_resources/NaCl_1x1_L2   -drive polymer-2_new.cons -perframe 100  -gopt 1000,1000 0.0,0.0    -T 1000 0.1  -Ftol 1e-6 -seed 654654  -stuck 300,0.2 -iParalel 1   -dt 0.05
#./$name -x common_resources/polymer-2_new   -g common_resources/NaCl_1x1_L2   -drive polymer-2_new.cons -perframe 100  -gopt 1000,1000 0.0,0.0    -T 1000 0.1  -Ftol 1e-6 -seed 654654  -stuck 500,400,0.2 -iParalel 1   -dt 0.01

#./$name -x common_resources/polymer-2_new -g common_resources/NaCl_1x1_L2  -iParalel 0 -T 100 0.01 -verb 2 -perframe 500
#./$name -x common_resources/polymer-2_new -g common_resources/NaCl_1x1_L2  -iParalel 0 -T 100 0.01 -verb 2 -perframe 500 -tricubic
#./$name -x common_resources/polymer-2_new -g common_resources/NaCl_1x1_L2  -iParalel 0 -verb 2 -perframe 500 -tricubic  -Ftol 0.0



#./$name -x common_resources/polymer-2_new   -g common_resources/NaCl_1x1_L2  -perframe 100  -gopt 1000,1000 0.0,0.0    -T 2000 0.1  -Ftol 1e-6 -seed 654654   -verb 0  -stuck 100,0.2 


#./$name -x BB.HNH-h.NHO-hh -b BB.HNH-h.NHO-hh.hbonds -perframe 1
#./$name -x BB.HNH-h.NHO-hh -b BB.HNH-h.NHO-hh.hbonds -perframe 50
#./$name -x BB.HNH-h.NHO-hh -b BB.HNH-h.NHO-hh.hbonds -dlvec -2.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0
#./$name -x BB.HNH-h.NHO-hh -b BB.HNH-h.NHO-hh.hbonds -dlvec -0.5,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0
#./$name  -x BB.HNH-h.NHO-hh -b BB.HNH-h.NHO-hh.hbonds -dlvec -0.5,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0 -latscan 10,20 0.0,0.5,0.0,0.0,0.0,0.0,0.0,0.0,0.0
#./$name  -x BB.HNH-h.NHO-hh -b BB.HNH-h.NHO-hh.hbonds -prelat 5,10000 -0.5,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0 -latscan 10,20 0.0,0.5,0.0,0.0,0.0,0.0,0.0,0.0,0.0


# ===== test Collision Damping

#./$name -x common_resources/polymer-2_new  -col_damp 10 -1.5 -1.6 0.9   0.0 0.5
#./$name -x common_resources/polymer-2_new  -col_damp 10 -1.0 -1.0 1.0   0.0 0.5

#./$name -x common_resources/O2 -perframe  1  -col_damp 5 -1.0 -1.0 0.1  0.0 0.5
#./$name -x common_resources/O2 -perframe  1  -col_damp 5 1.0 -1.0 0.0   0.0 0.5
#./$name -x common_resources/O2 -perframe  1  -col_damp 5 -1.0 1.0 0.0   0.0 0.5



#./$name -x common_resources/diamin_and_diether_C10 -perframe  1  -col_damp 1 0.5 0.3 0.0   0.0 0.5
#./$name -x common_resources/diamin_and_diether_C10- -perframe  10  -col_damp 1 -0.5 -0.3 0.0   0.0 0.5
#./$name -x common_resources/C8_diamin_diether-2 -perframe  10  -col_damp 1 -0.5 -0.3 0.0   0.0 0.5
#./$name -x common_resources/C8_diamin_diether -perframe  10  -col_damp 1 -0.5 -0.3 0.0   0.0 0.5
#./$name -x common_resources/nNonan_cross -perframe  10  -col_damp 1 -0.5 -0.3 0.0   0.0 0.5
#./$name -x common_resources/hydropentacene_cross -perframe  10  -col_damp 1 -0.5 -0.3 0.0   0.0 0.5

#./$name -x common_resources/nHexadecan      -perframe 100  -col_damp 10 -1.0 -1.0 0.1    0.0 0.5  
#./$name -x common_resources/nHexadecan_fold -perframe 100  -col_damp 10 -1.0 -1.0 0.1    0.0 0.5
#./$name -x common_resources/nHexadecan_fold -perframe 100  -col_damp 5  -1.0 -1.0 0.05   0.0 0.5
#./$name -x common_resources/nHexadecan_fold -perframe 100  -col_damp 5   1.0 -1.0 0.02   0.0 0.5
#./$name -x common_resources/nHexadecan_fold -perframe 100  -col_damp 5  -1.0 -1.0 0.02   0.0 0.5

#./$name -x common_resources/nHexadecan_fold -perframe  10  -col_damp 5  -1.0 -1.0 0.01   0.0 0.5
#./$name -x common_resources/nHexadecan_fold -perframe  10  -col_damp 5  -0.1 -0.01 0.005 0.0 0.5
#./$name -x common_resources/nHexadecan_fold -perframe  10  -col_damp 5  -0.5 0.01 0.005  0.0 0.5

#./$name -x common_resources/nHexadecan_fold -perframe  10  -col_damp 5  -1.0 -1.0 0.01   0.0 0.5
#./$name -x common_resources/nHexadecan_fold -perframe  10  -col_damp 5   1.0 -1.0 0.01   0.0 0.5
#./$name -x common_resources/nHexadecan_fold -perframe  10  -col_damp 5  -1.0  1.0 0.01   0.0 0.5
#./$name -x common_resources/nHexadecan_fold -perframe  10  -col_damp 2   1.0  0.1  0.001 0.0 0.5
#./$name -x common_resources/nHexadecan_fold -perframe  10  -col_damp 2   1.0  0.01 0.001 0.0 0.5


#./$name -x common_resources/polymer-2_new -perframe 100  -col_damp 2 -1.0 -1.0 0.020   0.0 0.5
#./$name -x common_resources/polymer-2_new -perframe 100  -col_damp 2 -1.0 -1.0 0.005   0.0 0.5
#./$name -x common_resources/polymer-2_new -perframe 100  -col_damp 2  1.0 -1.0 0.005   0.0 0.5
#./$name -x common_resources/polymer-2_new -perframe 100  -col_damp 2 -1.0  1.0 0.005   0.0 0.5

#./$name -x common_resources/polymer-2_new -perframe 100  -col_damp 4 -1.0 -1.0 0.020   0.0 0.5
#./$name -x common_resources/polymer-2_new -perframe 100  -col_damp 4 -1.0 -1.0 0.005   0.0 0.5
#./$name -x common_resources/polymer-2_new -perframe 100  -col_damp 4  1.0 -1.0 0.005   0.0 0.5
#./$name -x common_resources/polymer-2_new -perframe 100  -col_damp 4 -1.0  1.0 0.005   0.0 0.5

#./$name -x common_resources/polymer-2_new -perframe 100  -col_damp 10  1.0 -1.0 0.1   0.0 0.5
#./$name -x common_resources/polymer-2_new -perframe 100  -col_damp 10  1.0 -1.0 0.1   0.0 0.5
#./$name -x common_resources/polymer-2_new -perframe 100  -col_damp 10 -1.0  1.0 0.1   0.0 0.5
#./$name -x common_resources/polymer-2_new -perframe 100  -col_damp 10 -1.0 -1.0 0.1   0.0 0.5





#valgrind --log-file="valgrind.log" --leak-check=yes ./$name -x common_resources/H2O


