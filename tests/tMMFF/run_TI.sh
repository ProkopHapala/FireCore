#!/bin/bash

ln -s ../../cpp/common_resources data
ln -s ../../cpp/common_resources common_resources 

wd=`pwd`
cd ../../cpp/Build/libs/Molecular/
pwd
rm libMMFF_lib.so
make MMFF_lib
rm   libLattice2D_lib.so
make Lattice2D_lib
cd $wd

cd ../../cpp/Build/libs_SDL
rm libMolGUIlib.so
make MolGUIlib
cd $wd

# ---- Multiprocesing
ncpu=`nproc`
ncpu=$(($ncpu - 1))     # let one CPU free for user interaction
echo "compile using ncpu="$ncpu
OMP_NUM_THREADS=$ncpu
export OMP_NUM_THREADS

#rm *.bin



nMDsteps=1000000
nEQsteps=10000
t_damp=100
T=300
dt=0.01
nSteps=100
k=1.0
K=10.0
Ns=(5 10 20 30)


for N in "${Ns[@]}"
do
    python3 run.py --nMDsteps $nMDsteps \
                   --nEQsteps $nEQsteps \
                   --t_damp $t_damp \
                   --T $T \
                   --dt $dt \
                   --nSteps $nSteps \
                   --N $N \
                   --k $k \
                   --K $K
    ./results/TI_plot.sh
    cp results/TI_plot.png "results/entSpring_N${N}_MDsteps${nMDsteps}_EQsteps${nEQsteps}_tdamp${t_damp}_T${T}_dt${dt}_nSteps${nSteps}_k${k}_K${K}.png"
done
