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


nbSteps=100         # JE -> nProdSteps
nMDsteps=100000     # JE -> nrealization
nEQsteps=1000
t_damp=100
T=300
dt=0.01
K=10.0
Ns=(20)
lamda1=1.0
lamda2=3.0


for N in "${Ns[@]}"
do
    python3 run.py --nMDsteps $nMDsteps \
                   --nEQsteps $nEQsteps \
                   --t_damp $t_damp \
                   --T $T \
                   --dt $dt \
                   --nbSteps $nbSteps \
                   --N $N \
                   --K $K \
                   --lamda1 $lamda1 \
                   --lamda2 $lamda2 
    ./results/TI_plot.sh
    cp results/TI_plot.png "results/JE_entSpring_N${N}_MDsteps${nMDsteps}_EQsteps${nEQsteps}_tdamp${t_damp}_T${T}_dt${dt}_nSteps${nSteps}_k${k}_K${K}.png"
done
