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

# ====== Small Molecules @ Surface

#./$name -x common_resources/xyz/H2O       -g common_resources/xyz/NaCl_1x1_L3 -iParalel 1
#./$name -x common_resources/xyz/NH3       -g common_resources/xyz/NaCl_1x1_L3 -iParalel 1
#./$name -x common_resources/xyz/pyridine  -g common_resources/xyz/NaCl_1x1_L3 -iParalel 1
#./$name -x common_resources/xyz/HCOOH     -g common_resources/xyz/NaCl_1x1_L3 -iParalel 1

#./$name -x common_resources/xyz/nHexadecan_dicarboxylic -g common_resources/xyz/NaCl_1x1_L3 -iParalel 1  


#./$name -x common_resources/xyz/PTCDA           -g common_resources/xyz/NaCl_1x1_L3       -iParalel 1
#./$name -x common_resources/xyz/PTCDA           -g common_resources/xyz/NaCl_8x8_L3_step  -iParalel 1

#./$name -x common_resources/xyz/guanine          -g common_resources/xyz/NaCl_1x1_L3         -iParalel 0
#./$name -x common_resources/xyz/guanine          -g common_resources/xyz/NaCl_8x8_L3_step    -iParalel 1

./$name -x common_resources/xyz/uracil           -g common_resources/xyz/NaCl_1x1_L3         -iParalel 1
#./$name -x common_resources/xyz/uracil           -g common_resources/xyz/NaCl_8x8_L3_step    -iParalel 1
#./$name -x common_resources/xyz/uracil           -g common_resources/xyz/NaCl_8x8_L3_step    -iParalel 0;

#./$name -x common_resources/xyz/guanine-cytosine -g common_resources/xyz/NaCl_1x1_L3         -iParalel 1
#./$name -x common_resources/xyz/guanine-cytosine -g common_resources/xyz/NaCl_8x8_L3_step    -iParalel 1

# ====== Polymers @ Surfacef

#./$name  -x common_resources/xyz/polymer-2_new      -g common_resources/xyz/NaCl_1x1_L3   -iParalel 1 
#./$name  -x common_resources/xyz/polymer-2_new-OH   -g common_resources/xyz/NaCl_1x1_L3   -iParalel 1 
#./$name  -x common_resources/xyz/polymer-2_new      -g common_resources/xyz/NaCl_8x8_L3_step   -iParalel 1 
#./$name  -x common_resources/xyz/polymer-2_new-OH   -g common_resources/xyz/NaCl_8x8_L3_step   -iParalel 1 

#./$name -x tAttach/PNA_CG_poly-.mol2                -g common_resources/xyz/NaCl_1x1_L3                        -iParalel 1 
#./$name -x tAttach/PNA_CG_poly-.mol2   -nPBC 1,3,0  -g common_resources/xyz/NaCl_1x1_L3 -shift 0.0,0.0,2.0     -iParalel 1
#./$name -x tAttach/DANA_CG_poly-.mol2  -nPBC 1,3,0  -g common_resources/xyz/NaCl_1x1_L3 -shift 0.0,0.0,2.0    -iParalel 1
#./$name -x tAttach/TNA_CG_poly-.mol2   -nPBC 1,3,0  -g common_resources/xyz/NaCl_1x1_L3 -shift 0.0,0.0,2.0     -iParalel 1



# ==== Small Molecules  (without surface)

#./$name -x common_resources/mol/H2O.mol2                -iParalel 0 
#./$name -x common_resources/xyz/C2H4                    -iParalel 0 
#./$name -x common_resources/xyz/HCOOH                   -iParalel 0 
#./$name -x common_resources/xyz/butandiol               -iParalel 0 
#./$name -x common_resources/mol/porphirin.mol2          -iParalel 0 
#./$name -x common_resources/xyz/pyridine                -iParalel 0 
#./$name -x common_resources/xyz/formic_dimer            -iParalel 0 
#./$name -x common_resources/xyz/nHexadecan_dicarboxylic -iParalel 0 




