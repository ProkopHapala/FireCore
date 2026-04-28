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


./$name -m 100    -x common_resources/xyz/xylitol_centered     -g common_resources/xyz/surfaces_for_throughput/NaCl_3x3_Cl_hole         -iParalel 3 -T 300 0.2   -gopt 1000,100000 0.25,1.0 -verb 0 -perframe 100 
./$name -m 100    -x common_resources/xyz/xylitol_centered     -g common_resources/xyz/surfaces_for_throughput/NaCl_3x3_Cl_hole         -iParalel 3 -T 300 0.2   -gopt 1000,100000 0.25,1.0 -verb 0 -perframe 100 -nogridff
./$name -m 100    -x common_resources/xyz/xylitol_centered  -uff   -g common_resources/xyz/surfaces_for_throughput/NaCl_3x3_Cl_hole         -iParalel 3 -T 300 0.2   -gopt 1000,100000 0.25,1.0 -verb 0 -perframe 100 
./$name -m 100    -x common_resources/xyz/xylitol_centered  -uff   -g common_resources/xyz/surfaces_for_throughput/NaCl_3x3_Cl_hole         -iParalel 3 -T 300 0.2   -gopt 1000,100000 0.25,1.0 -verb 0 -perframe 100 -nogridff
