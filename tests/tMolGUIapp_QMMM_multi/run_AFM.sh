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

export OMP_NUM_THREADS=1

rm DEBUG*

./$name -m 1 -x common_resources/xyz/PTCDA          -g common_resources/xyz/NaCl_1x1_L2 -iParalel 1  -iMO 3 -substr_iso 0.05
#./$name -m 1 -x common_resources/xyz/PTCDA_SAM     -g common_resources/xyz/NaCl_1x1_L2 -iParalel 1
#./$name -m 1 -x common_resources/xyz/BPBA          -g common_resources/xyz/NaCl_1x1_L2 -iParalel 1
#./$name -m 1 -x common_resources/xyz/BPBA          -g common_resources/xyz/NaCl_1x1_L2 -iParalel 1 -e

