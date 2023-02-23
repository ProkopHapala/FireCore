name=MolGUIapp
dir=../../cpp/Build/apps/MolecularEditor
ln -s ../../cpp/common_resources data
ln -s ../../cpp/common_resources common_resources 

wd=`pwd`
cd $dir
pwd
rm $name
make $name
cd $wd

ln -s $dir/$name .

#rm *.bin *.xsf

#./$name -x common_resources/pyridine -n -g common_resources/NaCl_sym-center
#./$name -x common_resources/pyridine -n 110 -g common_resources/NaCl_sym-center
#./$name -x common_resources/polymer-2 -n 221 -g common_resources/NaCl_sym-center
#./$name -x common_resources/polydiacetylene -n 221 -g common_resources/NaCl_sym-center
#./$name -x common_resources/polydiacetylene -n 241 -g common_resources/NaCl_564
#./$name -x common_resources/polydiacetylene -n 241 -g common_resources/NaCl_564 -ng 4,-2,2,3
#./$name -x common_resources/polydiacetylene -n 141 -g common_resources/NaCl_564 -ng 2,-1,2,3
#./$name -x common_resources/polydiacetylene -n 141 -g common_resources/NaCl_1x1_L2 -ng 2,-1,3,4
#./$name -x common_resources/polydiacetylene_OH -n 141 -g common_resources/NaCl_1x1_L2 -ng 2,-1,3,4
#./$name -x common_resources/polydiacetylene_OH -g common_resources/NaCl_1x1_L2
#./$name -x common_resources/polydiacetylene -g common_resources/NaCl_1x1_L2 -subs 4,common_resources/-COOH.xyz
#./$name -x common_resources/polydiacetylene -g common_resources/NaCl_1x1_L2 -subs 4,common_resources/-COOH.xyz -n 141 -ng 2,-1,3,4 -q 0.05

./$name -x common_resources/formic_dimer

#./$name -x common_resources/oxalate -q 0.08
#./$name -x common_resources/oxalate          -n 211
#./$name -x common_resources/oxalate -q 0.08 -n 211
#./$name -x common_resources/oxalate -q 0.08 -n 221

#./$name -x common_resources/NaCl_mol -q 0.1
#./$name -x common_resources/H_ions

