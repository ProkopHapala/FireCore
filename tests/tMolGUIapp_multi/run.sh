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

#rm *.bin *.xsf

#./$name -x common_resources/formic_dimer
./$name -x common_resources/HCOOH -n 10