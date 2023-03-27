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

rm valgrind.log
#rm *.bin *.xsf

#./$name -x common_resources/formic_dimer
#./$name -x common_resources/HCOOH -m 10

./$name -x common_resources/HCOOH -m 10 -g common_resources/NaCl_1x1_L2


#valgrind --leak-check=yes --log-file="valgrind.log" ./$name -x common_resources/HCOOH -m 10
#valgrind --leak-check=yes --track-origins=yes  --log-file="valgrind.log" ./$name -x common_resources/HCOOH -m 10


