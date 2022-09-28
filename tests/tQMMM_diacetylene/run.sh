#!/bin/bash

fireball_cpp_2="../../cpp/Build/apps/MolecularEditor/FireCoreMMFFmini"
#FireCoreVisual="../../cpp/Build/apps/MolecularEditor/FireCoreVisual"
FireCoreVisual="../../cpp/Build/apps/MolecularEditor/FireCoreVisual"
FireCoreVisualOCL="../../cpp/Build/apps_OCL/MolecularEditorOCL/FireCoreVisualOCL"

MKL_PATH=$PATH:/home/prokop/SW/intel/mkl/lib/intel64
LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$MKL_PATH
export $LD_LIBRARY_PATH
#echo $LD_LIBRARY_PATH

rm answer.bas answer.xyz params.dat CHARGES *.out

#$fireball_cpp_2
$FireCoreVisual
#$FireCoreVisualOCL
