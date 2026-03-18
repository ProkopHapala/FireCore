#!/bin/bash -l
module load openbabel/3.1.1
set -e

# go to calculation directory
here=$(pwd)
cd energy_calc

# purge energy.out in case it contains "Clash!" so aims can run
echo -n > energy.out 

# input output files
infile="variables.in"
outfile="energy.out"
adsorbate_energy="adsorbate_opt_energy.txt"
substrate_energy="substrate_opt_energy.txt"

# free variables:
x=$(cat $infile | head -1 | tail -1)
y=$(cat $infile | head -2 | tail -1)
z=$(cat $infile | head -3 | tail -1)
Rx=$(cat $infile | head -4 | tail -1)
Ry=$(cat $infile | head -5 | tail -1)
Rz=$(cat $infile | head -6 | tail -1)

# run babel to convert aims input geometry file into xyz file
#obabel -ifhiaims adsorbate_opt.in -oxyz -Oadsorbate_opt.xyz &> .dump

# Using a script to build the interface with translational and rotational coordinates as variables
python interface_builder_boss_1.4.py -ia adsorbate_opt.in -is substrate_opt.in -t $x $y -z $z -r $Rx $Ry $Rz -ff aims -ss 6 8 4 -ca 2 -lc 3.632 
 
# run babel run convert aims geometry file into xyz file
#obabel -ifhiaims geometry.in -oxyz -Oacq.xyz &> .dump

# save the xyz structure
#cat acq.xyz >> $here/movie.xyz

# Check for the presence of keyword "Clash!" in the energy.out file, indicating that the interface_builder has created a configuration with adsorbate-substrate distance less than 0.5 Å, and if it does, skip the next steps
if grep -q "Clash!" energy.out ; then
  echo "Skipping energy calculation because energy file contains 'Clash'."
  cd $here
  exit 0
fi

# run FHI-aims with the control file in the directory
srun /scratch/project_2008059/FHI-aims/build_211214/aims.211214.x > acq.out 

# parse aims output file for total energy as well as input files for energies of the isolated substrate and adsorbate
#eline=($(grep "Total energy of the DFT / Hartree-Fock s.c.f. calculation"  acq.out))
#E=$(echo ${eline[11]})
#Eads=$(grep "TOTAL ENERGY (eV)" "$adsorbate_energy" | awk '{print $5}')
#Esub=$(grep "TOTAL ENERGY (eV)" "$substrate_energy" | awk '{print $5}')
#Efree=$(python -c "print(\"%.8f\" % ($Eads + ($Esub)))")
#Esys=$(python -c "print(\"%.8f\" % ($E - ($Efree)))")
#echo $Esys > $outfile

# extract energy values from files and convert to floats
Eads=$(grep "TOTAL ENERGY (eV)" "$adsorbate_energy" | awk '{print $4}')
Esub=$(grep "TOTAL ENERGY (eV)" "$substrate_energy" | awk '{print $4}')
E=$(grep "Total energy of the DFT / Hartree-Fock s.c.f. calculation" acq.out | awk '{print $12}')

# perform arithmetic operations on energy values
Efree=$(echo "$Eads + $Esub" | bc -l)
Esys=$(echo "$E - $Efree" | bc -l)

# write result to output file
echo "$Esys" > "$outfile"

# save aims output file
cat acq.out >> $here/aims.out

# clean
rm -rf acq*
rm -rf geometry.in

# return to original directory
cd $here

