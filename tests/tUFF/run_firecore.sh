#!/bin/bash

mols=("dioxygen" "dinitrogen" "water" "ethane" "formaldehyde" "H-p" "H-h_1" "H-h_2" "H-hh" "N-h" "N-hh" "O-p" "O-h" "HH-p_1" "HH-p_2" "HH-h_1" "HH-h_2" "HH-pp" "HH-hp" "HH-hh" "HH-h-p" "HH-hh-p" "HN-h" "HN-pp" "HN-hp_1" "HN-hp_2" "HN-hh" "HN-h-p" "HO-h" "NN-pp" "NN-hp" "NN-hh" "NO-h" "NO-p" "NO-h-p" "OO-h" "HHH-p" "HHH-h" "HHH-hp" "HHH-hh" "HHH-h-p" "HHN-hp" "HHN-hh" "HHO-p" "HHO-h" "HHO-hp" "HHO-hh" "HHO-h-p_1" "HHO-h-p_2" "HNH-p" "HNH-h" "HNH-hp" "HNH-hh" "HNH-h-p" "HNN-hp" "HNN-hh" "HNO-p" "HNO-h" "HNO-hp" "HNO-hh" "HNO-h-p" "NHO-hp" "NHO-hh" "OHO-p" "OHO-h_1" "OHO-h_2" "OHO-h-p" "NNN-hhh" "NNO-hp" "NNO-hh_1" "NNO-hh_2" "ONO-p" "ethylene" "isobutane" "benzene" "hexadecane" "naphthalene" "adenine" "cytosine" "thymine" "uracil" "guanine" "mod_adenine" "mod_guanine" "mod_uracyl" "penta_hb2_acceptor" "penta_hb2_donor" "penta_hb3_acceptor" "penta_hb3_donor" "penta_hb3_acceptor2" "hexa_hb3_acceptor" "hexa_hb3_donor" "naphta_hb2_acceptor" "naphta_hb2_donor") # "pentino" "gua_no_c2h4" "hb2_no_c2h4" "penta_no_c2h4" "penta_hb3_no_ch2" "penta_hb3_no_c2h4" "penta_para_no_ch2" "penta_para_no_c2h4" "napthn_ch2" "napthn_c2h4" "gua_no_ch2" "hb2_no_ch2" "hb3_no_ch2" "hb3_no_c2h4" "penta_no_ch2" 

cd ../../cpp/Build/libs/Molecular/
make -j4 MMFF_lib
cd - 1> /dev/null
#exit

#sed "s?MYbSimple?False?" run.py | sed "s?MYbConj?False?" > x.py

for (( i=0 ; i < ${#mols[@]} ; i++ ))
do
    mol=${mols[$i]}.log
    echo "run $(( i + 1 )) / ${#mols[@]} $mol"
    #python3 x.py "../data_UFF/xyz/$mol" 1> /dev/null
    python3 x.py "../data_UFF/xyz/$mol" 1> out
    mv out results/$mol
done

rm -f x.py mol.data.firecore

exit
