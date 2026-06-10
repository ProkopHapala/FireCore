# this function takes sysname as an input and set a global variable with all the xyz files needed
set_inputs() {
    local mysysname=$1
    if [ "$mysysname" == 'H2O-A1_H2O-D1-y' ] ; then
	files_inputs='H2O-A1_H2O-D1-y.xyz'
    elif [ "$mysysname" == 'CH2O-NH3' ] ; then
	files_inputs='CH2O-A1_NH3-D1-y.xyz CH2O-A1_NH3-D1-z.xyz NH3-D1_CH2O-A1-y.xyz NH3-D1_CH2O-A1-z.xyz NH3-A1_NH3-D1-y.xyz NH3-A1_NH3-D1-z.xyz NH3-D1_NH3-A1-y.xyz NH3-D1_NH3-A1-z.xyz'
    elif [ "$mysysname" == 'C4H3NO2-C4H3NO2' ] ; then
	files_inputs='C4H3NO2-A1_C4H3NO2-D1-y.xyz C4H3NO2-A1_C4H3NO2-D1-z.xyz C4H3NO2-D1_C4H3NO2-A1-y.xyz C4H3NO2-D1_C4H3NO2-A1-z.xyz'
    elif [ "$mysysname" == 'CH2NH-CH2NH' ] ; then
	files_inputs='CH2NH-A1_CH2NH-D1-y.xyz CH2NH-A1_CH2NH-D1-z.xyz CH2NH-D1_CH2NH-A1-y.xyz CH2NH-D1_CH2NH-A1-z.xyz'
    elif [ "$mysysname" == 'H2O-H2O' ] ; then
	files_inputs='H2O-A1_H2O-D1-y.xyz H2O-A1_H2O-D1-z.xyz H2O-D1_H2O-A1-y.xyz H2O-D1_H2O-A1-z.xyz'
    elif [ "$mysysname" == 'HBr-HBr' ] ; then
	files_inputs='HBr-A1_HBr-D1-y.xyz HBr-A1_HBr-D1-z.xyz HBr-D1_HBr-A1-y.xyz HBr-D1_HBr-A1-z.xyz'
    elif [ "$mysysname" == 'HCl-HCl' ] ; then
	files_inputs='HCl-A1_HCl-D1-y.xyz HCl-A1_HCl-D1-z.xyz HCl-D1_HCl-A1-y.xyz HCl-D1_HCl-A1-z.xyz'
    elif [ "$mysysname" == 'HCN-HCN' ] ; then
	files_inputs='HCN-A1_HCN-D1-y.xyz HCN-A1_HCN-D1-z.xyz HCN-D1_HCN-A1-y.xyz HCN-D1_HCN-A1-z.xyz'
    elif [ "$mysysname" == 'HCONH2-HCONH2' ] ; then
	files_inputs='HCONH2-A1_HCONH2-D1-y.xyz HCONH2-A1_HCONH2-D1-z.xyz HCONH2-D1_HCONH2-A1-y.xyz HCONH2-D1_HCONH2-A1-z.xyz'
    elif [ "$mysysname" == 'HCOOH-HCOOH' ] ; then
	files_inputs='HCOOH-A1_HCOOH-D1-y.xyz HCOOH-A1_HCOOH-D1-z.xyz HCOOH-A2_HCOOH-D1-y.xyz HCOOH-A2_HCOOH-D1-z.xyz HCOOH-D1_HCOOH-A1-y.xyz HCOOH-D1_HCOOH-A1-z.xyz HCOOH-D1_HCOOH-A2-y.xyz HCOOH-D1_HCOOH-A2-z.xyz'
    elif [ "$mysysname" == 'HF-HF' ] ; then
	files_inputs='HF-A1_HF-D1-y.xyz HF-A1_HF-D1-z.xyz HF-D1_HF-A1-y.xyz HF-D1_HF-A1-z.xyz'
    elif [ "$mysysname" == 'NH3-NH3' ] ; then
	files_inputs='NH3-A1_NH3-D1-y.xyz NH3-A1_NH3-D1-z.xyz NH3-D1_NH3-A1-y.xyz NH3-D1_NH3-A1-z.xyz'
    else
	echo "ERROR: $mysysname not implemented!"
	exit 1
    fi
}
