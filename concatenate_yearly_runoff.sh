#!/bin/bash

#module load  nco/4.7.9-gcc-4.8.5-hdf5-1.8.18-openmpi-2.0.4
module load  nco

CONFIG="eAMUXL025.L121"
RNF_DIR="../nemo_${CONFIG}/RNF"
YEARi=1979
YEARf=2018

#----------------------------------------------------------------------------

for YEAR in $(seq $YEARi $YEARf)
do

  ncrcat ${RNF_DIR}/runoff_${YEAR}_*_${CONFIG}.nc ${RNF_DIR}/runoff_y${YEAR}_${CONFIG}.nc
  if [ -f ${RNF_DIR}/runoff_y${YEAR}_${CONFIG}.nc ]; then
    rm -f ${RNF_DIR}/runoff_${YEAR}_*_${CONFIG}.nc
    echo "${RNF_DIR}/runoff_y${YEAR}_${CONFIG}.nc  [oK]"
  else
    echo "~!@#%^&* ERROR: ${RNF_DIR}/runoff_y${YEAR}_${CONFIG}.nc HAS NOT BEEN CREATED   >>>>> STOP !!"
    exit
  fi

done
