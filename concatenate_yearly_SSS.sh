#!/bin/bash

module load  nco/4.7.9-gcc-4.8.5-hdf5-1.8.18-openmpi-2.0.4

CONFIG="AMUXL025.L75"
SSS_DIR="$SHAREDELMER/input/nemo_${CONFIG}/SSS"
YEARi=1972
YEARf=1973

#----------------------------------------------------------------------------

for YEAR in $(seq $YEARi $YEARf)
do

  ncrcat ${SSS_DIR}/sss_${YEAR}_*_${CONFIG}.nc ${SSS_DIR}/sss_y${YEAR}_${CONFIG}.nc
  if [ -f ${SSS_DIR}/sss_y${YEAR}_${CONFIG}.nc ]; then
    rm -f ${SSS_DIR}/sss_${YEAR}_*_${CONFIG}.nc
    echo "${SSS_DIR}/sss_y${YEAR}_${CONFIG}.nc  [oK]"
  else
    echo "~!@#%^&* ERROR: ${SSS_DIR}/sss_y${YEAR}_${CONFIG}.nc HAS NOT BEEN CREATED   >>>>> STOP !!"
    exit
  fi

done
