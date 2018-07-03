#!/bin/bash

CONFIG="AMUXL12"
RNF_DIR="$SHAREDELMER/input/nemo_${CONFIG}/RNF"
YEARi=1972
YEARf=2016

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
