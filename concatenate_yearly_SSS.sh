#!/bin/bash

CONFIG="WED12"
SSS_DIR="$SHAREDELMER/input/nemo_WED12/SSS"
YEARi=1993
YEARf=2013

for YEAR in $(seq $YEARi $YEARf)
do

  ncrcat ${SSS_DIR}/sss_${YEAR}_??_${CONFIG}.nc ${SSS_DIR}/sss_y${YEAR}_${CONFIG}.nc
  if [ -f ${SSS_DIR}/sss_y${YEAR}_${CONFIG}.nc ]; then
    rm -f ${SSS_DIR}/sss_${YEAR}_??_${CONFIG}.nc
    echo "${SSS_DIR}/sss_y${YEAR}_${CONFIG}.nc  [oK]"
  else
    echo "~!@#%^&* ERROR: ${SSS_DIR}/sss_y${YEAR}_${CONFIG}.nc HAS NOT BEEN CREATED   >>>>> STOP !!"
    exit
  fi

done
