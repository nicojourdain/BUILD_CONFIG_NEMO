#!/bin/bash
##==============================================================================
## N. Jourdain, IGE-CNRS, Feb. 2017
##
## purpose : used to compile all fortran routines
##
##==============================================================================

##==============================================================================
## User's setting :

#module list
#module unload netcdf
#module load netcdf/4.1.1-intel

# Fortran compiler :
FC='ifort'

# Netcdf libraries :
#export NC_INC='-I /usr/local/netcdf/intel/4.1.1/include'
#export NC_LIB='-L /usr/local/netcdf/intel/4.1.1/lib -lnetcdf -lnetcdff'
export NC_INC="-I`nc-config --includedir`"
export NC_LIB="-L`nc-config --libs` -lnetcdff"

export GSW_DIR="./GSW-Fortran"

##==============================================================================

for file in src/*f90
do

 filobj=`basename $file | sed -e "s/\.f90/\.o/g"`
 namexe=`basename $file | sed -e "s/\.f90//g"`

 rm -f $namexe $filobj
 $FC -c $NC_INC -I${GSW_DIR}/modules $file
 $FC -o $namexe ${GSW_DIR}/modules/*.o ${GSW_DIR}/toolbox/*.o $filobj $NC_LIB
 
 if [ -f $namexe ]; then
  rm -f $filobj
  echo "$namexe   [oK]"
 else
  echo "~!@#%^&* ERROR : $namexe HAS NOT BEEN CREATED >>>>>>>> stop !!"
  exit
 fi

done
