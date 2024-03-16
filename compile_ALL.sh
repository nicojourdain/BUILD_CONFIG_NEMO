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

## ada :
## export NC_INC=" " ## empty on adapp with modules 1) hdf5/seq/1.8.9(default); 2) netcdf/seq/4.1.3(default); 3) intel/2013.1
## export NC_LIB=" " ## empty on adapp with modules 1) hdf5/seq/1.8.9(default); 2) netcdf/seq/4.1.3(default); 3) intel/2013.1
# other machines :
export NC_INC="-I`nc-config --includedir` `nc-config --fflags`"
export NC_LIB="`nc-config --flibs`"

export GSW_DIR="./GSW-Fortran"

##==============================================================================
echo "NC_INC=${NC_INC}"
echo "NC_LIB=${NC_LIB}"

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
