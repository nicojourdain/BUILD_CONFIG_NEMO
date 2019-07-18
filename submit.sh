#!/bin/bash
#################################################################################
# N. Jourdain, IGE-CNRS, Feb. 2017
#
# purpose: used to run an executable in a batch job
#
#################################################################################

date

# Random number :
RAND=$(( ( RANDOM % 10000 )  + 1 ))

rm -f tmptxp_${RAND}.sh

if [ ! -n "$1" ]; then
  echo "Usage: `basename $0` file_to_execute  [hh] [mem]                    "
  echo "       hh is the walltime in hours, with two digits (default is 02) "
  echo "       mem is the required memory in GB (default is 5GB)            "
  echo "       ex.: ./submit.sh extract_bathy_meter                         "
  echo "            ./submit.sh extract_bathy_meter 05 10                   "
  exit
fi

if [ $# -eq 1 ]; then
  walltime="02:00:00"
  mem=5
elif [ $# -eq 2 ]; then
  walltime="$2:00:00"
  mem=5
elif [ $# -eq 3 ]; then
  walltime="$2:00:00"
  mem=$3
elif [ $# -gt 3 ]; then
  echo "Usage: `basename $0` file_to_execute [hh]"
  echo "       (hh is the walltime in hours, with two digits, default is 02)"
  echo "       ex.: ./submit.sh extract_bathy_meter 05                      "
  exit
fi

echo "walltime=$walltime"
echo "mem=${mem}Gb"

#=====
if [ `hostname | cut -d"." -f2` == "occigen" ]; then

echo "host is occigen"
cat > tmptxp_${RAND}.sh << EOF
#!/bin/bash
#SBATCH -C HSW24
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --threads-per-core=1
#SBATCH -J submit_${1}
#SBATCH -e submit_${1}.e%j
#SBATCH -o submit_${1}.o%j
#SBATCH --mem=${mem}GB
#SBATCH --time=${walltime}
EOF

elif [ `hostname | cut -c 1-3` == "ada" ]; then

echo "host is adapp"
cat > tmptxp_${RAND}.sh << EOF
# @ job_type = serial
# @ requirements = (Feature == "prepost")
# @ wall_clock_limit = ${walltime}
# @ job_name = submit_${1}
# @ output = \$(job_name).\$(jobid)
# @ error = \$(job_name).\$(jobid)
# @ wall_clock_limit = ${walltime}
# @ as_limit = ${mem}.0gb
# @ queue
EOF

else

echo "default host"
echo '#!/bin/bash' > tmptxp_${RAND}.sh
echo " "
echo "WARNING: You may need to add a specific header in submit.sh if `hostname` enables batch jobs"
echo " "

fi
#=====

echo "./$1" >> tmptxp_${RAND}.sh

chmod +x tmptxp_${RAND}.sh

echo "Launching $1 on  `hostname`"

if [ `hostname | cut -d"." -f2` == "occigen" ]; then
  sbatch ./tmptxp_${RAND}.sh
elif [ `hostname | cut -c 1-3` == "ada" ]; then
  llsubmit ./tmptxp_${RAND}.sh
else
  ./tmptxp_${RAND}.sh
fi

rm -f tmptxp_${RAND}.sh

date
