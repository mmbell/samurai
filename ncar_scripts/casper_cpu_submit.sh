#!/bin/bash -l
#PBS -N SAMURAI
#PBS -A NTDD0004
#PBS -l select=1:ncpus=36:ompthreads=36:mem=350GB
#PBS -q casper
#PBS -l walltime=00:30:00
#PBS -j oe
#PBS -k eod

cd $PBS_O_WORKDIR
cd ..
export SAMURAI_ROOT=$(pwd)

##################
# Build the code #
##################
cd $SAMURAI_ROOT/ncar_scripts
./ncar_build.sh

##############
# Run a case #
##############
suffix="casper_cpu"
for i in beltrami # supercell hurricane # hurricane_4panel
do
  ./ncar_run.sh /glade/campaign/cisl/asap/samurai/cases/preprocessed/${i}_preprocessed.xml >& log_${i}_$suffix
  if [ ! -d  ${i}_${suffix} ]; then
     mkdir ${i}_${suffix}
  fi
  mv samurai* log* timing* ${i}_${suffix}
done
