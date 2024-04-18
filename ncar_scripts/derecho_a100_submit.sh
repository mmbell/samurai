#!/bin/bash -l
#PBS -N SAMURAI
#PBS -A NTDD0004
#PBS -l select=1:ncpus=64:ompthreads=1:mem=100GB:ngpus=1
#PBS -q main
#PBS -l walltime=02:30:00
#PBS -j oe
#PBS -k eod
 
cd $PBS_O_WORKDIR
cd ..
export SAMURAI_ROOT=$(pwd)

##################
# Build the code #
##################
cd ncar_scripts 
./ncar_build.sh gpu

##############
# Run a case #
##############
suffix="derecho_gpu"
for i in beltrami # supercell hurricane # hurricane_4panel
do
  ./ncar_run.sh /glade/campaign/cisl/asap/samurai/cases/preprocessed/${i}_preprocessed.xml >& log_${i}_$suffix
  if [ ! -d  ${i}_${suffix} ]; then
     mkdir ${i}_${suffix}
  fi
  mv $SAMURAI_ROOT/run/samurai* $SAMURAI_ROOT/run/timing* log* ${i}_${suffix}
done
