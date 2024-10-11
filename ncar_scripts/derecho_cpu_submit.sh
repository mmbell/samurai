#!/bin/bash -l
#PBS -N SAMURAI
#PBS -A NEOL0013
#PBS -l select=1:ncpus=128:ompthreads=128:mem=230GB
#PBS -q main
#PBS -l walltime=02:30:00
#PBS -j oe
#PBS -k eod
 
cd $PBS_O_WORKDIR
cd ..
export SAMURAI_ROOT=$(pwd)

echo $SAMURAI_ROOT

ID=`date '+%Y%m%d%H%M'`
##################
# Build the code #
##################
cd ncar_scripts 
./ncar_build.sh

##############
# Run a case #
##############
suffix="derecho_cpu"
for i in beltrami #supercell hurricane typhoonChanthu2020 # hurricane_4panel
do
  ./ncar_run.sh $SAMURAI_ROOT/ncar_scripts/TDRP/${i}.tdrp >& log_${i}_$suffix.$ID
  echo "log_${i}_$suffix.$ID is done"
  if [ ! -d  ${i}_${suffix} ]; then
     mkdir ${i}_${suffix}
  fi
  mv $SAMURAI_ROOT/run/timing.0 ${i}_${suffix}/timing.$ID
  mv $SAMURAI_ROOT/run/samurai* log* ${i}_${suffix}
done
