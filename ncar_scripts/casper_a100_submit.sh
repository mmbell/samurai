#!/bin/bash -l
#PBS -N SAMURAI
#PBS -A NEOL0013
#PBS -l select=1:ncpus=36:ompthreads=1:mem=700GB:ngpus=1
#PBS -l gpu_type=a100
#PBS -q casper
#PBS -l walltime=00:30:00
#PBS -j oe
#PBS -k eod
 
cd $PBS_O_WORKDIR
cd ..
export SAMURAI_ROOT=$(pwd)

ID=`date '+%Y%m%d%H%M'`
##################
# Build the code #
##################
sed -i 's/cc70/cc80/g' CMakeLists.txt
sed -i 's/cc90/cc80/g' CMakeLists.txt


cd ncar_scripts 
./ncar_build.sh gpu

##############
# Run a case #
##############
suffix="casper_gpu"
for i in beltrami supercell hurricane typhoonChanthu2020  # hurricane_4panel
do
  ./ncar_run.sh $SAMURAI_ROOT/ncar_scripts/TDRP/${i}.tdrp >& log_${i}_$suffix.$ID
  if [ ! -d  ${i}_${suffix} ]; then
     mkdir ${i}_${suffix}
  fi
  mv $SAMURAI_ROOT/run/timing.0 ${i}_${suffix}/timing.$ID
  mv $SAMURAI_ROOT/run/samurai* log* ${i}_${suffix}
done
