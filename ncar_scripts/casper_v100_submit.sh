#!/bin/bash -l
#PBS -N SAMURAI
#PBS -A NTDD0004
#PBS -l select=1:ncpus=36:ompthreads=1:mem=700GB:ngpus=1
#PBS -l gpu_type=v100
#PBS -q casper
#PBS -l walltime=00:30:00
#PBS -j oe
#PBS -k eod
 
cd $PBS_O_WORKDIR
cd ..
export SAMURAI_ROOT=$(pwd)

ID=`date '+%Y%m%d%H%M'
##################
# Build the code #
##################
# Use cc70 for V100 GPU on Casper
sed -i 's/cc80/cc70/g' CMakeLists.txt

cd ncar_scripts 
./ncar_build.sh gpu

##############
# Run a case #
##############
suffix="casper_gpu"
for i in beltrami supercell hurricane # hurricane_4panel
do
  ./ncar_run.sh /glade/campaign/cisl/asap/samurai/cases/preprocessed/${i}_preprocessed.xml >& log_${i}_$suffix.$ID
  if [ ! -d  ${i}_${suffix} ]; then
     mkdir ${i}_${suffix}
  fi
  mv $SAMURAI_ROOT/run/timing.0 ${i}_${suffix}/timing.$ID
  mv $SAMURAI_ROOT/run/samurai* log* ${i}_${suffix}
done
