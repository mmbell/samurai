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

##################
# Build the code #
##################
cd ..

# Use cc70 for V100 GPU on Casper
sed -i 's/cc80/cc70/g' CMakeLists.txt

./bin/ncar_build.sh gpu

##############
# Run a case #
##############
# Generate a "run" folder
if [ ! -d run ]; then
  mkdir run
fi
cd run
suffix="casper_gpu"
for i in beltrami # supercell hurricane # hurricane_4panel
do
  ../bin/ncar_run.sh /glade/campaign/cisl/asap/samurai/cases/preprocessed/${i}_preprocessed.xml >& log_${i}_$suffix
  mv timing.0 timing_${i}_$suffix
  if [ ! -d  ${i}_${suffix} ]; then
     mkdir ${i}_${suffix}
  fi
  mv samurai* log* timing* ${i}_${suffix}
done
