#!/bin/bash

EXE=$(readlink -f ./build/release/bin/samurai )

# Determine if we're running on the GPU or not by checking the libraries we've linked:
GPU=$(grep -c libacc ${EXE})

module list

if [ "${GPU}" == "0" ]; then
  echo "CPU run; using 36 threads"
  export OMP_NUM_THREADS=36

  # Minor annoyance - Casper doesn't have libcpuset on it, so we can't use omplace if we're on Casper but doing a CPU run:
  SYSTEM=$(hostname | cut -c1-6)
  if [ "${SYSTEM}" == "casper" ]; then
    unset OMPLACE
  else    
    OMPLACE=omplace
  fi
else
  module load cuda/10.1
  echo "GPU run; using 1 thread"
  export OMP_NUM_THREADS=1
  unset OMPLACE
  # OMPLACE=cuda-gdb
  # OMPLACE=nvvp
  #printenv 
fi

${OMPLACE} ${EXE} $*

