#!/bin/bash

# This script should be executed under /path_to_samurai_root/run directory

EXE="../build/release/bin/samurai"

# Determine if we're running on the GPU or not by checking the libraries we've linked:
GPU=$(grep -c libacc ${EXE})

module list

if [ "${GPU}" == "0" ]; then
   if [ "$NCAR_HOST" == "casper" ]; then
      echo "CPU run; using 36 threads"
      export OMP_NUM_THREADS=36
   else
      # Assume it is Derecho
      echo "CPU run; using 128 threads"
      export OMP_NUM_THREADS=128
   fi
else
  echo "GPU run; using 1 thread and 1 GPU"
  export OMP_NUM_THREADS=1
fi

${EXE} $*
