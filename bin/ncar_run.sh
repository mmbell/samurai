#!/bin/bash

# This script should be executed under /path_to_samurai_root/run directory

EXE="../build/release/bin/samurai"

# Determine if we're running on the GPU or not by checking the libraries we've linked:
GPU=$(grep -c libacc ${EXE})

module list

if [ "${GPU}" == "0" ]; then
  echo "CPU run; using one cpu core"
else
  echo "GPU run; using one GPU"
fi

${EXE} $*
