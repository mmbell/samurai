#!/bin/bash

if [[ ! -v SAMURAI_ROOT ]]; then
    echo "Please set the SAMURAI_ROOT environment variable first."
    echo "Use export SAMURAI_ROOT=/path_to_samurai_root_directory"
    exit 911
fi
cd $SAMURAI_ROOT

# Check the number of arguments
if [ "$#" -ne 2 ]; then
    echo "Error: Incorrect number of arguments provided."
    echo "Usage: $0 [ID] [INPUT_NAMELIST]"
    exit 911
fi
ID=${1}
INPUT_NAMELIST=${2}

# Generate a "run" folder
if [ -d "run_${ID}" ]; then
    rm -rf run_${ID}
fi
mkdir run_${ID}
cd run_${ID}

EXE="../bin/samurai"

echo ${EXE}
# Determine if we're running on the GPU or not by checking the libraries we've linked:
GPU=$(grep -c libacc ${EXE})
echo ${GPU}
MODE="GPU"
module --force purge
if [ "$NCAR_HOST" == "casper" ]; then
    module load ncarenv/23.10
else
    # Assume it is Derecho
    module load ncarenv/23.09
fi
module reset
if [ "$MODE" == "GPU" ]; then
   if [ "$NCAR_HOST" == "casper" ]; then 
      module load nvhpc/23.7
   else
      # Assume it is Derecho
      module load nvhpc/24.3
   fi
   module load cuda/12.2.1
else
    module load intel/2023.2.1
fi
module load fftw/3.3.10
module load netcdf/4.9.2
module load cmake/3.26.3
module load ncarcompilers/1.0.0
module list
# Export the path to the pre-built LROSE library, which is currently pointed to Jian's directory
if [ "$NCAR_HOST" == "casper" ]; then
    if [ "$MODE" == "GPU" ]; then
        export LROSE_INSTALL_DIR=/glade/work/sunjian/lrose/casper/nvhpc/23.7/lrose-core-20240525
    else
        export LROSE_INSTALL_DIR=/glade/work/sunjian/lrose/casper/intel/2023.2.1/lrose-core-20240525
    fi
else
    # Assume it is Derecho
    if [ "$MODE" == "GPU" ]; then
        export LROSE_INSTALL_DIR=/glade/work/sunjian/lrose/derecho/nvhpc/24.3/lrose-core-20240525
    else
        export LROSE_INSTALL_DIR=/glade/work/sunjian/lrose/derecho/intel/2023.2.1/lrose-core-20240525
    fi
fi

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
  nvidia-smi
  echo "GPU run; using 1 thread and 1 GPU"
  export OMP_NUM_THREADS=1
fi

${EXE} -params $INPUT_NAMELIST
