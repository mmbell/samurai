#!/bin/bash

# If an argument is supplied, if it's 'GPU' (case insensitive), assign that to the variable in CMakeLists.txt 
# (note that this assumes we're running from the source root.)  Default is CPU mode.
MODE=${1^^}
if [ "$MODE" == "GPU" ]; then
  if [ "$COMPILER" != "pgi" ]; then
    echo "NOTE: You currently need the PGI compiler to build in GPU mode.  You have: ${COMPILER}.  Exiting..."
    exit 1
  fi
  sed -i 's/MODE CPU/MODE GPU/g' ./CMakeLists.txt
else
  sed -i 's/MODE GPU/MODE CPU/g' ./CMakeLists.txt     # This will switch it to CPU if it was GPU
fi


# Normally we'd clear our old modules -- this hopefully isn't necessary, as we should be able to build/link with
# multiple compilers.  We'll try this out, and if it fails, add an explicit purge/load sequence again:
# module purge


# One thing we DO need is to load up the LROSE module, which is currently in Brian's directory, not the system-wide
# ones.  So we need to modify MODULEPATH first.  We should have LROSE modules for Intel 18+, GNU 8+ and PGI 19 (and 20?). 
# You should be able to use whatever compiler is in your current environment, and it'll pick up the right LROSE.

MODULEPATH=${MODULEPATH}:/glade/u/home/bdobbins/Software/Modules
module load LROSE

if [[ -z "$CUDA_HOME" ]]; then
  module load fftw
fi

module load cmake

rm -rf CMakeFiles/
rm CMakeCache.txt

cmake .

make -j 4 
