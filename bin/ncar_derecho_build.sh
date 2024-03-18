#!/bin/bash

# Assume that this script is run under the current directory

# If an argument is supplied and if it's 'GPU' (case insensitive), assign that to the variable 
#    in CMakeLists.txt file under the root directory. Default is CPU mode.

MODE=${1^^}
if [ "$MODE" == "GPU" ]; then
  sed -i 's/MODE CPU/MODE GPU/g' CMakeLists.txt
  # Somehow the Release mode will break the GPU code; change back to Debug mode
  sed -i 's/CMAKE_BUILD_TYPE Release/CMAKE_BUILD_TYPE Debug/g' CMakeLists.txt
else
  sed -i 's/MODE GPU/MODE CPU/g' CMakeLists.txt     # This will switch it to CPU if it was GPU
fi

# Load desired modules on Derecho

module purge 
module load ncarenv/23.09
if [ "$MODE" == "GPU" ]; then
    module load nvhpc/23.7
    module load cuda/12.2.1
else
    module load intel/2023.2.1
fi
module load fftw/3.3.10
module load netcdf/4.9.2
module load cmake/3.26.3
module load ncarcompilers/1.0.0

# Export the path to the pre-built LROSE library, which is currently pointed to Jian's directory on Derecho
# Currently the available LROSE library is built by intel/2023.2.1 and nvhpc/23.7

if [ "$MODE" == "GPU" ]; then
    export LROSE_INSTALL_DIR=/glade/work/sunjian/lrose/derecho/nvhpc23.7
else
    export LROSE_INSTALL_DIR=/glade/work/sunjian/lrose/derecho/intel2023.2.1
fi

# Generate an empty folder

if [ ! -d "build" ]; then
    mkdir build
fi
cd build
rm -rf CMakeFiles/
rm CMakeCache.txt

cmake ..

make -j 8 VERBOSE=1
