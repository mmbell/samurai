#!/bin/bash

if [[ ! -v SAMURAI_ROOT ]]; then
    echo "Please set the SAMURAI_ROOT environment variable first."
    echo "Use export SAMURAI=/path_to_samurai_root_directory"
    exit 911
fi
cd $SAMURAI_ROOT

# If an argument is supplied and if it's 'GPU' (case insensitive), assign that to the variable 
#    in CMakeLists.txt file under the root directory. Default is CPU mode.
MODE=${1^^}

# Load desired modules
module purge 
if [ "$NCAR_HOST" == "casper" ]; then
    module load ncarenv/23.10
else
    # Assume it is Derecho
    module load ncarenv/23.09
fi
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

# Export the path to the pre-built LROSE library, which is currently pointed to Jian's directory
# Currently the available LROSE library is built by intel/2023.2.1 and nvhpc/23.7

if [ "$NCAR_HOST" == "casper" ]; then
    if [ "$MODE" == "GPU" ]; then
        export LROSE_INSTALL_DIR=/glade/work/sunjian/lrose/casper/nvhpc23.7
    else
        export LROSE_INSTALL_DIR=/glade/work/sunjian/lrose/casper/intel2023.2.1
    fi
else
    # Assume it is Derecho
    if [ "$MODE" == "GPU" ]; then
        export LROSE_INSTALL_DIR=/glade/work/sunjian/lrose/derecho/nvhpc23.7
    else
        export LROSE_INSTALL_DIR=/glade/work/sunjian/lrose/derecho/intel2023.2.1
    fi
fi

# Generate a "build" folder
if [ -d "build" ]; then
    rm -rf build
fi
mkdir build
cd build

if [ "$MODE" == "GPU" ]; then
   cmake -DUSE_GPU=true -DDEBUG_COMPILE=false ..
else
   cmake -DDEBUG_COMPILE=false ..
fi

make -j 8 VERBOSE=1

ln -sf $(pwd)/release/bin/samurai ../bin
