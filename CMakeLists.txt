# CMakeLists for samurai

cmake_minimum_required(VERSION 3.0.2)
project(samurai)

#set(CMAKE_VERBOSE_MAKEFILE ON)

# NCAR: Options are 'CPU' or 'GPU' - CPU enables OpenMP, GPU enables OpenACC(default: CPU)
# When invoking cmake, add -DUSE_GPU=true to the args. 
string(TOLOWER "${USE_GPU}" use_gpu_lower )
if(use_gpu_lower STREQUAL "true")
  set(MODE GPU)
else()
  set(MODE CPU)
endif()

# Standard: Options are 'Debug' or 'Release'
string(TOLOWER "${DEBUG_COMPILE}" debug_compile_lower )
if(debug_compile_lower STREQUAL "true")
  set(CMAKE_BUILD_TYPE Debug)
else()
  set(CMAKE_BUILD_TYPE Release)
endif()

# max iterations in the solver(src/CostFunction.cpp) default=2000
set(SOLVER_MAXITER 2000)

# max interations for inner loop with TN solver default=4000
set(SOLVER_INNER_MAXITER 4000)

# 1= truncated Newton solver, 2 = Samurai CG(orig)
set(SOLVER_SAMURAI 1)

# solver convergence tolerance
set(SOLVER_CONV_TOL 1.0e-4)

# solver convergence tolerance for the inner loop with TN(usually = conv_tol for best results)
set(SOLVER_INNER_CONV_TOL 1.0e-4)

# 1 = print out extra solver converegnce info.(0 =  no extra printing)
set(SOLVER_VERBOSE 1)

# epsilon for the Samurai NCG solver
set(SOLVER_CG_EPSILON 1.0e-18)

# specify beta for Samurai NCG solver(1 = PR (orig), 2 = PRP+(recommended), 3 = FR, 4 = DY, 5 = HZ)
set(SOLVER_CG_BETA_TYPE 2)

# type of convergenve for Samurai CG(1 = ~step size(orig), 2 = ||g(X)||/||g(X0)||)
set(SOLVER_CG_CONV_TYPE 1)

# Controls whether to just do the I/O portion of the run; this is just for testing
set(IO_BENCHMARK 0)

# Controls whether we write out the observations after reading them all in(probably safe to be 0)
set(IO_WRITEOBS 0)

# ----- Users shouldn't need to change anything below this ----- 

# Add our own CMake Modules to the search path:
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/modules/")
set(CMAKE_PREFIX_PATH ${CMAKE_MODULE_PATH})

# Finding pkg-config

find_package(PkgConfig REQUIRED)

# Find packages

if(MODE STREQUAL "GPU")
  # For GPU versions
  find_package(OpenACC REQUIRED)
  # NCAR: Target options for GPU version - unused in CPU version  
  set(OpenACC_ACCEL_TARGET "cc80,managed,lineinfo")
  # NCAR: If you want to add custom C++ flags, add them here
  set(ADDITIONAL_FLAGS "-cudalib=cufft")
  # NCAR: If you want to add custom C++ flags, add them here  
  #set(ADDITIONAL_FLAGS "-Minfo=accel")
  # Optionally find CUDA - if found, we'll use CuFFTW, if not, we'll need FFTW
  find_package(CUDA QUIET)
  if(NOT CUDA_CUFFT_LIBRARIES)
    find_package(FFTW REQUIRED)
  endif()
else()
  # For CPU versions(required; run with 1 thread if we want serial times)  
  find_package(OpenMP REQUIRED)
  find_package(FFTW REQUIRED)
  # set(ADDITIONAL_FLAGS "-pg") 
endif()

# LROSE

find_package(LROSE REQUIRED)

#  hdf5 and netcdf

find_package (HDF5)
find_package (NETCDF)
if (DEFINED HDF5_hdf5_LIBRARY_RELEASE)
  get_filename_component(HDF5_INSTALL_PREFIX ${HDF5_hdf5_LIBRARY_RELEASE} DIRECTORY)
endif()
message("netCDF_INSTALL_PREFIX: ${netCDF_INSTALL_PREFIX}")
message("HDF5_INSTALL_PREFIX: ${HDF5_INSTALL_PREFIX}")
message("HDF5_C_INCLUDE_DIR: ${HDF5_C_INCLUDE_DIR}")

# recurse into src directory for the build

add_subdirectory(src)  

