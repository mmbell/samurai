cmake_minimum_required ( VERSION 3.0.2 )
project(samurai)

# set(CMAKE_VERBOSE_MAKEFILE ON)

set ( MODE GPU )                                   # NCAR: Options are 'CPU' or 'GPU' - CPU enables OpenMP, GPU enables OpenACC (default: CPU)
set ( CMAKE_BUILD_TYPE Release )                     # Standard: Options are 'Debug' or 'Release'


set ( SOLVER_MAXITER 2000 )        # max iterations in the solver (src/CostFunction.cpp) default=2000
set ( SOLVER_INNER_MAXITER 4000 )  # max interations for inner loop with TN solver default=4000
set ( SOLVER_SAMURAI 1)            # 1= truncated Newton solver, 2 = Samurai CG (orig)
set ( SOLVER_CONV_TOL 1.0e-4 )     # solver convergence tolerance
set ( SOLVER_INNER_CONV_TOL 1.0e-4)  # solver convergence tolerance for the inner loop with TN (usually = conv_tol for best results)
set ( SOLVER_VERBOSE 1)            # 1 = print out extra solver converegnce info. (0 =  no extra printing)
set ( SOLVER_CG_EPSILON 1.0e-18 )  # epsilon for the Samurai NCG solver
set ( SOLVER_CG_BETA_TYPE 2 )      # specify beta for Samurai NCG solver (1 = PR  (orig), 2 = PRP+ (recommended), 3 = FR, 4 = DY, 5 = HZ)
set ( SOLVER_CG_CONV_TYPE 1)       # type of convergenve for Samurai CG (1 = ~step size (orig), 2 = ||g(X)||/||g(X0)||)
set ( IO_BENCHMARK 0 )  # Controls whether to just do the I/O portion of the run; this is just for testing
set ( IO_WRITEOBS 0 )   # Controls whether we write out the observations after reading them all in (probably safe to be 0)

# ----- Users shouldn't need to change anything below this, hopefully ----- 


# Code standards:
set ( CMAKE_CXX_STANDARD 11 )     # C++11

# Configure the solver include file:
configure_file ( src/solver.inc.in src/solver.inc )

# Set RPATH so that executable finds lrose libraries without needing LD_LIBRARY_PATH
SET(CMAKE_SKIP_BUILD_RPATH  FALSE)
SET(CMAKE_BUILD_WITH_INSTALL_RPATH TRUE) 
SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

set ( CMAKE_INSTALL_RPATH ${LROSE_LIB_DIR} )
set ( CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE )

# Add our own CMake Modules to the search path:
set ( CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/" )
set ( CMAKE_PREFIX_PATH ${CMAKE_MODULE_PATH} )

# List of packages:
if (MODE STREQUAL "GPU")
  find_package ( OpenACC REQUIRED )                           # For GPU versions
  set ( OpenACC_ACCEL_TARGET "tesla:cc70,managed,lineinfo" )  # NCAR: Target options for GPU version - unused in CPU version
  set ( ADDITIONAL_FLAGS "-Mcudalib=cufft" )    # NCAR: If you want to add custom C++ flags, add them here
  #set ( ADDITIONAL_FLAGS "-Minfo=accel" )    # NCAR: If you want to add custom C++ flags, add them here
  # Optionally find CUDA - if found, we'll use CuFFTW, if not, we'll need FFTW
  find_package ( CUDA QUIET )
	if ( NOT CUDA_CUFFT_LIBRARIES )
    find_package( FFTW REQUIRED )
  endif()
else()
  find_package ( OpenMP REQUIRED )              # For CPU versions (required; run with 1 thread if we want serial times)
  find_package ( FFTW REQUIRED )
  # set ( ADDITIONAL_FLAGS "-pg" ) 
endif()
find_package ( Lrose COMPONENTS Radx kd tdrp netcdf REQUIRED )
find_package ( GeographicLib REQUIRED ) 

# Set up the include and link paths:

include_directories ( ${FFTW_INCLUDES} ${NETCDF_INCLUDES} include ${GeographicLib_INCLUDE_DIRS} ${PETSC_INC} ${LROSE_INCLUDE_DIRS} )
link_directories    ( ${GeographicLib_LIBRARY_DIRS} ${LROSE_LIB_DIR} )

# Add our additional flags:
add_compile_options( ${ADDITIONAL_FLAGS} )

# Add the OpenMP / OpenACC flags to certain files:
set_source_files_properties ( src/VarDriver3D.cpp src/CostFunction3D.cpp src/CostFunction.cpp src/CostFunction3D.h src/CostFunction.h src/RecursiveFilter.cpp PROPERTIES COMPILE_FLAGS "${OpenMP_CXX_FLAGS} ${OpenACC_CXX_FLAGS}" )

# Set up the GPTL macros:
set_source_files_properties ( ${gptl_SRCS} PROPERTIES COMPILE_FLAGS "-DTHREADED_OMP -DHAVE_NANOTIME " )

# We release an executable and both static and shared libraries

#add_library ( samLibShared SHARED ${APP_TYPE} ${lib_SRCS} ${UIS} ${RSCS} ${TRS} ${MOCS})
#add_library ( samLibStatic STATIC ${APP_TYPE} ${lib_SRCS} ${UIS} ${RSCS} ${TRS} ${MOCS})



if(NOT TARGET OpenMP::OpenMP_CXX)
    add_library(OpenMP_TARGET INTERFACE)
    add_library(OpenMP::OpenMP_CXX ALIAS OpenMP_TARGET)
    target_compile_options(OpenMP_TARGET INTERFACE ${OpenMP_CXX_FLAGS})
    find_package(Threads REQUIRED)
    target_link_libraries(OpenMP_TARGET INTERFACE Threads::Threads)
    target_link_libraries(OpenMP_TARGET INTERFACE ${OpenMP_CXX_FLAGS})
endif()

set ( EXECUTABLE_OUTPUT_PATH build/release/bin )
set ( LIBRARY_OUTPUT_PATH build/release/lib )

set ( samurai_HDRS
	src/BandedMatrix.h 
	src/BkgdAdapter.h
	src/BkgdObsLoaders.h
  src/BSpline.h 
  src/CostFunction.h 
  src/CostFunction3D.h	
	src/CostFunctionXYZ.h
	src/CostFunctionXYP.h
	src/CostFunctionRTZ.h
	src/CostFunctionCOAMPS.cpp
  src/Dorade.h
	src/ErrorData.h
  src/FrameCenter.h
	src/MetObs.h 
  src/Observation.h 
  src/precision.h 
	src/Projection.h
  src/RecursiveFilter.h 
	src/ReferenceState.h	
  src/read_dorade.h 
  src/VarDriver.h
	src/timers.h
  src/VarDriver3D.h
	src/Xml.h
	)

set ( library_HDRS
	src/samurai.h
	)

set ( common_SRCS
	src/BkgdArr.cpp
	src/BkgdStream.cpp
	src/BkgdObsLoaders.cpp
	src/BSpline.cpp
	src/BSplineD.cpp
	src/BSplineF.cpp
	src/CostFunction.cpp
	src/CostFunction3D.cpp
	src/CostFunctionXYZ.cpp
	src/CostFunctionXYP.cpp
	src/CostFunctionRTZ.cpp
	src/CostFunctionCOAMPS.cpp
	src/Dorade.cpp
	src/ErrorData.cpp
  src/FrameCenter.cpp
	src/MetObs.cpp
	src/mac_debug.xcconfig
	src/mac_release.xcconfig
	src/Observation.cpp
	src/Projection.cpp
	src/RecursiveFilter.cpp
	src/ReferenceState.cpp
	src/timers.cpp
	src/VarDriver.cpp
	src/VarDriver3D.cpp
	src/Xml.cpp
	src/Args.cpp
	src/Params.cc	# generated file
	src/XML/tinyxml2.cpp
	src/LineSplit.cpp
	src/datetime.cpp
	src/FileList.cpp
  src/HashMap.cpp
	)

set (gptl_SRCS 
  src/timing/gptl.c
  src/timing/GPTLget_memusage.c
  src/timing/GPTLprint_memusage.c
  src/timing/GPTLutil.c
  )

set (samurai_SRCS
  ${common_SRCS}
  ${gptl_SRCS}
	src/main.cpp
	)

set (lib_SRCS
	${common_SRCS}		
  ${gptl_SRCS}
	src/samurai.cpp
	)

add_executable ( samurai ${APP_TYPE} ${samurai_SRCS} ${UIS} ${RSCS} ${TRS} ${MOCS})
target_compile_definitions( samurai PRIVATE IO_BENCHMARK=${IO_BENCHMARK} IO_WRITEOBS=${IO_WRITEOBS} )

# NCAR: Unfortunately the PGI compiler doesn't work out of the box - these are to fix various issues (maybe there's a better way to do this?):
# Note that because of how PGI treats command-line arguments, this is still imperfect - eg, adding '-O1' to a release build doesn't completely
# override the earlier, default '-O3' - rather, the '-O3' sets various things on (like vectorization), that are NOT turned off by the subsequent
# '-O1'.  This is annoying, and I don't currently have a workaround.
if (CMAKE_CXX_COMPILER_ID STREQUAL "PGI")
  if (CMAKE_BUILD_TYPE STREQUAL "Release")
    set ( CMAKE_CXX_FLAGS "" )
    set ( CMAKE_CXX_FLAGS_RELEASE "" )
    set_source_files_properties ( ${samurai_SRCS} PROPERTIES COMPILE_FLAGS "-O3 -fast" )
    if (MODE STREQUAL "GPU")
      set_source_files_properties ( src/CostFunction3D.cpp src/CostFunction.cpp src/RecursiveFilter.cpp PROPERTIES COMPILE_FLAGS "-O3 -fast ${OpenACC_CXX_FLAGS} -ta=${OpenACC_ACCEL_TARGET}" )
      target_link_libraries ( samurai -ta=${OpenACC_ACCEL_TARGET} )
    else()
       set_source_files_properties ( src/CostFunction3D.cpp src/RecursiveFilter.cpp PROPERTIES COMPILE_FLAGS "-O3 -fast ${OpenMP_CXX_FLAGS}" )
    endif()
    set_source_files_properties ( src/BSplineD.cpp PROPERTIES COMPILE_FLAGS "-O1" ) # BSpineD.cpp has issues above -O1 with PGI 19.3 & 19.9 at least
  else()
    if (MODE STREQUAL "GPU")
      set_source_files_properties ( src/CostFunction3D.cpp src/CostFunction.cpp src/RecursiveFilter.cpp PROPERTIES COMPILE_FLAGS "${OpenACC_CXX_FLAGS} -ta=${OpenACC_ACCEL_TARGET}" )
      target_link_libraries ( samurai -ta=${OpenACC_ACCEL_TARGET} )
    endif()
  endif()
else() # Intel, GNU:
  target_compile_options(samurai PRIVATE -Wall -Wno-unknown-pragmas )
  set_source_files_properties ( src/CostFunction3D.cpp src/BkgdObsLoaders.cpp src/RecursiveFilter.cpp PROPERTIES COMPILE_FLAGS "-O3 -Ofast ${OpenMP_CXX_FLAGS} " )
endif()

# If we're in GPU mode, check if we have cuFFT 
if ( CUDA_CUFFT_LIBRARIES )
  add_definitions( -DUSE_CUFFTW )
  message( "<< Linking with cuFFTW libraries : ${CUDA_TOOLKIT_ROOT_DIR}/lib64/libcufftw.so >>")
  target_link_libraries( samurai ${CUDA_TOOLKIT_ROOT_DIR}/lib64/libcufftw.so )
  # add_definitions( -DUSE_CUFFT )
  # message( "<< Linking with cuFFT libraries : ${CUDA_TOOLKIT_ROOT_DIR}/lib64/libcufft.so >>")
  # target_link_libraries( samurai ${CUDA_TOOLKIT_ROOT_DIR}/lib64/libcufft.so )
endif()

    

# Generate the argument files

add_custom_command (
  OUTPUT src/Params.hh src/Params.cc
  COMMAND ${TDRP_EXECUTABLE} -c++ -f src/paramdef.samurai && rsync -av Params.* src/
  DEPENDS src/paramdef.samurai
  )


target_link_libraries(samurai
  Geographic ${FFTW_LIBRARIES}
  ${LROSE_LIBRARIES}
  OpenMP::OpenMP_CXX
  bz2 z curl
  )
#target_link_libraries(samLibStatic
#  Geographic fftw3
#  ${LROSE_LIBRARIES}   
#  OpenMP::OpenMP_CXX
#  #bz2 z
#  )
#target_link_libraries(samLibShared
#  Geographic fftw3
#  ${LROSE_LIBRARIES}
#  OpenMP::OpenMP_CXX
#  #  bz2 z
#)  

#set_property ( TARGET samLibShared PROPERTY OUTPUT_NAME "samurai" )
#set_property ( TARGET samLibStatic PROPERTY OUTPUT_NAME "samurai" )

#target_link_libraries ( samurai samLibStatic )
target_link_libraries ( samurai ${OpenACC_CXX_FLAGS}  ${OpenMP_CXX_LIBRARIES} )

#target_link_libraries ( samLibShared
#	${HDF5_LIBRARY} ${NETCDF_LIBRARY} ${FFTW_LIBRARY} ${LROSE_LIBRARIES}
#	${NETCDF_CPP_LIBRARY} ${GEOGRAPHIC_LIBRARY} ${OPENMP_LINK})
      
#target_link_libraries ( samLibStatic
#	${HDF5_LIBRARY} ${NETCDF_LIBRARY} ${FFTW_LIBRARY} ${LROSE_LIBRARIES}
#	${NETCDF_CPP_LIBRARY} ${GEOGRAPHIC_LIBRARY} ${OPENMP_LINK})

 # Where to install files
      
install ( TARGETS samurai DESTINATION bin )
install ( FILES src/samurai.h DESTINATION include )
#install ( TARGETS samLibShared DESTINATION lib )
#install ( TARGETS samLibStatic DESTINATION lib )


message( "<< INFO: ${CMAKE_BUILD_TYPE} build for ${MODE} - CG using ${SOLVER_MAXITER} iterations and epsilon of ${SOLVER_CONV_TOL} >>")
