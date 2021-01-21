# CMakeLists for src-level of samurai

project(samurai)

# Compile flags

set(CMAKE_CXX_STANDARD 11) # C++11
add_compile_options(${ADDITIONAL_FLAGS})

# If needed, generate TDRP Params.cc and Params.hh files
# from their associated paramdef.<app> file

makeTdrpParams()

# Set up the include and link paths:

include_directories(${LROSE_PREFIX}/include)
include_directories(${FFTW_INCLUDES})
include_directories(${NETCDF_INCLUDES})
#include_directories(${GeographicLib_INCLUDE_DIRS})
include_directories(${PETSC_INC})

# Configure the solver include file:

configure_file(solver.inc.in ${CMAKE_CURRENT_SOURCE_DIR}/solver.inc)

# Add the OpenMP / OpenACC flags to certain files:

set(CPROPS "${OpenMP_CXX_FLAGS} ${OpenACC_CXX_FLAGS}")
set_source_files_properties(VarDriver3D.cpp PROPERTIES COMPILE_FLAGS ${CPROPS})
set_source_files_properties(CostFunction3D.cpp PROPERTIES COMPILE_FLAGS ${CPROPS})
set_source_files_properties(CostFunction.cpp PROPERTIES COMPILE_FLAGS ${CPROPS})
set_source_files_properties(CostFunction.h PROPERTIES COMPILE_FLAGS ${CPROPS})
set_source_files_properties(RecursiveFilter.cpp PROPERTIES COMPILE_FLAGS ${CPROPS})

# Set up the GPTL macros:
set_source_files_properties(${gptl_SRCS} PROPERTIES COMPILE_FLAGS "-DTHREADED_OMP -DHAVE_NANOTIME ")

# We release an executable and both static and shared libraries

#add_library(samLibShared SHARED ${APP_TYPE} ${lib_SRCS} ${UIS} ${RSCS} ${TRS} ${MOCS})
#add_library(samLibStatic STATIC ${APP_TYPE} ${lib_SRCS} ${UIS} ${RSCS} ${TRS} ${MOCS})

if(NOT TARGET OpenMP::OpenMP_CXX)
  add_library(OpenMP_TARGET INTERFACE)
  add_library(OpenMP::OpenMP_CXX ALIAS OpenMP_TARGET)
  target_compile_options(OpenMP_TARGET INTERFACE ${OpenMP_CXX_FLAGS})
  find_package(Threads REQUIRED)
  target_link_libraries(OpenMP_TARGET INTERFACE Threads::Threads)
  target_link_libraries(OpenMP_TARGET INTERFACE ${OpenMP_CXX_FLAGS})
endif()

set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR}/release/bin)
set(LIBRARY_OUTPUT_PATH ${CMAKE_BINARY_DIR}/release/lib)

# source lists

set(samurai_HDRS
  BandedMatrix.h 
  BkgdAdapter.h
  BkgdObsLoaders.h
  BSpline.h 
  CostFunction.h 
  CostFunction3D.h	
  CostFunctionXYZ.h
  CostFunctionXYP.h
  CostFunctionRTZ.h
  CostFunctionCOAMPS.cpp
  Dorade.h
  ErrorData.h
  FrameCenter.h
  MetObs.h 
  Observation.h 
  precision.h 
  Projection.h
  RecursiveFilter.h 
  ReferenceState.h	
  read_dorade.h 
  VarDriver.h
  timers.h
  VarDriver3D.h
  Xml.h
 )

set(library_HDRS
  samurai.h
 )

set(common_SRCS
  BkgdArr.cpp
  BkgdStream.cpp
  BkgdObsLoaders.cpp
  BSpline.cpp
  BSplineD.cpp
  BSplineF.cpp
  CostFunction.cpp
  CostFunction3D.cpp
  CostFunctionXYZ.cpp
  CostFunctionXYP.cpp
  CostFunctionRTZ.cpp
  CostFunctionCOAMPS.cpp
  Dorade.cpp
  ErrorData.cpp
  FrameCenter.cpp
  MetObs.cpp
  mac_debug.xcconfig
  mac_release.xcconfig
  Observation.cpp
  Projection.cpp
  RecursiveFilter.cpp
  ReferenceState.cpp
  timers.cpp
  VarDriver.cpp
  VarDriver3D.cpp
  Xml.cpp
  Args.cpp
  Params.cc	# generated file
  XML/tinyxml2.cpp
  LineSplit.cpp
  datetime.cpp
  FileList.cpp
  HashMap.cpp
 )

set(gptl_SRCS 
  timing/gptl.c
  timing/GPTLget_memusage.c
  timing/GPTLprint_memusage.c
  timing/GPTLutil.c
 )

set(samurai_SRCS
  ${common_SRCS}
  ${gptl_SRCS}
  main.cpp
 )

set(lib_SRCS
  ${common_SRCS}		
  ${gptl_SRCS}
  samurai.cpp
 )

# link

#link_directories (${GeographicLib_LIBRARY_DIRS})
link_directories (${LROSE_PREFIX}/lib)

# Where to install files
      
set(INSTALL_PREFIX $ENV{SAMURAI_INSTALL_DIR})

if(NOT INSTALL_PREFIX)
  set(INSTALL_PREFIX ${LROSE_PREFIX})
endif(NOT INSTALL_PREFIX)

# Set RPATH so that executable finds lrose libraries without needing LD_LIBRARY_PATH
SET(CMAKE_SKIP_BUILD_RPATH FALSE)
SET(CMAKE_BUILD_WITH_INSTALL_RPATH TRUE) 
SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
set(CMAKE_INSTALL_RPATH ${INSTALL_PREFIX}/lib)

# do the compile

add_executable(${PROJECT_NAME} ${APP_TYPE} ${samurai_SRCS} ${UIS} ${RSCS} ${TRS} ${MOCS})
target_compile_definitions(${PROJECT_NAME} PRIVATE IO_BENCHMARK=${IO_BENCHMARK} IO_WRITEOBS=${IO_WRITEOBS})

# NCAR: Unfortunately the PGI compiler doesn't work out of the box.
# These are to fix various issues(maybe there's a better way to do this?):
# Note that because of how PGI treats command-line arguments,
# this is still imperfect - eg, adding '-O1' to a release build doesn't completely
# override the earlier, default '-O3' - rather,
# the '-O3' sets various things on(like vectorization), that are NOT turned
# off by the subsequent '-O1'.
# This is annoying, and I don't currently have a workaround.

if(CMAKE_CXX_COMPILER_ID STREQUAL "PGI")
  if(CMAKE_BUILD_TYPE STREQUAL "Release")
    set(CMAKE_CXX_FLAGS "")
    set(CMAKE_CXX_FLAGS_RELEASE "")
    set_source_files_properties(${samurai_SRCS} PROPERTIES COMPILE_FLAGS "-O3 -fast")
    if(MODE STREQUAL "GPU")
      set_source_files_properties(CostFunction3D.cpp CostFunction.cpp RecursiveFilter.cpp PROPERTIES COMPILE_FLAGS "-O3 -fast ${OpenACC_CXX_FLAGS} -ta=${OpenACC_ACCEL_TARGET}")
      target_link_libraries(${PROJECT_NAME} -ta=${OpenACC_ACCEL_TARGET})
    else()
      set_source_files_properties(CostFunction3D.cpp RecursiveFilter.cpp PROPERTIES COMPILE_FLAGS "-O3 -fast ${OpenMP_CXX_FLAGS}")
    endif()
    set_source_files_properties(BSplineD.cpp PROPERTIES COMPILE_FLAGS "-O1") # BSpineD.cpp has issues above -O1 with PGI 19.3 & 19.9 at least
  else()
    if(MODE STREQUAL "GPU")
      set_source_files_properties(CostFunction3D.cpp CostFunction.cpp RecursiveFilter.cpp PROPERTIES COMPILE_FLAGS "${OpenACC_CXX_FLAGS} -ta=${OpenACC_ACCEL_TARGET}")
      target_link_libraries(${PROJECT_NAME} -ta=${OpenACC_ACCEL_TARGET})
    endif()
  endif()
else() # Intel, GNU:
  target_compile_options(${PROJECT_NAME} PRIVATE -Wall -Wno-unknown-pragmas)
  set_source_files_properties(CostFunction3D.cpp BkgdObsLoaders.cpp RecursiveFilter.cpp PROPERTIES COMPILE_FLAGS "-O3 -Ofast ${OpenMP_CXX_FLAGS} ")
endif()

# If we're in GPU mode, check if we have cuFFT 
if(CUDA_CUFFT_LIBRARIES)
  add_definitions(-DUSE_CUFFTW)
  message("<< Linking with cuFFTW libraries : ${CUDA_TOOLKIT_ROOT_DIR}/lib64/libcufftw.so >>")
  target_link_libraries(${PROJECT_NAME} ${CUDA_TOOLKIT_ROOT_DIR}/lib64/libcufftw.so)
  # add_definitions(-DUSE_CUFFT)
  # message("<< Linking with cuFFT libraries : ${CUDA_TOOLKIT_ROOT_DIR}/lib64/libcufft.so >>")
  # target_link_libraries(${PROJECT_NAME} ${CUDA_TOOLKIT_ROOT_DIR}/lib64/libcufft.so)
endif()

#target_link_libraries(${PROJECT_NAME} Geographic)
target_link_libraries(${PROJECT_NAME} ${FFTW_LIBRARIES})
target_link_libraries(${PROJECT_NAME} ${LROSE_LIBRARIES})
target_link_libraries(${PROJECT_NAME} OpenMP::OpenMP_CXX)
target_link_libraries(${PROJECT_NAME} bz2)
target_link_libraries(${PROJECT_NAME} z)
target_link_libraries(${PROJECT_NAME} curl)

#set_property(TARGET samLibShared PROPERTY OUTPUT_NAME "samurai")
#set_property(TARGET samLibStatic PROPERTY OUTPUT_NAME "samurai")

#target_link_libraries(${PROJECT_NAME} samLibStatic)

target_link_libraries(${PROJECT_NAME} ${OpenACC_CXX_FLAGS})
target_link_libraries(${PROJECT_NAME} ${OpenMP_CXX_LIBRARIES})

#target_link_libraries(samLibShared
#	${HDF5_LIBRARY} ${NETCDF_LIBRARY} ${FFTW_LIBRARY} ${LROSE_LIBRARIES}
#	${NETCDF_CPP_LIBRARY} ${GEOGRAPHIC_LIBRARY} ${OPENMP_LINK})
      
#target_link_libraries(samLibStatic
#	${HDF5_LIBRARY} ${NETCDF_LIBRARY} ${FFTW_LIBRARY} ${LROSE_LIBRARIES}
#	${NETCDF_CPP_LIBRARY} ${GEOGRAPHIC_LIBRARY} ${OPENMP_LINK})

# perform the install

install(TARGETS ${PROJECT_NAME} DESTINATION ${INSTALL_PREFIX}/bin)
install(FILES samurai.h DESTINATION ${INSTALL_PREFIX}/include)

#install(TARGETS samLibShared DESTINATION lib)
#install(TARGETS samLibStatic DESTINATION lib)

message("<< INFO: ${CMAKE_BUILD_TYPE} build for ${MODE} - CG using ${SOLVER_MAXITER} iterations and epsilon of ${SOLVER_CONV_TOL} >>")

