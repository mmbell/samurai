cmake_minimum_required ( VERSION 3.5 )
project(samurai)

set (CMAKE_CXX_STANDARD 11)
add_definitions (-Wall)

set ( CMAKE_BUILD_TYPE Debug )
set ( CMAKE_CXX_FLAGS_DEBUG "-v -gdwarf-4 -g -std=c++11" )

set (NETCDF_CXX)

# Where to find *Find*.cmake and *-config.cmake

set ( CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/" )
set ( CMAKE_PREFIX_PATH ${CMAKE_MODULE_PATH} )

# Tell CMake to run moc automatically when needed
set(CMAKS_AUTOMOC ON)

find_package(Qt5 COMPONENTS Core Widgets Xml REQUIRED)

find_package ( Lrose COMPONENTS Radx kd tdrp netcdf REQUIRED )
find_package ( LibZip REQUIRED )
find_package ( GeographicLib REQUIRED )
find_package ( FFTW REQUIRED )
find_package ( OpenMP REQUIRED )

# Set RPATH so that executable finds lrose libraries without needing LD_LIBRARY_PATH

SET(CMAKE_SKIP_BUILD_RPATH  FALSE)
SET(CMAKE_BUILD_WITH_INSTALL_RPATH TRUE) 

SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

set ( CMAKE_INSTALL_RPATH ${LROSE_LIB_DIR} )
set ( CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE )

# Use RUNPATH, not RPATH

if (${CMAKE_SYSTEM_NAME} MATCHES "Linux")
  set ( CMAKE_EXE_LINKER_FLAGS "-Wl,--enable-new-dtags" )
endif()

set ( OPENMP_FLAG ${OpenMP_CXX_FLAGS} )

#set(OpenMP_C_FLAGS "-I/usr/local/Cellar/open-mpi/4.0.1_2/include")
#set(OpenMP_C_LIB_NAMES "libmpi")
#set(OpenMP_C_LIB_NAMES "-L/usr/local/Cellar/open-mpi/4.0.1_2/lib -L/usr/local/opt/libevent/lib -lmpi")

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
	)

set (samurai_SRCS
	src/main.cpp
	)

set (lib_SRCS
	${common_SRCS}		
	src/samurai.cpp
	)

# Some soures contain OpenMP pragmas
      
set_source_files_properties ( src/CostFunction3D.cpp src/VarDriver3D.cpp src/BkgdObsLoaders.cpp
	PROPERTIES COMPILE_FLAGS ${OPENMP_FLAG}
	)

include_directories ( /usr/local/include ${LROSE_INCLUDE_DIRS} )
link_directories( /usr/local/lib ${LROSE_LIB_DIR}  /usr/lib/x86_64-linux-gnu/hdf5/serial )

# Generate the argument files

add_custom_command (
  OUTPUT src/Params.hh src/Params.cc
  COMMAND ${TDRP_EXECUTABLE} -c++ -f src/paramdef.samurai && rsync -av Params.* src/
  DEPENDS src/paramdef.samurai
  )

# We release an executable and both static and shared libraries

add_library ( samLibShared SHARED ${APP_TYPE} ${lib_SRCS} ${UIS} ${RSCS} ${TRS} ${MOCS})
add_library ( samLibStatic STATIC ${APP_TYPE} ${lib_SRCS} ${UIS} ${RSCS} ${TRS} ${MOCS})
add_executable ( samurai ${APP_TYPE} ${samurai_SRCS} ${UIS} ${RSCS} ${TRS} ${MOCS})

target_link_libraries(samurai
  Qt5::Widgets Qt5::Xml
  Geographic fftw3
  ${LROSE_LIBRARIES}
  OpenMP::OpenMP_CXX
  bz2 z
  )
target_link_libraries(samLibStatic
  Qt5::Widgets Qt5::Xml
  Geographic fftw3
  ${LROSE_LIBRARIES}   
  OpenMP::OpenMP_CXX
  bz2 z
  )
target_link_libraries(samLibShared
  Qt5::Widgets Qt5::Xml
  Geographic fftw3
  ${LROSE_LIBRARIES}
  OpenMP::OpenMP_CXX
  bz2 z
)  

set_property ( TARGET samLibShared PROPERTY OUTPUT_NAME "samurai" )
set_property ( TARGET samLibStatic PROPERTY OUTPUT_NAME "samurai" )

target_link_libraries ( samurai samLibStatic )

target_link_libraries ( samLibShared
	${HDF5_LIBRARY} ${NETCDF_LIBRARY} ${FFTW_LIBRARY} ${LROSE_LIBRARIES}
	${NETCDF_CPP_LIBRARY} ${GEOGRAPHIC_LIBRARY} ${OPENMP_LINK})
      
target_link_libraries ( samLibStatic
	${HDF5_LIBRARY} ${NETCDF_LIBRARY} ${FFTW_LIBRARY} ${LROSE_LIBRARIES}
	${NETCDF_CPP_LIBRARY} ${GEOGRAPHIC_LIBRARY} ${OPENMP_LINK})

 # Where to install files
      
install ( TARGETS samurai DESTINATION bin )
install ( FILES src/samurai.h DESTINATION include )
install ( TARGETS samLibShared DESTINATION lib )
install ( TARGETS samLibStatic DESTINATION lib )
