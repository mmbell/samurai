project(samurai)
cmake_minimum_required ( VERSION 3.0 )
set ( CMAKE_BUILD_TYPE Debug )
set ( CMAKE_CXX_FLAGS_DEBUG "-gdwarf-4 -Og" )
# set ( CMAKE_BUILD_TYPE Release )
#set ( CMAKE_CXX_COMPILER "icpc" )
#set ( CMAKE_CXX_COMPILER "g++-4.2" )
#set ( CMAKE_CXX_COMPILER "/home/mbell/bin/g++" )

# RPATH handling
# TODO Plug qt5 header and library finder
# link_directories( /home/bpmelli/.linuxbrew/opt/qt/lib )

# use, i.e. don't skip the full RPATH for the build tree
SET(CMAKE_SKIP_BUILD_RPATH  FALSE)

SET(CMAKE_BUILD_WITH_INSTALL_RPATH TRUE) 

# add the automatically determined parts of the RPATH
# which point to directories outside the build tree to the install RPATH
SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

set ( CMAKE_CXX_FLAGS "-g -std=c++11  -fPIC" )
# set (CMAKE_SHARED_LINKER_FLAGS "-Wl,-rpath,/home/mbell/lib64 -Wl,-rpath,/home/mbell/lib" )
# set (CMAKE_EXE_LINKER_FLAGS "-Wl,-rpath,/home/mbell/lib64 -Wl,-rpath,/home/mbell/lib -Wl,-rpath,/usr/local/lrose/lib")

# set (CMAKE_SHARED_LINKER_FLAGS "-Wl,-rpath,/home/mbell/lib64 -Wl,-rpath,/home/mbell/lib" )
# set (CMAKE_EXE_LINKER_FLAGS "-Wl,-rpath,/home/mbell/lib66 -Wl,-rpath,/home/mbell/lib -Wl,-rpath,/usr/local/lrose/lib")
set (CMAKE_EXE_LINKER_FLAGS "-Wl,-rpath,/usr/local/lrose/lib")

add_definitions ( -Wall )

if (${CMAKE_SYSTEM_NAME} MATCHES "Linux")
	if(${CMAKE_CXX_COMPILER} MATCHES "icpc")
		set ( OPENMP_FLAG "-openmp" )
		set ( OPENMP_LINK "-openmp" )
        else()
                set ( OPENMP_FLAG "-fopenmp" )
                set ( OPENMP_LINK "-lgomp" )
        endif(${CMAKE_CXX_COMPILER} MATCHES "icpc")
elseif (${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
	if(${CMAKE_GENERATOR} MATCHES "Makefile")
		find_package ( OpenMP REQUIRED )
		set ( OPENMP_FLAG "-fopenmp" )
		set ( OPENMP_LINK "-lgomp" )
		#set ( APP_TYPE MACOSX_BUNDLE )
	else()
		set ( OPENMP_FLAG "-fopenmp" )
	        set ( OPENMP_LINK "" )
	endif(${CMAKE_GENERATOR} MATCHES "Makefile")
else()
	set ( OPENMP_FLAG "" )
	set ( OPENMP_LINK "" )
endif()
set ( CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")
set ( CMAKE_PREFIX_PATH 
	${CMAKE_PREFIX_PATH}
	/sw
	/opt
	/opt/local
	/Users/$ENV{USER}/Development
	/home/bpmelli/release/lib
#	/home/mbell/include
#	/home/mbell/lib
	./include
	./lib
	)

# Tell CMake to run moc automatically when needed
set(CMAKS_AUTOMOC ON)
# Find the QT Widget library
# find_package(Qt5Widgets)
find_package(Qt5 COMPONENTS Core Widgets Xml REQUIRED)

find_package ( HDF5 REQUIRED )
find_package ( NetCDF REQUIRED )
find_package ( GeographicLib REQUIRED )
find_package ( FFTW REQUIRED )
# include ( ${QT_USE_FILE} )

set ( LROSE_INCLUDE_DIR	    "${LROSE_PREFIX}/include" )
# set ( QT_INCLUDE_DIR        "${QT5_PREFIX}/include" )
# set ( QT_QTCORE_INCLUDE_DIR "${QT_INCLUDE_DIR}/QtCore" )
# set ( QT_QTGUI_INCLUDE_DIR  "${QT_INCLUDE_DIR}/QtGui" )
# set ( QT_QTXML_INCLUDE_DIR  "${QT_INCLUDE_DIR}/QtXml" )

set ( LROSE_LIBRARY     "-L ${LROSE_PREFIX}/lib -lkd" )
# set ( QT_QTCORE_LIBRARY "-L ${QT5_PREFIX}/lib -lQt5Core" )
# set ( QT_QTCORE_LIBRARY -lQt5Core )
# set ( QT_QTGUI_LIBRARY -lQt5Gui )
# set ( QT_QTXML_LIBRARY -lQt5Xml )


include_directories (
	${CMAKE_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}
	${HDF5_INCLUDE_DIR} ${NETCDF_INCLUDE_DIR}
	${GEOGRAPHIC_INCLUDE_DIR} ${FFTW_INCLUDE_DIR} ${LROSE_INCLUDE_DIR} )
      
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
	)

set (samurai_SRCS
#	${common_SRCS}	
	src/main.cpp
	)
      
set (lib_SRCS
	${common_SRCS}		
	src/samurai.cpp
	)

set_source_files_properties ( src/CostFunction3D.cpp src/VarDriver3D.cpp src/BkgdObsLoaders.cpp
	PROPERTIES COMPILE_FLAGS ${OPENMP_FLAG}
	)

add_library ( samLibShared SHARED ${APP_TYPE} ${lib_SRCS} ${UIS} ${RSCS} ${TRS} ${MOCS})
      
add_library ( samLibStatic STATIC ${APP_TYPE} ${lib_SRCS} ${UIS} ${RSCS} ${TRS} ${MOCS})
      
add_executable ( samurai ${APP_TYPE} ${samurai_SRCS} ${UIS} ${RSCS} ${TRS} ${MOCS})


target_link_libraries(samurai      Qt5::Widgets Qt5::Xml netcdf netcdf_c++ Geographic fftw3)
target_link_libraries(samLibStatic Qt5::Widgets Qt5::Xml netcdf netcdf_c++ Geographic fftw3)
target_link_libraries(samLibShared Qt5::Widgets Qt5::Xml netcdf netcdf_c++ Geographic fftw3)

set_property ( TARGET samLibShared PROPERTY OUTPUT_NAME "samurai" )
set_property ( TARGET samLibStatic PROPERTY OUTPUT_NAME "samurai" )

target_link_libraries ( samurai samLibStatic )

target_link_libraries ( samLibShared
	${HDF5_LIBRARY} ${NETCDF_LIBRARY} ${FFTW_LIBRARY} ${LROSE_LIBRARY}
	${NETCDF_CPP_LIBRARY} ${GEOGRAPHIC_LIBRARY} ${OPENMP_LINK})
      
target_link_libraries ( samLibStatic
	${HDF5_LIBRARY} ${NETCDF_LIBRARY} ${FFTW_LIBRARY} ${LROSE_LIBRARY}
	${NETCDF_CPP_LIBRARY} ${GEOGRAPHIC_LIBRARY} ${OPENMP_LINK})
install ( TARGETS samurai DESTINATION bin )
install ( FILES samurai.h DESTINATION include )
install ( TARGETS samLibShared DESTINATION lib )
install ( TARGETS samLibStatic DESTINATION lib )
