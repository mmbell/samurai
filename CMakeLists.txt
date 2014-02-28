project(samurai)
cmake_minimum_required ( VERSION 2.6 )
#set ( CMAKE_BUILD_TYPE Debug )
set ( CMAKE_BUILD_TYPE Release )
#set ( CMAKE_CXX_COMPILER "icpc" )
#set ( CMAKE_CXX_COMPILER "g++-4.2" )
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
	./include
	./lib
	)
find_package ( Qt4 REQUIRED ) 
find_package ( HDF5 REQUIRED )
find_package ( NetCDF REQUIRED )
find_package ( Geographic REQUIRED )
find_package ( FFTW REQUIRED )
include ( ${QT_USE_FILE} )
include_directories (
	${CMAKE_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}
	${QT_QTCORE_INCLUDE_DIR} ${QT_QTGUI_INCLUDE_DIR} ${QT_QTXML_INCLUDE_DIR}
	${HDF5_INCLUDE_DIR} ${NETCDF_INCLUDE_DIR}
	${GEOGRAPHIC_INCLUDE_DIR} ${FFTW_INCLUDE_DIR} )
set ( EXECUTABLE_OUTPUT_PATH build/release )
set ( samurai_HDRS
	src/BandedMatrix.h 
        src/BSpline.h 
        src/CostFunction.h 
        src/CostFunction3D.h	
	src/CostFunctionXYZ.h
	src/CostFunctionXYP.h
	src/CostFunctionRTZ.h
	src/MetObs.h 
        src/Observation.h 
        src/RecursiveFilter.h 
        src/VarDriver.h 
        src/VarDriver3D.h
	src/ReferenceState.h	
        src/Dorade.h 
        src/read_dorade.h 
        src/precision.h 
        src/FrameCenter.h
	)

set ( samurai_SRCS
	src/BSpline.cpp
	src/BSplineD.cpp
	src/BSplineF.cpp
	src/CostFunction.cpp
	src/CostFunction3D.cpp
	src/CostFunctionXYZ.cpp
	src/CostFunctionXYP.cpp
	src/CostFunctionRTZ.cpp
	src/main.cpp
	src/MetObs.cpp
	src/Observation.cpp
	src/RecursiveFilter.cpp
	src/VarDriver.cpp
	src/VarDriver3D.cpp
	src/ReferenceState.cpp
	src/Dorade.cpp
        src/FrameCenter.cpp
	
	src/mac_debug.xcconfig
	src/mac_release.xcconfig
	)
set_source_files_properties ( src/CostFunction3D.cpp src/VarDriver3D.cpp
	PROPERTIES COMPILE_FLAGS ${OPENMP_FLAG}
	)
add_executable ( samurai ${APP_TYPE} ${samurai_SRCS} ${UIS} ${RSCS} ${TRS} ${MOCS})
target_link_libraries ( samurai  
	${QT_QTCORE_LIBRARY} ${QT_QTGUI_LIBRARY} ${QT_QTXML_LIBRARY} 
	${HDF5_LIBRARY} ${NETCDF_LIBRARY} ${FFTW_LIBRARY}
	${NETCDF_CPP_LIBRARY} ${GEOGRAPHIC_LIBRARY} ${OPENMP_LINK})
install ( TARGETS samurai DESTINATION bin )
