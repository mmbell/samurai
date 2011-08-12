project(samurai)
cmake_minimum_required ( VERSION 2.6 )
set ( CMAKE_BUILD_TYPE Release )
add_definitions ( -Wall )
if (${CMAKE_SYSTEM_NAME} MATCHES "Linux")
	set ( CMAKE_CXX_COMPILER "icpc" )
	set ( OPENMP_FLAG "-openmp" )
	set ( OPENMP_LINK "-openmp" )
endif (${CMAKE_SYSTEM_NAME} MATCHES "Linux")
if (${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
	find_package ( OpenMP REQUIRED )
	set ( OPENMP_FLAG "-fopenmp" )
	set ( OPENMP_LINK "-lgomp" )
endif (${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
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
find_package ( CURL REQUIRED )
find_package ( HDF5 REQUIRED )
find_package ( NetCDF REQUIRED )
find_package ( Geographic REQUIRED )
include ( ${QT_USE_FILE} )
include_directories (
	${CMAKE_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}
	${QT_QTCORE_INCLUDE_DIR} ${QT_QTGUI_INCLUDE_DIR} ${QT_QTXML_INCLUDE_DIR}
	${CURL_INCLUDE_DIR} ${HDF5_INCLUDE_DIR} ${NETCDF_INCLUDE_DIR}
	${GEOGRAPHIC_INCLUDE_DIR} )
set ( EXECUTABLE_OUTPUT_PATH build/release )
set ( samurai_HDRS
	BandedMatrix.h 
        BSpline.h 
        CostFunction.h 
        CostFunctionXYZ.h  
        MetObs.h 
        Observation.h 
        RecursiveFilter.h 
        VarDriver.h 
        VarDriverXYZ.h 
        Dorade.h 
        read_dorade.h 
        precision.h 
        FrameCenter.h
	)

set ( samurai_SRCS
	BSpline.cpp
	BSplineD.cpp
	BSplineF.cpp
	CostFunction.cpp
	CostFunctionXYZ.cpp
	main.cpp
	MetObs.cpp
	Observation.cpp
	RecursiveFilter.cpp
	VarDriver.cpp
	VarDriverXYZ.cpp
	Dorade.cpp
        FrameCenter.cpp
	
	mac_debug.xcconfig
	mac_release.xcconfig
	)
set_source_files_properties ( CostFunctionXYZ.cpp
				PROPERTIES COMPILE_FLAGS ${OPENMP_FLAG}
				)
add_executable ( samurai ${samurai_SRCS} ${UIS} ${RSCS} ${TRS} ${MOCS} )
target_link_libraries ( samurai  
	${QT_QTCORE_LIBRARY} ${QT_QTGUI_LIBRARY} ${QT_QTXML_LIBRARY} 
	${CURL_LIBRARY} ${HDF5_LIBRARY} ${NETCDF_LIBRARY}
	${NETCDF_CPP_LIBRARY} ${GEOGRAPHIC_LIBRARY} ${OPENMP_LINK})
