# User can tell us where to look

set ( LROSE_PREFIX $ENV{LROSE_ROOT_DIR})

# If not, use the default location
if (NOT LROSE_PREFIX)
  set (LROSE_PREFIX "/usr/local/lrose")
endif ( NOT LROSE_PREFIX )

set ( LROSE_INCLUDE_DIRS ${LROSE_PREFIX}/include )
set ( LROSE_LIB_DIR ${LROSE_PREFIX}/lib )
set ( LROSE_BIN_DIR ${LROSE_PREFIX}/bin )
set ( LROSE_DEFINITIONS -L ${LROSE_LIB_DIR} )
set ( LROSE_LIBRARIES -lkd -ltdrp -lRadx -lNcxx -lnetcdf -lhdf5_cpp -lhdf5)

find_program ( TDRP_EXECUTABLE tdrp_gen PATHS ${LROSE_BIN_DIR} /usr/local/bin )

message ( STATUS "lrose-config: CMAKE_INSTALL_PREFIX: " ${CMAKE_INSTALL_PREFIX} )
message ( STATUS "lrose-config: LROSE_ROOT_DIR: " ${LROSE_ROOT_DIR} )
message ( STATUS "lrose-config: LROSE_BIN_DIR: " ${LROSE_BIN_DIR} )
message ( STATUS "lrose-config: TDRP_EXECUTABLE: " ${TDRP_EXECUTABLE} )

