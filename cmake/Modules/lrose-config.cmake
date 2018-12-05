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
set ( LROSE_LIBRARIES -lkd -ltdrp -lRadx -lhdf5_cpp -lhdf5 -lnetcdf -lNcxx -lnetcdf_c++ )

set ( TDRP_EXECUTABLE ${LROSE_BIN_DIR}/tdrp_gen )
