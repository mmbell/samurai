# - Try to find NETCDF
# Once done, this will define
#
#  NETCDF_FOUND - system has NETCDF
#  NETCDF_INCLUDE_DIRS - the NETCDF include directories
#  NETCDF_LIBRARIES - link these to use NETCDF

include(LibFindMacros)

# Dependencies
#libfind_package(NETCDF)

# Use pkg-config to get hints about paths
libfind_pkg_check_modules(NETCDF_PKGCONF netcdf)

# Include dir
find_path(NETCDF_INCLUDE_DIR
  NAMES netcdfcpp.h
  PATHS ${NETCDF_PKGCONF_INCLUDE_DIRS}
)

# Finally the libraries
find_library(NETCDF_LIBRARY
  NAMES netcdf
  PATHS ${NETCDF_PKGCONF_LIBRARY_DIRS}
)

# Finally the libraries
find_library(NETCDF_CPP_LIBRARY
  NAMES netcdf_c++
  PATHS ${NETCDF_PKGCONF_LIBRARY_DIRS}
)

# Set the include dir variables and the libraries and let libfind_process do the rest.
# NOTE: Singular variables for this library, plural for libraries this this lib depends on.
set(NETCDF_PROCESS_INCLUDES NETCDF_INCLUDE_DIR NETCDF_INCLUDE_DIRS)
set(NETCDF_PROCESS_LIBS NETCDF_LIBRARY NETCDF_CPP_LIBRARY NETCDF_LIBRARIES)
libfind_process(NETCDF)
