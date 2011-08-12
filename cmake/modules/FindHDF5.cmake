# - Try to find HDF5
# Once done, this will define
#
#  HDF5_FOUND - system has HDF5
#  HDF5_INCLUDE_DIRS - the HDF5 include directories
#  HDF5_LIBRARIES - link these to use HDF5

include(LibFindMacros)

# Dependencies
#libfind_package(HDF5 HDF5)

# Use pkg-config to get hints about paths
#libfind_pkg_check_modules(HDF5_PKGCONF hdf5)

# Include dir
find_path(HDF5_INCLUDE_DIR
  NAMES hdf5.h
  PATHS ${HDF5_PKGCONF_INCLUDE_DIRS}
)

# Finally the library itself
find_library(HDF5_LIBRARY
  NAMES hdf5
  PATHS ${HDF5_PKGCONF_LIBRARY_DIRS}
)

# Set the include dir variables and the libraries and let libfind_process do the rest.
# NOTE: Singular variables for this library, plural for libraries this this lib depends on.
set(HDF5_PROCESS_INCLUDES HDF5_INCLUDE_DIR HDF5_INCLUDE_DIRS)
set(HDF5_PROCESS_LIBS HDF5_LIBRARY HDF5_LIBRARIES)
libfind_process(HDF5)
