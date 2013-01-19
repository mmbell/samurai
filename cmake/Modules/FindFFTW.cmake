# - Try to find FFTW
# Once done, this will define
#
#  FFTW_FOUND - system has FFTW
#  FFTW_INCLUDE_DIRS - the FFTW include directories
#  FFTW_LIBRARIES - link these to use FFTW

include(LibFindMacros)

# Dependencies
#libfind_package(FFTW FFTW)

# Use pkg-config to get hints about paths
#libfind_pkg_check_modules(FFTW_PKGCONF hdf5)

# Include dir
find_path(FFTW_INCLUDE_DIR
  NAMES fftw3.h
  PATHS ${FFTW_PKGCONF_INCLUDE_DIRS}
)

# Finally the library itself
find_library(FFTW_LIBRARY
  NAMES fftw3
  PATHS ${FFTW_PKGCONF_LIBRARY_DIRS}
)

# Set the include dir variables and the libraries and let libfind_process do the rest.
# NOTE: Singular variables for this library, plural for libraries this this lib depends on.
set(FFTW_PROCESS_INCLUDES FFTW_INCLUDE_DIR FFTW_INCLUDE_DIRS)
set(FFTW_PROCESS_LIBS FFTW_LIBRARY FFTW_LIBRARIES)
libfind_process(FFTW)
