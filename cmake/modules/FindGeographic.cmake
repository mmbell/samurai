# - Try to find GEOGRAPHIC
# Once done, this will define
#
#  GEOGRAPHIC_FOUND - system has GEOGRAPHIC
#  GEOGRAPHIC_INCLUDE_DIRS - the GEOGRAPHIC include directories
#  GEOGRAPHIC_LIBRARIES - link these to use GEOGRAPHIC

include(LibFindMacros)

# Dependencies
#libfind_package(GEOGRAPHIC)

# Use pkg-config to get hints about paths
#libfind_pkg_check_modules(GEOGRAPHIC_PKGCONF geographic)

# Include dir
find_path(GEOGRAPHIC_INCLUDE_DIR
  NAMES GeographicLib
  PATHS ${GEOGRAPHIC_PKGCONF_INCLUDE_DIRS}
)

# Finally the library itself
find_library(GEOGRAPHIC_LIBRARY
  NAMES Geographic
  PATHS ${GEOGRAPHIC_PKGCONF_LIBRARY_DIRS}
)

# Set the include dir variables and the libraries and let libfind_process do the rest.
# NOTE: Singular variables for this library, plural for libraries this this lib depends on.
set(GEOGRAPHIC_PROCESS_INCLUDES GEOGRAPHIC_INCLUDE_DIR GEOGRAPHIC_INCLUDE_DIRS)
set(GEOGRAPHIC_PROCESS_LIBS GEOGRAPHIC_LIBRARY GEOGRAPHIC_LIBRARIES)
libfind_process(GEOGRAPHIC)
