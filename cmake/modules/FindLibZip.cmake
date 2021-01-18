# Try to find libzip. Once done, this will define:
#
#   LIBZIP_FOUND - variable which returns the result of the search
#   LIBZIP_INCLUDE_DIRS - list of include directories
#   LIBZIP_LIBRARIES - options for the linker

#=============================================================================
# Copyright 2012 Benjamin Eikel
#
# Distributed under the OSI-approved BSD License (the "License");
# see accompanying file Copyright.txt for details.
#
# This software is distributed WITHOUT ANY WARRANTY; without even the
# implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
# See the License for more information.
#=============================================================================
# (To distribute this file outside of CMake, substitute the full
#  License text for the above reference.)

find_package(PkgConfig)
pkg_check_modules(PC_LIBZIP QUIET libzip)

find_path(LIBZIP_INCLUDE_DIR
	zip.h
	HINTS ${PC_LIBZIP_INCLUDEDIR} ${PC_LIBZIP_INCLUDE_DIRS}
)
find_library(LIBZIP_LIBRARY
	zip
	HINTS ${PC_LIBZIP_LIBDIR} ${PC_LIBZIP_LIBRARY_DIRS}
)

set(LIBZIP_INCLUDE_DIRS ${LIBZIP_INCLUDE_DIR})
set(LIBZIP_LIBRARIES ${LIBZIP_LIBRARY})

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(libzip DEFAULT_MSG
	LIBZIP_INCLUDE_DIR
	LIBZIP_LIBRARY
)

mark_as_advanced(
	LIBZIP_INCLUDE_DIR
	LIBZIP_LIBRARY
	)
      
