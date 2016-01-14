# SAMURAI

## Spline Analysis at Mesoscale Utilizing Radar and Aircraft Instrumentation

SAMURAI is a variational analysis technique developed based primarily on the work of Ooyama (1987), Ooyama (2002) and Gao et al.(2004). The SAMURAI analysis yields a maximum likelihood estimate of the atmospheric state for a given set of observations and error estimates by minimizing a variational cost function. 

The technique has several advantages over traditional objective analysis techniques, including:

 + observational error specifications for different instrumentation
 + more complex observation operators for remote sensing data
 + the addition of  balance constraints such as mass continuity
 + including a priori background estimates of the atmospheric state when available. 

A distinguishing characteristic of the SAMURAI technique compared to other variational solvers is that the analysis can be performed directly in an axisymmetric cylindrical coordinate system or in a 3D Cartesian coordinate system. The two-dimensional solver improves the computational efficiency and minimizes potential errors in mass conservation that arise when interpolating from a three-dimensional domain. 

Another distinguishing characteristic from other variational solvers is the use of a Galerkin approach, which is similar to the Fourier spectral transform but uses the cubic B-spline as a basis Ooyama (2002). The disadvantage of the B-spline basis is that it is not orthogonal and therefore requires an extra matrix  to obtain the spline coefficients, but this is a fair trade-off with its other desirable characteristics. The basis is computationally efficient and continuously differentiable to second order, allowing for efficient, accurate interpolation to observation locations, flexible incorporation of boundary conditions, and high numerical accuracy of kinematic derivatives. 

The analysis is performed in a manner similar to the spectral transform method Machenhauer (1979), transforming to and from the spline coefficients and physical space at each step of the cost function minimization. 

A more technical description of SAMURAI is given in the documentation and in Bell et al. (2012).

Please reference Bell et al. (2012) if you use SAMURAI in your research. A description of the 3D version including analytic tests is currently in preparation, and will serve as a useful citation for v1.x shortly.

## Compilation

To compile, use [CMake] (http://www.cmake.org) from the top-level directory:

     $ cmake .

to create a Makefile or Xcode project for your machine. Run `make` or build via Xcode to create the `samurai` binary.
Run `make install` to install the binary to /usr/local/bin

To install samurai in an alternate location:

     $ cmake -DCMAKE_INSTALL_PREFIX=/some/other/location .

To use an alternate compiler (for example the Intel C++ compiler):

     $ cmake -DCMAKE_CXX_COMPILER="icpc" . 

The program is known to work with GCC and Intel compilers, but has not been tested with other compilers.

A few external libraries are required:

[Geographiclib] (http://geographiclib.sourceforge.net) is used for geolocation of data and map projection.

[Qt] (http://qt.nokia.com/products) is used for helper classes including XML parsing. A graphical user interface based on Qt will be available in a future release.

[NetCDF] (http://www.unidata.ucar.edu/software/netcdf) is used for output of the gridded analysis results.

[cURL] (http://curl.haxx.se/libcurl) and [HDF5] (http://www.hdfgroup.org/HDF5) are prerequisites for NetCDF4.

[FFTW] (http://www.fftw.org/) is used for Fourier filtering in periodic domains

### Notes for Mac OSX Yosemite and El Capitan

[Homebrew] (http://brew.sh) can provide all of the libraries listed above to run SAMURAI on the Mac, but a few tweaks are required.

NetCDF introduced new C++ bindings in later versions, and these are not compiled by default. To install a compatible netCDF version use:

     $ brew install --with-cxx-compat --with-fortran netcdf

The default Apple compiler (clang) does not support openMP and the CMake will fail. SAMURAI can run in serial mode but it is much slower. 
A better option is to install GCC via Homebrew, but some recent versions (ie., 5.2) actually have problems with openMP as well. To fix this problem use:

     $ brew install gcc --without-multilib

You then need to force CMake to use the new compiler, either with the CMakelists.txt file, the command line option given above, or

     $ export CXX="g++-5.2" 

in the terminal will do the same thing.

## Running SAMURAI

The program takes a single argument in the form of an XML configuration file describing the run. Observational data and reference frame information should be placed in a subdirectory specified in the configuration. 

The current interface is available for 3D Cartesian and 3D cylindrical analysis. The 2D Cylindrical interface will be added soon. 

Several utility scripts are included in the `util` subdirectory, and more will be added as available.

See the Wiki page for the User's Manual with more details

## Contributing to SAMURAI

* Check out the latest master to make sure the feature hasn't been implemented or the bug hasn't been fixed yet
* Check out the [issue tracker](http://github.com/mmbell/samurai/issues) to make sure someone already hasn't requested it and/or contributed it
* Fork the project
* Start a feature/bugfix branch
* Commit and push until you are happy with your contribution, then open a Pull Request.
* Make sure to add tests for the feature/bugfix. This is important so I don't break it in a future version unintentionally.

## Copyright

Copyright (c) 2012 Michael Bell

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

See LICENSE for a copy of the GNU General Public License, or [gnu.org] (http://www.gnu.org/licenses/).

### References
Bell, M. M., M. T. Montgomery, and K. A. Emanuel, 2012: Air-sea enthalpy and momentum exchange at major hurricane wind speeds observed during CBLAST. *J. Atmos. Sci.*, **69**, 3197-3122.

Gao, J., M. Xue, K. Brewster, and K. K. Droegemeier, 2004: A three-dimensional variational data analysis method with recursive filter for Doppler radars. *J. Atmos. Oceanic Technol.*, **21**, 457–469.

Machenhauer, B., 1979: The spectral method. *Numerical methods used in atmospheric models*, A. Kasahara, Ed., GARP Publications Series No 17, WMO and ICSU, Vol. 2, pp. 121–275.

Ooyama, K. V., 1987: Scale controlled objective analysis. *Mon. Wea. Rev.*, **115**, 2479–2506.

Ooyama, K. V., 2002: The cubic-spline transform method: Basic definitions and tests in a 1d single domain. *Mon. Wea. Rev.*, **130**, 2392–2415.

