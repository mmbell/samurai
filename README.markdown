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

A more technical description of SAMURAI is given in the documentation, and in

Bell, M. M., M. T. Montgomery, and K. A. Emanuel, 2012: Air-sea enthalpy and momentum exchange at major hurricane wind speeds observed during CBLAST. *J. Atmos. Sci.*, *accepted with minor revisions*.

Please reference this publication if you use SAMURAI in your research. A description of the 3D version including analytic tests is currently in preparation.

## Compilation

To compile, use [CMake] (http://www.cmake.org) from the top-level directory:

     $ cmake .

to create a Makefile or Xcode project for your machine. Run `make` or build via Xcode to create the `samurai` binary.

A few external libraries are required:

[Geographiclib] (http://geographiclib.sourceforge.net) is used for geolocation of data and map projection.

[Qt] (http://qt.nokia.com/products) is used for helper classes including XML parsing. A graphical user interface based on Qt will be available in a future release.

[NetCDF] (http://www.unidata.ucar.edu/software/netcdf) is used for output of the gridded analysis results.

[cURL] (http://curl.haxx.se/libcurl) and [HDF5] (http://www.hdfgroup.org/HDF5) are prerequisites for NetCDF4.

## Running SAMURAI

The program takes a single argument in the form of an XML configuration file describing the run. Observational data and reference frame information should be placed in a `vardata` subdirectory. 

The current interface is only available for the 3D Cartesian analysis. The 2D Cylindrical interface will be added soon. More complete documentation on how to run the software will be also be added shortly.

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

Copyright (c) 2011 Michael Bell

See LICENSE for details.

### References
Gao, J., M. Xue, K. Brewster, and K. K. Droegemeier, 2004: A three-dimensional variational data analysis method with recursive filter for Doppler radars. *J. Atmos. Oceanic Technol.*, **21**, 457–469.

Machenhauer, B., 1979: The spectral method. *Numerical methods used in atmospheric models*, A. Kasahara, Ed., GARP Publications Series No 17, WMO and ICSU, Vol. 2, pp. 121–275.

Ooyama, K. V., 1987: Scale controlled objective analysis. *Mon. Wea. Rev.*, **115**, 2479–2506.

Ooyama, K. V., 2002: The cubic-spline transform method: Basic definitions and tests in a 1d single domain. *Mon. Wea. Rev.*, **130**, 2392–2415.

