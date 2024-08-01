/* 
 *  NetCDF.cpp
 *  samurai
 *
 *  Created by Annette Foerster on 8/15/13. 
 *  Adapted from $Id: pres_temp_4D_rd.cpp,v 1.11 2006/08/22 19:22:06 ed Exp $ 
 *  (http://www.unidata.ucar.edu/software/netcdf/docs/netcdf-tutorial)
 *  Copyright 2013 Michael Bell. All rights reserved.
 *
 */


#include "NetCDF.h"
#include <iostream>
#include <QString>

NetCDF::NetCDF(const int& metFile_idim, const int& metFile_jdim, const int& metFile_kdim) :c_p(1005.7), g(9.81), f(5.85e-05), pi(3.14159265358979323846)
 {
	NDIMS = 4;
	NALT = metFile_kdim;
	NX = metFile_idim;      // NRADIUS or NLON
	NY = metFile_jdim;      // NTHETA or NLAT
	NREC = 0;
	NC_ERR = 2;
   
  std::cout << "Constructor \n";
}

NetCDF::~NetCDF() {
}




