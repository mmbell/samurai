/* Adapted from $Id: pres_temp_4D_rd.cpp,v 1.11 2006/08/22 19:22:06 ed Exp $ 
(http://www.unidata.ucar.edu/software/netcdf/docs/netcdf-tutorial)

*/

#include "NetCDF.h"
#include <iostream>
#include <cmath>

NetCDF::NetCDF() {
	NDIMS = 4;
	NALT = 33;
	NRADIUS = 73;
	NTHETA = 73;
	NREC = 0;
	NC_ERR = 2;
	
	radius = new float [NRADIUS];
	theta = new float [NTHETA];
	altitude = new float [NALT];
	u = new float[NALT*NRADIUS*NTHETA];
	
}

NetCDF::~NetCDF() {

	delete[] radius;
	delete[] theta;
	delete[] altitude;
	delete[] u;

}

int NetCDF::readNetCDF(const char* filename) {

	NcError err(NcError::verbose_nonfatal);
    NcFile dataFile(filename, NcFile::ReadOnly);
     
    if(!dataFile.is_valid())
    	return NC_ERR;
     
	// Get pointers to the latitude and longitude variables.
    NcVar *radiusVar, *thetaVar, *altVar;
    if (!(radiusVar = dataFile.get_var("radius")))
    	return NC_ERR;
    if (!(thetaVar = dataFile.get_var("theta")))
        return NC_ERR;
    if (!(altVar = dataFile.get_var("altitude")))
        return NC_ERR;    
     
    // Get the lat/lon data from the file.
    if (!radiusVar->get(radius, NRADIUS))
    	return NC_ERR;
    if (!thetaVar->get(theta, NTHETA))
    	return NC_ERR;
    if (!altVar->get(altitude, NALT))
    	return NC_ERR;

     // Get pointers to the pressure and temperature variables.
    NcVar *uVar;

    if (!(uVar = dataFile.get_var("U")))
    	return NC_ERR;
   

    if (!uVar->set_cur(NREC, 0, 0, 0))
		return NC_ERR;
	 
	if (!uVar->get(u, 1, NALT, NRADIUS, NTHETA))
		return NC_ERR;
      
   return 0;

}


bool NetCDF::getPoint(const int &i,const int &j,const int &k, double &radius_out, double &theta_out, double &alt_out, double &u_out)
{
	if ((i<0) or (i > NRADIUS-1)) return false;
	if ((j<0) or (j > NTHETA-1)) return false;
	if ((k<0) or (k > NALT-1)) return false;
	
	radius_out = radius[i];
	theta_out = theta[j];
	alt_out = altitude[k];
	u_out= u[k*NTHETA*NRADIUS + j*NTHETA + i];

	return true;	
}

