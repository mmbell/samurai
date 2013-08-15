/*
 *  NetCDF.h
 *
 *  Created by Michael Bell on 5/29/13.
 *  Based on example code from netCDF library
 *  Copyright 2013. All rights reserved.
 *
 */

#ifndef NETCDF_H
#define NETCDF_H

#include <netcdfcpp.h>
	 
class NetCDF  
{

public:
	NetCDF();
	~NetCDF();
	
	int readNetCDF(const char* filename);
	bool getPoint(const int &i,const int &j,const int &k, double &radius_out, double &theta_out, double &alt_out, double &u_out);	
	
private:
	int NDIMS, NALT, NRADIUS, NTHETA, NREC, NC_ERR;

    float* radius;
	float* theta;
	float* altitude;	
	float* u;
	
};


#endif
