 /*
 *  NetCDF.h
 *  samurai
 *
 *  Created by Annette Foerster.
 *  Based on example code from netCDF library
 *  Copyright 2013 Michael Bell. All rights reserved.
 *
 */

#ifndef NETCDF_XYZ_H
#define NETCDF_XYZ_H

#include "NetCDF.h"
//#include <netcdfcpp.h>
#include <Ncxx/Nc3File.hh>
#include <string>
	 
class NetCDF_XYZ : public NetCDF  
{

public:
    NetCDF_XYZ(const int& metFile_idim, const int& metFile_jdim, const int& metFile_kdim);

 	~NetCDF_XYZ();
  
    int readNetCDF(std::string filename);
	double getValue(const int &i,const int &j,const int &k,const std::string& varName);	
	int getValue(const std::string& varName);
	double getDerivative(const int &i,const int &j,const int &k, const std::string &var, const int &der);
	double calc_A(const int &i,const int &j,const int &k);
	double calc_B(const int &i,const int &j,const int &k);
	double calc_C(const int &i,const int &j,const int &k);
	double calc_D(const int &i,const int &j,const int &k);
    double calc_E(const int &i,const int &j,const int &k);
  
protected:
  	int NLON, NLAT;

  	float* longitude;
	float* latitude;
	float* altitude;	
	int    obtime;
	float* u;
	float* v;
	float* w;
	float* dudx;
	float* dvdx;
	float* dwdx;
	float* dudy;
	float* dvdy;
	float* dwdy;
	float* dudz;
	float* dvdz;
	float* dwdz;
  	float* thetarhobar;
	float* dpibardx;
	float* dpibardy;
	  
	std::string varName;  
	
};


#endif
