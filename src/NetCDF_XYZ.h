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
#include <netcdfcpp.h>
#include <QString>
	 
class NetCDF_XYZ : public NetCDF  
{

public:
	NetCDF_XYZ();
  ~NetCDF_XYZ();
  
  int readNetCDF(const char* filename);
	double getValue(const int &i,const int &j,const int &k,const QString &varName);	
	double getDerivative(const int &i,const int &j,const int &k, const QString &var, const int &der);
	double calc_A(const int &i,const int &j,const int &k);
	double calc_B(const int &i,const int &j,const int &k);
	double calc_C(const int &i,const int &j,const int &k);
	double calc_D(const int &i,const int &j,const int &k);
  
  
protected:
  int NLON, NLAT;

  float* longitude;
	float* latitude;
	float* altitude;	
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
	float* azimuth;
	float* radius;
	float* vtbar;
  float* uprime;
  float* vprime;
  
  float* dpiprimedx;
  float* dpiprimedy; 
  float* dpiprimedz; 
  float* thetarhoprime;
	  
	QString varName;  
	
};


#endif
