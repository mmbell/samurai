 /*
 *  NetCDF.h
 *  samurai
 *
 *  Created by Annette Foerster on 8/15/13.
 *  Based on example code from netCDF library
 *  Copyright 2013 Michael Bell. All rights reserved.
 *
 */

#ifndef NETCDF_H
#define NETCDF_H

#include <netcdfcpp.h>
#include <QString>
	 
class NetCDF  
{

public:
	NetCDF();
	~NetCDF();
	
	int readNetCDF(const char* filename);
	double getValue(const int &i,const int &j,const int &k,const QString &varName);	
	double calc_A(const int &i,const int &j,const int &k);
	double calc_B(const int &i,const int &j,const int &k);
	double calc_C(const int &i,const int &j,const int &k);
	double calc_D(const int &i,const int &j,const int &k);
	double calcDerivative(const int &i,const int &j,const int &k, const QString &var, const int &der);

	
private:
	int NDIMS, NALT, NRADIUS, NTHETA, NREC, NC_ERR;

    float* radius;
	float* theta;
	float* altitude;	
	float* u;
	float* v;
	float* w;
	float* dudr;
	float* dvdr;
	float* dwdr;
	float* dudt;
	float* dvdt;
	float* dwdt;
	float* dudz;
	float* dvdz;
	float* dwdz;
	float* rhoa;
	float* pibar;
	float* thetarhobar;
	float* vbar;
	float* vp;
	
	QString varName;
		
	const float c_p;
	const float g;
	const double f; //this needs to be made dynamic eventually, here lat assumed to be 22 deg north
	const double pi;
};


#endif
