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
	virtual ~NetCDF();
	
	virtual int readNetCDF(const char* filename)=0;
	virtual double getValue(const int &i,const int &j,const int &k,const QString &varName)=0;	
	virtual double getDerivative(const int &i,const int &j,const int &k, const QString &var, const int &der)=0;
	virtual double calc_A(const int &i,const int &j,const int &k)=0;
	virtual double calc_B(const int &i,const int &j,const int &k)=0;
	virtual double calc_C(const int &i,const int &j,const int &k)=0;
	virtual double calc_D(const int &i,const int &j,const int &k)=0;

	

protected:
	int NDIMS, NALT, NX, NY, NREC, NC_ERR;
		
	const float c_p;
	const float g;
	const double f; //this needs to be made dynamic eventually, here lat assumed to be 22 deg north
	const double pi;
};


#endif
