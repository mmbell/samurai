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

NetCDF::NetCDF() :c_p(1005.7), g(9.81), f(0.448432045656147e-05)
 {
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
	v = new float[NALT*NRADIUS*NTHETA];
	w = new float[NALT*NRADIUS*NTHETA];
	dudr = new float[NALT*NRADIUS*NTHETA];
	dvdr = new float[NALT*NRADIUS*NTHETA];
	dwdr = new float[NALT*NRADIUS*NTHETA];
	dudt = new float[NALT*NRADIUS*NTHETA];
	dvdt = new float[NALT*NRADIUS*NTHETA];
	dwdt = new float[NALT*NRADIUS*NTHETA];
	dudz = new float[NALT*NRADIUS*NTHETA];
	dvdz = new float[NALT*NRADIUS*NTHETA];
	dwdz = new float[NALT*NRADIUS*NTHETA];	
	rhoa = new float[NALT*NRADIUS*NTHETA];
	pibar = new float[NALT*NRADIUS*NTHETA];
	thetarhobar = new float[NALT*NRADIUS*NTHETA];	
	vbar = new float[NALT*NRADIUS*NTHETA];
	vp = new float[NALT*NRADIUS*NTHETA];	
}

NetCDF::~NetCDF() {

	delete[] radius;
	delete[] theta;
	delete[] altitude;
	delete[] u;
	delete[] v;
	delete[] w;
	delete[] dudr;
	delete[] dvdr;
	delete[] dwdr;
	delete[] dudt;
	delete[] dvdt;
	delete[] dwdt;
	delete[] dudz;
	delete[] dvdz;
	delete[] dwdz;	
	delete[] rhoa;
	delete[] pibar;
	delete[] thetarhobar;	
	delete[] vbar;
	delete[] vp;	
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
    NcVar *uVar,*vVar,*wVar,*dudrVar,*dvdrVar,*dwdrVar,*dudtVar,*dvdtVar,*dwdtVar,*dudzVar,*dvdzVar,*dwdzVar,*rhoaVar,*pibarVar,*thetarhobarVar,*vbarVar,*vpVar;

    if (!(uVar = dataFile.get_var("U")))
    	return NC_ERR;
    if (!(vVar = dataFile.get_var("V")))
    	return NC_ERR;
    if (!(wVar = dataFile.get_var("W")))
    	return NC_ERR;    	   
    if (!(dudrVar = dataFile.get_var("DUDR")))
    	return NC_ERR;   
    if (!(dvdrVar = dataFile.get_var("DVDR")))
    	return NC_ERR;   
    if (!(dwdrVar = dataFile.get_var("DWDR")))
    	return NC_ERR;   
    if (!(dudtVar = dataFile.get_var("DUDT")))
    	return NC_ERR;   
    if (!(dvdtVar = dataFile.get_var("DVDT")))
    	return NC_ERR;   
    if (!(dwdtVar = dataFile.get_var("DWDT")))
    	return NC_ERR;   
    if (!(dudzVar = dataFile.get_var("DUDZ")))
    	return NC_ERR;   
    if (!(dvdzVar = dataFile.get_var("DVDZ")))
    	return NC_ERR;   
    if (!(dwdzVar = dataFile.get_var("DWDZ")))
    	return NC_ERR;       	
    if (!(rhoaVar = dataFile.get_var("RHOA")))
    	return NC_ERR;   
    if (!(pibarVar = dataFile.get_var("PIBAR")))
    	return NC_ERR;   
    if (!(thetarhobarVar = dataFile.get_var("THETARHOBAR")))
    	return NC_ERR;    
    if (!(vbarVar = dataFile.get_var("VBAR")))
    	return NC_ERR;   
    if (!(vpVar = dataFile.get_var("VP")))
    	return NC_ERR;   
    	
    if (!uVar->set_cur(NREC, 0, 0, 0))
		return NC_ERR;
    if (!vVar->set_cur(NREC, 0, 0, 0))
		return NC_ERR;
    if (!wVar->set_cur(NREC, 0, 0, 0))
		return NC_ERR;
    if (!dudrVar->set_cur(NREC, 0, 0, 0))
		return NC_ERR;
    if (!dvdrVar->set_cur(NREC, 0, 0, 0))
		return NC_ERR;
    if (!dwdrVar->set_cur(NREC, 0, 0, 0))
		return NC_ERR;
    if (!dudtVar->set_cur(NREC, 0, 0, 0))
		return NC_ERR;
    if (!dvdtVar->set_cur(NREC, 0, 0, 0))
		return NC_ERR;
    if (!dwdtVar->set_cur(NREC, 0, 0, 0))
		return NC_ERR;
    if (!dudzVar->set_cur(NREC, 0, 0, 0))
		return NC_ERR;
    if (!dvdzVar->set_cur(NREC, 0, 0, 0))
		return NC_ERR;
    if (!dwdzVar->set_cur(NREC, 0, 0, 0))
		return NC_ERR;		
    if (!rhoaVar->set_cur(NREC, 0, 0, 0))
		return NC_ERR;
    if (!pibarVar->set_cur(NREC, 0, 0, 0))
		return NC_ERR;
    if (!thetarhobarVar->set_cur(NREC, 0, 0, 0))
		return NC_ERR;	
    if (!vbarVar->set_cur(NREC, 0, 0, 0))
		return NC_ERR;
    if (!vpVar->set_cur(NREC, 0, 0, 0))
		return NC_ERR;	
											 
	if (!uVar->get(u, 1, NALT, NRADIUS, NTHETA))
		return NC_ERR;
  	if (!vVar->get(v, 1, NALT, NRADIUS, NTHETA))
		return NC_ERR;
 	if (!wVar->get(w, 1, NALT, NRADIUS, NTHETA))
		return NC_ERR; 
	if (!dudrVar->get(dudr, 1, NALT, NRADIUS, NTHETA))
		return NC_ERR;
  	if (!dvdrVar->get(dvdr, 1, NALT, NRADIUS, NTHETA))
		return NC_ERR;
 	if (!dwdrVar->get(dwdr, 1, NALT, NRADIUS, NTHETA))
		return NC_ERR; 
	if (!dudtVar->get(dudt, 1, NALT, NRADIUS, NTHETA))
		return NC_ERR;
  	if (!dvdtVar->get(dvdt, 1, NALT, NRADIUS, NTHETA))
		return NC_ERR;
 	if (!dwdtVar->get(dwdt, 1, NALT, NRADIUS, NTHETA))
		return NC_ERR; 
	if (!dudzVar->get(dudz, 1, NALT, NRADIUS, NTHETA))
		return NC_ERR;
  	if (!dvdzVar->get(dvdz, 1, NALT, NRADIUS, NTHETA))
		return NC_ERR;
 	if (!dwdzVar->get(dwdz, 1, NALT, NRADIUS, NTHETA))
		return NC_ERR; 
	if (!rhoaVar->get(rhoa, 1, NALT, NRADIUS, NTHETA))
		return NC_ERR;
  	if (!pibarVar->get(pibar, 1, NALT, NRADIUS, NTHETA))
		return NC_ERR;
 	if (!thetarhobarVar->get(thetarhobar, 1, NALT, NRADIUS, NTHETA))
		return NC_ERR; 	
  	if (!vbarVar->get(vbar, 1, NALT, NRADIUS, NTHETA))
		return NC_ERR;
 	if (!vpVar->get(vp, 1, NALT, NRADIUS, NTHETA))
		return NC_ERR; 	
									      
   return 0;

}


double NetCDF::getValue(const int &i,const int &j,const int &k, const QString &varName)
{
	//Returned values are all in SI units
	double value_out;

	if ((i<0) or (i > NRADIUS-1)) return false;
	if ((j<0) or (j > NTHETA-1)) return false;
	if ((k<0) or (k > NALT-1)) return false;
	
	if (varName=="R") {
	value_out= radius[i]*1000.0;
	} else if (varName=="Z") {
	value_out= altitude[k]*1000.0;
	} else if (varName=="U") {
	value_out= u[k*NTHETA*NRADIUS + j*NTHETA + i];
	} else if (varName=="V") {
	value_out= v[k*NTHETA*NRADIUS + j*NTHETA + i];
	} else if (varName=="W") {
	value_out= w[k*NTHETA*NRADIUS + j*NTHETA + i];	
	} else if (varName=="DUDR") {
	value_out= dudr[k*NTHETA*NRADIUS + j*NTHETA + i]/1.0E5;	
	} else if (varName=="DVDR") {
	value_out= dvdr[k*NTHETA*NRADIUS + j*NTHETA + i]/1.0E5;		
	} else if (varName=="DWDR") {
	value_out= dwdr[k*NTHETA*NRADIUS + j*NTHETA + i]/1.0E5;	
	} else if (varName=="DUDT") {
	value_out= dudt[k*NTHETA*NRADIUS + j*NTHETA + i]/1.0E5;	
	} else if (varName=="DVDT") {
	value_out= dvdt[k*NTHETA*NRADIUS + j*NTHETA + i]/1.0E5;		
	} else if (varName=="DWDT") {
	value_out= dwdt[k*NTHETA*NRADIUS + j*NTHETA + i]/1.0E5;	
	} else if (varName=="DUDZ") {
	value_out= dudz[k*NTHETA*NRADIUS + j*NTHETA + i]/1.0E5;	
	} else if (varName=="DVDZ") {
	value_out= dvdz[k*NTHETA*NRADIUS + j*NTHETA + i]/1.0E5;		
	} else if (varName=="DWDZ") {
	value_out= dwdz[k*NTHETA*NRADIUS + j*NTHETA + i]/1.0E5;	
	} else if (varName=="RHOA") {
	value_out= rhoa[k*NTHETA*NRADIUS + j*NTHETA + i];	
	} else if (varName=="PIBAR") {
	value_out= pibar[k*NTHETA*NRADIUS + j*NTHETA + i];		
	} else if (varName=="THETARHOBAR") {
	value_out= thetarhobar[k*NTHETA*NRADIUS + j*NTHETA + i];	
	} else if (varName=="VBAR") {
	value_out= vbar[k*NTHETA*NRADIUS + j*NTHETA + i];		
	} else if (varName=="VP") {
	value_out= vp[k*NTHETA*NRADIUS + j*NTHETA + i];
	} else {
	std::cout << "Requested Variable unknown. \n";
	return false;
	}
	
	return value_out;	
}


double NetCDF::calc_A(const int &i,const int &j,const int &k)
{
	QString var;
	var = "THETARHOBAR";
	double thetarhobar = this->getValue(i,j,k,var);
	var = "U";
	double u = this->getValue(i,j,k,var);
	var = "DUDR";
	double dudr = this->getValue(i,j,k,var);
	var = "V";
	double v = this->getValue(i,j,k,var);	
	var = "DUDT";
	double dudlambda = this->getValue(i,j,k,var);
	var = "W";
	double w = this->getValue(i,j,k,var);	
	var = "DUDZ";
	double dudz = this->getValue(i,j,k,var);
	var = "VBAR";
	double vbar = this->getValue(i,j,k,var);	
	var = "VP";
	double vprime = this->getValue(i,j,k,var);
	var = "R";
	double r = this->getValue(i,j,k,var);

	double a =1.0/c_p/thetarhobar*(u*dudr+v*dudlambda+w*dudz-2.0*vbar*vprime/r-1.0*vprime*vprime/r-f*vprime);
	return a;	
}


double NetCDF::calc_B(const int &i,const int &j,const int &k)
{
	QString var;
	var = "THETARHOBAR";
	double thetarhobar = this->getValue(i,j,k,var);
	var = "U";
	double u = this->getValue(i,j,k,var);
	var = "DVDR";
	double dvdr = this->getValue(i,j,k,var);
	var = "V";
	double v = this->getValue(i,j,k,var);	
	var = "DVDT";
	double dvdlambda = this->getValue(i,j,k,var);
	var = "W";
	double w = this->getValue(i,j,k,var);	
	var = "DVDZ";
	double dvdz = this->getValue(i,j,k,var);
	var = "R";
	double r = this->getValue(i,j,k,var);

	double b =1.0*r/c_p/thetarhobar*(u*dvdr+v*dvdlambda+w*dvdz+1.0*u*v/r+f*u);
	return b;	
}

double NetCDF::calc_C(const int &i,const int &j,const int &k)
{
	QString var;
	var = "THETARHOBAR";
	double thetarhobar = this->getValue(i,j,k,var);
	var = "U";
	double u = this->getValue(i,j,k,var);
	var = "DWDR";
	double dwdr = this->getValue(i,j,k,var);
	var = "V";
	double v = this->getValue(i,j,k,var);	
	var = "DWDT";
	double dwdlambda = this->getValue(i,j,k,var);
	var = "W";
	double w = this->getValue(i,j,k,var);	
	var = "DWDZ";
	double dwdz = this->getValue(i,j,k,var);

	double c =1.0/c_p/thetarhobar* (u*dwdr+v*dwdlambda+w*dwdz);
	return c;	
}