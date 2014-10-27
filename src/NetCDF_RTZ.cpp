/* 
 *  NetCDF_RTZ.cpp
 *  samurai
 *
 *  Created by Annette Foerster. 
 *  Adapted from $Id: pres_temp_4D_rd.cpp,v 1.11 2006/08/22 19:22:06 ed Exp $ 
 *  (http://www.unidata.ucar.edu/software/netcdf/docs/netcdf-tutorial)
 *  Copyright 2013 Michael Bell. All rights reserved.
 *
 */

#include "NetCDF.h"
#include "NetCDF_RTZ.h"
#include <iostream>
#include <QString>

NetCDF_RTZ::NetCDF_RTZ()
:NetCDF::NetCDF()
 {
	NRADIUS = NX;
	NTHETA = NY;
   
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
	pip = new float[NALT*NRADIUS*NTHETA];
	thetarhop = new float[NALT*NRADIUS*NTHETA];	
	vbar = new float[NALT*NRADIUS*NTHETA];
	vp = new float[NALT*NRADIUS*NTHETA];	
 }
 
NetCDF_RTZ::~NetCDF_RTZ()
{
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
	delete[] pip;
	delete[] thetarhop;	
	delete[] vbar;
	delete[] vp;	
}

int NetCDF_RTZ::readNetCDF(const char* filename) {

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
	NcVar *pipVar, *thetarhopVar;

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
    if (!(pipVar = dataFile.get_var("PIP")))
    	return NC_ERR;   
    if (!(thetarhopVar = dataFile.get_var("THETARHOP")))
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
    if (!pipVar->set_cur(NREC, 0, 0, 0))
		return NC_ERR;
    if (!thetarhopVar->set_cur(NREC, 0, 0, 0))
		return NC_ERR;	
    if (!vbarVar->set_cur(NREC, 0, 0, 0))
		return NC_ERR;
    if (!vpVar->set_cur(NREC, 0, 0, 0))
		return NC_ERR;	
											 
	if (!uVar->get(u, 1, NALT, NTHETA, NRADIUS))
		return NC_ERR;
  	if (!vVar->get(v, 1, NALT, NTHETA, NRADIUS))
		return NC_ERR;
 	if (!wVar->get(w, 1, NALT, NTHETA, NRADIUS))
		return NC_ERR; 
	if (!dudrVar->get(dudr, 1, NALT, NTHETA, NRADIUS))
		return NC_ERR;
  	if (!dvdrVar->get(dvdr, 1, NALT, NTHETA, NRADIUS))
		return NC_ERR;
 	if (!dwdrVar->get(dwdr, 1, NALT, NTHETA, NRADIUS))
		return NC_ERR; 
	if (!dudtVar->get(dudt, 1, NALT, NTHETA, NRADIUS))
		return NC_ERR;
  	if (!dvdtVar->get(dvdt, 1, NALT, NTHETA, NRADIUS))
		return NC_ERR;
 	if (!dwdtVar->get(dwdt, 1, NALT, NTHETA, NRADIUS))
		return NC_ERR; 
	if (!dudzVar->get(dudz, 1, NALT, NTHETA, NRADIUS))
		return NC_ERR;
  	if (!dvdzVar->get(dvdz, 1, NALT, NTHETA, NRADIUS))
		return NC_ERR;
 	if (!dwdzVar->get(dwdz, 1, NALT, NTHETA, NRADIUS))
		return NC_ERR; 
	if (!rhoaVar->get(rhoa, 1, NALT, NTHETA, NRADIUS))
		return NC_ERR;
  	if (!pibarVar->get(pibar, 1, NALT, NTHETA, NRADIUS))
		return NC_ERR;
 	if (!thetarhobarVar->get(thetarhobar, 1, NALT, NTHETA, NRADIUS))
		return NC_ERR; 	
  	if (!pipVar->get(pip, 1, NALT, NTHETA, NRADIUS))
		return NC_ERR;
 	if (!thetarhopVar->get(thetarhop, 1, NALT, NTHETA, NRADIUS))
		return NC_ERR; 	
  	if (!vbarVar->get(vbar, 1, NALT, NTHETA, NRADIUS))
		return NC_ERR;
 	if (!vpVar->get(vp, 1, NALT, NTHETA, NRADIUS))
		return NC_ERR; 	
									      
   return 0;

}


double NetCDF_RTZ::getValue(const int &i,const int &j,const int &k, const QString &varName)
{
	//Returned values are all in SI units
	double value_out;

	if ((i<0) or (i > NRADIUS-1)) return false;
	if ((j<0) or (j > NTHETA-1)) return false;
	if ((k<0) or (k > NALT-1)) return false;
	
	if (varName=="R") {
	value_out= radius[i]*1000.0;
	} else if (varName=="LAMBDA") {
	value_out= theta[j];	
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
	} else if (varName=="PIP") {
	value_out= pip[k*NTHETA*NRADIUS + j*NTHETA + i];		
	} else if (varName=="THETARHOP") {
	value_out= thetarhop[k*NTHETA*NRADIUS + j*NTHETA + i];	
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


double NetCDF_RTZ::calc_A(const int &i,const int &j,const int &k)
{
	double thetarhobar = this->getValue(i,j,k,(QString)"THETARHOBAR");
	double u = this->getValue(i,j,k,(QString)"U");
	double dudr = this->getValue(i,j,k,(QString)"DUDR");
	double v = this->getValue(i,j,k,(QString)"V");	
	double dudlambda = this->getValue(i,j,k,(QString)"DUDT");
	double w = this->getValue(i,j,k,(QString)"W");	
	double dudz = this->getValue(i,j,k,(QString)"DUDZ");
	double vbar = this->getValue(i,j,k,(QString)"VBAR");	
	double vprime = this->getValue(i,j,k,(QString)"VP");
	double r = this->getValue(i,j,k,(QString)"R");

	double a =(u*dudr+v*dudlambda+w*dudz-2.0*vbar*vprime/r-vprime*vprime/r-f*vprime); //(c_p*thetarhobar);
	return a;	
}


double NetCDF_RTZ::calc_B(const int &i,const int &j,const int &k)
{
	double thetarhobar = this->getValue(i,j,k,(QString)"THETARHOBAR");
	double u = this->getValue(i,j,k,(QString)"U");
	double dvdr = this->getValue(i,j,k,(QString)"DVDR");
	double v = this->getValue(i,j,k,(QString)"V");
	double dvdlambda = this->getValue(i,j,k,(QString)"DVDT");
	double w = this->getValue(i,j,k,(QString)"W");	
	double dvdz = this->getValue(i,j,k,(QString)"DVDZ");
	double r = this->getValue(i,j,k,(QString)"R");

	double b =(u*dvdr+v*dvdlambda+w*dvdz+u*v/r+f*u); //(c_p*thetarhobar);
	double dpipdlambda = this->getDerivative(i,j,k,(QString)"PIP",2);
	double scaled = (c_p*thetarhobar)*dpipdlambda;
	return b;	
}

double NetCDF_RTZ::calc_C(const int &i,const int &j,const int &k)
{
	double thetarhobar = this->getValue(i,j,k,(QString)"THETARHOBAR");
	double u = this->getValue(i,j,k,(QString)"U");
	double dwdr = this->getValue(i,j,k,(QString)"DWDR");
	double v = this->getValue(i,j,k,(QString)"V");	
	double dwdlambda = this->getValue(i,j,k,(QString)"DWDT");
	double w = this->getValue(i,j,k,(QString)"W");	
	double dwdz = this->getValue(i,j,k,(QString)"DWDZ");

	double c =(u*dwdr+v*dwdlambda+w*dwdz); //(c_p*thetarhobar);
	return c;	
}

double NetCDF_RTZ::calc_D(const int &i,const int &j,const int &k)
{
	double u = this->getValue(i,j,k,(QString)"U");
	double w = this->getValue(i,j,k,(QString)"W");	
	double dthetarhobardr = this->getDerivative(i,j,k,(QString)"THETARHOBAR",1);
	double dthetarhobardz = this->getDerivative(i,j,k,(QString)"THETARHOBAR",3);	
	
	double d = -u*dthetarhobardr-w*dthetarhobardz;
	return d;	
}


double NetCDF_RTZ::getDerivative(const int &i,const int &j,const int &k, const QString &var, const int &der)
{
	// input variable "der" specifies direction of derivation: 1=dr, 2=dlambda, 3=dz
	double derivative;
	QString derDir;
	switch ( der ) {
	  case 1:
	    derDir = "R";
	    if (i==0){
	      derivative = (this->getValue(i+1,j,k,var)-this->getValue(i,j,k,var))/(this->getValue(i+1,j,k,derDir)-this->getValue(i,j,k,derDir));
	    } else if (i== NRADIUS-1) {
	      derivative = (this->getValue(i,j,k,var)-this->getValue(i-1,j,k,var))/(this->getValue(i,j,k,derDir)-this->getValue(i-1,j,k,derDir));	    
	    } else {
	      derivative = (this->getValue(i+1,j,k,var)-this->getValue(i-1,j,k,var))/(this->getValue(i+1,j,k,derDir)-this->getValue(i-1,j,k,derDir));	    
	    }
	    break;
	  case 2:
	    derDir = "LAMBDA";
	    if (j==0){
	      derivative = (this->getValue(i,j+1,k,var)-this->getValue(i,j,k,var))/(this->getValue(i,j+1,k,derDir)-this->getValue(i,j,k,derDir));
	    } else if (j== NTHETA-1) {
	      derivative = (this->getValue(i,j,k,var)-this->getValue(i,j-1,k,var))/(this->getValue(i,j,k,derDir)-this->getValue(i,j-1,k,derDir));	    
	    } else {
	      derivative = (this->getValue(i,j+1,k,var)-this->getValue(i,j-1,k,var))/(this->getValue(i,j+1,k,derDir)-this->getValue(i,j-1,k,derDir));	    
	    }
	    derivative = derivative*180/pi;
	    break;
	  case 3:
	  	derDir = "Z";
	    if (k==0){
	      derivative = (this->getValue(i,j,k+1,var)-this->getValue(i,j,k,var))/(this->getValue(i,j,k+1,derDir)-this->getValue(i,j,k,derDir));
	    } else if (k== NALT-1) {
	      derivative = (this->getValue(i,j,k,var)-this->getValue(i,j,k-1,var))/(this->getValue(i,j,k,derDir)-this->getValue(i,j,k-1,derDir));	    
	    } else {
	      derivative = (this->getValue(i,j,k+1,var)-this->getValue(i,j,k-1,var))/(this->getValue(i,j,k+1,derDir)-this->getValue(i,j,k-1,derDir));	    
	    }	  	   
	    break; 
	  default:
		std::cout << "Unknown value for calculating derivative. Valid options are 0, 1 and 2.\n";
		exit(1);
	}
	return derivative;
}
