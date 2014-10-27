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
#include "NetCDF_XYZ.h"
#include <iostream>
#include <QString>
#include <GeographicLib/TransverseMercatorExact.hpp>


NetCDF_XYZ::NetCDF_XYZ()
:NetCDF::NetCDF()
 {
  NLON = NX;
  NLAT = NY; 
   
  longitude = new float [NLON];
	latitude = new float [NLAT];
	altitude = new float [NALT];
	u = new float[NALT*NLON*NLAT];
	v = new float[NALT*NLON*NLAT];
	w = new float[NALT*NLON*NLAT];
	dudx = new float[NALT*NLON*NLAT];
	dvdx = new float[NALT*NLON*NLAT];
	dwdx = new float[NALT*NLON*NLAT];
	dudy = new float[NALT*NLON*NLAT];
	dvdy = new float[NALT*NLON*NLAT];
	dwdy = new float[NALT*NLON*NLAT];
	dudz = new float[NALT*NLON*NLAT];
	dvdz = new float[NALT*NLON*NLAT];
	dwdz = new float[NALT*NLON*NLAT];	
	thetarhobar = new float[NALT*NLON*NLAT];	
  dpibardx = new float[NALT*NLON*NLAT];	
  dpibardy = new float[NALT*NLON*NLAT];	
  azimuth = new float[NALT*NLON*NLAT];	
  radius = new float[NALT*NLON*NLAT];	
  vtbar = new float[NALT*NLON*NLAT];	
  uprime = new float[NALT*NLON*NLAT];	
  vprime = new float[NALT*NLON*NLAT];	

  dpiprimedx = new float[NALT*NLON*NLAT];	
  dpiprimedy = new float[NALT*NLON*NLAT];	   
  dpiprimedz = new float[NALT*NLON*NLAT];	   
	thetarhoprime = new float[NALT*NLON*NLAT];	 
   
 }
 
NetCDF_XYZ::~NetCDF_XYZ()
{
	delete[] longitude;
	delete[] latitude;
	delete[] altitude;
	delete[] u;
	delete[] v;
	delete[] w;
	delete[] dudx;
	delete[] dvdx;
	delete[] dwdx;
	delete[] dudy;
	delete[] dvdy;
	delete[] dwdy;
	delete[] dudz;
	delete[] dvdz;
	delete[] dwdz;	
  delete[] thetarhobar;
  delete[] dpibardx;
  delete[] dpibardy;
  delete[] azimuth;
  delete[] radius;
  delete[] vtbar;
  delete[] uprime;
  delete[] vprime;
 
  delete[] dpiprimedx;
  delete[] dpiprimedy; 
  delete[] dpiprimedz; 
  delete[] thetarhoprime;

}

int NetCDF_XYZ::readNetCDF(const char* filename) {
	NcError err(NcError::verbose_nonfatal);
    NcFile dataFile(filename, NcFile::ReadOnly);
     
    if(!dataFile.is_valid())
    	return NC_ERR;
     
	// Get pointers to the latitude and longitude variables.
    NcVar *lonVar, *latVar, *altVar;
    if (!(lonVar = dataFile.get_var("lon")))
    	return NC_ERR;
    if (!(latVar = dataFile.get_var("lat")))
        return NC_ERR;
    if (!(altVar = dataFile.get_var("z")))
        return NC_ERR;  
     
    // Get the lat/lon data from the file.
    if (!lonVar->get(longitude, NLON))
    	return NC_ERR;
    if (!latVar->get(latitude, NLAT))
    	return NC_ERR;
    if (!altVar->get(altitude, NALT))
    	return NC_ERR;

     // Get pointers to the pressure and temperature variables.
    NcVar *uVar,*vVar,*wVar,*dudxVar,*dvdxVar,*dwdxVar,*dudyVar,*dvdyVar,*dwdyVar,*dudzVar,*dvdzVar,*dwdzVar, *trbVar, *dpibdxVar, *dpibdyVar, *azVar, *rVar, *vtbVar, *upVar, *vpVar, *trpVar, *dpipdxVar, *dpipdyVar, *dpipdzVar;

    if (!(uVar = dataFile.get_var("u")))
    	return NC_ERR;
    if (!(vVar = dataFile.get_var("v")))
    	return NC_ERR;
    if (!(wVar = dataFile.get_var("w")))
    	return NC_ERR;    	     
    if (!(dudxVar = dataFile.get_var("dudx")))
    	return NC_ERR;   
    if (!(dvdxVar = dataFile.get_var("dvdx")))
    	return NC_ERR;   
    if (!(dwdxVar = dataFile.get_var("dwdx")))
    	return NC_ERR;   
    if (!(dudyVar = dataFile.get_var("dudy")))
    	return NC_ERR;   
    if (!(dvdyVar = dataFile.get_var("dvdy")))
    	return NC_ERR;   
    if (!(dwdyVar = dataFile.get_var("dwdy")))
    	return NC_ERR;   
    if (!(dudzVar = dataFile.get_var("dudz")))
    	return NC_ERR;   
    if (!(dvdzVar = dataFile.get_var("dvdz")))
    	return NC_ERR;   
    if (!(dwdzVar = dataFile.get_var("dwdz")))
    	return NC_ERR;       	  
    if (!(trbVar = dataFile.get_var("trb")))
    	return NC_ERR; 
    if (!(dpibdxVar = dataFile.get_var("dpibdx")))
    	return NC_ERR; 
    if (!(dpibdyVar = dataFile.get_var("dpibdy")))
    	return NC_ERR; 
    if (!(azVar = dataFile.get_var("az")))
    	return NC_ERR; 
    if (!(rVar = dataFile.get_var("r")))
    	return NC_ERR;       
    if (!(vtbVar = dataFile.get_var("vtb")))
    	return NC_ERR;       
    if (!(upVar = dataFile.get_var("up")))
    	return NC_ERR;       
    if (!(vpVar = dataFile.get_var("vp")))
    	return NC_ERR; 
      
    if (!(dpipdxVar = dataFile.get_var("dpipdx")))
    	return NC_ERR; 
    if (!(dpipdyVar = dataFile.get_var("dpipdy")))
    	return NC_ERR; 
    if (!(dpipdzVar = dataFile.get_var("dpipdz")))
    	return NC_ERR; 
    if (!(trpVar = dataFile.get_var("trp")))
    	return NC_ERR; 
                	
    if (!uVar->set_cur(NREC, 0, 0, 0))
		return NC_ERR;
    if (!vVar->set_cur(NREC, 0, 0, 0))
		return NC_ERR;
    if (!wVar->set_cur(NREC, 0, 0, 0))
		return NC_ERR;
    if (!dudxVar->set_cur(NREC, 0, 0, 0))
		return NC_ERR;
    if (!dvdxVar->set_cur(NREC, 0, 0, 0))
		return NC_ERR;
    if (!dwdxVar->set_cur(NREC, 0, 0, 0))
		return NC_ERR;
    if (!dudyVar->set_cur(NREC, 0, 0, 0))
		return NC_ERR;
    if (!dvdyVar->set_cur(NREC, 0, 0, 0))
		return NC_ERR;
    if (!dwdyVar->set_cur(NREC, 0, 0, 0))
		return NC_ERR;
    if (!dudzVar->set_cur(NREC, 0, 0, 0))
		return NC_ERR;
    if (!dvdzVar->set_cur(NREC, 0, 0, 0))
		return NC_ERR;
    if (!dwdzVar->set_cur(NREC, 0, 0, 0))
		return NC_ERR;			
    if (!trbVar->set_cur(NREC, 0, 0, 0))
    return NC_ERR; 
    if (!dpibdxVar->set_cur(NREC, 0, 0, 0))
    return NC_ERR; 
    if (!dpibdyVar->set_cur(NREC, 0, 0, 0))
    return NC_ERR; 
    if (!azVar->set_cur(NREC, 0, 0, 0))
    return NC_ERR; 
    if (!rVar->set_cur(NREC, 0, 0, 0))
    return NC_ERR;       
    if (!vtbVar->set_cur(NREC, 0, 0, 0))
    return NC_ERR;       
    if (!upVar->set_cur(NREC, 0, 0, 0))
    return NC_ERR;       
    if (!vpVar->set_cur(NREC, 0, 0, 0))
    	return NC_ERR; 

    if (!trpVar->set_cur(NREC, 0, 0, 0))
    return NC_ERR; 
    if (!dpipdxVar->set_cur(NREC, 0, 0, 0))
    return NC_ERR; 
    if (!dpipdyVar->set_cur(NREC, 0, 0, 0))
    return NC_ERR; 
    if (!dpipdzVar->set_cur(NREC, 0, 0, 0))
    return NC_ERR; 
											 
	if (!uVar->get(u, 1, NALT, NLAT, NLON))
		return NC_ERR;
  	if (!vVar->get(v, 1, NALT, NLAT, NLON))
		return NC_ERR;
 	if (!wVar->get(w, 1, NALT, NLAT, NLON))
		return NC_ERR; 
	if (!dudxVar->get(dudx, 1, NALT, NLAT, NLON))
		return NC_ERR;
  	if (!dvdxVar->get(dvdx, 1, NALT, NLAT, NLON))
		return NC_ERR;
 	if (!dwdxVar->get(dwdx, 1, NALT, NLAT, NLON))
		return NC_ERR; 
	if (!dudyVar->get(dudy, 1, NALT, NLAT, NLON))
		return NC_ERR;
  	if (!dvdyVar->get(dvdy, 1, NALT, NLAT, NLON))
		return NC_ERR;
 	if (!dwdyVar->get(dwdy, 1, NALT, NLAT, NLON))
		return NC_ERR; 
	if (!dudzVar->get(dudz, 1, NALT, NLAT, NLON))
		return NC_ERR;
  	if (!dvdzVar->get(dvdz, 1, NALT, NLAT, NLON))
		return NC_ERR;
 	if (!dwdzVar->get(dwdz, 1, NALT, NLAT, NLON))
		return NC_ERR; 
  if (!trbVar->get(thetarhobar, 1, NALT, NLAT, NLON))
    return NC_ERR; 
  if (!dpibdxVar->get(dpibardx, 1, NALT, NLAT, NLON))
    return NC_ERR; 
  if (!dpibdyVar->get(dpibardy, 1, NALT, NLAT, NLON))
    return NC_ERR; 
  if (!azVar->get(azimuth, 1, NALT, NLAT, NLON))
    return NC_ERR; 
  if (!rVar->get(radius, 1, NALT, NLAT, NLON))
    return NC_ERR;       
  if (!vtbVar->get(vtbar, 1, NALT, NLAT, NLON))
    return NC_ERR;       
  if (!upVar->get(uprime, 1, NALT, NLAT, NLON))
    return NC_ERR;       
  if (!vpVar->get(vprime, 1, NALT, NLAT, NLON))
    return NC_ERR; 	

  if (!trpVar->get(thetarhoprime, 1, NALT, NLAT, NLON))
    return NC_ERR; 
  if (!dpipdxVar->get(dpiprimedx, 1, NALT, NLAT, NLON))
    return NC_ERR; 
  if (!dpipdyVar->get(dpiprimedy, 1, NALT, NLAT, NLON))
    return NC_ERR; 
  if (!dpipdzVar->get(dpiprimedz, 1, NALT, NLAT, NLON))
    return NC_ERR; 
									      
   return 0;
}


double NetCDF_XYZ::getValue(const int &i,const int &j,const int &k, const QString &varName)
{
	//Returned values are all in SI units
	double value_out;

	if ((i<0) or (i > NLON-1)) return false;
	if ((j<0) or (j > NLAT-1)) return false;
	if ((k<0) or (k > NALT-1)) return false;
	
	if (varName=="lon") {
	value_out= longitude[i];
	} else if (varName=="lat") {
	value_out= latitude[j];	
	} else if (varName=="z") {
	value_out= altitude[k]*1000.0;
	} else if (varName=="u") {
	value_out= u[k*NLAT*NLON + j*NLAT + i]*1000.0;
	} else if (varName=="v") {
	value_out= v[k*NLAT*NLON + j*NLAT + i]*1000.0;
	} else if (varName=="w") {
	value_out= w[k*NLAT*NLON + j*NLAT + i]*1000.0;	
	} else if (varName=="dudx") {
	value_out= dudx[k*NLAT*NLON + j*NLAT + i]/1.0E5;	
	} else if (varName=="dvdx") {
	value_out= dvdx[k*NLAT*NLON + j*NLAT + i]/1.0E5;		
	} else if (varName=="dwdx") {
	value_out= dwdx[k*NLAT*NLON + j*NLAT + i]/1.0E5;	
	} else if (varName=="dudy") {
	value_out= dudy[k*NLAT*NLON + j*NLAT + i]/1.0E5;	
	} else if (varName=="dvdy") {
	value_out= dvdy[k*NLAT*NLON + j*NLAT + i]/1.0E5;		
	} else if (varName=="dwdy") {
	value_out= dwdy[k*NLAT*NLON + j*NLAT + i]/1.0E5;	
	} else if (varName=="dudz") {
	value_out= dudz[k*NLAT*NLON + j*NLAT + i]/1.0E5;	
	} else if (varName=="dvdz") {
	value_out= dvdz[k*NLAT*NLON + j*NLAT + i]/1.0E5;		
	} else if (varName=="dwdz") {
	value_out= dwdz[k*NLAT*NLON + j*NLAT + i]/1.0E5;	
	} else if (varName=="trb") {
	value_out= thetarhobar[k*NLAT*NLON + j*NLAT + i];	
	} else if (varName=="dpibdx") {
	value_out= dpibardx[k*NLAT*NLON + j*NLAT + i];	
	} else if (varName=="dpibdy") {
	value_out= dpibardy[k*NLAT*NLON + j*NLAT + i];	 
	} else if (varName=="az") {
	value_out= azimuth[k*NLAT*NLON + j*NLAT + i];	
	} else if (varName=="r") {
	value_out= radius[k*NLAT*NLON + j*NLAT + i];	
	} else if (varName=="vtb") {
	value_out= vtbar[k*NLAT*NLON + j*NLAT + i];	
	} else if (varName=="up") {
	value_out= uprime[k*NLAT*NLON + j*NLAT + i]*1000.0;	
	} else if (varName=="vp") {
	value_out= vprime[k*NLAT*NLON + j*NLAT + i]*1000.0;	
 
	} else if (varName=="dpipdx") {
	value_out= dpiprimedx[k*NLAT*NLON + j*NLAT + i];	
	} else if (varName=="dpipdy") {
	value_out= dpiprimedy[k*NLAT*NLON + j*NLAT + i];	  
	} else if (varName=="dpipdz") {
	value_out= dpiprimedz[k*NLAT*NLON + j*NLAT + i];
	} else if (varName=="trp") {
	value_out= thetarhoprime[k*NLAT*NLON + j*NLAT + i];	 
   
	} else {
	std::cout << "Requested Variable unknown. \n";
  std::cout << varName.toStdString() << "\n";
	return false;
	}
	
	return value_out;	
}


double NetCDF_XYZ::calc_A(const int &i,const int &j,const int &k)
{
	double thetarhobar = this->getValue(i,j,k,(QString)"trb");
	double u = this->getValue(i,j,k,(QString)"u");
	double dudx = this->getValue(i,j,k,(QString)"dudx");
	double v = this->getValue(i,j,k,(QString)"v");	
	double dudy = this->getValue(i,j,k,(QString)"dudy");
	double w = this->getValue(i,j,k,(QString)"w");	
	double dudz = this->getValue(i,j,k,(QString)"dudz");
	double azimuth = this->getValue(i,j,k,(QString)"az");	
	double vtbar = this->getValue(i,j,k,(QString)"vtb");
	double radius = this->getValue(i,j,k,(QString)"r");	
	double vprime = this->getValue(i,j,k,(QString)"vp");  
  //double dpibdx = this->getValue(i,j,k,(QString)"dpibdx");
  //double trp = this->getValue(i,j,k,(QString)"trp");
  double dpipdx = this->getValue(i,j,k,(QString)"dpipdx");
  float c_p = 1005.7;

  if (thetarhobar==-999 or u==-999 or dudx==-999 or v==-999 or dudy==-999 or w==-999 or dudz==-999 or azimuth==-999 or vtbar==-999 or radius==-999 or vprime==-999){
    return -999;}
  
  //double a = (u*dudx+v*dudy+w*dudz+vtbar*vtbar/radius*cos(azimuth)-f*vprime);
  //double a = (u*dudx+v*dudy+w*dudz+vtbar*vtbar/radius*cos(azimuth)-f*vprime)+c_p*dpibdx*trp;		//If pip only
  double a = (u*dudx+v*dudy+w*dudz+vtbar*vtbar/radius*cos(azimuth)-f*vprime)+c_p*thetarhobar*dpipdx;   //If trp only  
  
	return a;	
}


double NetCDF_XYZ::calc_B(const int &i,const int &j,const int &k)
{
	double thetarhobar = this->getValue(i,j,k,(QString)"trb");
	double u = this->getValue(i,j,k,(QString)"u");
	double dvdx = this->getValue(i,j,k,(QString)"dvdx");
	double v = this->getValue(i,j,k,(QString)"v");	
	double dvdy = this->getValue(i,j,k,(QString)"dvdy");
	double w = this->getValue(i,j,k,(QString)"w");	
	double dvdz = this->getValue(i,j,k,(QString)"dvdz");
	double azimuth = this->getValue(i,j,k,(QString)"az");	
	double vtbar = this->getValue(i,j,k,(QString)"vtb");
	double radius = this->getValue(i,j,k,(QString)"r");	
	double uprime = this->getValue(i,j,k,(QString)"up"); 
  //double dpibdy = this->getValue(i,j,k,(QString)"dpibdy");
	//double trp = this->getValue(i,j,k,(QString)"trp");
  double dpipdy = this->getValue(i,j,k,(QString)"dpipdy");
  float c_p = 1005.7;

  if (thetarhobar==-999 or u==-999 or dvdx==-999 or v==-999 or dvdy==-999 or w==-999 or dvdz==-999 or azimuth==-999 or vtbar==-999 or radius==-999 or uprime==-999){
    return -999;}
    
  //double b = (u*dvdx+v*dvdy+w*dvdz+vtbar*vtbar/radius*sin(azimuth)+f*uprime);
  //double b = (u*dvdx+v*dvdy+w*dvdz+vtbar*vtbar/radius*sin(azimuth)+f*uprime)+c_p*dpibdy*trp;   //If pip only
  double b = (u*dvdx+v*dvdy+w*dvdz+vtbar*vtbar/radius*sin(azimuth)+f*uprime)+c_p*thetarhobar*dpipdy;   //If trp only
  return b;	
}

double NetCDF_XYZ::calc_C(const int &i,const int &j,const int &k)
{
	double thetarhobar = this->getValue(i,j,k,(QString)"trb");
	double u = this->getValue(i,j,k,(QString)"u");
	double dwdx = this->getValue(i,j,k,(QString)"dwdx");
	double v = this->getValue(i,j,k,(QString)"v");	
	double dwdy = this->getValue(i,j,k,(QString)"dwdy");
	double w = this->getValue(i,j,k,(QString)"w");	
	double dwdz = this->getValue(i,j,k,(QString)"dwdz");
  double dpipdz = this->getValue(i,j,k,(QString)"dpipdz");
	//double trp = this->getValue(i,j,k,(QString)"trp");
  float g = 9.81*1000.0;
  
  if (thetarhobar==-999 or u==-999 or dwdx==-999 or v==-999 or dwdy==-999 or w==-999 or dwdz==-999){
    return -999;}

  //double c = (u*dwdx+v*dwdy+w*dwdz);
  //double c = (u*dwdx+v*dwdy+w*dwdz)-g/thetarhobar*trp;    //If pip only
  double c = (u*dwdx+v*dwdy+w*dwdz)+c_p*thetarhobar*dpipdz;  //If trp only

  return c;	
}

double NetCDF_XYZ::calc_D(const int &i,const int &j,const int &k)
{
  double u = this->getValue(i,j,k,(QString)"u");
  double dtrbdx = this->getDerivative(i,j,k,(QString)"trb",1);
  double v = this->getValue(i,j,k,(QString)"v");
  double dtrbdy = this->getDerivative(i,j,k,(QString)"trb",2);
  double w = this->getValue(i,j,k,(QString)"w");
  double dtrbdz = this->getDerivative(i,j,k,(QString)"trb",3);

 if (u==-999 or dtrbdx==-999 or v==-999 or dtrbdy==-999 or w==-999 or dtrbdz==-999){
    return -999;}

  double d = u*dtrbdx+v*dtrbdy+w*dtrbdz;

	return d;	
}


double NetCDF_XYZ::getDerivative(const int &i,const int &j,const int &k, const QString &var, const int &der)
{
	// input variable "der" specifies direction of derivation: 1=dx, 2=dy, 3=dz
	double derivative;
	QString derDir;
  // Geographic functions
  GeographicLib::TransverseMercatorExact tm = GeographicLib::TransverseMercatorExact::UTM;
  double referenceLon = -90.0;  //arbitrary
  double x1,x2,y1,y2;      

	switch ( der ) {
	  case 1:
	    derDir = "X";  
	    if (i==0){
        //need to convert lat/lon differences to kms
        tm.Forward(referenceLon,this->getValue(i+1,j,k,"lat"),this->getValue(i+1,j,k,"lon"),x1,y1);
        tm.Forward(referenceLon,this->getValue(i,j,k,"lat"),this->getValue(i,j,k,"lon"),x2,y2);
        double distance = sqrt((x1-x2)*(x1-x2) + (y1-y2)*(y1-y2));
	      derivative = (this->getValue(i+1,j,k,var)-this->getValue(i,j,k,var))/distance;
	    } else if (i== NLON-1) {
        tm.Forward(referenceLon,this->getValue(i,j,k,"lat"),this->getValue(i,j,k,"lon"),x1,y1);
        tm.Forward(referenceLon,this->getValue(i-1,j,k,"lat"),this->getValue(i-1,j,k,"lon"),x2,y2);
        double distance = sqrt((x1-x2)*(x1-x2) + (y1-y2)*(y1-y2));
	      derivative = (this->getValue(i,j,k,var)-this->getValue(i-1,j,k,var))/distance;	    
	    } else {
        tm.Forward(referenceLon,this->getValue(i+1,j,k,"lat"),this->getValue(i+1,j,k,"lon"),x1,y1);
        tm.Forward(referenceLon,this->getValue(i-1,j,k,"lat"),this->getValue(i-1,j,k,"lon"),x2,y2);
        double distance = sqrt((x1-x2)*(x1-x2) + (y1-y2)*(y1-y2));
	      derivative = (this->getValue(i+1,j,k,var)-this->getValue(i-1,j,k,var))/distance;	    
	    }
	    break;
	  case 2:
	    derDir = "Y";
	    if (j==0){
        tm.Forward(referenceLon,this->getValue(i,j+1,k,"lat"),this->getValue(i,j+1,k,"lon"),x1,y1);
        tm.Forward(referenceLon,this->getValue(i,j,k,"lat"),this->getValue(i,j,k,"lon"),x2,y2);
        double distance = sqrt((x1-x2)*(x1-x2) + (y1-y2)*(y1-y2));       
	      derivative = (this->getValue(i,j+1,k,var)-this->getValue(i,j,k,var))/distance;
	    } else if (j== NLAT-1) {
        tm.Forward(referenceLon,this->getValue(i,j,k,"lat"),this->getValue(i,j,k,"lon"),x1,y1);
        tm.Forward(referenceLon,this->getValue(i,j-1,k,"lat"),this->getValue(i,j-1,k,"lon"),x2,y2);
        double distance = sqrt((x1-x2)*(x1-x2) + (y1-y2)*(y1-y2)); 
	      derivative = (this->getValue(i,j,k,var)-this->getValue(i,j-1,k,var))/distance;	    
	    } else {
        tm.Forward(referenceLon,this->getValue(i,j+1,k,"lat"),this->getValue(i,j+1,k,"lon"),x2,y2);
        tm.Forward(referenceLon,this->getValue(i,j-1,k,"lat"),this->getValue(i,j-1,k,"lon"),x2,y2);
        double distance = sqrt((x1-x2)*(x1-x2) + (y1-y2)*(y1-y2)); 
	      derivative = (this->getValue(i,j+1,k,var)-this->getValue(i,j-1,k,var))/distance;	    
	    }
	    break;
	  case 3:
	  	derDir = "Z";
	    if (k==0){
	      derivative = (this->getValue(i,j,k+1,var)-this->getValue(i,j,k,var))/(this->getValue(i,j,k+1,"z")-this->getValue(i,j,k,"z"));
	    } else if (k== NALT-1) {
	      derivative = (this->getValue(i,j,k,var)-this->getValue(i,j,k-1,var))/(this->getValue(i,j,k,"z")-this->getValue(i,j,k-1,"z"));	    
	    } else {
	      derivative = (this->getValue(i,j,k+1,var)-this->getValue(i,j,k-1,var))/(this->getValue(i,j,k+1,"z")-this->getValue(i,j,k-1,"z"));	    
	    }	  	   
	    break; 
	  default:
		std::cout << "Unknown value for calculating derivative. Valid options are 1, 2 and 3.\n";
		exit(1);
	}
	return derivative;
}
