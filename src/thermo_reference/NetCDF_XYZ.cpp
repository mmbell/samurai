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


NetCDF_XYZ::NetCDF_XYZ(const int& metFile_idim, const int& metFile_jdim, const int& metFile_kdim)
:NetCDF::NetCDF(metFile_idim, metFile_jdim, metFile_kdim)
 {
  NLON = metFile_idim;
  NLAT = metFile_jdim;
  NALT = metFile_kdim; 
   
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
    NcVar *uVar,*vVar,*wVar,*dudxVar,*dvdxVar,*dwdxVar,*dudyVar,*dvdyVar,*dwdyVar,*dudzVar,*dvdzVar,*dwdzVar, *trbVar, *dpibdxVar, *dpibdyVar;

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
	value_out= u[k*NLAT*NLON + j*NLAT + i];
	} else if (varName=="v") {
	value_out= v[k*NLAT*NLON + j*NLAT + i];
	} else if (varName=="w") {
	value_out= w[k*NLAT*NLON + j*NLAT + i];	
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
	value_out= dpibardx[k*NLAT*NLON + j*NLAT + i]/1000.0;	
	} else if (varName=="dpibdy") {
	value_out= dpibardy[k*NLAT*NLON + j*NLAT + i]/1000.0;	 
        } else if (varName=="A") {
        value_out= this->calc_A(i,j,k); 
        } else if (varName=="B") {
       	value_out= this->calc_B(i,j,k);
       	} else if (varName=="C") {
       	value_out= this->calc_C(i,j,k);

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
        double dpibdx = this->getValue(i,j,k,(QString)"dpibdx");
        float c_p = 1005.7;

  if (thetarhobar==-999 or u==-999 or dudx*1.0E5==-999 or v==-999 or dudy*1.0E5==-999 or w==-999 or dudz*1.0E5==-999 or dpibdx*1000.0==-999){
    return -999;}
  
        double a = 1.0/(c_p*thetarhobar)* (u*dudx+v*dudy+w*dudz-f*v)+dpibdx;   
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
        double dpibdy = this->getValue(i,j,k,(QString)"dpibdy");
        float c_p = 1005.7;

  if (thetarhobar==-999 or u==-999 or dvdx*1.0E5==-999 or v==-999 or dvdy*1.0E5==-999 or w==-999 or dvdz*1.0E5==-999 or dpibdy*1000.0==-999){
    return -999;}
    
    double b = 1.0/(c_p*thetarhobar)*(u*dvdx+v*dvdy+w*dvdz+f*u)+dpibdy;   // trp neglected
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
        float c_p = 1005.7;
        float g = 9.81;
  
  if (thetarhobar==-999 or u==-999 or dwdx*1.0E5==-999 or v==-999 or dwdy*1.0E5==-999 or w==-999 or dwdz*1.0E5==-999){
    return -999;}

    double c = 1.0/(c_p*thetarhobar)*(u*dwdx+v*dwdy+w*dwdz);
    return c;	
}

double NetCDF_XYZ::calc_D(const int &i,const int &j,const int &k)
{

  double thetarhobar = this->getValue(i,j,k,(QString)"trb");
  double dAdz = this->getDerivative(i,j,k,(QString)"A",3)*1000.0;   // per km
  double dCdx = this->getDerivative(i,j,k,(QString)"C",1)*1000.0;   // per km
  float c_p = 1005.7;
  float g = 9.81;  

 if (thetarhobar==-999 or dAdz==-999000 or dCdx==-999000){
    return -999;}

  double d = (dAdz-dCdx)*thetarhobar*thetarhobar*(-c_p/g);

	return d;	
}

double NetCDF_XYZ::calc_E(const int &i,const int &j,const int &k)
{
  double thetarhobar = this->getValue(i,j,k,(QString)"trb");
  double dBdz = this->getDerivative(i,j,k,(QString)"B",3)*1000.0;  // per km
  double dCdy = this->getDerivative(i,j,k,(QString)"C",2)*1000.0;  // per km
  float	c_p = 1005.7;
  float	g = 9.81;
  
 if (thetarhobar==-999 or dBdz==-999000 or dCdy==-999000){
    return -999;}

  double e = (dBdz-dCdy)*thetarhobar*thetarhobar*(-c_p/g);

        return e;
}


double NetCDF_XYZ::getDerivative(const int &i,const int &j,const int &k, const QString &var, const int &der)
{
	// input variable "der" specifies direction of derivation: 1=dx, 2=dy, 3=dz
	double derivative;
	QString derDir;
  // Geographic functions
  GeographicLib::TransverseMercatorExact tm = GeographicLib::TransverseMercatorExact::UTM();
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
              if (this->getValue(i+1,j,k,var)==-999 || this->getValue(i,j,k,var)==-999){
                 derivative = -999;
              } else {
                 derivative = (this->getValue(i+1,j,k,var)-this->getValue(i,j,k,var))/distance;
              }
            } else if (i== NLON-1) {
        tm.Forward(referenceLon,this->getValue(i,j,k,"lat"),this->getValue(i,j,k,"lon"),x1,y1);
        tm.Forward(referenceLon,this->getValue(i-1,j,k,"lat"),this->getValue(i-1,j,k,"lon"),x2,y2);
        double distance = sqrt((x1-x2)*(x1-x2) + (y1-y2)*(y1-y2));
              if (this->getValue(i,j,k,var)==-999 || this->getValue(i-1,j,k,var)==-999){
                 derivative = -999;
              } else {
              derivative = (this->getValue(i,j,k,var)-this->getValue(i-1,j,k,var))/distance;
              }
            } else {
        tm.Forward(referenceLon,this->getValue(i+1,j,k,"lat"),this->getValue(i+1,j,k,"lon"),x1,y1);
        tm.Forward(referenceLon,this->getValue(i-1,j,k,"lat"),this->getValue(i-1,j,k,"lon"),x2,y2);
        double distance = sqrt((x1-x2)*(x1-x2) + (y1-y2)*(y1-y2));
              if (this->getValue(i+1,j,k,var)==-999 || this->getValue(i-1,j,k,var)==-999){
                 derivative = -999;
              } else {
                 derivative = (this->getValue(i+1,j,k,var)-this->getValue(i-1,j,k,var))/distance;
              }
            }
            break;
          case 2:
             derDir = "Y";
            if (j==0){
        tm.Forward(referenceLon,this->getValue(i,j+1,k,"lat"),this->getValue(i,j+1,k,"lon"),x1,y1);
        tm.Forward(referenceLon,this->getValue(i,j,k,"lat"),this->getValue(i,j,k,"lon"),x2,y2);
        double distance = sqrt((x1-x2)*(x1-x2) + (y1-y2)*(y1-y2));
              if (this->getValue(i,j+1,k,var)==-999 || this->getValue(i,j,k,var)==-999){
                 derivative = -999;
              } else {
                 derivative = (this->getValue(i,j+1,k,var)-this->getValue(i,j,k,var))/distance;
              }
            } else if (j== NLAT-1) {
        tm.Forward(referenceLon,this->getValue(i,j,k,"lat"),this->getValue(i,j,k,"lon"),x1,y1);
        tm.Forward(referenceLon,this->getValue(i,j-1,k,"lat"),this->getValue(i,j-1,k,"lon"),x2,y2);
        double distance = sqrt((x1-x2)*(x1-x2) + (y1-y2)*(y1-y2));
              if (this->getValue(i,j,k,var)==-999 || this->getValue(i,j-1,k,var)==-999){
                 derivative = -999;
              } else {
              derivative = (this->getValue(i,j,k,var)-this->getValue(i,j-1,k,var))/distance;
              }
            } else {
        tm.Forward(referenceLon,this->getValue(i,j+1,k,"lat"),this->getValue(i,j+1,k,"lon"),x2,y2);
        tm.Forward(referenceLon,this->getValue(i,j-1,k,"lat"),this->getValue(i,j-1,k,"lon"),x2,y2);
        double distance = sqrt((x1-x2)*(x1-x2) + (y1-y2)*(y1-y2));
              if (this->getValue(i,j+1,k,var)==-999 || this->getValue(i,j-1,k,var)==-999){
                 derivative = -999;
              } else {
              derivative = (this->getValue(i,j+1,k,var)-this->getValue(i,j-1,k,var))/distance;
              }
            }
            break;
          case 3:
                derDir = "Z";
            if (k==0){
              if (this->getValue(i,j,k+1,var)==-999 || this->getValue(i,j,k,var)==-999){
                 derivative = -999;
              } else {
                derivative = (this->getValue(i,j,k+1,var)-this->getValue(i,j,k,var))/(this->getValue(i,j,k+1,"z")-this->getValue(i,j,k,"z"));
              }
            } else if (k== NALT-1) {
              if (this->getValue(i+1,j,k,var)==-999 || this->getValue(i,j,k,var)==-999){
                 derivative = -999;
              } else {
                 derivative = (this->getValue(i,j,k,var)-this->getValue(i,j,k-1,var))/(this->getValue(i,j,k,"z")-this->getValue(i,j,k-1,"z"));
              }   
            } else {
              if (this->getValue(i,j,k+1,var)==-999 || this->getValue(i,j,k-1,var)==-999){
                 derivative = -999;
              } else {
                 derivative = (this->getValue(i,j,k+1,var)-this->getValue(i,j,k-1,var))/(this->getValue(i,j,k+1,"z")-this->getValue(i,j,k-1,"z"));
              }
            }
            break;
          default:
                std::cout << "Unknown value for calculating derivative. Valid options are 1, 2 and 3.\n";
                exit(1);
        }
	return derivative;
}
