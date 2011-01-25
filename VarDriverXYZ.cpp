/*
 *  VarDriverXYZ.cpp
 *  samurai
 *
 *  Copyright 2008 Michael Bell. All rights reserved.
 *
 */

#include "VarDriverXYZ.h"
#include "Dorade.h"
#include <iterator>
#include <fstream>
#include <cmath>
#include <QTextStream>
#include <QFile>
#include <QVector>
#include <iomanip>
#include "RecursiveFilter.h"

VarDriverXYZ::VarDriverXYZ()
	: VarDriver()
{
	numVars = 7;
	maxIter = 1.;
}

VarDriverXYZ::~VarDriverXYZ()
{
	for (unsigned int yi = 0; yi < maxJdim; yi++) {
		delete[] BG[yi];
		delete[] BGsave[yi];
	}
	delete[] BG;
	delete[] BGsave;
	delete[] obs;
	delete[] bgObs;
	delete obCostXYZ;
}


void VarDriverXYZ::preProcessMetObs()
{
	
	vector<real> rhoP;
	
	// Cressman for the reflectivity
	double ROI = 1.25*(configHash.value("xincr").toFloat());
	double RSquare = ROI*ROI;
	double ROIsquare2 = ROI*sqrt(2.);
	
	// Check the data directory for files
	QDir dataPath("./vardata");
	dataPath.setFilter(QDir::Files);
	dataPath.setSorting(QDir::Name);
	QStringList filenames = dataPath.entryList();
	
	int processedFiles = 0;
	QList<MetObs>* metData = new QList<MetObs>;
	cout << "Found " << filenames.size() << " data files to read..." << endl;
	for (int i = 0; i < filenames.size(); ++i) {
		metData->clear();
		QString file = filenames.at(i);
		QStringList fileparts = file.split(".");
		if (fileparts.isEmpty()) {
			cout << "Unknown file! " << file.toAscii().data() << endl;
			continue;
		}
		QString suffix = fileparts.last();
		QString prefix = fileparts.first();
		if (prefix == "swp") {
			// Switch it to suffix
			suffix = "swp";
		}
		cout << "Processing " << file.toAscii().data() << " of type " << suffix.toAscii().data() << endl;
		QFile metFile(dataPath.filePath(file));
		
		// Read different types of files
		switch (dataSuffix.value(suffix)) {
			case (frd):
				if (!read_frd(metFile, metData))
					cout << "Error reading frd file" << endl;
				break;
			case (cls):
				if (!read_cls(metFile, metData))
					cout << "Error reading cls file" << endl;
				break;
			case (sec):
				if (!read_sec(metFile, metData))
					cout << "Error reading sec file" << endl;
				break;
			case (ten):
				if (!read_ten(metFile, metData))
					cout << "Error reading ten file" << endl;
				break;
			case (swp):
				if (!read_dorade(metFile, metData))
					cout << "Error reading swp file" << endl;
				break;
			case (sfmr):
				if (!read_sfmr(metFile, metData))
					cout << "Error reading sfmr file" << endl;
				break;
			case (wwind):
				if (!read_wwind(metFile, metData))
					cout << "Error reading wwind file" << endl;
				break;
			case (eol):
				if (!read_eol(metFile, metData))
					cout << "Error reading eol file" << endl;
				break;
			case (qscat):
				if (!read_qscat(metFile, metData))
					cout << "Error reading wwind file" << endl;
				break;
			case (ascat):
				if (!read_ascat(metFile, metData))
					cout << "Error reading wwind file" << endl;
				break;
			case (nopp):
				if (!read_nopp(metFile, metData))
					cout << "Error reading wwind file" << endl;
				break;
			case (cimss):
				if (!read_cimss(metFile, metData))
					cout << "Error reading cimss file" << endl;
				break;
			case (dwl):
				if (!read_dwl(metFile, metData))
					cout << "Error reading dwl file" << endl;
				break;
			case (cen):
				continue;				
			default:
				cout << "Unknown data type, skipping..." << endl;
				continue;
		}
		
		processedFiles++;
		
		// Process the metObs into Observations
		QDateTime startTime = tcVector.front().getTime();
		QDateTime endTime = tcVector.back().getTime();
		for (int i = 0; i < metData->size(); ++i) {
			
			// Make sure the ob is within the time limits
			MetObs metOb = metData->at(i);
			QDateTime obTime = metOb.getTime();
			QString obstring = obTime.toString(Qt::ISODate);
			QString tcstart = startTime.toString(Qt::ISODate);
			QString tcend = endTime.toString(Qt::ISODate);		
			if ((obTime < startTime) or (obTime > endTime)) continue;
			int tci = startTime.secsTo(obTime);
			if ((tci < 0) or (tci > (int)tcVector.size())) {
				cout << "Time problem with observation " << tci << endl;
				continue;
			}
			
			// Our generic observation
			Observation varOb;
			
			// Get the X, Y & Z
			double latrad = tcVector[tci].getLat() * Pi/180.0;
			double fac_lat = 111.13209 - 0.56605 * cos(2.0 * latrad)
			+ 0.00012 * cos(4.0 * latrad) - 0.000002 * cos(6.0 * latrad);
			double fac_lon = 111.41513 * cos(latrad)
			- 0.09455 * cos(3.0 * latrad) + 0.00012 * cos(5.0 * latrad);
			double obY = (metOb.getLat() - tcVector[tci].getLat())*fac_lat;
			double obX = (metOb.getLon() - tcVector[tci].getLon())*fac_lon;
			double heightm = metOb.getAltitude();
			double obZ = heightm/1000.;
			// Make sure the ob is in the domain
			if ((obX < imin) or (obX > imax) or
				(obY < jmin) or (obY > jmax) or
				(obZ < kmin) or (obZ > kmax))
				continue;
			
			varOb.setCartesianX(obX);
			varOb.setCartesianY(obY);
			varOb.setAltitude(obZ);
			varOb.setTime(obTime.toTime_t());
			
			real Um = tcVector[tci].getUmean();
			real Vm = tcVector[tci].getVmean();

			// Reference states			
			real rhoBar = getReferenceVariable(rhoaref, heightm);
			real qBar = getReferenceVariable(qvbhypref, heightm);
			real tBar = getReferenceVariable(tempref, heightm);
			
			// Initialize the weights
			varOb.setWeight(0., 0);
			varOb.setWeight(0., 1);
			varOb.setWeight(0., 2);
			varOb.setWeight(0., 3);
			varOb.setWeight(0., 4);
			varOb.setWeight(0., 5);
			varOb.setWeight(0., 6);
			double u, v, w, rho, rhoa, qv, tempk, rhov, rhou, rhow, wspd; 
			switch (metOb.getObType()) {
				case (MetObs::dropsonde):
					varOb.setType(MetObs::dropsonde);
					u = metOb.getCartesianUwind();
					v = metOb.getCartesianVwind();
					w = metOb.getVerticalVelocity();
					rho = metOb.getMoistDensity();
					rhoa = metOb.getAirDensity();
					qv = metOb.getQv();
					tempk = metOb.getTemperature();
					
					// Separate obs for each measurement
					// rho v 1 m/s error
					if ((u != -999) and (rho != -999)) {
						// rho u 1 m/s error
						varOb.setWeight(1., 0);
						rhou = rho*(u - Um);
						//cout << "RhoU: " << rhou << endl;
						varOb.setOb(rhou);
						varOb.setError(1.0);
						obVector.push_back(varOb);
						varOb.setWeight(0., 0);
						
						varOb.setWeight(1., 1);
						rhov = rho*(v - Vm);
						varOb.setOb(rhov);
						varOb.setError(1.0);
						obVector.push_back(varOb);
						varOb.setWeight(0., 1);
						
					}
					if ((w != -999) and (rho != -999)) {
						// rho w 1.5 m/s error
						varOb.setWeight(1., 2);
						rhow = rho*w;
						varOb.setOb(rhow);
						varOb.setError(1.5);
						obVector.push_back(varOb);
						varOb.setWeight(0., 2);
					}
					if (tempk != -999) {
						// temperature 1 K error
						varOb.setWeight(1., 3);
						varOb.setOb(tempk - tBar);
						varOb.setError(1.0);
						obVector.push_back(varOb);
						varOb.setWeight(0., 3);
					}
					if (qv != -999) {
						// Qv 0.5 g/kg error
						varOb.setWeight(1., 4);
						qv = bhypTransform(qv);
						varOb.setOb(qv-qBar);
						varOb.setError(0.5);
						obVector.push_back(varOb);
						varOb.setWeight(0., 4);
					}
					if (rhoa != -999) {
						// Rho prime .1 kg/m^3 error
						varOb.setWeight(1., 5);
						varOb.setOb((rhoa-rhoBar)*100);
						varOb.setError(1.0);
						obVector.push_back(varOb);
						varOb.setWeight(0., 5);
					}
					
					break;
					
				case (MetObs::flightlevel):
					varOb.setType(MetObs::flightlevel);
					u = metOb.getCartesianUwind();
					v = metOb.getCartesianVwind();
					w = metOb.getVerticalVelocity();
					rho = metOb.getMoistDensity();
					rhoa = metOb.getAirDensity();
					qv = metOb.getQv();
					tempk = metOb.getTemperature();
					
					// Separate obs for each measurement
					// rho v 1 m/s error
					if ((u != -999) and (rho != -999)) {						
						// rho u 1 m/s error
						varOb.setWeight(1., 0);
						rhou = rho*(u - Um);
						varOb.setOb(rhou);
						varOb.setError(1.0);
						obVector.push_back(varOb);
						varOb.setWeight(0., 0);
						
						varOb.setWeight(1., 1);
						rhov = rho*(v - Vm);
						varOb.setOb(rhov);
						varOb.setError(1.0);
						obVector.push_back(varOb);
						varOb.setWeight(0., 1);
					}
					if ((w != -999) and (rho != -999)) {
						// rho w 1 dm/s error
						varOb.setWeight(1., 2);
						rhow = rho*w;
						varOb.setOb(rhow);
						varOb.setError(0.25);
						obVector.push_back(varOb);
						varOb.setWeight(0., 2);
					}
					if (tempk != -999) {
						// temperature 1 K error
						varOb.setWeight(1., 3);
						varOb.setOb(tempk - tBar);
						varOb.setError(1.0);
						obVector.push_back(varOb);
						varOb.setWeight(0., 3);
					}
					if (qv != -999) {
						// Qv 0.5 g/kg error
						varOb.setWeight(1., 4);
						qv = bhypTransform(qv);
						varOb.setOb(qv-qBar);
						varOb.setError(0.5);
						obVector.push_back(varOb);
						varOb.setWeight(0., 4);
					}
					if (rhoa != -999) {
						// Rho prime .1 kg/m^3 error
						varOb.setWeight(1., 5);
						varOb.setOb((rhoa-rhoBar)*100);
						varOb.setError(1.0);
						obVector.push_back(varOb);
						varOb.setWeight(0., 5);
					}
					
					break;

				case (MetObs::sfmr):
					varOb.setType(MetObs::sfmr);
					wspd = metOb.getWindSpeed();
					// This needs to be redone for the Cartesian case
					//vBG = 1.e3*bilinearField(obX, obY, 0);
					//uBG = -1.e5*bilinearField(obX, 20., 1)/(rad*20.);
					varOb.setWeight(1., 0);
					//varOb.setWeight(1., 1);
					varOb.setOb(wspd);
					varOb.setError(10.0);
					obVector.push_back(varOb);
					break;
				
				case (MetObs::qscat):
					varOb.setType(MetObs::qscat);
					u = metOb.getCartesianUwind();
					v = metOb.getCartesianVwind();
					if (u != -999) {
						// rho u 1 m/s error
						varOb.setWeight(1., 0);
						rhou = (u - Um);
						//cout << "RhoU: " << rhou << endl;
						varOb.setOb(rhou);
						varOb.setError(2.5);
						obVector.push_back(varOb);
						varOb.setWeight(0., 0);
						
						varOb.setWeight(1., 1);
						// Multiply by rho later from grid values
						rhov = (v - Vm);
						varOb.setOb(rhov);
						varOb.setError(2.5);
						obVector.push_back(varOb);
						varOb.setWeight(0., 1);						
					}
					break;
					
				case (MetObs::ascat):
					varOb.setType(MetObs::ascat);
					u = metOb.getCartesianUwind();
					v = metOb.getCartesianVwind();
					if (u != -999) {						
						// rho u 1 m/s error
						varOb.setWeight(1., 0);
						rhou = (u - Um);
						//cout << "RhoU: " << rhou << endl;
						varOb.setOb(rhou);
						varOb.setError(2.5);
						obVector.push_back(varOb);
						varOb.setWeight(0., 0);
						
						varOb.setWeight(1., 1);
						// Multiply by rho later from grid values
						rhov = (v - Vm);
						varOb.setOb(rhov);
						varOb.setError(2.5);
						obVector.push_back(varOb);
						varOb.setWeight(0., 1);						
					}
					break;
					
				case (MetObs::AMV):
					varOb.setType(MetObs::AMV);
					u = metOb.getCartesianUwind();
					v = metOb.getCartesianVwind();
					if (u != -999) {
						// rho u 10 m/s error
						varOb.setWeight(1., 0);
						rhou = (u - Um);
						//cout << "RhoU: " << rhou << endl;
						varOb.setOb(rhou);
						varOb.setError(10.);
						obVector.push_back(varOb);
						varOb.setWeight(0., 0);
						
						varOb.setWeight(1., 1);
						// Multiply by rho later from grid values
						rhov = (v - Vm);
						varOb.setOb(rhov);
						varOb.setError(10.);
						obVector.push_back(varOb);
						varOb.setWeight(0., 1);						
					}
					break;
					
				case (MetObs::lidar):
				{
					varOb.setType(MetObs::lidar);
					// Geometry terms
					double az = metOb.getAzimuth()*Pi/180.;
					double el = metOb.getElevation()*Pi/180.;
					double uWgt = sin(az)*cos(el);
					double vWgt = cos(az)*cos(el);
					double wWgt = sin(el);
					
					// Fall speed is assumed zero for now since we are dealing with aerosols, perhaps this could be done better
					// Reflectivity observations only used for QC -- need to move this into a pre-processing step
					double db = metOb.getReflectivity();
					double vr = metOb.getRadialVelocity();
					//if ((db < 22.0) or (fabs(vr) > 10.0)) continue;
					if (fabs(vr) > 10.0) continue;
					double w_term = 0.0;  
					double Vdopp = vr - w_term*sin(el) - Um*sin(az)*cos(el) - Vm*cos(az)*cos(el);
					
					varOb.setWeight(uWgt, 0);
					varOb.setWeight(vWgt, 1);
					varOb.setWeight(wWgt, 2);
					
					// Theoretically, rhoPrime could be included as a prognostic variable here...
					// However, adding another unknown without an extra equation makes the problem even more underdetermined
					// so assume it is small and ignore it
					// double rhopWgt = -Vdopp;
					//varOb.setWeight(rhopWgt, 5);
					
					// Set the error according to the spectrum width and potential fall speed error (assume 1 m/s here)
					double DopplerError = metOb.getSpectrumWidth() + fabs(wWgt)*1.;
					if (DopplerError < 1.0) DopplerError = 1.0;
					varOb.setError(DopplerError);
					varOb.setOb(Vdopp);
					obVector.push_back(varOb);
					varOb.setWeight(0., 0);	
					varOb.setWeight(0., 1);	
					varOb.setWeight(0., 2);
					
					break;
				}	
				case (MetObs::radar):
				{
					varOb.setType(MetObs::radar);
					// Geometry terms
					double az = metOb.getAzimuth()*Pi/180.;
					double el = metOb.getElevation()*Pi/180.;
					double uWgt = sin(az)*cos(el);
					double vWgt = cos(az)*cos(el);
					double wWgt = sin(el);
					
					// Fall speed
					double Z = metOb.getReflectivity();
					double H = metOb.getAltitude();
					double ZZ=pow(10.0,(Z*0.1));
					double zeroC = 4800.;
					double hlow= zeroC; 
					double hhi= hlow + 1000;
					
					/* density correction term (rhoo/rho)*0.45 [rho(Z)=rho_o exp-(z/H), where 
					 H is the scale height = 9.58125 from Gray's inner 2 deg composite] 
					 0.45 density correction from Beard (1985, JOAT pp 468-471) 
					 Adjusted to use Jordan hydrostatic scale height -MB */
					//double DCOR=exp(0.45*metOb.getAltitude()*0.0001068);
					real rho = getReferenceVariable(rhoref, H);
					real rhosfc = getReferenceVariable(rhoref, 0.);
					real DCOR = pow((rhosfc/rho),(double)0.45);
					
					// The snow relationship (Atlas et al., 1973) --- VT=0.817*Z**0.063  (m/s) 
					double VTS=-DCOR * (0.817*pow(ZZ,(double)0.063));
					
					// The rain relationship (Joss and Waldvogel,1971) --- VT=2.6*Z**.107 (m/s) */
					double VTR=-DCOR * (2.6*pow(ZZ,(double).107));
					/* Test if height is in the transition region between SNOW and RAIN
					   defined as hlow in km < H < hhi in km
					   if in the transition region do a linear weight of VTR and VTS */
					if ((Z > 20) and (Z <= 30)) {
						double WEIGHTR=(Z-20)/(10);
						double WEIGHTS=1.-WEIGHTR;
						VTS=(VTR*WEIGHTR+VTS*WEIGHTS)/(WEIGHTR+WEIGHTS);
					} else if (Z > 30) {
						VTS=VTR;
					}
					double w_term=VTR*(hhi-H)/1000 + VTS*(H-hlow)/1000;  
					if (H < hlow) w_term=VTR; 
					if (H > hhi) w_term=VTS;
					double Vdopp = metOb.getRadialVelocity() - w_term*sin(el) - Um*sin(az)*cos(el) - Vm*cos(az)*cos(el);
					
					varOb.setWeight(uWgt, 0);
					varOb.setWeight(vWgt, 1);
					varOb.setWeight(wWgt, 2);
					
					// Theoretically, rhoPrime could be included as a prognostic variable here...
					// However, adding another unknown without an extra equation makes the problem even more underdetermined
					// so assume it is small and ignore it
					// double rhopWgt = -Vdopp;
					//varOb.setWeight(rhopWgt, 5);
					
					// Set the error according to the spectrum width and potential fall speed error (assume 2 m/s?)
					double DopplerError = metOb.getSpectrumWidth() + fabs(wWgt)*2.;
					if (DopplerError < 1.0) DopplerError = 1.0;
					varOb.setError(DopplerError);
					varOb.setOb(Vdopp);
					obVector.push_back(varOb);
					varOb.setWeight(0., 0);	
					varOb.setWeight(0., 1);	
					varOb.setWeight(0., 2);
					
					// Reflectivity observations
					QString gridref = configHash.value("gridreflectivity");
					double qr = 0.;
					if (gridref == "qr") {
						// Do the gridding as part of the variational synthesis using Z-M relationships
						// Z-M relationships from Gamache et al (1993) JAS
						double rainmass = pow(ZZ/14630.,(double)0.6905);
						double icemass = pow(ZZ/670.,(double)0.5587);
						if ((Z > 20) and (Z <= 30)) {
							double WEIGHTR=(Z-20)/(10);
							double WEIGHTS=1.-WEIGHTR;
							icemass=(rainmass*WEIGHTR+icemass*WEIGHTS)/(WEIGHTR+WEIGHTS);
						} else if (Z > 30) {
							icemass=rainmass;
						}
						
						double precipmass = rainmass*(hhi-H)/1000 + icemass*(H-hlow)/1000;
						if (H < hlow) precipmass = rainmass;
						if (H > hhi) precipmass = icemass;
						qr = bhypTransform(precipmass/rhoBar);
						
						/* Include an observation of this quantity in the variational synthesis
						varOb.setOb(qr);
						varOb.setWeight(1., 6);
						varOb.setError(1.0);
						obVector.push_back(varOb); */
						
					} else if (gridref == "dbz") {
						//qr =  bhypTransform(ZZ);
						qr = bhypTransform(Z + 35.);
						/* Include an observation of this quantity in the variational synthesis
						 varOb.setOb(qr);
						 varOb.setWeight(1., 6);
						 varOb.setError(1.0);
						 obVector.push_back(varOb); */
						
					}

					// Do a Cressman interpolation of the maximum reflectivity/qr in a grid box
					for (int zi = 0; zi < (kdim-1); zi++) {	
						for (int zmu = -1; zmu <= 1; zmu += 2) {
							real zPos = kmin + kincr * (zi + (0.5*sqrt(1./3.) * zmu + 0.5));
							if (fabs(zPos-obZ) > ROIsquare2) continue;
							for (int xi = 0; xi < (idim-1); xi++) {
								for (int xmu = -1; xmu <= 1; xmu += 2) {
									real xPos = imin + iincr * (xi + (0.5*sqrt(1./3.) * xmu + 0.5));
									if (fabs(xPos-obX) > ROIsquare2) continue;
									
									for (int yi = 0; yi < (jdim-1); yi++) {
										for (int ymu = -1; ymu <= 1; ymu += 2) {
											real yPos = jmin + jincr * (yi + (0.5*sqrt(1./3.) * ymu + 0.5));
											if (fabs(yPos-obY) > ROIsquare2) continue;
											real rSquare = (obX-xPos)*(obX-xPos) + (obY-yPos)*(obY-yPos) + (obZ-zPos)*(obZ-zPos); 
											int bgI = xi*2 + (xmu+1)/2;
											int bgJ = yi*2 + (ymu+1)/2;
											int bgK = zi*2 + (zmu+1)/2;
											int bIndex = numVars*(idim-1)*2*(jdim-1)*2*bgK + numVars*(idim-1)*2*bgJ +numVars*bgI;
											if (rSquare < RSquare) {
												real weight = (RSquare - rSquare)/(RSquare + rSquare);
												//if (qr > bgU[bIndex +6]) bgU[bIndex +6] = qr;
												bgU[bIndex +6] += weight*qr;
												bgWeights[bIndex] += weight;
											}
										}
									}
								}
							}
						}
					}
					
					break;
				}	
							
			}

		} 
		cout << obVector.size() << " total observations." << endl;
	}
	
	delete metData;
	
	// Finish Cressman for reflectivity
	Observation varOb;
	varOb.setTime(configHash.value("reftime").toFloat());	
	varOb.setWeight(0., 0);
	varOb.setWeight(0., 1);
	varOb.setWeight(1., 2);
	varOb.setWeight(0., 3);
	varOb.setWeight(0., 4);
	varOb.setWeight(0., 5);
	varOb.setWeight(0., 6);	
	varOb.setError(1.);
	varOb.setOb(0.);
	for (int xi = 0; xi < (idim-1); xi++) {
		for (int xmu = -1; xmu <= 1; xmu += 2) {
			real xPos = imin + iincr * (xi + (0.5*sqrt(1./3.) * xmu + 0.5));
			
			for (int yi = 0; yi < (jdim-1); yi++) {
				for (int ymu = -1; ymu <= 1; ymu += 2) {
					real yPos = jmin + jincr * (yi + (0.5*sqrt(1./3.) * ymu + 0.5));
					real maxrefHeight = -1;
					for (int zi = 0; zi < (kdim-1); zi++) {	
						for (int zmu = -1; zmu <= 1; zmu += 2) {
							real zPos = kmin + kincr * (zi + (0.5*sqrt(1./3.) * zmu + 0.5));
							
							int bgI = xi*2 + (xmu+1)/2;
							int bgJ = yi*2 + (ymu+1)/2;
							int bgK = zi*2 + (zmu+1)/2;
							int bIndex = numVars*(idim-1)*2*(jdim-1)*2*bgK + numVars*(idim-1)*2*bgJ +numVars*bgI;
							if (bgWeights[bIndex] != 0) {
								bgU[bIndex +6] /= bgWeights[bIndex];
							}
							if (bgU[bIndex +6] > 0) {
								maxrefHeight = zPos;
							}
						}
					}
					if (maxrefHeight > 0) {
						varOb.setCartesianX(xPos);
						varOb.setCartesianY(yPos);
						varOb.setAltitude(maxrefHeight);					
						obVector.push_back(varOb);
					}
				}
			}
		}
	}	
	cout << obVector.size() << " total observations including pseudo W obs" << endl;
	
	// Write the Obs to a summary text file
	ofstream obstream("Observations.in");
	// Header messes up reload
	/*ostream_iterator<string> os(obstream, "\t ");
	*os++ = "Type";
	*os++ = "r";
	*os++ = "z";
	*os++ = "NULL";
	*os++ = "Observation";
	*os++ = "Inverse Error";
	*os++ = "Weight 1";
	*os++ = "Weight 2";
	*os++ = "Weight 3";
	*os++ = "Weight 4";
	*os++ = "Weight 5";
	*os++ = "Weight 6";
	obstream << endl; */

	ostream_iterator<double> od(obstream, "\t ");
	for (unsigned int i=0; i < obVector.size(); i++) {
		Observation ob = obVector.at(i);
		*od++ = ob.getOb();
		*od++ = ob.getInverseError();
		*od++ = ob.getCartesianX();
		*od++ = ob.getCartesianY();
		*od++ = ob.getAltitude();		
		*od++ = ob.getType();
		*od++ = ob.getTime();
		for (unsigned int var = 0; var < numVars; var++)
			*od++ = ob.getWeight(var);

		obstream << endl;	
	}
	
	// Load the observations into a vector
	obs = new real[obVector.size()*14];
	for (unsigned int m=0; m < obVector.size(); m++) {
		int n = m*14;
		Observation ob = obVector.at(m);
		obs[n] = ob.getOb();
		obs[n+1] = ob.getInverseError();
		obs[n+2] = ob.getCartesianX();
		obs[n+3] = ob.getCartesianY();
		obs[n+4] = ob.getAltitude();
		obs[n+5] = ob.getType();
		obs[n+6] = ob.getTime();
		for (unsigned int var = 0; var < numVars; var++) {
			obs[n+7+var] = ob.getWeight(var);
		}
	}	
	
	// All done preprocessing
	if (!processedFiles) {
		cout << "No files processed, nothing to do :(" << endl;
		// return 0;
	} else {
		cout << "Finished preprocessing " << processedFiles << " files." << endl;
	}
	
}


bool VarDriverXYZ::loadMetObs()
{
	
	// Our generic observation
	Observation varOb;
	double wgt[numVars];
	double xPos, yPos, zPos, ob, error;
	int type, time;
	ifstream obstream("./Observations.in");
	while (obstream >> ob >> error >> xPos >> yPos >> zPos >> type >> time
		   >> wgt[0] >> wgt[1] >> wgt[2] >> wgt[3] >> wgt[4] >> wgt[5] >> wgt[6])
	{
		varOb.setOb(ob);
		varOb.setCartesianX(xPos);
		varOb.setCartesianY(yPos);
		varOb.setAltitude(zPos);
		varOb.setType(type);
		varOb.setTime(time);
		varOb.setError(1./error);
		for (unsigned int var = 0; var < numVars; var++)
			varOb.setWeight(wgt[var],var);
		obVector.push_back(varOb);
	}
	
	// Load the observations into a vector
	obs = new real[obVector.size()*14];
	for (unsigned int m=0; m < obVector.size(); m++) {
		int n = m*14;
		Observation ob = obVector.at(m);
		obs[n] = ob.getOb();
		obs[n+1] = ob.getInverseError();
		obs[n+2] = ob.getCartesianX();
		obs[n+3] = ob.getCartesianY();
		obs[n+4] = ob.getAltitude();
		obs[n+5] = ob.getType();
		obs[n+6] = ob.getTime();
		for (unsigned int var = 0; var < numVars; var++) {
			obs[n+7+var] = ob.getWeight(var);
		}
		
	}	
	
	return true;
}

int VarDriverXYZ::loadBackgroundObs()
{
	//SplineD::Debug(1);
	QList<real> bgIn;
	QVector<real> logheights, uBG, vBG, wBG, tBG, qBG, rBG;
	SplineD* bgSpline;
	int time;
	double lat, lon, alt, u, v, w, t, qv, rhoa;
	double bgX, bgY, bgZ;
	float ROI = configHash.value("backgroundroi").toFloat();
	double RSquare = ROI*ROI;
	double ROIsquare2 = ROI*sqrt(2.);
	ifstream bgstream("./Background.in");
	cout << "Loading background onto Gaussian mish with " << ROI << " km radius of influence" << endl;
	
	while (bgstream >> time >> lat >> lon >> alt >> u >> v >> w >> t >> qv >> rhoa)
	{
	
		// Process the metObs into Observations
		QDateTime startTime = tcVector.front().getTime();
		QDateTime endTime = tcVector.back().getTime();
		
		// Make sure the bg is within the time limits
		QDateTime bgTime;
		bgTime.setTimeSpec(Qt::UTC);
		bgTime.setTime_t(time);
		QString obstring = bgTime.toString(Qt::ISODate);
		QString tcstart = startTime.toString(Qt::ISODate);
		QString tcend = endTime.toString(Qt::ISODate);		
		if ((bgTime < startTime) or (bgTime > endTime)) continue;
		int tci = startTime.secsTo(bgTime);
		if ((tci < 0) or (tci > (int)tcVector.size())) {
			cout << "Time problem with observation " << tci << "secs more than center entries" << endl;
			continue;
		}
		
		// Get the X, Y & Z
		double latrad = tcVector[tci].getLat() * Pi/180.0;
		double fac_lat = 111.13209 - 0.56605 * cos(2.0 * latrad)
		+ 0.00012 * cos(4.0 * latrad) - 0.000002 * cos(6.0 * latrad);
		double fac_lon = 111.41513 * cos(latrad)
		- 0.09455 * cos(3.0 * latrad) + 0.00012 * cos(5.0 * latrad);
		bgY = (lat - tcVector[tci].getLat())*fac_lat;
		bgX = (lon - tcVector[tci].getLon())*fac_lon;
		double heightm = alt;
		bgZ = heightm/1000.;
		// Make sure the ob is in the Cressman domain
		if ((bgX < (imin-ROIsquare2)) or (bgX > (imax+ROIsquare2)) or
			(bgY < (jmin-ROIsquare2)) or (bgY > (jmax+ROIsquare2))
			or (bgZ < kmin)) //Allow for higher values for interpolation purposes
			continue;
				
		real Um = tcVector[tci].getUmean();
		real Vm = tcVector[tci].getVmean();
		
		// Reference states			
		real rhoBar = getReferenceVariable(rhoaref, heightm);
		real qBar = getReferenceVariable(qvbhypref, heightm);
		real tBar = getReferenceVariable(tempref, heightm);

		real rho = rhoa + rhoa*qv/1000.;
		real rhou = rho*(u - Um);
		real rhov = rho*(v - Vm);
		real rhow = rho*w;
		real tprime = t - tBar;
		qv = bhypTransform(qv);
		real qvprime = qv-qBar;
		real rhoprime = (rhoa-rhoBar)*100;
		real logZ = log(bgZ);
		// We assume here that the background precipitation field is always zero
		real qr = 0.;
		bgIn << bgX << bgY << logZ << time << rhou << rhov << rhow << tprime << qvprime << rhoprime << qr ;
		if (logheights.size() == 0) {
			// First column
			logheights.push_back(logZ);
			uBG.push_back(rhou);
			vBG.push_back(rhov);
			wBG.push_back(rhow);
			tBG.push_back(tprime);
			qBG.push_back(qvprime);
			rBG.push_back(rhoprime);
		} else if (logZ > logheights.back()) {
			// Same column
			logheights.push_back(logZ);
			uBG.push_back(rhou);
			vBG.push_back(rhov);
			wBG.push_back(rhow);
			tBG.push_back(tprime);
			qBG.push_back(qvprime);
			rBG.push_back(rhoprime);
		} else {
			// Solve for the spline
			bgSpline = new SplineD(&logheights.front(), logheights.size(), uBG.data(), 0, SplineBase::BC_ZERO_SECOND);
			if (!bgSpline->ok())
			{
				cerr << "bgSpline setup failed." << endl;
				return -1;
			}
			/* for (int zi = 0; zi < logheights.size(); zi++) {	
				real logzPos = logheights.at(zi);
				real z = exp(logzPos);
				real s = bgSpline->evaluate(logzPos);
				real u = uBG.at(zi);
				cout << setprecision(4) << z << "\t" << logzPos << "\t" << s << "\t" << u << endl;
			} */
			// Cressman interpolation in horizontal, b-Spline interpolation on log height in vertical
			for (int zi = 0; zi < (kdim-1); zi++) {	
				for (int zmu = -1; zmu <= 1; zmu += 2) {
					real zPos = kmin + kincr * (zi + (0.5*sqrt(1./3.) * zmu + 0.5));
					real logzPos = log(zPos);
					
					for (int xi = 0; xi < (idim-1); xi++) {
						for (int xmu = -1; xmu <= 1; xmu += 2) {
							real xPos = imin + iincr * (xi + (0.5*sqrt(1./3.) * xmu + 0.5));
							if (fabs(xPos-bgX) > ROIsquare2) continue;
							
							for (int yi = 0; yi < (jdim-1); yi++) {
								for (int ymu = -1; ymu <= 1; ymu += 2) {
									real yPos = jmin + jincr * (yi + (0.5*sqrt(1./3.) * ymu + 0.5));
									if (fabs(yPos-bgY) > ROIsquare2) continue;
									
									real rSquare = (bgX-xPos)*(bgX-xPos) + (bgY-yPos)*(bgY-yPos);
									int bgI = xi*2 + (xmu+1)/2;
									int bgJ = yi*2 + (ymu+1)/2;
									int bgK = zi*2 + (zmu+1)/2;
									int bIndex = numVars*(idim-1)*2*(jdim-1)*2*bgK + numVars*(idim-1)*2*bgJ +numVars*bgI;
									if (rSquare < RSquare) {
										real weight = (RSquare - rSquare)/(RSquare + rSquare);
										if (logzPos > logheights.front()) {
											bgSpline->solve(uBG.data());
											bgU[bIndex] += weight*(bgSpline->evaluate(logzPos));
											bgSpline->solve(vBG.data());
											bgU[bIndex +1] += weight*(bgSpline->evaluate(logzPos));
											bgSpline->solve(wBG.data());
											bgU[bIndex +2] += weight*(bgSpline->evaluate(logzPos));
											bgSpline->solve(tBG.data());
											bgU[bIndex +3] += weight*(bgSpline->evaluate(logzPos));
											bgSpline->solve(qBG.data());
											bgU[bIndex +4] += weight*(bgSpline->evaluate(logzPos));
											bgSpline->solve(rBG.data());
											bgU[bIndex +5] += weight*(bgSpline->evaluate(logzPos));
											bgWeights[bIndex] += weight;
										} else {
											// Below the spline interpolation
											bgU[bIndex] += weight*uBG.front();
											bgU[bIndex +1] += weight*vBG.front();
											bgU[bIndex +2] += weight*wBG.front();
											bgU[bIndex +3] += weight*tBG.front();
											bgU[bIndex +4] += weight*qBG.front();
											bgU[bIndex +5] += weight*rBG.front();
											bgWeights[bIndex] += weight;
										}											
									}
								}
							}
						}
					}
				}
			}
			
			delete bgSpline;
			logheights.clear();
			uBG.clear();
			vBG.clear();
			wBG.clear();
			tBG.clear();
			qBG.clear();
			rBG.clear();
			
			logheights.push_back(log(bgZ));
			uBG.push_back(rhou);
			vBG.push_back(rhov);
			wBG.push_back(rhow);
			tBG.push_back(tprime);
			qBG.push_back(qvprime);
			rBG.push_back(rhoprime);			
		}
	}				

	// Solve for the last spline
	bgSpline = new SplineD(&logheights.front(), logheights.size(), &uBG[0], 2, SplineBase::BC_ZERO_SECOND);
	if (!bgSpline->ok())
	{
		cerr << "bgSpline setup failed." << endl;
		return -1;
	}
	
	// Cressman interpolation in horizontal, b-Spline interpolation on log height in vertical
	for (int zi = 0; zi < (kdim-1); zi++) {	
		for (int zmu = -1; zmu <= 1; zmu += 2) {
			real zPos = kmin + kincr * (zi + (0.5*sqrt(1./3.) * zmu + 0.5));
			real logzPos = log(zPos);
			
			for (int xi = 0; xi < (idim-1); xi++) {
				for (int xmu = -1; xmu <= 1; xmu += 2) {
					real xPos = imin + iincr * (xi + (0.5*sqrt(1./3.) * xmu + 0.5));
					if (fabs(xPos-bgX) > ROIsquare2) continue;
					
					for (int yi = 0; yi < (jdim-1); yi++) {
						for (int ymu = -1; ymu <= 1; ymu += 2) {
							real yPos = jmin + jincr * (yi + (0.5*sqrt(1./3.) * ymu + 0.5));
							if (fabs(yPos-bgY) > ROIsquare2) continue;
							
							real rSquare = (bgX-xPos)*(bgX-xPos) + (bgY-yPos)*(bgY-yPos);
							int bgI = xi*2 + (xmu+1)/2;
							int bgJ = yi*2 + (ymu+1)/2;
							int bgK = zi*2 + (zmu+1)/2;
							int bIndex = numVars*(idim-1)*2*(jdim-1)*2*bgK + numVars*(idim-1)*2*bgJ +numVars*bgI;
							if (rSquare < RSquare) {
								real weight = (RSquare - rSquare)/(RSquare + rSquare);
								if (logzPos > logheights.front()) {
									bgSpline->solve(uBG.data());
									bgU[bIndex] += weight*(bgSpline->evaluate(logzPos));
									bgSpline->solve(vBG.data());
									bgU[bIndex +1] += weight*(bgSpline->evaluate(logzPos));
									bgSpline->solve(wBG.data());
									bgU[bIndex +2] += weight*(bgSpline->evaluate(logzPos));
									bgSpline->solve(tBG.data());
									bgU[bIndex +3] += weight*(bgSpline->evaluate(logzPos));
									bgSpline->solve(qBG.data());
									bgU[bIndex +4] += weight*(bgSpline->evaluate(logzPos));
									bgSpline->solve(rBG.data());
									bgU[bIndex +5] += weight*(bgSpline->evaluate(logzPos));
									bgWeights[bIndex] += weight;
								} else {
									// Below the spline interpolation
									bgU[bIndex] += weight*uBG.front();
									bgU[bIndex +1] += weight*vBG.front();
									bgU[bIndex +2] += weight*wBG.front();
									bgU[bIndex +3] += weight*tBG.front();
									bgU[bIndex +4] += weight*qBG.front();
									bgU[bIndex +5] += weight*rBG.front();
									bgWeights[bIndex] += weight;
								}								
							}
						}
					}
				}
			}
		}
	}
	
	delete bgSpline;
	logheights.clear();
	uBG.clear();
	vBG.clear();
	wBG.clear();
	tBG.clear();
	qBG.clear();
	rBG.clear();
	
	
	int numbgObs = bgIn.size()*7/11;
	if (numbgObs > 0) {
		// Check interpolation
		for (int zi = 0; zi < (kdim-1); zi++) {	
			for (int zmu = -1; zmu <= 1; zmu += 2) {
				real zPos = kmin + kincr * (zi + (0.5*sqrt(1./3.) * zmu + 0.5));
				
				for (int xi = 0; xi < (idim-1); xi++) {
					for (int xmu = -1; xmu <= 1; xmu += 2) {
						real xPos = imin + iincr * (xi + (0.5*sqrt(1./3.) * xmu + 0.5));
						
						for (int yi = 0; yi < (jdim-1); yi++) {
							for (int ymu = -1; ymu <= 1; ymu += 2) {
								real yPos = jmin + jincr * (yi + (0.5*sqrt(1./3.) * ymu + 0.5));
								int bgI = xi*2 + (xmu+1)/2;
								int bgJ = yi*2 + (ymu+1)/2;
								int bgK = zi*2 + (zmu+1)/2;
								int bIndex = numVars*(idim-1)*2*(jdim-1)*2*bgK + numVars*(idim-1)*2*bgJ +numVars*bgI;
								for (unsigned int var = 0; var < numVars; var++) {
									if (bgWeights[bIndex] != 0) {
										bgU[bIndex +var] /= bgWeights[bIndex];
									} else {
										cout << "Empty background mish at " << xPos << ", " << yPos << ", " << zPos << endl;
									}
								}				
								bgWeights[bIndex] = 0.;
							}
						}
					}
				}
			}
		}	
	} else {
		cout << "No background observations loaded" << endl;
	}
	
	
	// Load the observations into a vector
	bgObs = new real[numbgObs*14];
	for (int m=0; m < numbgObs*14; m++) bgObs[m] = 0.;
	int p = 0;
	for (int m=0; m < bgIn.size(); m+=11) {
		real bgX = bgIn[m];
		real bgY = bgIn[m+1];
		real bgZ = exp(bgIn[m+2]);
		real bgTime = bgIn[m+3];
		if ((bgX < imin) or (bgX > imax) or
			(bgY < jmin) or (bgY > jmax) or
			(bgZ < kmin) or (bgZ > kmax)) {
			numbgObs -= 7;
			continue;
		}
		for (unsigned int n = 0; n < numVars; n++) {
			bgObs[p] = bgIn[m+4+n];
			// Error of background = 1
			bgObs[p+1] = 1.;
			bgObs[p+2] = bgX;
			bgObs[p+3] = bgY;
			bgObs[p+4] = bgZ;
			// Null type
			bgObs[p+5] = -1;
			bgObs[p+6] = bgTime;
			bgObs[p+7+n] = 1.;
			p += 14;
		}
	}	

	cout << numbgObs << " background observations loaded" << endl;	
	return numbgObs;
}

bool VarDriverXYZ::initialize(const QString& xmlfile)
{
	// Run a XYZ vortex background field
	cout << "Initializing SAMURAI XYZ" << endl;
	
	// Read XML configuration
	if (!readXMLconfig(xmlfile)) {
		cout << "Error reading XML configuration, quitting...\n";
		exit(-1);
	}
	
	// Set the initial background to zero
	imin = configHash.value("xmin").toFloat();
	imax = configHash.value("xmax").toFloat();
	iincr = configHash.value("xincr").toFloat();
	idim = (int)((imax - imin)/iincr) + 1;
	jmin = configHash.value("ymin").toFloat();
	jmax = configHash.value("ymax").toFloat();
	jincr = configHash.value("yincr").toFloat();
	jdim = (int)((jmax - jmin)/jincr) + 1;		
	kmin = configHash.value("zmin").toFloat();
	kmax = configHash.value("zmax").toFloat();
	kincr = configHash.value("zincr").toFloat();
	kdim = (int)((kmax - kmin)/kincr) + 1;
	
	// Define the sizes of the arrays we are passing to the cost function
	int uStateSize = 8*(idim-1)*(jdim-1)*(kdim-1)*(numVars);
	int bStateSize = idim*jdim*kdim*numVars;
	cout << "Physical (mish) State size = " << uStateSize << "\n";
	cout << "Nodal State size = " << bStateSize << ", Grid dimensions:\n";
	cout << "xMin\txMax\txIncr\tyMin\tyMax\tyIncr\tzMin\tzMax\tzIncr\n";
	cout << imin << "\t" <<  imax << "\t" <<  iincr << "\t";
	cout << jmin << "\t" <<  jmax << "\t" <<  jincr << "\t";
	cout << kmin << "\t" <<  kmax << "\t" <<  kincr << "\n\n";
	
	// Load the BG into a empty vector
	bgU = new real[uStateSize];
	bgWeights = new real[uStateSize];
	for (int i=0; i < uStateSize; i++) {
		bgU[i] = 0.;
		bgWeights[i] = 0.;
	}
		
	
	// Define the Reference state
	if (configHash.value("refstate") == "jordan") {
		referenceState = jordan;
	} else {
		cout << "Reference state not defined!\n";
		exit(-1);
	}
	
	cout << "Reference profile: Z\t\tQv\tRhoa\tRho\tH\tTemp\tPressure\n";
	for (real k = kmin; k < kmax+kincr; k+= kincr) {
		cout << "                   " << k << "\t";
		for (int i = 0; i < 6; i++) {
			real var = getReferenceVariable(i, k*1000);
			if (i == 0) var = bhypInvTransform(var);
			cout << setw(9) << setprecision(4)  << var << "\t";
		}
		cout << "\n";
	}
	cout << setprecision(9);
	
	// Read in the TC centers
	// Ideally, create a time-based spline from limited center fixes here
	// but just load 1 second centers into vector for now
	readTCcenters();
	
	// Get the reference center
	QTime reftime = QTime::fromString(configHash.value("reftime"), "hh:mm:ss");
	QString refstring = reftime.toString();
	bool foundref = false;
	for (unsigned int tci = 0; tci < tcVector.size(); tci++) {
		QDateTime tctime = tcVector[tci].getTime();
		if (reftime == tctime.time()) {
			QString tempstring;
			QDate refdate = tctime.date();
			QDateTime unixtime(refdate, reftime, Qt::UTC);
			configHash.insert("reflat", tempstring.setNum(tcVector[tci].getLat()));
			configHash.insert("reflon", tempstring.setNum(tcVector[tci].getLon()));
			configHash.insert("reftime", tempstring.setNum(unixtime.toTime_t()));
			cout << "Found matching reference time " << refstring.toStdString()
			<< " at " << tcVector[tci].getLat() << ", " << tcVector[tci].getLon() << "\n";
			foundref = true;
			break;
		}
	}
	if (!foundref) {
		cout << "Error finding reference time, please check date and time in XML file\n";
		exit(-1);
	}
	
	bool loadBG = configHash.value("load_background").toInt();
	int numbgObs = 0;
	if (loadBG) {
		// Set up the Gaussian Grid by a previous samurai analysis
		numbgObs = loadBackgroundObs();
	}
	
	bool adjustBG = configHash.value("adjust_background").toInt();
	if (adjustBG) {
		/* Set the minimum filter length to the background resolution, not the analysis resolution
		 to avoid artifacts when running interpolating to small mesoscale grids */
		float hfilter = configHash.value("xfilter").toFloat();
		float ares = hfilter*iincr;
		float bgres = configHash.value("backgroundroi").toFloat();
		if (ares < bgres) {
			QString bgfilter;
			bgfilter.setNum(bgres/iincr);
			configHash.insert("xfilter", bgfilter);
			configHash.insert("yfilter", bgfilter);
		}
			
		// Adjust the background field to the spline mish
		bgCostXYZ = new CostFunctionXYZ(numbgObs, bStateSize);
		
		bgCostXYZ->initialize(configHash, bgU, bgObs);
		bgCostXYZ->initState();
		bgCostXYZ->minimize();
		// Increment the variables
		bgCostXYZ->updateBG();
		
		delete bgCostXYZ;
		
		// Reset the horizontal filter
		if (ares < bgres) {
			QString afilter;
			afilter.setNum(hfilter);
			configHash.insert("xfilter", afilter);
			configHash.insert("yfilter", afilter);
		}
	}
	
	// Read in the observations, process them into weights and positions
	// Either preprocess from raw observations or load an already processed Observations.out file
	bool preprocess = true;
	if (preprocess) {
		preProcessMetObs();
	} else {
		loadMetObs();
	}
	cout << "Number of New Observations: " << obVector.size() << endl;		
	
	obCostXYZ = new CostFunctionXYZ(obVector.size(), bStateSize);
	obCostXYZ->initialize(configHash, bgU, obs);
	
	return true;
}


bool VarDriverXYZ::run()
{
	// CQRMS not used currently
	double CQRMS = 999;
	int iter=0;
	while ((CQRMS > CQTOL) and (iter < maxIter)) {
		iter++;
		cout << "Outer Loop Iteration: " << iter << endl;
		obCostXYZ->initState();
		obCostXYZ->minimize();
		// Increment the variables
		obCostXYZ->updateBG();
	}	
/*	cout << "Increment RMS Tolerance of " << CQTOL << " reached in "
		<< iter << " iterations. Writing analysis results..." << endl; */

	return true;

}


