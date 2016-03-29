/*
 *  CostFunctionRTZ.cpp
 *  samurai
 *
 *  Copyright 2008 Michael Bell. All rights reserved.
 *
 */

#include "CostFunctionRTZ.h"
#include <cmath>
#include <QTextStream>
#include <QDir>
#include <QDateTime>
#include <netcdfcpp.h>
#include <GeographicLib/TransverseMercatorExact.hpp>

CostFunctionRTZ::CostFunctionRTZ(const int& numObs, const int& stateSize)
	: CostFunction3D(numObs, stateSize)
{
}

CostFunctionRTZ::~CostFunctionRTZ()
{
}
bool CostFunctionRTZ::outputAnalysis(const QString& suffix, real* Astate)
{

	cout << "Outputting " << suffix.toStdString() << "...\n";
	// H --> to Mish for output
	QString samuraiout = "samurai_RTZ_" + suffix + ".out";
	ofstream samuraistream;
	if (configHash->value("output_txt") == "true") {
		samuraistream.open(outputPath.absoluteFilePath(samuraiout).toAscii().data());
		samuraistream << "R\tT\tZ\tu\tv\tw\tVorticity\tDivergence\tqv\trho\tT\tP\tTheta\tTheta_e\tTheta_es\t";
		samuraistream << "udr\tudt\tudz\tvdr\tvdt\tvdz\twdr\twdt\twdz\trhowdz\tMC residual\tdBZ\n";
		samuraistream.precision(10);
	}

	int analysisDim = 51;
	int analysisSize = (iDim-2)*(jDim-2)*(kDim-2);
	finalAnalysis = new real[analysisSize*analysisDim];
	real Pi = acos(-1.0);

	for (int iIndex = 1; iIndex < iDim-1; iIndex++) {
		for (int imu = 0; imu <= 1; imu++) {
			real i = iMin + DI * (iIndex + (0.5 * imu));
			real r = i*1000;
			if (i > ((iDim-2)*DI + iMin)) continue;

			for (int jIndex = 1; jIndex < jDim-1; jIndex++) {
				for (int jmu = 0; jmu <= 1; jmu++) {
					real j = jMin + DJ * (jIndex + (0.5 * jmu));
					if (j > ((jDim-2)*DJ + jMin)) continue;

					real tpw = 0;

					for (int kIndex = 1; kIndex < kDim-1; kIndex++) {
						for (int kmu = 0; kmu <= 1; kmu++) {
							real k = kMin + DK * (kIndex + (0.5 * kmu));
							if (k > ((kDim-2)*DK + kMin)) continue;

							real heightm = 1000*k;
							real rhoBar = refstate->getReferenceVariable(ReferenceVariable::rhoaref, heightm);
							real qBar = refstate->getReferenceVariable(ReferenceVariable::qvbhypref, heightm);
							real tBar = refstate->getReferenceVariable(ReferenceVariable::tempref, heightm);

							int ii = (int)((i - iMin)*DIrecip);
							int jj = (int)((j - jMin)*DJrecip);
							int kk = (int)((k - kMin)*DKrecip);
							real ibasis = 0.;
							real jbasis = 0.;
							real kbasis = 0.;
							real idbasis = 0.;
							real jdbasis = 0.;
							real kdbasis = 0.;
							real rhov = 0.;
							real rhou = 0.;
							real rhow = 0.;
							real rhovdr = 0.; real rhoudr = 0.; real rhowdr = 0.;
							real rhovdt = 0.; real rhoudt = 0.; real rhowdt = 0.;
							real rhovdz = 0.; real rhoudz = 0.; real rhowdz = 0.;
							real tprime = 0.; real tdr = 0.; real tdt = 0.; real tdz = 0.;
							real rhoadr = 0.; real rhoadt = 0.; real rhoadz = 0.;
							real qvdr = 0.; real qvdt = 0.; real qvdz = 0.;
							real pdr = 0.; real pdt = 0.; real pdz = 0.;
							real qvprime = 0.;
							real rhoprime = 0.;
							real qrprime = 0.;
							for (int var = 0; var < varDim; var++) {
								for (int kkNode = (kk-1); kkNode <= (kk+2); ++kkNode) {
									int kNode = kkNode;
									if ((kNode < 0) or (kNode >= kDim)) continue;
									for (int iiNode = (ii-1); iiNode <= (ii+2); ++iiNode) {
										int iNode = iiNode;
										if ((iNode < 0) or (iNode >= iDim)) continue;
										for (int jjNode = (jj-1); jjNode <= (jj+2); ++jjNode) {
											int jNode = jjNode;
											if ((jNode < 0) or (jNode >= jDim)) continue;
											ibasis = Basis(iNode, i, iDim-1, iMin, DI, DIrecip, 0, iBCL[var], iBCR[var]);
											jbasis = Basis(jNode, j, jDim-1, jMin, DJ, DJrecip, 0, jBCL[var], jBCR[var]);
											kbasis = Basis(kNode, k, kDim-1, kMin, DK, DKrecip, 0, kBCL[var], kBCR[var]);
											idbasis = Basis(iNode, i, iDim-1, iMin, DI, DIrecip, 1, iBCL[var], iBCR[var]);
											jdbasis = Basis(jNode, j, jDim-1, jMin, DJ, DJrecip, 1, jBCL[var], jBCR[var]);
											kdbasis = Basis(kNode, k, kDim-1, kMin, DK, DKrecip, 1, kBCL[var], kBCR[var]);
											real basis3x = ibasis*jbasis*kbasis;
											int aIndex = varDim*iDim*jDim*kNode + varDim*iDim*jNode +varDim*iNode;
											switch (var) {
												case 0:
												rhou +=  Astate[aIndex] * basis3x;
												rhoudr += Astate[aIndex] * idbasis * jbasis * kbasis;
												rhoudt += 180 * Astate[aIndex] * ibasis * jdbasis * kbasis / (i * Pi);
												rhoudz += Astate[aIndex] * ibasis * jbasis * kdbasis;
												break;
												case 1:
												rhov +=  Astate[aIndex + 1] * basis3x;
												rhovdr += Astate[aIndex + 1] * idbasis * jbasis * kbasis;
												rhovdt += 180 * Astate[aIndex + 1] * ibasis * jdbasis * kbasis / (i * Pi);
												rhovdz += Astate[aIndex + 1] * ibasis * jbasis * kdbasis;
												break;
												case 2:
												rhow += Astate[aIndex + 2] * basis3x;
												rhowdr += Astate[aIndex + 2] * idbasis * jbasis * kbasis;
												rhowdt += 180 * Astate[aIndex + 2] * ibasis * jdbasis * kbasis / (i * Pi);
												rhowdz += Astate[aIndex + 2] * ibasis * jbasis * kdbasis;
												break;
												case 3:
												tprime += Astate[aIndex + 3] * basis3x;
												tdr += Astate[aIndex + 3] * idbasis * jbasis * kbasis;
												tdt += 180 * Astate[aIndex + 3] * ibasis * jdbasis * kbasis / (i * Pi);
												tdz += Astate[aIndex + 3] * ibasis * jbasis * kdbasis;
												break;
												case 4:
												qvprime += Astate[aIndex + 4] * basis3x;
												qvdr += Astate[aIndex + 4] * idbasis * jbasis * kbasis;
												qvdt += 180 * Astate[aIndex + 4] * ibasis * jdbasis * kbasis / (i * Pi);
												qvdz += Astate[aIndex + 4] * ibasis * jbasis * kdbasis;
												break;
												case 5:
												rhoprime += Astate[aIndex + 5] * basis3x;
												rhoadr += Astate[aIndex + 5] * idbasis * jbasis * kbasis;
												rhoadt += 180 * Astate[aIndex + 5] * ibasis * jdbasis * kbasis / (i * Pi);
												rhoadz += Astate[aIndex + 5] * ibasis * jbasis * kdbasis;
												break;
												case 6:
												qrprime += Astate[aIndex + 6] * basis3x;
												break;
											}
										}
									}
								}
							}

							// Save mish values for future iterations
							int uJ = (jIndex-1)*2 + jmu;
							int uI = (iIndex-1)*2 + imu;
							int uK = (kIndex-1)*2 + kmu;
							int uIndex = varDim*uiDim*ujDim*uK +varDim*uiDim*uJ +varDim*uI;

							bgFields[uIndex] = rhou;
							bgFields[uIndex + 1] = rhov;
							bgFields[uIndex + 2] = rhow;
							bgFields[uIndex + 3] = tprime;
							bgFields[uIndex + 4] = qvprime;
							bgFields[uIndex + 5] = rhoprime;
							bgFields[uIndex + 6] = qrprime;

							if ((configHash->value("output_mish") == "false")
							and (imu or jmu or kmu)) continue;

							// Output it
							real rhoa = rhoBar + rhoprime / 100;
							real qv = refstate->bhypInvTransform(qBar + qvprime);
							real qbardz = 1000. * refstate->getReferenceVariable(ReferenceVariable::qvbhypref, heightm, 1);
							// qv derivatives multipled by 2 to account for hyperbolic transform
							qvdr = 2.0*qvdr;
							qvdt = 2.0*qvdt;
							qvdz = 2.0*(qbardz + qvdz);

							real qr;
							QString gridref = configHash->value("qr_variable");
							if (gridref == "dbz") {
								qr = qrprime*10. - 35.;
								if (qr < -35.) {
									qr = -999.;
								}
							} else {
								qr = refstate->bhypInvTransform(qrprime);
							}
							real rhoq = qv * rhoa / 1000.;
							real rho = rhoa + rhoq;
							real v = rhov / rho;
							real u = rhou / rho;
							real w = rhow / rho;
							real wspd = sqrt(u*u + v*v);
							real temp = tBar + tprime;
							real tbardz = 1000. * refstate->getReferenceVariable(ReferenceVariable::tempref, heightm, 1);
							tdz = tbardz + tdz;

							real h = 1005.7*temp + 2.501e3*qv + 9.81*heightm;
							real airpress = temp*rhoa*287./100.;
							real satvp =  exp(-6096.9385 / temp + 16.635794 - 2.711193e-2 * temp
														+ 1.673952e-5 * temp*temp + 2.433502 * log(temp));
							real vp = temp*rhoq*461./100.;
							//real vp = airpress * qv / (622 + qv);
							real press = airpress + vp;

							real pprime = press - refstate->getReferenceVariable(ReferenceVariable::pressref, heightm)/100.;
							real hprime = h - refstate->getReferenceVariable(ReferenceVariable::href, heightm);

							real RoverCp = 0.2854*(1 - 0.00028*qv);
							real theta = temp * pow((1000/press), RoverCp);
							real lcl = 2840/(3.5*log(temp) - log(vp) - 4.805) + 55.0;
							real thetae = theta * exp(((3.376/lcl) - 0.00254) * qv * (1 + 0.00081 * qv));
							real qvsat = 622 * satvp / airpress;
							real relhum = -999.;
							real thetaes = -999.;
							if (satvp != 0) {
								relhum = 100*vp/satvp;
								lcl = 2840/(3.5*log(temp) - log(satvp) - 4.805) + 55.0;
								thetaes = theta * exp(((3.376/lcl) - 0.00254) * qvsat * (1 + 0.00081 * qvsat));
							} else {
								relhum = -999.;
								thetaes = -999.;
							}
							if (relhum > 100.) {
								relhum = 100.0;
								vp = satvp;
								qv = qvsat;
							}
							real dewp = -999.0;
							if (vp != 0) {
								dewp = 237.3 * log(vp/6.1078) / (17.2694 - log(vp/6.1078)) + 273.15;
							}

							// Calculate the kinematic derivatives
							// rhoa derivatives divided by 100
							// qv derivatives multipled by 2 to account for hyperbolic transform, not exact but close enough
							rhoadr /= 100.;
							rhoadt /= 100.;
							rhoadz /= 100.;
							real rhodr = rhoadr * (1. + qv/1000.) + rhoa * qvdr/1000.;
							real rhodt = rhoadt * (1. + qv/1000.) + rhoa * qvdt/1000.;
							real rhodz = rhoadz * (1. + qv/1000.) + rhoa * qvdz/1000.;
							real rhobardz = 1000 * refstate->getReferenceVariable(ReferenceVariable::rhoref, heightm, 1);
							rhodz += rhobardz;
							real rhoabardz = 1000 * refstate->getReferenceVariable(ReferenceVariable::rhoaref, heightm, 1);
							rhoadz += rhoabardz;

							// Units 10-5
							real udr = 100. * (rhoudr - u*rhodr) / rho;
							real udt = 100. * (rhoudt - u*rhodt) / rho;
							real udz = 100. * (rhoudz - u*rhodz) / rho;

							real vdr = 100. * (rhovdr - v*rhodr) / rho;
							real vdt = 100. * (rhovdt - v*rhodt) / rho;
							real vdz = 100. * (rhovdz - v*rhodz) / rho;

							real wdr = 100. * (rhowdr - w*rhodr) / rho;
							real wdt = 100. * (rhowdt - w*rhodt) / rho;
							real wdz = 100. * (rhowdz - w*rhodz) / rho;

							// Vorticity units are 10-5
							real vorticity = 1.0e5 * (vdr * 1.0e-5 + v/r - udt * 1.0e-5);
							real divergence = 1.0e5 * (udr * 1.0e-5 + u/r + vdt * 1.0e-5);
							real s1 = 1.0e5 * (udr * 1.0e-5 + u/r - vdt * 1.0e-5);
							real s2 = 1.0e5 * (vdr * 1.0e-5 + v/r + udt * 1.0e-5);
							real strain = sqrt(s1*s1 + s2*s2);
							real okuboweiss = vorticity*vorticity - s1*s1 -s2*s2;
							real mcresidual = 1.0e5 * (rhoudr * 1.0e-5 + rhou / r + rhovdt * 1.0e-5 + rhowdz * 1.0e-5);

							// Add Coriolis parameter to relative vorticity
							real latReference = configHash->value("ref_lat").toFloat();
							real Coriolisf = 2 * 7.2921 * sin(latReference*Pi/180); // Units 10^-5 s-1
							real absVorticity = vorticity + Coriolisf;

							// Thermodynamic derivatives
							pdr = (tdr*rhoa + rhoadr*temp)*287./100. + (tdr*rhoq + (rhoadr*qv + qvdr*rhoa)*temp/1000.0)*461./100.;
							pdt = (tdt*rhoa + rhoadt*temp)*287./100. + (tdt*rhoq + (rhoadt*qv + qvdt*rhoa)*temp/1000.0)*461./100.;
							pdz = (tdz*rhoa + rhoadz*temp)*287./100. + (tdz*rhoq + (rhoadz*qv + qvdz*rhoa)*temp/1000.0)*461./100.;

							QString refmask = configHash->value("mask_reflectivity");
							if (refmask != "None") {
								real refthreshold = refmask.toFloat();
								if (qr < refthreshold) {
									u = -999.;
									v = -999.;
									w = -999.;
									wspd = -999.;
									relhum = -999.;
									hprime = -999.;
									qvprime = -999.;
									rhoprime = -999.;
									tprime = -999.;
									pprime = -999.;
									vorticity = -999.;
									absVorticity = -999.;
									divergence = -999.;
									okuboweiss = -999.;
									strain = -999.;
									tpw = -999.;
									rhou = -999.;
									rhov = -999.;
									rhow = -999.;
									rho = -999.;
									press = -999.;
									temp = -999.;
									qv = -999.;
									h = -999.;
									qr = -999.;
									udr = -999.; udt = -999.; udz = -999.;
									vdr = -999.; vdt = -999.; vdz = -999.;
									wdr = -999.; wdt = -999.; wdz = -999.;
									tdr = -999.; tdt = -999.; tdz = -999.;
									qvdr = -999.; qvdt = -999.; qvdz = -999.;
									pdr = -999.; pdt = -999.; pdz = -999.;
									rhodr = -999.; rhodt = -999.; rhodz = -999.;
									dewp = -999.;
									theta = -999.; thetae = -999.; thetaes = -999.;
								}
							}

							// Avoid singularity at origin
							if (r==0) {
								u = -999.;
								v = -999.;
								w = -999.;
								wspd = -999.;
								relhum = -999.;
								hprime = -999.;
								qvprime = -999.;
								rhoprime = -999.;
								tprime = -999.;
								pprime = -999.;
								vorticity = -999.;
								absVorticity = -999.;
								divergence = -999.;
								okuboweiss = -999.;
								strain = -999.;
								tpw = -999.;
								rhou = -999.;
								rhov = -999.;
								rhow = -999.;
								rho = -999.;
								press = -999.;
								temp = -999.;
								qv = -999.;
								h = -999.;
								qr = -999.;
								udr = -999.; udt = -999.; udz = -999.;
								vdr = -999.; vdt = -999.; vdz = -999.;
								wdr = -999.; wdt = -999.; wdz = -999.;
								tdr = -999.; tdt = -999.; tdz = -999.;
								qvdr = -999.; qvdt = -999.; qvdz = -999.;
								pdr = -999.; pdt = -999.; pdz = -999.;
								rhodr = -999.; rhodt = -999.; rhodz = -999.;
								dewp = -999.;
								theta = -999.; thetae = -999.; thetaes = -999.;
							}


							if (configHash->value("output_txt") == "true") {
								samuraistream << scientific << i << "\t" << j << "\t"  << k
								<< "\t" << u << "\t" << v << "\t" << w << "\t" << vorticity << "\t" << divergence
								<< "\t" << qv << "\t" << rho << "\t" << temp << "\t" << press
								<< "\t" << theta << "\t" << thetae << "\t" << thetaes << "\t"
								<< udr << "\t" << udt << "\t" << udz << "\t"
								<< vdr << "\t" << vdt << "\t" << vdz << "\t"
								<< wdr << "\t" << wdt << "\t" << wdz << "\t"
								<< rhowdz * 100. << "\t" << mcresidual << "\t" << qr << "\n";
							}

							// Sum up the TPW in the vertical, top level is tpw
							tpw += qv * rhoa * DK;

							// On the nodes
							if (!imu and !jmu and !kmu){
								int fIndex = (iDim-2)*(jDim-2)*(kDim-2);
								int posIndex = (iDim-2)*(jDim-2)*(kIndex-1) + (iDim-2)*(jIndex-1) + (iIndex-1);
								finalAnalysis[fIndex * 0 + posIndex] = u;
								finalAnalysis[fIndex * 1 + posIndex] = v;
								finalAnalysis[fIndex * 2 + posIndex] = w;
								finalAnalysis[fIndex * 3 + posIndex] = wspd;
								finalAnalysis[fIndex * 4 + posIndex] = relhum;
								finalAnalysis[fIndex * 5 + posIndex] = hprime;
								if (qvprime != -999) {
									finalAnalysis[fIndex * 6 + posIndex] = 2*qvprime;
								} else {
									finalAnalysis[fIndex * 6 + posIndex] = -999;
								}
								finalAnalysis[fIndex * 7 + posIndex] = rhoprime;
								finalAnalysis[fIndex * 8 + posIndex] = tprime;
								finalAnalysis[fIndex * 9 + posIndex] = pprime;
								finalAnalysis[fIndex * 10 + posIndex] = vorticity;
								finalAnalysis[fIndex * 11 + posIndex] = divergence;
								finalAnalysis[fIndex * 12 + posIndex] = okuboweiss;
								finalAnalysis[fIndex * 13 + posIndex] = strain;
								finalAnalysis[fIndex * 14 + posIndex] = tpw;
								finalAnalysis[fIndex * 15 + posIndex] = rhou;
								finalAnalysis[fIndex * 16 + posIndex] = rhov;
								finalAnalysis[fIndex * 17 + posIndex] = rhow;
								finalAnalysis[fIndex * 18 + posIndex] = rho;
								finalAnalysis[fIndex * 19 + posIndex] = press;
								finalAnalysis[fIndex * 20 + posIndex] = temp;
								finalAnalysis[fIndex * 21 + posIndex] = qv;
								finalAnalysis[fIndex * 22 + posIndex] = h;
								finalAnalysis[fIndex * 23 + posIndex] = qr;
								finalAnalysis[fIndex * 24 + posIndex] = absVorticity;
								finalAnalysis[fIndex * 25 + posIndex] = dewp;
								finalAnalysis[fIndex * 26 + posIndex] = theta;
								finalAnalysis[fIndex * 27 + posIndex] = thetae;
								finalAnalysis[fIndex * 28 + posIndex] = thetaes;
								finalAnalysis[fIndex * 29 + posIndex] = udr;
								finalAnalysis[fIndex * 30 + posIndex] = vdr;
								finalAnalysis[fIndex * 31 + posIndex] = wdr;
								finalAnalysis[fIndex * 32 + posIndex] = udt;
								finalAnalysis[fIndex * 33 + posIndex] = vdt;
								finalAnalysis[fIndex * 34 + posIndex] = wdt;
								finalAnalysis[fIndex * 35 + posIndex] = udz;
								finalAnalysis[fIndex * 36 + posIndex] = vdz;
								finalAnalysis[fIndex * 37 + posIndex] = wdz;
								finalAnalysis[fIndex * 38 + posIndex] = tdr;
								finalAnalysis[fIndex * 39 + posIndex] = tdt;
								finalAnalysis[fIndex * 40 + posIndex] = tdz;
								finalAnalysis[fIndex * 41 + posIndex] = qvdr;
								finalAnalysis[fIndex * 42 + posIndex] = qvdt;
								finalAnalysis[fIndex * 43 + posIndex] = qvdz;
								finalAnalysis[fIndex * 44 + posIndex] = pdr;
								finalAnalysis[fIndex * 45 + posIndex] = pdt;
								finalAnalysis[fIndex * 46 + posIndex] = pdz;
								finalAnalysis[fIndex * 47 + posIndex] = rhodr;
								finalAnalysis[fIndex * 48 + posIndex] = rhodt;
								finalAnalysis[fIndex * 49 + posIndex] = rhodz;
								finalAnalysis[fIndex * 50 + posIndex] = mcresidual;
							}
						}
					}
				}
			}
		}
	}

	QString fileName = "samurai_RTZ_" + suffix;
	QString outFileName = outputPath.absoluteFilePath(fileName);

	// Write the Obs to a summary text file
	if (configHash->value("output_qc") == "true") {
		QString qcout = "samurai_QC_" + suffix + ".out";
		QString qcFileName = outputPath.absoluteFilePath(qcout);
		ofstream qcstream(qcFileName.toAscii().data());
		ostream_iterator<string> os(qcstream, "\t ");
		*os++ = "Observation";
		*os++ = "Inverse Error";
		*os++ = "R";
		*os++ = "T";
		*os++ = "Z";
		*os++ = "Type";
		*os++ = "Time";
		*os++ = "rhou";
		*os++ = "rhov";
		*os++ = "rhow";
		*os++ = "T'";
		*os++ = "qv'";
		*os++ = "rhoa'";
		*os++ = "qr";
		*os++ = "Analysis";
		*os++ = "Background";
		qcstream << endl;
		qcstream.precision(10);

		ostream_iterator<real> od(qcstream, "\t ");
		for (int m = 0; m < mObs; m++) {
			int mi = m*(7+varDim*derivDim);
			real i = obsVector[mi+2];
			real j = obsVector[mi+3];
			real k = obsVector[mi+4];
			real tempsum = 0;
			int ii = (int)((i - iMin)*DIrecip);
			int jj = (int)((j - jMin)*DJrecip);
			int kk = (int)((k - kMin)*DKrecip);
			real ibasis = 0;
			real jbasis = 0;
			real kbasis = 0;
			for (int var = 0; var < varDim; var++) {
				for (int d = 0; d < derivDim; d++) {
					int wgt_index = mi + (7*(d+1)) + var;
					if (!obsVector[wgt_index]) continue;
					for (int kkNode = (kk-1); kkNode <= (kk+2); ++kkNode) {
						int kNode = kkNode;
						if ((kNode < 0) or (kNode >= kDim)) continue;
						kbasis = Basis(kNode, k, kDim-1, kMin, DK, DKrecip, derivative[d][2], kBCL[var], kBCR[var]);
						for (int iiNode = (ii-1); iiNode <= (ii+2); ++iiNode) {
							int iNode = iiNode;
							if ((iNode < 0) or (iNode >= iDim)) continue;
							ibasis = Basis(iNode, i, iDim-1, iMin, DI, DIrecip, derivative[d][1], iBCL[var], iBCR[var]);
							for (int jjNode = (jj-1); jjNode <= (jj+2); ++jjNode) {
								int jNode = jjNode;
								if ((jNode < 0) or (jNode >= jDim)) continue;
								int aIndex = varDim*iDim*jDim*kNode + varDim*iDim*jNode +varDim*iNode;
								jbasis = Basis(jNode, j, jDim-1, jMin, DJ, DJrecip, derivative[d][0], jBCL[var], jBCR[var]);
								tempsum += Astate[aIndex + var] * ibasis * jbasis * kbasis * obsVector[wgt_index];
							}
						}
					}
				}
			}
			for (int t=0; t < 6; t++) {
				*od++ = obsVector[mi+t];
			}
			int unixtime = (int)obsVector[mi+6];
			QDateTime obtime;
			obtime.setTime_t(unixtime);
			obtime.setTimeSpec(Qt::UTC);
			QString timestring = obtime.toString("hh:mm:ss.zzz");
			qcstream << timestring.toStdString() << "\t";

			// Multiply the weight by the ob -- Observations.in has individual weights alreadt
			// Only non-derivative for now
			for (int t=7; t<14; t++) {
				*od++ = obsVector[mi+t] * obsVector[mi];
			}

			*od++ = tempsum;
			*od++ = obsVector[mi]-innovation[m];
			qcstream << endl;

		}
	}

	adjustInternalDomain(-1);

	// Write out to a netCDF file
	if (configHash->value("output_netcdf") == "true") {
		QString cdfFileName = outFileName + ".nc";
		if (!writeNetCDF(outputPath.absoluteFilePath(cdfFileName)))
		cout << "Error writing netcdf file " << cdfFileName.toStdString() << endl;
	}
	// Write out to an asi file
	if (configHash->value("output_asi") == "true") {
		QString asiFileName = outFileName + ".asi";
		if (!writeAsi(outputPath.absoluteFilePath(asiFileName)))
		cout << "Error writing asi file " << asiFileName.toStdString() << endl;
	}
	// Set the domain back
	adjustInternalDomain(1);

	// Free the memory for the analysis variables
	delete[] finalAnalysis;

	return true;

}

bool CostFunctionRTZ::writeNetCDF(const QString& netcdfFileName)
{
	NcError err(NcError::verbose_nonfatal);
	int NC_ERR = 0;

	// Create the file.
	NcFile dataFile(netcdfFileName.toAscii(), NcFile::Replace);

	// Check to see if the file was created.
	if(!dataFile.is_valid())
		return NC_ERR;

	// Define the dimensions. NetCDF will hand back an ncDim object for
	// each.
	NcDim *lvlDim, *radDim, *thetaDim, *timeDim;
	if (!(radDim = dataFile.add_dim("radius", iDim)))
		return NC_ERR;
	if (!(thetaDim = dataFile.add_dim("theta", jDim)))
		return NC_ERR;
	if (!(lvlDim = dataFile.add_dim("altitude", kDim)))
		return NC_ERR;
	// Add an unlimited dimension...
	if (!(timeDim = dataFile.add_dim("time")))
		return NC_ERR;

	// Define the coordinate variables.
    NcVar *lvlVar, *timeVar, *radVar, *thetaVar;
    /* NcVar *latVar, *lonVar;
	if (!(lonVar = dataFile.add_var("longitude", ncFloat, radDim)))
		return NC_ERR;
	if (!(latVar = dataFile.add_var("latitude", ncFloat, thetaDim)))
		return NC_ERR; */
  if (!(radVar = dataFile.add_var("radius", ncFloat, radDim)))
		return NC_ERR;
	if (!(thetaVar = dataFile.add_var("theta", ncFloat, thetaDim)))
		return NC_ERR;
	if (!(lvlVar = dataFile.add_var("altitude", ncFloat, lvlDim)))
		return NC_ERR;
	if (!(timeVar = dataFile.add_var("time", ncInt, timeDim)))
		return NC_ERR;

	// Define units attributes for coordinate vars. This attaches a
	// text attribute to each of the coordinate variables, containing
	// the units.
	/* if (!latVar->add_att("units", "degrees_north"))
		return NC_ERR;
	if (!lonVar->add_att("units", "degrees_east"))
		return NC_ERR; */
  if (!radVar->add_att("units", "km"))
		return NC_ERR;
	if (!thetaVar->add_att("units", "degrees"))
		return NC_ERR;
	if (!lvlVar->add_att("units", "km"))
		return NC_ERR;
	if (!timeVar->add_att("units", "seconds since 1970-01-01 00:00:00 +0000"))
		return NC_ERR;

	// Define the netCDF variables
	NcVar *u, *v, *w, *wspd, *relhum, *hprime, *qvprime, *rhoprime, *tprime, *pprime;
	NcVar *vorticity, *divergence, *okuboweiss, *strain, *tpw, *rhou, *rhov, *rhow;
	NcVar *rho, *press, *temp, *qv, *h, *qr, *absVorticity;
  NcVar *dudr, *dvdr, *dwdr, *dudt, *dvdt, *dwdt, *dudz, *dvdz, *dwdz;
	NcVar *dtdr, *dqdr, *dpdr, *dtdt, *dqdt, *dpdt, *dtdz, *dqdz, *dpdz;
  NcVar *drhodr, *drhodt, *drhodz;
	NcVar *dewp, *theta, *thetae, *thetaes, *mcresidual;

	if (!(u = dataFile.add_var("U", ncFloat, timeDim,
                               lvlDim, thetaDim, radDim)))
		return NC_ERR;
	if (!(v = dataFile.add_var("V", ncFloat, timeDim,
							   lvlDim, thetaDim, radDim)))
		return NC_ERR;
	if (!(w = dataFile.add_var("W", ncFloat, timeDim,
							   lvlDim, thetaDim, radDim)))
		return NC_ERR;
	if (!(wspd = dataFile.add_var("WSPD", ncFloat, timeDim,
                                  lvlDim, thetaDim, radDim)))
		return NC_ERR;
	if (!(relhum = dataFile.add_var("RH", ncFloat, timeDim,
                                    lvlDim, thetaDim, radDim)))
		return NC_ERR;
	if (!(hprime = dataFile.add_var("HP", ncFloat, timeDim,
                                    lvlDim, thetaDim, radDim)))
		return NC_ERR;
	if (!(qvprime = dataFile.add_var("QVP", ncFloat, timeDim,
                                     lvlDim, thetaDim, radDim)))
		return NC_ERR;
	if (!(rhoprime = dataFile.add_var("RHOAP", ncFloat, timeDim,
                                      lvlDim, thetaDim, radDim)))
		return NC_ERR;
	if (!(tprime = dataFile.add_var("TP", ncFloat, timeDim,
									lvlDim, thetaDim, radDim)))
		return NC_ERR;
	if (!(pprime = dataFile.add_var("PP", ncFloat, timeDim,
									lvlDim, thetaDim, radDim)))
		return NC_ERR;
	if (!(vorticity = dataFile.add_var("VORT", ncFloat, timeDim,
                                       lvlDim, thetaDim, radDim)))
		return NC_ERR;
	if (!(divergence = dataFile.add_var("DIV", ncFloat, timeDim,
                                        lvlDim, thetaDim, radDim)))
		return NC_ERR;
	if (!(okuboweiss = dataFile.add_var("OW", ncFloat, timeDim,
                                        lvlDim, thetaDim, radDim)))
		return NC_ERR;
	if (!(strain = dataFile.add_var("STRAIN", ncFloat, timeDim,
                                    lvlDim, thetaDim, radDim)))
		return NC_ERR;
	if (!(tpw = dataFile.add_var("TPW", ncFloat, timeDim,
                                 lvlDim, thetaDim, radDim)))
		return NC_ERR;
	if (!(rhou = dataFile.add_var("RHOU", ncFloat, timeDim,
                                  lvlDim, thetaDim, radDim)))
		return NC_ERR;
	if (!(rhov = dataFile.add_var("RHOV", ncFloat, timeDim,
                                  lvlDim, thetaDim, radDim)))
		return NC_ERR;
	if (!(rhow = dataFile.add_var("RHOW", ncFloat, timeDim,
                                  lvlDim, thetaDim, radDim)))
		return NC_ERR;
	if (!(rho = dataFile.add_var("RHOA", ncFloat, timeDim,
                                 lvlDim, thetaDim, radDim)))
		return NC_ERR;
	if (!(press = dataFile.add_var("P", ncFloat, timeDim,
                                   lvlDim, thetaDim, radDim)))
		return NC_ERR;
	if (!(temp = dataFile.add_var("T", ncFloat, timeDim,
								  lvlDim, thetaDim, radDim)))
		return NC_ERR;
	if (!(qv = dataFile.add_var("QV", ncFloat, timeDim,
								lvlDim, thetaDim, radDim)))
		return NC_ERR;
	if (!(h = dataFile.add_var("H", ncFloat, timeDim,
                               lvlDim, thetaDim, radDim)))
		return NC_ERR;
  if (configHash->value("qr_variable") == "dbz") {
    	if (!(qr = dataFile.add_var("DBZ", ncFloat, timeDim,
                                    lvlDim, thetaDim, radDim)))
            return NC_ERR;
	} else {
      if (!(qr = dataFile.add_var("QR", ncFloat, timeDim,
                                    lvlDim, thetaDim, radDim)))
            return NC_ERR;
    }
  if (!(absVorticity = dataFile.add_var("ABSVORT", ncFloat, timeDim,
                                       lvlDim, thetaDim, radDim)))
		return NC_ERR;
	if (!(dewp = dataFile.add_var("DEWPOINT", ncFloat, timeDim,
								  lvlDim, thetaDim, radDim)))
		return NC_ERR;
	if (!(theta = dataFile.add_var("THETA", ncFloat, timeDim,
								   lvlDim, thetaDim, radDim)))
		return NC_ERR;
	if (!(thetae = dataFile.add_var("THETAE", ncFloat, timeDim,
									lvlDim, thetaDim, radDim)))
		return NC_ERR;
	if (!(thetaes = dataFile.add_var("THETAES", ncFloat, timeDim,
									 lvlDim, thetaDim, radDim)))
		return NC_ERR;
  if (!(dudr = dataFile.add_var("DUDR", ncFloat, timeDim,
                                  lvlDim, thetaDim, radDim)))
		return NC_ERR;
	if (!(dvdr = dataFile.add_var("DVDR", ncFloat, timeDim,
                                  lvlDim, thetaDim, radDim)))
		return NC_ERR;
	if (!(dwdr = dataFile.add_var("DWDR", ncFloat, timeDim,
                                  lvlDim, thetaDim, radDim)))
		return NC_ERR;
	if (!(dudt = dataFile.add_var("DUDT", ncFloat, timeDim,
                                  lvlDim, thetaDim, radDim)))
		return NC_ERR;
	if (!(dvdt = dataFile.add_var("DVDT", ncFloat, timeDim,
                                  lvlDim, thetaDim, radDim)))
		return NC_ERR;
	if (!(dwdt = dataFile.add_var("DWDT", ncFloat, timeDim,
                                  lvlDim, thetaDim, radDim)))
		return NC_ERR;
	if (!(dudz = dataFile.add_var("DUDZ", ncFloat, timeDim,
                                  lvlDim, thetaDim, radDim)))
		return NC_ERR;
	if (!(dvdz = dataFile.add_var("DVDZ", ncFloat, timeDim,
                                  lvlDim, thetaDim, radDim)))
		return NC_ERR;
	if (!(dwdz = dataFile.add_var("DWDZ", ncFloat, timeDim,
                                  lvlDim, thetaDim, radDim)))
		return NC_ERR;
	if (!(dtdr = dataFile.add_var("DTDR", ncFloat, timeDim,
                                  lvlDim, thetaDim, radDim)))
		return NC_ERR;
	if (!(dqdr = dataFile.add_var("DQVDR", ncFloat, timeDim,
                                  lvlDim, thetaDim, radDim)))
		return NC_ERR;
	if (!(dpdr = dataFile.add_var("DPDR", ncFloat, timeDim,
                                  lvlDim, thetaDim, radDim)))
		return NC_ERR;
	if (!(dtdt = dataFile.add_var("DTDT", ncFloat, timeDim,
                                  lvlDim, thetaDim, radDim)))
		return NC_ERR;
	if (!(dqdt = dataFile.add_var("DQVDT", ncFloat, timeDim,
                                  lvlDim, thetaDim, radDim)))
		return NC_ERR;
	if (!(dpdt = dataFile.add_var("DPDT", ncFloat, timeDim,
                                  lvlDim, thetaDim, radDim)))
		return NC_ERR;
	if (!(dtdz = dataFile.add_var("DTDZ", ncFloat, timeDim,
                                  lvlDim, thetaDim, radDim)))
		return NC_ERR;
	if (!(dqdz = dataFile.add_var("DQVDZ", ncFloat, timeDim,
                                  lvlDim, thetaDim, radDim)))
		return NC_ERR;
	if (!(dpdz = dataFile.add_var("DPDZ", ncFloat, timeDim,
                                  lvlDim, thetaDim, radDim)))
		return NC_ERR;
  if (!(drhodr = dataFile.add_var("DRHODR", ncFloat, timeDim,
                                    lvlDim, thetaDim, radDim)))
		return NC_ERR;
	if (!(drhodt = dataFile.add_var("DRHODT", ncFloat, timeDim,
                                    lvlDim, thetaDim, radDim)))
		return NC_ERR;
	if (!(drhodz = dataFile.add_var("DRHODZ", ncFloat, timeDim,
                                    lvlDim, thetaDim, radDim)))
		return NC_ERR;
	if (!(mcresidual = dataFile.add_var("MCRESIDUAL", ncFloat, timeDim,
										lvlDim, thetaDim, radDim)))
		return NC_ERR;

	// Define units attributes for data variables.
	if (!u->add_att("units", "m s-1"))
		return NC_ERR;
	if (!v->add_att("units", "m s-1"))
		return NC_ERR;
	if (!w->add_att("units", "m s-1"))
		return NC_ERR;
	if (!wspd->add_att("units", "m s-1"))
		return NC_ERR;
	if (!relhum->add_att("units", "percent"))
		return NC_ERR;
	if (!hprime->add_att("units", "kJ"))
		return NC_ERR;
	if (!qvprime->add_att("units", "g kg-1"))
		return NC_ERR;
	if (!rhoprime->add_att("units", "kg m-3"))
		return NC_ERR;
	if (!tprime->add_att("units", "K"))
		return NC_ERR;
	if (!pprime->add_att("units", "hPa"))
		return NC_ERR;
	if (!vorticity->add_att("units", "10-5s-1"))
		return NC_ERR;
	if (!divergence->add_att("units", "10-5s-1"))
		return NC_ERR;
	if (!okuboweiss->add_att("units", "10-10s-1"))
		return NC_ERR;
	if (!strain->add_att("units", "10-5s-1"))
		return NC_ERR;
	if (!tpw->add_att("units", "mm"))
		return NC_ERR;
	if (!rhou->add_att("units", "kg m-2s-1"))
		return NC_ERR;
	if (!rhov->add_att("units", "kg m-2s-1"))
		return NC_ERR;
	if (!rhow->add_att("units", "kg m-2s-1"))
		return NC_ERR;
	if (!rho->add_att("units", "kg m-3"))
		return NC_ERR;
	if (!press->add_att("units", "hPa"))
		return NC_ERR;
	if (!temp->add_att("units", "K"))
		return NC_ERR;
	if (!qv->add_att("units", "g kg-1"))
		return NC_ERR;
	if (!h->add_att("units", "kJ"))
		return NC_ERR;
  if (configHash->value("qr_variable") == "dbz") {
    if (!qr->add_att("units", "dBZ"))
        return NC_ERR;
    } else {
        if (!qr->add_att("units", "g kg-1"))
          return NC_ERR;
    }
  if (!absVorticity->add_att("units", "10-5s-1"))
		return NC_ERR;
	if (!dewp->add_att("units", "K"))
		return NC_ERR;
	if (!theta->add_att("units", "K"))
		return NC_ERR;
	if (!thetae->add_att("units", "K"))
		return NC_ERR;
	if (!thetaes->add_att("units", "K"))
		return NC_ERR;
	if (!dudr->add_att("units", "10-5s-1"))
		return NC_ERR;
	if (!dvdr->add_att("units", "10-5s-1"))
		return NC_ERR;
	if (!dwdr->add_att("units", "10-5s-1"))
		return NC_ERR;
	if (!dudt->add_att("units", "10-5s-1"))
		return NC_ERR;
	if (!dvdt->add_att("units", "10-5s-1"))
		return NC_ERR;
	if (!dwdt->add_att("units", "10-5s-1"))
		return NC_ERR;
	if (!dudz->add_att("units", "10-5s-1"))
		return NC_ERR;
	if (!dvdz->add_att("units", "10-5s-1"))
		return NC_ERR;
	if (!dwdz->add_att("units", "10-5s-1"))
		return NC_ERR;
  if (!dtdr->add_att("units", "K km-1"))
		return NC_ERR;
	if (!dqdr->add_att("units", "g kg-1 km-1"))
		return NC_ERR;
	if (!dpdr->add_att("units", "hPa km-1"))
		return NC_ERR;
	if (!dtdt->add_att("units", "K km-1"))
		return NC_ERR;
	if (!dqdt->add_att("units", "g kg-1 km-1"))
		return NC_ERR;
	if (!dpdt->add_att("units", "hPa km-1"))
		return NC_ERR;
	if (!dtdz->add_att("units", "K km-1"))
		return NC_ERR;
	if (!dqdz->add_att("units", "g kg-1 km-1"))
		return NC_ERR;
	if (!dpdz->add_att("units", "hPa km-1"))
		return NC_ERR;
	if (!drhodr->add_att("units", "kg m-3 km-1"))
		return NC_ERR;
	if (!drhodt->add_att("units", "kg m-3 km-1"))
		return NC_ERR;
	if (!drhodz->add_att("units", "kg m-3 km-1"))
		return NC_ERR;
	if (!mcresidual->add_att("units", "10-5s-1"))
		return NC_ERR;

	// Define long names for data variables.
	if (!u->add_att("long_name", "u wind component"))
		return NC_ERR;
	if (!v->add_att("long_name", "v wind component"))
		return NC_ERR;
	if (!w->add_att("long_name", "w wind component"))
		return NC_ERR;
	if (!wspd->add_att("long_name", "wind speed"))
		return NC_ERR;
	if (!relhum->add_att("long_name", "relative humidity"))
		return NC_ERR;
	if (!hprime->add_att("long_name", "moist static energy perturbation"))
		return NC_ERR;
	if (!qvprime->add_att("long_name", "water vapor mixing ratio perturbation"))
		return NC_ERR;
	if (!rhoprime->add_att("long_name", "air density perturbation"))
		return NC_ERR;
	if (!tprime->add_att("long_name", "temperature perturbation"))
		return NC_ERR;
	if (!pprime->add_att("long_name", "pressure perturbation"))
		return NC_ERR;
	if (!vorticity->add_att("long_name", "vertical vorticity"))
		return NC_ERR;
	if (!divergence->add_att("long_name", "horizontal divergence"))
		return NC_ERR;
	if (!okuboweiss->add_att("long_name", "Okubo-Weiss parameter"))
		return NC_ERR;
	if (!strain->add_att("long_name", "horizontal strain"))
		return NC_ERR;
	if (!tpw->add_att("long_name", "total precipitable water"))
		return NC_ERR;
	if (!rhou->add_att("long_name", "mass-weighted u wind component"))
		return NC_ERR;
	if (!rhov->add_att("long_name", "mass-weighted v wind component"))
		return NC_ERR;
	if (!rhow->add_att("long_name", "mass-weighted w wind component"))
		return NC_ERR;
	if (!rho->add_att("long_name", "density"))
		return NC_ERR;
	if (!press->add_att("long_name", "pressure"))
		return NC_ERR;
	if (!temp->add_att("long_name", "temperature"))
		return NC_ERR;
	if (!qv->add_att("long_name", "water vapor mixing ratio"))
		return NC_ERR;
	if (!h->add_att("long_name", "moist static energy"))
		return NC_ERR;
  if (configHash->value("qr_variable") == "dbz") {
    if (!qr->add_att("long_name", "radar reflectivity"))
        return NC_ERR;
    } else {
        if (!qr->add_att("long_name", "precipitation mixing ratio"))
          return NC_ERR;
    }
  if (!absVorticity->add_att("long_name", "absolute vertical vorticity"))
		return NC_ERR;
	if (!dewp->add_att("long_name", "dewpoint temperature"))
		return NC_ERR;
	if (!theta->add_att("long_name", "potential temperature"))
		return NC_ERR;
	if (!thetae->add_att("long_name", "equivalent potential temperature"))
		return NC_ERR;
	if (!thetaes->add_att("long_name", "saturation equivalent potential temperature"))
		return NC_ERR;
  if (!dudr->add_att("long_name", "wind gradient"))
		return NC_ERR;
	if (!dvdr->add_att("long_name", "wind gradient"))
		return NC_ERR;
	if (!dwdr->add_att("long_name", "wind gradient"))
		return NC_ERR;
	if (!dudt->add_att("long_name", "wind gradient"))
		return NC_ERR;
	if (!dvdt->add_att("long_name", "wind gradient"))
		return NC_ERR;
	if (!dwdt->add_att("long_name", "wind gradient"))
		return NC_ERR;
	if (!dudz->add_att("long_name", "wind gradient"))
		return NC_ERR;
	if (!dvdz->add_att("long_name", "wind gradient"))
		return NC_ERR;
	if (!dwdz->add_att("long_name", "wind gradient"))
		return NC_ERR;
	if (!dtdr->add_att("long_name", "temperature gradient"))
		return NC_ERR;
	if (!dqdr->add_att("long_name", "moisture gradient"))
		return NC_ERR;
	if (!dpdr->add_att("long_name", "pressure gradient"))
		return NC_ERR;
	if (!dtdt->add_att("long_name", "temperature gradient"))
		return NC_ERR;
	if (!dqdt->add_att("long_name", "moisture gradient"))
		return NC_ERR;
	if (!dpdt->add_att("long_name", "pressure gradient"))
		return NC_ERR;
	if (!dtdz->add_att("long_name", "temperature gradient"))
		return NC_ERR;
	if (!dqdz->add_att("long_name", "moisture gradient"))
		return NC_ERR;
	if (!dpdz->add_att("long_name", "pressure gradient"))
		return NC_ERR;
  if (!drhodr->add_att("long_name", "density gradient"))
		return NC_ERR;
	if (!drhodt->add_att("long_name", "density gradient"))
		return NC_ERR;
	if (!drhodz->add_att("long_name", "density gradient"))
		return NC_ERR;
	if (!mcresidual->add_att("long_name", "residual from mass continuity equation"))
		return NC_ERR;

	// Define missing data
	if (!u->add_att("missing_value", -999.f))
		return NC_ERR;
	if (!v->add_att("missing_value", -999.f))
		return NC_ERR;
	if (!w->add_att("missing_value", -999.f))
		return NC_ERR;
	if (!wspd->add_att("missing_value", -999.f))
		return NC_ERR;
	if (!relhum->add_att("missing_value", -999.f))
		return NC_ERR;
	if (!hprime->add_att("missing_value", -999.f))
		return NC_ERR;
	if (!qvprime->add_att("missing_value", -999.f))
		return NC_ERR;
	if (!rhoprime->add_att("missing_value", -999.f))
		return NC_ERR;
	if (!tprime->add_att("missing_value", -999.f))
		return NC_ERR;
	if (!pprime->add_att("missing_value", -999.f))
		return NC_ERR;
	if (!vorticity->add_att("missing_value", -999.f))
		return NC_ERR;
	if (!divergence->add_att("missing_value", -999.f))
		return NC_ERR;
	if (!okuboweiss->add_att("missing_value", -999.f))
		return NC_ERR;
	if (!strain->add_att("missing_value", -999.f))
		return NC_ERR;
	if (!tpw->add_att("missing_value", -999.f))
		return NC_ERR;
	if (!rhou->add_att("missing_value", -999.f))
		return NC_ERR;
	if (!rhov->add_att("missing_value", -999.f))
		return NC_ERR;
	if (!rhow->add_att("missing_value", -999.f))
		return NC_ERR;
	if (!rho->add_att("missing_value", -999.f))
		return NC_ERR;
	if (!press->add_att("missing_value", -999.f))
		return NC_ERR;
	if (!temp->add_att("missing_value", -999.f))
		return NC_ERR;
	if (!qv->add_att("missing_value", -999.f))
		return NC_ERR;
	if (!h->add_att("missing_value", -999.f))
		return NC_ERR;
	if (!qr->add_att("missing_value", -999.f))
		return NC_ERR;
	if (!absVorticity->add_att("missing_value", -999.f))
		return NC_ERR;
	if (!dewp->add_att("missing_value", -999.f))
		return NC_ERR;
	if (!theta->add_att("missing_value", -999.f))
		return NC_ERR;
	if (!thetae->add_att("missing_value", -999.f))
		return NC_ERR;
	if (!thetaes->add_att("missing_value", -999.f))
		return NC_ERR;
  if (!dudr->add_att("missing_value", -999.f))
		return NC_ERR;
	if (!dvdr->add_att("missing_value", -999.f))
		return NC_ERR;
	if (!dwdr->add_att("missing_value", -999.f))
		return NC_ERR;
	if (!dudt->add_att("missing_value", -999.f))
		return NC_ERR;
	if (!dvdt->add_att("missing_value", -999.f))
		return NC_ERR;
	if (!dwdt->add_att("missing_value", -999.f))
		return NC_ERR;
	if (!dudz->add_att("missing_value", -999.f))
		return NC_ERR;
	if (!dvdz->add_att("missing_value", -999.f))
		return NC_ERR;
	if (!dwdz->add_att("missing_value", -999.f))
		return NC_ERR;
  if (!dtdr->add_att("missing_value", -999.f))
		return NC_ERR;
	if (!dqdr->add_att("missing_value", -999.f))
		return NC_ERR;
	if (!dpdr->add_att("missing_value", -999.f))
		return NC_ERR;
	if (!dtdt->add_att("missing_value", -999.f))
		return NC_ERR;
	if (!dqdt->add_att("missing_value", -999.f))
		return NC_ERR;
	if (!dpdt->add_att("missing_value", -999.f))
		return NC_ERR;
	if (!dtdz->add_att("missing_value", -999.f))
		return NC_ERR;
	if (!dqdz->add_att("missing_value", -999.f))
		return NC_ERR;
	if (!dpdz->add_att("missing_value", -999.f))
		return NC_ERR;
	if (!drhodr->add_att("missing_value", -999.f))
		return NC_ERR;
	if (!drhodt->add_att("missing_value", -999.f))
		return NC_ERR;
	if (!drhodz->add_att("missing_value", -999.f))
		return NC_ERR;
	if (!mcresidual->add_att("missing_value", -999.f))
		return NC_ERR;

	// Define _Fill_Value for NCL users
	if (!u->add_att("_FillValue", -999.f))
		return NC_ERR;
	if (!v->add_att("_FillValue", -999.f))
		return NC_ERR;
	if (!w->add_att("_FillValue", -999.f))
		return NC_ERR;
	if (!wspd->add_att("_FillValue", -999.f))
		return NC_ERR;
	if (!relhum->add_att("_FillValue", -999.f))
		return NC_ERR;
	if (!hprime->add_att("_FillValue", -999.f))
		return NC_ERR;
	if (!qvprime->add_att("_FillValue", -999.f))
		return NC_ERR;
	if (!rhoprime->add_att("_FillValue", -999.f))
		return NC_ERR;
	if (!tprime->add_att("_FillValue", -999.f))
		return NC_ERR;
	if (!pprime->add_att("_FillValue", -999.f))
		return NC_ERR;
	if (!vorticity->add_att("_FillValue", -999.f))
		return NC_ERR;
	if (!divergence->add_att("_FillValue", -999.f))
		return NC_ERR;
	if (!okuboweiss->add_att("_FillValue", -999.f))
		return NC_ERR;
	if (!strain->add_att("_FillValue", -999.f))
		return NC_ERR;
	if (!tpw->add_att("_FillValue", -999.f))
		return NC_ERR;
	if (!rhou->add_att("_FillValue", -999.f))
		return NC_ERR;
	if (!rhov->add_att("_FillValue", -999.f))
		return NC_ERR;
	if (!rhow->add_att("_FillValue", -999.f))
		return NC_ERR;
	if (!rho->add_att("_FillValue", -999.f))
		return NC_ERR;
	if (!press->add_att("_FillValue", -999.f))
		return NC_ERR;
	if (!temp->add_att("_FillValue", -999.f))
		return NC_ERR;
	if (!qv->add_att("_FillValue", -999.f))
		return NC_ERR;
	if (!h->add_att("_FillValue", -999.f))
		return NC_ERR;
	if (!qr->add_att("_FillValue", -999.f))
		return NC_ERR;
	if (!absVorticity->add_att("_FillValue", -999.f))
		return NC_ERR;
	if (!dewp->add_att("_FillValue", -999.f))
		return NC_ERR;
	if (!theta->add_att("_FillValue", -999.f))
		return NC_ERR;
	if (!thetae->add_att("_FillValue", -999.f))
		return NC_ERR;
	if (!thetaes->add_att("_FillValue", -999.f))
		return NC_ERR;
  if (!dudr->add_att("_FillValue", -999.f))
		return NC_ERR;
	if (!dvdr->add_att("_FillValue", -999.f))
		return NC_ERR;
	if (!dwdr->add_att("_FillValue", -999.f))
		return NC_ERR;
	if (!dudt->add_att("_FillValue", -999.f))
		return NC_ERR;
	if (!dvdt->add_att("_FillValue", -999.f))
		return NC_ERR;
	if (!dwdt->add_att("_FillValue", -999.f))
		return NC_ERR;
	if (!dudz->add_att("_FillValue", -999.f))
		return NC_ERR;
	if (!dvdz->add_att("_FillValue", -999.f))
		return NC_ERR;
	if (!dwdz->add_att("_FillValue", -999.f))
		return NC_ERR;
	if (!dtdr->add_att("_FillValue", -999.f))
		return NC_ERR;
	if (!dqdr->add_att("_FillValue", -999.f))
		return NC_ERR;
	if (!dpdr->add_att("_FillValue", -999.f))
		return NC_ERR;
	if (!dtdt->add_att("_FillValue", -999.f))
		return NC_ERR;
	if (!dqdt->add_att("_FillValue", -999.f))
		return NC_ERR;
	if (!dpdt->add_att("_FillValue", -999.f))
		return NC_ERR;
	if (!dtdz->add_att("_FillValue", -999.f))
		return NC_ERR;
	if (!dqdz->add_att("_FillValue", -999.f))
		return NC_ERR;
	if (!dpdz->add_att("_FillValue", -999.f))
		return NC_ERR;
	if (!drhodr->add_att("_FillValue", -999.f))
		return NC_ERR;
	if (!drhodt->add_att("_FillValue", -999.f))
		return NC_ERR;
	if (!drhodz->add_att("_FillValue", -999.f))
		return NC_ERR;
	if (!mcresidual->add_att("_FillValue", -999.f))
		return NC_ERR;

	// Write the coordinate variable data to the file.
	/* real *lons = new real[iDim];
	real *lats = new real[jDim]; */
	real *levs = new real[kDim];
  real *radius = new real[iDim];
  real *thetadeg = new real[jDim];
	int time[2];

	// Reference time and position from center file
	time[0] = configHash->value("ref_time").toInt();

	/* real latReference = configHash->value("ref_lat").toFloat();
	real lonReference = configHash->value("ref_lon").toFloat();
	real refX, refY;
	GeographicLib::TransverseMercatorExact tm = GeographicLib::TransverseMercatorExact::UTM;
	tm.Forward(lonReference, latReference, lonReference, refX, refY); */

	for (int iIndex = 0; iIndex < iDim; iIndex++) {
		real i = (iMin + DI * iIndex);
        radius[iIndex] = i;
		/* real j = (jMin + DJ * (jDim/2))*1000;
		real latnull = 0;
		tm.Reverse(lonReference,refX + i, refY + j, latnull, lons[iIndex]);
        x[iIndex] = i/1000; */
	}

	for (int jIndex = 0; jIndex < jDim; jIndex++) {
		real j = (jMin + DJ * jIndex);
        thetadeg[jIndex] = j;
        /* real i = (iMin + DI * (iDim/2))*1000;
		real lonnull = 0;
		tm.Reverse(lonReference,refX + i, refY + j, lats[jIndex], lonnull);
        y[jIndex] = j/1000; */
	}

	/* if (!lonVar->put(lons, iDim))
		return NC_ERR;

	if (!latVar->put(lats, jDim))
		return NC_ERR; */

    if (!radVar->put(radius, iDim))
		return NC_ERR;

	if (!thetaVar->put(thetadeg, jDim))
		return NC_ERR;

	for (int kIndex = 0; kIndex < kDim; kIndex++) {
		real k = kMin + DK * kIndex;
		levs[kIndex] = k;
	}
	if (!lvlVar->put(levs, kDim))
		return NC_ERR;

	if (!timeVar->put(time, 1))
		return NC_ERR;

	// Write the data.
	for (int rec = 0; rec < 1; rec++)
	{
		if (!u->put_rec(&finalAnalysis[0], rec))
			return NC_ERR;
		if (!v->put_rec(&finalAnalysis[iDim*jDim*kDim*1], rec))
			return NC_ERR;
		if (!w->put_rec(&finalAnalysis[iDim*jDim*kDim*2], rec))
			return NC_ERR;
		if (!wspd->put_rec(&finalAnalysis[iDim*jDim*kDim*3], rec))
			return NC_ERR;
		if (!relhum->put_rec(&finalAnalysis[iDim*jDim*kDim*4], rec))
			return NC_ERR;
		if (!hprime->put_rec(&finalAnalysis[iDim*jDim*kDim*5], rec))
			return NC_ERR;
		if (!qvprime->put_rec(&finalAnalysis[iDim*jDim*kDim*6], rec))
			return NC_ERR;
		if (!rhoprime->put_rec(&finalAnalysis[iDim*jDim*kDim*7], rec))
			return NC_ERR;
		if (!tprime->put_rec(&finalAnalysis[iDim*jDim*kDim*8], rec))
			return NC_ERR;
		if (!pprime->put_rec(&finalAnalysis[iDim*jDim*kDim*9], rec))
			return NC_ERR;
		if (!vorticity->put_rec(&finalAnalysis[iDim*jDim*kDim*10], rec))
			return NC_ERR;
		if (!divergence->put_rec(&finalAnalysis[iDim*jDim*kDim*11], rec))
			return NC_ERR;
		if (!okuboweiss->put_rec(&finalAnalysis[iDim*jDim*kDim*12], rec))
			return NC_ERR;
		if (!strain->put_rec(&finalAnalysis[iDim*jDim*kDim*13], rec))
			return NC_ERR;
		if (!tpw->put_rec(&finalAnalysis[iDim*jDim*kDim*14], rec))
			return NC_ERR;
		if (!rhou->put_rec(&finalAnalysis[iDim*jDim*kDim*15], rec))
			return NC_ERR;
		if (!rhov->put_rec(&finalAnalysis[iDim*jDim*kDim*16], rec))
			return NC_ERR;
		if (!rhow->put_rec(&finalAnalysis[iDim*jDim*kDim*17], rec))
			return NC_ERR;
		if (!rho->put_rec(&finalAnalysis[iDim*jDim*kDim*18], rec))
			return NC_ERR;
		if (!press->put_rec(&finalAnalysis[iDim*jDim*kDim*19], rec))
			return NC_ERR;
		if (!temp->put_rec(&finalAnalysis[iDim*jDim*kDim*20], rec))
			return NC_ERR;
		if (!qv->put_rec(&finalAnalysis[iDim*jDim*kDim*21], rec))
			return NC_ERR;
		if (!h->put_rec(&finalAnalysis[iDim*jDim*kDim*22], rec))
			return NC_ERR;
		if (!qr->put_rec(&finalAnalysis[iDim*jDim*kDim*23], rec))
			return NC_ERR;
    if (!absVorticity->put_rec(&finalAnalysis[iDim*jDim*kDim*24], rec))
			return NC_ERR;
		if (!dewp->put_rec(&finalAnalysis[iDim*jDim*kDim*25], rec))
			return NC_ERR;
		if (!theta->put_rec(&finalAnalysis[iDim*jDim*kDim*26], rec))
			return NC_ERR;
		if (!thetae->put_rec(&finalAnalysis[iDim*jDim*kDim*27], rec))
			return NC_ERR;
		if (!thetaes->put_rec(&finalAnalysis[iDim*jDim*kDim*28], rec))
			return NC_ERR;
    if (!dudr->put_rec(&finalAnalysis[iDim*jDim*kDim*29], rec))
			return NC_ERR;
		if (!dvdr->put_rec(&finalAnalysis[iDim*jDim*kDim*30], rec))
			return NC_ERR;
		if (!dwdr->put_rec(&finalAnalysis[iDim*jDim*kDim*31], rec))
			return NC_ERR;
		if (!dudt->put_rec(&finalAnalysis[iDim*jDim*kDim*32], rec))
			return NC_ERR;
		if (!dvdt->put_rec(&finalAnalysis[iDim*jDim*kDim*33], rec))
			return NC_ERR;
		if (!dwdt->put_rec(&finalAnalysis[iDim*jDim*kDim*34], rec))
			return NC_ERR;
		if (!dudz->put_rec(&finalAnalysis[iDim*jDim*kDim*35], rec))
			return NC_ERR;
		if (!dvdz->put_rec(&finalAnalysis[iDim*jDim*kDim*36], rec))
			return NC_ERR;
		if (!dwdz->put_rec(&finalAnalysis[iDim*jDim*kDim*37], rec))
			return NC_ERR;
		if (!dtdr->put_rec(&finalAnalysis[iDim*jDim*kDim*38], rec))
			return NC_ERR;
		if (!dtdt->put_rec(&finalAnalysis[iDim*jDim*kDim*39], rec))
			return NC_ERR;
		if (!dtdz->put_rec(&finalAnalysis[iDim*jDim*kDim*40], rec))
			return NC_ERR;
		if (!dqdr->put_rec(&finalAnalysis[iDim*jDim*kDim*41], rec))
			return NC_ERR;
		if (!dqdt->put_rec(&finalAnalysis[iDim*jDim*kDim*42], rec))
			return NC_ERR;
		if (!dqdz->put_rec(&finalAnalysis[iDim*jDim*kDim*43], rec))
			return NC_ERR;
		if (!dpdr->put_rec(&finalAnalysis[iDim*jDim*kDim*44], rec))
			return NC_ERR;
		if (!dpdt->put_rec(&finalAnalysis[iDim*jDim*kDim*45], rec))
			return NC_ERR;
		if (!dpdz->put_rec(&finalAnalysis[iDim*jDim*kDim*46], rec))
			return NC_ERR;
		if (!drhodr->put_rec(&finalAnalysis[iDim*jDim*kDim*47], rec))
			return NC_ERR;
		if (!drhodt->put_rec(&finalAnalysis[iDim*jDim*kDim*48], rec))
			return NC_ERR;
		if (!drhodz->put_rec(&finalAnalysis[iDim*jDim*kDim*49], rec))
			return NC_ERR;
		if (!mcresidual->put_rec(&finalAnalysis[iDim*jDim*kDim*50], rec))
			return NC_ERR;
	}

	// The file is automatically closed by the destructor. This frees
	// up any internal netCDF resources associated with the file, and
	// flushes any buffers.
	/* delete[] lats;
	delete[] lons; */
	delete[] levs;
	delete[] radius;
    delete[] thetadeg;
	return true;

}

bool CostFunctionRTZ::writeAsi(const QString& asiFileName)
{
	// Initialize header
	int id[511];
	for (int n = 1; n <= 510; n++) {
		id[n]=-999;
	}

	// Calculate headers
	QStringList fieldNames;
	fieldNames  << "U" << "V" << "W" << "WS" << "RH"<< "HP" << "QP" << "RP" << "TP" << "PP" << "VO" << "DV" << "OW" << "S" << "PW"
	<< "MU" << "MV" << "MW" << "RO" << "PS" << "TK" << "QV" << "HH" << "DZ" << "AV" << "DP" << "TH" << "TE" << "TS";
	id[175] = fieldNames.size();
	for(int n = 0; n < id[175]; n++) {
		QString name_1 = fieldNames.at(n).left(1);
		QString name_2 = fieldNames.at(n).mid(1,1);
		int int_1 = *name_1.toAscii().data();
		int int_2 = *name_2.toAscii().data();
		id[176 + (5 * n)] = (int_1 * 256) + int_2;
		id[177 + (5 * n)] = 8224;
		id[178 + (5 * n)] = 8224;
		id[179 + (5 * n)] = 8224;
		id[180 + (5 * n)] = 1;
	}

	// Polar file
	id[16] = 20559;
	id[17] = 19521;

	/* Lat and Lon
	 id[33] = (int)latReference;
	 id[34] = (int)((latReference - (float)id[33]) * 60.);
	 id[35] = (int)((((latReference - (float)id[33]) * 60.) - (float)id[34]) * 60.) * 100;
	 if (lonReference < 0) {
	 lonReference += 360.;
	 }
	 id[36] = (int)lonReference;
	 id[37] = (int)((lonReference - (float)id[36]) * 60.);
	 id[38] = (int)((((lonReference - (float)id[36]) * 60.) - (float)id[37]) * 60.) * 100; */

	id[33] = 0;
	id[34] = 0;
	id[35] = 0;
	id[36] = 0;
	id[37] = 0;
	id[38] = 0;
	id[40] = 90;

	// Scale factors
	id[68] = 100;
	id[69] = 64;

	// X Header
	id[160] = (int)iMin*100;
	id[161] = (int)iMax*100;
	id[162] = (int)iDim;
	id[163] = (int)(DI * 1000);
	id[164] = 1;

	// Y Header
	id[165] = (int)jMin*100;
	id[166] = (int)(jMax-DJ)*64;
	id[167] = (int)jDim-1;
	id[168] = (int)(DJ * 64);
	id[169] = 2;

	// Z Header
	id[170] = (int)kMin*1000;
	id[171] = (int)kMax*1000;
	id[172] = (int)kDim;
	id[173] = int(DK * 1000);
	id[174] = 3;

	// Number of radars
	id[303] = 1;

	// Index of center
	id[309] = (int)(1);
	id[310] = (int)(1);
	id[311] = 0;

	// Write ascii file for grid2ps
	//Message::toScreen("Trying to write cappi to "+outFileName);
	QFile asiFile(asiFileName);
	if(!asiFile.open(QIODevice::WriteOnly)) {
		cout << "Can't open CAPPI file for writing" << endl;
		return false;
	}

	QTextStream out(&asiFile);

	// Write header
    int line = 0;
	for (int n = 1; n <= 510; n++) {
		line++;
		out << qSetFieldWidth(8) << id[n];
		if (line == 10) {
			out << endl;
            line = 0;
		}
	}

	// Write data
	for(int k = 0; k < kDim; k++) {
		out << reset << "level" << qSetFieldWidth(2) << k+1 << endl;
		for(int j = 0; j < jDim-1; j++) {
			out << reset << "azimuth" << qSetFieldWidth(3) << j+1 << endl;
			for(int n = 0; n < fieldNames.size(); n++) {
				out << reset << left << fieldNames.at(n) << endl;
				int line = 0;
				for (int i = 0; i < iDim;  i++){
					out << reset << qSetRealNumberPrecision(3) << scientific << qSetFieldWidth(10) <<
					finalAnalysis[iDim*jDim*kDim*n +iDim*jDim*k + iDim*j + i];
					line++;
					if (line == 8) {
						out << endl;
						line = 0;
					}
				}
				if (line != 0) {
					out << endl;
				}
			}
		}
	}

	return true;
}
