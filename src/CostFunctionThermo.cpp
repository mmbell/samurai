/*
 *  CostFunctionThermo.cpp
 *  samurai
 *
 *  Copyright 2008 Michael Bell. All rights reserved.
 *
 */

#include "CostFunctionThermo.h"
#include <cmath>
#include <QString>
#include <QStringList>
#include <QTextStream>
#include <QFile>
#include <QDir>
#include <QDateTime>
#include <netcdfcpp.h>
#include <GeographicLib/TransverseMercatorExact.hpp>

CostFunctionThermo::CostFunctionThermo(const int& numObs, const int& stateSize)
	: CostFunction3D(numObs, stateSize)
{   
}

CostFunctionThermo::~CostFunctionThermo()
{
}

void CostFunctionThermo::finalize()
{	    
}

void CostFunctionThermo::initialize(const QHash<QString, QString>* config, real* bgU, real* obs, ReferenceState* ref)
{
	// Initialize number of variables
	configHash = config;

	/* Set the output path */
	outputPath.setPath(configHash->value("output_directory"));
	
	// Horizontal boundary conditions
	iBCL[0] = bcHash.value(configHash->value("i_pip_bcL"));
    iBCR[0] = bcHash.value(configHash->value("i_pip_bcR"));
	iBCL[1] = bcHash.value(configHash->value("i_thetarhop_bcL"));
    iBCR[1] = bcHash.value(configHash->value("i_thetarhop_bcR"));
	iBCL[2] = bcHash.value(configHash->value("i_ftheta_bcL"));
    iBCR[2] = bcHash.value(configHash->value("i_ftheta_bcR"));

	jBCL[0] = bcHash.value(configHash->value("j_pip_bcL"));
    jBCR[0] = bcHash.value(configHash->value("j_pip_bcR"));
	jBCL[1] = bcHash.value(configHash->value("j_thetarhop_bcL"));
    jBCR[1] = bcHash.value(configHash->value("j_thetarhop_bcR"));
	jBCL[2] = bcHash.value(configHash->value("j_ftheta_bcL"));
    jBCR[2] = bcHash.value(configHash->value("j_ftheta_bcR"));

	kBCL[0] = bcHash.value(configHash->value("k_pip_bcL"));
    kBCR[0] = bcHash.value(configHash->value("k_pip_bcR"));
	kBCL[1] = bcHash.value(configHash->value("k_thetarhop_bcL"));
    kBCR[1] = bcHash.value(configHash->value("k_thetarhop_bcR"));
	kBCL[2] = bcHash.value(configHash->value("k_ftheta_bcL"));
    kBCR[2] = bcHash.value(configHash->value("k_ftheta_bcR"));

	// Define the Reference state
    refstate = ref;
	
	// Assign local object pointers
	bgFields = bgU;
	rawObs = obs;

	iMin = configHash->value("i_min").toFloat();
	iMax = configHash->value("i_max").toFloat();
	DI = configHash->value("i_incr").toFloat();
	iDim = (int)((iMax - iMin)/DI) + 1;
	jMin = configHash->value("j_min").toFloat();
	jMax = configHash->value("j_max").toFloat();
	DJ = configHash->value("j_incr").toFloat();
	jDim = (int)((jMax - jMin)/DJ) + 1;
	kMin = configHash->value("k_min").toFloat();
	kMax = configHash->value("k_max").toFloat();
	DK = configHash->value("k_incr").toFloat();
	kDim = (int)((kMax - kMin)/DK) + 1;

	DIrecip = 1./DI;
	DJrecip = 1./DJ;
	DKrecip = 1./DK;
		
	// Adjust the internal, variable domain to include boundaries
	adjustInternalDomain(1);
    
	// Define nodes with internal domain
	int nodes = iDim*jDim*kDim;
    
	// Set up the initial recursive filter
	iFilter = new RecursiveFilter(4);
	jFilter = new RecursiveFilter(4);
	kFilter = new RecursiveFilter(4);
	
	// Allocate memory for the needed arrays
	// These are common to all CostFunctions
	currState = new real[nState];
	currGradient = new real[nState];
	tempState = new real[nState];
	tempGradient = new real[nState];
	xt = new real[nState];
	df = new real[nState];

	// These are local to this one
	CTHTd = new real[nState];
	stateU = new real[nState];
	obsVector = new real[mObs*(obMetaSize+varDim*derivDim)];
	HCq = new real[mObs+nodes];
	innovation = new real[mObs+nodes];	
	bgState = new real[nState];
	bgStdDev = new real[nState];
	stateA = new real[nState];
	stateB = new real[nState];
	stateC = new real[nState];
    
    if (iBCL[0] == PERIODIC) {
        iLDim = iDim-2;   
    } else {
        iLDim = 4;
    }
    if (jBCL[0] == PERIODIC) {
        jLDim = jDim-2;   
    } else {
        jLDim = 4;
    }
    kLDim = 4;
    for (int var = 0; var < varDim; ++var) {
        iRank[var] = iDim - rankHash[iBCL[var]] - rankHash[iBCR[var]];
        if (iBCL[var] == PERIODIC) iRank[var]--;
        jRank[var] = jDim - rankHash[jBCL[var]] - rankHash[jBCR[var]];
        if (jBCL[var] == PERIODIC) jRank[var]--;
        kRank[var] = kDim - rankHash[kBCL[var]] - rankHash[kBCR[var]];
        // Need to verify that Rank is sufficient (>2?)
        iL[var] = new real[iRank[var]*iLDim];
        jL[var] = new real[jRank[var]*jLDim];
        kL[var] = new real[kRank[var]*kLDim];
        iGamma[var] = new real[iRank[var]*iDim];
        jGamma[var] = new real[jRank[var]*jDim];
        kGamma[var] = new real[kRank[var]*kDim];

    }
	// Precalculate the basis functions for lookup table option
	basisappx = configHash->value("spline_approximation").toInt();
	if (basisappx > 0) {
		basis0 = new real[2000000];
		basis1 = new real[2000000];
		fillBasisLookup();
	}
    
    // Initialize the Fourier transforms
    iFFTin = (double*) fftw_malloc(sizeof(double) * iDim);
    iFFTout = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * iDim);
    iForward = fftw_plan_dft_r2c_1d(iDim, iFFTin, iFFTout, FFTW_MEASURE);
    iBackward = fftw_plan_dft_c2r_1d(iDim, iFFTout, iFFTin, FFTW_MEASURE);
    jFFTin = (double*) fftw_malloc(sizeof(double) * jDim);
    jFFTout = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * jDim);
    jForward = fftw_plan_dft_r2c_1d(jDim, jFFTin, jFFTout, FFTW_MEASURE);
    jBackward = fftw_plan_dft_c2r_1d(jDim, jFFTout, jFFTin, FFTW_MEASURE);
    iMaxWavenumber = configHash->value("i_max_wavenumber").toFloat();
    jMaxWavenumber = configHash->value("j_max_wavenumber").toFloat();
}

void CostFunctionThermo::initState(const int iteration)
{

	// Clear the state vector
	cout << "Initializing state vector..." << endl;
	for (int n = 0; n < nState; n++) {
		currState[n] = 0.0;
		currGradient[n] = 0.0;
		tempState[n] = 0.0;
		tempGradient[n] = 0.0;
		xt[n] = 0.0;
		df[n] = 0.0;
		stateA[n] = 0.0;
		stateB[n] = 0.0;
		stateC[n] = 0.0;
	}
	
	
	// Initialize background errors and filter scales
	bgError[0] = configHash->value("bg_pip_error").toFloat();
	bgError[1] = configHash->value("bg_thetarhop_error").toFloat();
	bgError[2] = configHash->value("bg_ftheta_error").toFloat();
	bgError[3] = 0;
	bgError[4] = 0;
	bgError[5] = 0;
	bgError[6] = 0;
	
	// Initialize filter scales
	// Set up the recursive filter
	iFilterScale = configHash->value("i_filter_length").toFloat();
	jFilterScale = configHash->value("j_filter_length").toFloat();
	kFilterScale = configHash->value("k_filter_length").toFloat();		
	iFilter->setFilterLengthScale(iFilterScale);
	jFilter->setFilterLengthScale(jFilterScale);
	kFilter->setFilterLengthScale(kFilterScale);
	
	// Set up the spline matrices
	setupSplines();

	// Flag whether or not to print the subgrid information
	outputMish = configHash->value("output_mish").toInt();
	
	// Mass continuity weight
	mcWeight = configHash->value("mc_weight").toFloat();
	cout << "Mass continuity weight set to " << mcWeight << endl;
	
	if (iteration == 1) {
		cout << "Initializing background..." << endl;
		// Set up the background state
		for (int n = 0; n < nState; n++) {
			bgState[n] = 0.0;
			bgStdDev[n] = 0.0;
		}
		
		// SB Transform on the original bg fields
		SBtransform(bgFields, stateB);
		
		// SA transform = bg B's -> bg A's
		SAtransform(stateB, bgState);
	}
	
	// Using a constant bg error variance for now, but this could be variable across the nodes
	for (int var = 0; var < varDim; var++) {
		for (int iIndex = 0; iIndex < iDim; iIndex++) {
			for (int jIndex = 0; jIndex < jDim; jIndex++) {
				for (int kIndex = 0; kIndex < kDim; kIndex++) {
					int bIndex = varDim*iDim*jDim*kIndex + varDim*iDim*jIndex +varDim*iIndex + var;
					bgStdDev[bIndex] = bgError[var];
				}
			}
		}
	}
	
	// Compute and display the variable BG errors and RMS of values
	for (int var = 0; var < varDim; var++) {
		real varScale = 0;
		for (int iIndex = 0; iIndex < iDim; iIndex++) {
			for (int jIndex = 0; jIndex < jDim; jIndex++) {
				for (int kIndex = 0; kIndex < kDim; kIndex++) {
					int bIndex = varDim*iDim*jDim*kIndex + varDim*iDim*jIndex +varDim*iIndex + var;
					varScale += bgState[bIndex] * bgState[bIndex];
				}
			}
		}
		varScale = sqrt(varScale/(iDim*jDim*kDim));
		if (varScale) {
			real errPct = 100*bgError[var]/varScale;
			cout << "Variable " << var << " RMS = " << varScale << "\t BG Error = " << bgError[var] 
			<< " ( " << errPct << " %)" << endl;
		} else {
			cout << "Variable " << var << " RMS = " << varScale << "\t BG Error = " << bgError[var] 
			<< " ( Infinite! %)" << endl;
		}
	}
	
	// Load the obs locally and weight the nonlinear observation operators by interpolated bg fields
	//obAdjustments();
	obsVector= rawObs;

	// d = y - HXb
	calcInnovation();
	
	// Output the original background field
	outputAnalysis("background", bgState);

	cout << "Beginning analysis...\n";
		
	// HTd
	calcHTranspose(innovation, stateC);
	
	SCtranspose(stateC, stateA);
	
	// S^T (Inverse SA transform) yield B's, put it in the tempState
	SAtranspose(stateA, CTHTd);
				
}	

bool CostFunctionThermo::outputAnalysis(const QString& suffix, real* Astate)
{
	
	cout << "Outputting " << suffix.toStdString() << "...\n";
	// H --> to Mish for output
    QString samuraiout = "samurai_RTZ_" + suffix + ".out";
    ofstream samuraistream;
    if (configHash->value("output_txt") == "true") {
        samuraistream.open(outputPath.absoluteFilePath(samuraiout).toAscii().data());
        samuraistream << "R\tT\tZ\tpip\tthetarhop\tftheta\n";
        samuraistream.precision(10);
    }
    
    int analysisDim = 12;
    int analysisSize = (iDim-2)*(jDim-2)*(kDim-2);
	finalAnalysis = new real[analysisSize*analysisDim];
	
    real Pi = acos(-1.0);
	for (int iIndex = 1; iIndex < iDim-1; iIndex++) {
		for (int ihalf = 0; ihalf <= outputMish; ihalf++) {
			for (int imu = -ihalf; imu <= ihalf; imu++) {
				real i = iMin + DI * (iIndex + (0.5*sqrt(1./3.) * imu + 0.5*ihalf));
                real r = i*1000;
				if (i > ((iDim-1)*DI + iMin)) continue;
				
				for (int jIndex = 1; jIndex < jDim-1; jIndex++) {
					for (int jhalf =0; jhalf <=outputMish; jhalf++) {
						for (int jmu = -jhalf; jmu <= jhalf; jmu++) {
							real j = jMin + DJ * (jIndex + (0.5*sqrt(1./3.) * jmu + 0.5*jhalf));
							if (j > ((jDim-1)*DJ + jMin)) continue;	
														
							for (int kIndex = 1; kIndex < kDim-1; kIndex++) {
								for (int khalf =0; khalf <=outputMish; khalf++) {
									for (int kmu = -khalf; kmu <= khalf; kmu++) {
										real k = kMin + DK * (kIndex + (0.5*sqrt(1./3.) * kmu + 0.5*khalf));
										if (k > ((kDim-1)*DK + kMin)) continue;	
										                                        
										int ii = (int)((i - iMin)*DIrecip);
										int jj = (int)((j - jMin)*DJrecip);
										int kk = (int)((k - kMin)*DKrecip);
										real ibasis = 0.;
										real jbasis = 0.;
										real kbasis = 0.;
										real idbasis = 0.;
										real jdbasis = 0.;
										real kdbasis = 0.;
										real pip = 0.;
										real thetarhop = 0.;
										real ftheta = 0.;
										real pipdr = 0.; real thetarhopdr = 0.; real fthetadr = 0.;
										real pipdt = 0.; real thetarhopdt = 0.; real fthetadt = 0.;
										real pipdz = 0.; real thetarhopdz = 0.; real fthetadz = 0.;
			
                                        for (int var = 0; var < varDim; var++) {
                                            for (int kNode = max(kk-1,0); kNode <= min(kk+2,kDim-1); ++kNode) {
                                                for (int iiNode = (ii-1); iiNode <= (ii+2); ++iiNode) {
                                                    int iNode = iiNode;
                                                    if ((iBCL[var] == PERIODIC) and (iNode < 1)) iNode = iDim-3;
                                                    if ((iBCR[var] == PERIODIC) and (iNode > (iDim-3))) iNode = iiNode - (iDim-3);
                                                    if ((iNode < 0) or (iNode >= iDim)) continue;
                                                    for (int jjNode = (jj-1); jjNode <= (jj+2); ++jjNode) {
                                                        int jNode = jjNode;
                                                        if ((jBCL[var] == PERIODIC) and (jNode < 1)) jNode = jDim-3;
                                                        if ((jBCR[var] == PERIODIC) and (jNode > (jDim-3))) jNode = jjNode - (jDim-3);
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
																pip +=  Astate[aIndex] * basis3x;
																pipdr += Astate[aIndex] * idbasis * jbasis * kbasis;
																pipdt += 180 * Astate[aIndex] * ibasis * jdbasis * kbasis / (i * Pi);
																pipdz += Astate[aIndex] * ibasis * jbasis * kdbasis;
																break;
															case 1:
																thetarhop +=  Astate[aIndex + 1] * basis3x;
																thetarhopdr += Astate[aIndex + 1] * idbasis * jbasis * kbasis;
																thetarhopdt += 180 * Astate[aIndex + 1] * ibasis * jdbasis * kbasis / (i * Pi);
																thetarhopdz += Astate[aIndex + 1] * ibasis * jbasis * kdbasis;
																break;
															case 2:
																ftheta += Astate[aIndex + 2] * basis3x;
																fthetadr += Astate[aIndex + 2] * idbasis * jbasis * kbasis;
																fthetadt += 180 * Astate[aIndex + 2] * ibasis * jdbasis * kbasis / (i * Pi);
																fthetadz += Astate[aIndex + 2] * ibasis * jbasis * kdbasis;
																break;
														}
													}
												}
											}
										}
 										                                        
                                        
										// Avoid singularity at origin
										if (r==0) {
											pip = -999.;
											thetarhop = -999.;
											ftheta = -999.;
										}
										
										
                                        if (configHash->value("output_txt") == "true") {
                                            samuraistream << scientific << i << "\t" << j << "\t"  << k 
                                            << "\t" << pip << "\t" << thetarhop << "\t" << ftheta << "\n";
										}
                                        							
										// On the nodes	
										if (!ihalf and !jhalf and !khalf){
											int fIndex = (iDim-2)*(jDim-2)*(kDim-2); 
											int posIndex = (iDim-2)*(jDim-2)*(kIndex-1) + (iDim-2)*(jIndex-1) + (iIndex-1);
											finalAnalysis[fIndex * 0 + posIndex] = pip;
											finalAnalysis[fIndex * 1 + posIndex] = thetarhop;
											finalAnalysis[fIndex * 2 + posIndex] = ftheta;											
											finalAnalysis[fIndex * 3 + posIndex] = pipdr;
											finalAnalysis[fIndex * 4 + posIndex] = thetarhopdr;
											finalAnalysis[fIndex * 5 + posIndex] = fthetadr;
											finalAnalysis[fIndex * 6 + posIndex] = pipdt;
											finalAnalysis[fIndex * 7 + posIndex] = thetarhopdt;
											finalAnalysis[fIndex * 8 + posIndex] = fthetadt;
											finalAnalysis[fIndex * 9 + posIndex] = pipdz;
											finalAnalysis[fIndex * 10 + posIndex] = thetarhopdz;
											finalAnalysis[fIndex * 11 + posIndex] = fthetadz;
										}
									}
								}
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
        *os++ = "pip";
        *os++ = "thetarhop";
        *os++ = "ftheta";
        *os++ = "Analysis";
        *os++ = "Background";
        qcstream << endl;
        qcstream.precision(10);
        
        ostream_iterator<real> od(qcstream, "\t ");
        for (int m = 0; m < mObs; m++) {
            int mi = m*(obMetaSize+varDim*derivDim);
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
                    int wgt_index = mi + obMetaSize +d*varDim + var;
                    if (!obsVector[wgt_index]) continue;
                    for (int kNode = max(kk-1,0); kNode <= min(kk+2,kDim-1); ++kNode) {
                        kbasis = Basis(kNode, k, kDim-1, kMin, DK, DKrecip, derivative[d][2], kBCL[var], kBCR[var]);
                        for (int iiNode = (ii-1); iiNode <= (ii+2); ++iiNode) {
                            int iNode = iiNode;
                            if ((iBCL[var] == PERIODIC) and (iNode < 1)) iNode = iDim-3;
                            if ((iBCR[var] == PERIODIC) and (iNode > (iDim-3))) iNode = iiNode - (iDim-3);
                            if ((iNode < 0) or (iNode >= iDim)) continue;
                            ibasis = Basis(iNode, i, iDim-1, iMin, DI, DIrecip, derivative[d][1], iBCL[var], iBCR[var]);
                            for (int jjNode = (jj-1); jjNode <= (jj+2); ++jjNode) {
                                int jNode = jjNode;
                                if ((jBCL[var] == PERIODIC) and (jNode < 1)) jNode = jDim-3;
                                if ((jBCR[var] == PERIODIC) and (jNode > (jDim-3))) jNode = jjNode - (jDim-3);
                                if ((jNode < 0) or (jNode >= jDim)) continue;
                                int aIndex = varDim*iDim*jDim*kNode + varDim*iDim*jNode +varDim*iNode;
                                jbasis = Basis(jNode, j, jDim-1, jMin, DJ, DJrecip, derivative[d][0], jBCL[var], jBCR[var]);
                                tempsum += Astate[aIndex + var] * ibasis * jbasis * kbasis * obsVector[wgt_index];
                            }
                        }
                    }
                }
            }				
            for (int t=0; t < obMetaSize-1; t++) {
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
            for (int t=obMetaSize; t<obMetaSize+varDim; t++) {
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

bool CostFunctionThermo::writeNetCDF(const QString& netcdfFileName)
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
	NcVar *pip, *thetarhop, *ftheta;
    NcVar *dpipdr, *dthetarhopdr, *dfthetadr, *dpipdt, *dthetarhopdt, *dfthetadt, *dpipdz, *dthetarhopdz, *dfthetadz;
	
	if (!(pip = dataFile.add_var("PIP", ncFloat, timeDim, 
                               lvlDim, thetaDim, radDim)))
		return NC_ERR;
	if (!(thetarhop = dataFile.add_var("THETARHOP", ncFloat, timeDim, 
							   lvlDim, thetaDim, radDim)))
		return NC_ERR;
	if (!(ftheta = dataFile.add_var("FTHETA", ncFloat, timeDim, 
							   lvlDim, thetaDim, radDim)))

		return NC_ERR;
    if (!(dpipdr = dataFile.add_var("DPIPDR", ncFloat, timeDim, 
                                  lvlDim, thetaDim, radDim)))
		return NC_ERR;
	if (!(dthetarhopdr = dataFile.add_var("DTHETARHOPDR", ncFloat, timeDim, 
                                  lvlDim, thetaDim, radDim)))
		return NC_ERR;
	if (!(dfthetadr = dataFile.add_var("DFTHETADR", ncFloat, timeDim, 
                                  lvlDim, thetaDim, radDim)))
		return NC_ERR;
	if (!(dpipdt = dataFile.add_var("DPIPDT", ncFloat, timeDim, 
                                  lvlDim, thetaDim, radDim)))
		return NC_ERR;
	if (!(dthetarhopdt = dataFile.add_var("DTHETARHOPDT", ncFloat, timeDim, 
                                  lvlDim, thetaDim, radDim)))
		return NC_ERR;
	if (!(dfthetadt = dataFile.add_var("DFTHETADT", ncFloat, timeDim, 
                                  lvlDim, thetaDim, radDim)))
		return NC_ERR;
	if (!(dpipdz = dataFile.add_var("DPIPDZ", ncFloat, timeDim, 
                                  lvlDim, thetaDim, radDim)))
		return NC_ERR;
	if (!(dthetarhopdz = dataFile.add_var("DTHETARHOPDZ", ncFloat, timeDim, 
                                  lvlDim, thetaDim, radDim)))
		return NC_ERR;
	if (!(dfthetadz = dataFile.add_var("DFTHETADZ", ncFloat, timeDim, 
                                  lvlDim, thetaDim, radDim)))
		return NC_ERR;

	
	// Define units attributes for data variables.
	if (!pip->add_att("units", "???"))
		return NC_ERR;
	if (!thetarhop->add_att("units", "???"))
		return NC_ERR;
	if (!ftheta->add_att("units", "???"))
		return NC_ERR;
	if (!dpipdr->add_att("units", "10-5s-1")) 
		return NC_ERR;
	if (!dthetarhopdr->add_att("units", "10-5s-1")) 
		return NC_ERR;
	if (!dfthetadr->add_att("units", "10-5s-1")) 
		return NC_ERR;
	if (!dpipdt->add_att("units", "10-5s-1")) 
		return NC_ERR;
	if (!dthetarhopdt->add_att("units", "10-5s-1")) 
		return NC_ERR;
	if (!dfthetadt->add_att("units", "10-5s-1")) 
		return NC_ERR;
	if (!dpipdz->add_att("units", "10-5s-1")) 
		return NC_ERR;
	if (!dthetarhopdz->add_att("units", "10-5s-1")) 
		return NC_ERR;
	if (!dfthetadz->add_att("units", "10-5s-1")) 
		return NC_ERR;	
	
	// Define long names for data variables.
	if (!pip->add_att("long_name", "pi prime"))
		return NC_ERR;
	if (!thetarhop->add_att("long_name", "theta rho prime"))
		return NC_ERR;
	if (!ftheta->add_att("long_name", "f theta"))
		return NC_ERR;
    if (!dpipdr->add_att("long_name", "pi prime gradient")) 
		return NC_ERR;	
	if (!dthetarhopdr->add_att("long_name", "theta rho prime gradient")) 
		return NC_ERR;	
	if (!dfthetadr->add_att("long_name", "f theta gradient")) 
		return NC_ERR;
	if (!dpipdt->add_att("long_name", "pi prime gradient")) 
		return NC_ERR;	
	if (!dthetarhopdt->add_att("long_name", "theta rho prime gradient")) 
		return NC_ERR;	
	if (!dfthetadt->add_att("long_name", "f theta gradient")) 
		return NC_ERR;
	if (!dpipdz->add_att("long_name", "pi prime gradient")) 
		return NC_ERR;	
	if (!dthetarhopdz->add_att("long_name", "theta rho prime gradient")) 
		return NC_ERR;	
	if (!dfthetadz->add_att("long_name", "f theta gradient")) 
		return NC_ERR;
	
	// Define missing data
	if (!pip->add_att("missing_value", -999.f))
		return NC_ERR;
	if (!thetarhop->add_att("missing_value", -999.f))
		return NC_ERR;
	if (!ftheta->add_att("missing_value", -999.f))
		return NC_ERR;
    if (!dpipdr->add_att("missing_value", -999.f))
		return NC_ERR;	
	if (!dthetarhopdr->add_att("missing_value", -999.f))
		return NC_ERR;	
	if (!dfthetadr->add_att("missing_value", -999.f))
		return NC_ERR;
	if (!dpipdt->add_att("missing_value", -999.f))
		return NC_ERR;	
	if (!dthetarhopdt->add_att("missing_value", -999.f))
		return NC_ERR;	
	if (!dfthetadt->add_att("missing_value", -999.f))
		return NC_ERR;
	if (!dpipdz->add_att("missing_value", -999.f))
		return NC_ERR;	
	if (!dthetarhopdz->add_att("missing_value", -999.f))
		return NC_ERR;	
	if (!dfthetadz->add_att("missing_value", -999.f))
		return NC_ERR;
	
	// Define _Fill_Value for NCL users
	if (!pip->add_att("_FillValue", -999.f))
		return NC_ERR;
	if (!thetarhop->add_att("_FillValue", -999.f))
		return NC_ERR;
	if (!ftheta->add_att("_FillValue", -999.f))
		return NC_ERR;
    if (!dpipdr->add_att("_FillValue", -999.f))
		return NC_ERR;	
	if (!dthetarhopdr->add_att("_FillValue", -999.f))
		return NC_ERR;	
	if (!dfthetadr->add_att("_FillValue", -999.f))
		return NC_ERR;
	if (!dpipdt->add_att("_FillValue", -999.f))
		return NC_ERR;	
	if (!dthetarhopdt->add_att("_FillValue", -999.f))
		return NC_ERR;	
	if (!dfthetadt->add_att("_FillValue", -999.f))
		return NC_ERR;
	if (!dpipdz->add_att("_FillValue", -999.f))
		return NC_ERR;	
	if (!dthetarhopdz->add_att("_FillValue", -999.f))
		return NC_ERR;	
	if (!dfthetadz->add_att("_FillValue", -999.f))
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
		if (!pip->put_rec(&finalAnalysis[0], rec))
			return NC_ERR;
		if (!thetarhop->put_rec(&finalAnalysis[iDim*jDim*kDim*1], rec))
			return NC_ERR;
		if (!ftheta->put_rec(&finalAnalysis[iDim*jDim*kDim*2], rec))
			return NC_ERR;
        if (!dpipdr->put_rec(&finalAnalysis[iDim*jDim*kDim*3], rec)) 
			return NC_ERR;
		if (!dthetarhopdr->put_rec(&finalAnalysis[iDim*jDim*kDim*4], rec)) 
			return NC_ERR;
		if (!dfthetadr->put_rec(&finalAnalysis[iDim*jDim*kDim*5], rec)) 
			return NC_ERR;
		if (!dpipdt->put_rec(&finalAnalysis[iDim*jDim*kDim*6], rec)) 
			return NC_ERR;
		if (!dthetarhopdt->put_rec(&finalAnalysis[iDim*jDim*kDim*7], rec)) 
			return NC_ERR;
		if (!dfthetadt->put_rec(&finalAnalysis[iDim*jDim*kDim*8], rec)) 
			return NC_ERR;
		if (!dpipdz->put_rec(&finalAnalysis[iDim*jDim*kDim*9], rec)) 
			return NC_ERR;
		if (!dthetarhopdz->put_rec(&finalAnalysis[iDim*jDim*kDim*10], rec)) 
			return NC_ERR;
		if (!dfthetadz->put_rec(&finalAnalysis[iDim*jDim*kDim*11], rec)) 
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

bool CostFunctionThermo::writeAsi(const QString& asiFileName)
{
	cout << "Writing Asi not implemented yet" << endl;
	
	return true;
}	

