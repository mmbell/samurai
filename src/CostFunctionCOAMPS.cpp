/*
 *  CostFunctionCOAMPS.cpp
 *  samurai
 *
 *  Copyright 2015 Michael Bell. All rights reserved.
 *
 */

#include "CostFunctionCOAMPS.h"
#include <cmath>
#include "datetime.h"
//#include <netcdfcpp.h>
#include <Ncxx/Nc3File.hh>

// #include <euclid/GeographicLib/LambertConformalConic.hpp>

CostFunctionCOAMPS::CostFunctionCOAMPS(const Projection& proj, const int& numObs, const int& stateSize)
  : CostFunction3D(proj, numObs, stateSize)
{
}

CostFunctionCOAMPS::~CostFunctionCOAMPS()
{
}

void CostFunctionCOAMPS::setSigmas(float *sigmas, int size) {
  sigmaTable = sigmas;
  sDim = size;
}

bool CostFunctionCOAMPS::outputAnalysis_thermo(const std::string& suffix, real* Astate)
{
}
bool CostFunctionCOAMPS::outputAnalysis(const std::string& suffix, real* Astate)
{
/* NCAR FIXME
  cout << "Outputting " << suffix << "...\n";
  // H --> to Mish for output
  std::string samuraiout = "samurai_COAMPS_" + suffix + ".out";
  ofstream samuraistream;
  if ((*configHash)["output_txt"] == "true") {
    samuraistream.open(outputPath + "/" + samuraiout);
    samuraistream << "X\tY\tZ\tu\tv\tw\tVorticity\tDivergence\tqv\trho\tT\tP\tTheta\tTheta_e\tTheta_es\t";
    samuraistream << "udx\tudy\tudz\tvdx\tvdy\tvdz\twdx\twdy\twdz\trhowdz\tMC residual\tdBZ\n";
    samuraistream.precision(10);
  }

  // sigmaTable and sDim might have been set with a call to setSigmas
  // if not, use default values

  float defaultSigmas[] { 31885.0, 29385.0, 25085.0, 21985.0, 19685.0, 17885.0,
      16385.0, 15072.5, 13910.0, 12860.0, 11885.0, 10955.0, 10065.0,
      9215.0, 8405.0, 7635.0, 6905.0, 6215.0, 5565.0, 4955.0, 4385.0,
      3855.0, 3365.0, 2915.0, 2505.0, 2135.0, 1805.0, 1515.0, 1265.0,
      1050.0, 860.0, 690.0, 540.0, 410.0, 300.0, 210.0, 140.0, 90.0,
      55.0, 30.0, 10.0, 0.0 };
  
  float *sigma = sigmaTable;
  
  if (sigma == NULL) {
    sDim = 42;
    sigma = defaultSigmas;
  }
  
  int analysisDim = 51;

  int analysisSize = (iDim - 2) * (jDim - 2) * (sDim - 2);
  finalAnalysis = new real[analysisSize * analysisDim];
  if (finalAnalysis == NULL) {
    std::cerr << "CostFunctionCOAMPS::outputAnalysis: Failed to allocate finalAnalysis["
	      << analysisSize * analysisDim << "]. " << std::endl;
    return false;
  }
  real gausspoint = 0.5 * sqrt(1. / 3.);


  for (int iIndex = 1; iIndex < iDim - 1; iIndex++) {
    for (int ihalf = 0; ihalf <= mishFlag; ihalf++) {
      for (int imu = -ihalf; imu <= ihalf; imu++) {
	real i = iMin + DI * (iIndex + (gausspoint * imu + 0.5 * ihalf));
	if (i > ((iDim - 1) * DI + iMin))
	  continue;

	for (int jIndex = 1; jIndex < jDim - 1; jIndex++) {
	  for (int jhalf =0; jhalf <= mishFlag; jhalf++) {
	    for (int jmu = -jhalf; jmu <= jhalf; jmu++) {
	      real j = jMin + DJ * (jIndex + (gausspoint * jmu + 0.5 * jhalf));
	      if (j > ((jDim - 1) * DJ + jMin))
		continue;

	      real tpw = 0;
	      for (int kIndex = 1; kIndex < sDim - 1; kIndex++) {
		real k = (sigma[kIndex] * (sigma[0] - sigma[41])/sigma[0] + sigma[41]) / 1000.0;
		real imu, jmu, kmu, ihalf, jhalf, khalf;
		imu = jmu = kmu = ihalf = jhalf = khalf = 0.0;
		if (k < kMin) continue;
		if (k > ((kDim - 1) * DK + kMin)) continue;

		real heightm = 1000 * k;
		real rhoBar = refstate->getReferenceVariable(ReferenceVariable::rhoaref, heightm);
		real qBar = refstate->getReferenceVariable(ReferenceVariable::qvbhypref, heightm);
		real tBar = refstate->getReferenceVariable(ReferenceVariable::tempref, heightm);

		int ii = (int)((i - iMin) * DIrecip);
		int jj = (int)((j - jMin) * DJrecip);
		int kk = (int)((k - kMin) * DKrecip);
		real ibasis = 0.;
		real jbasis = 0.;
		real kbasis = 0.;
		real idbasis = 0.;
		real jdbasis = 0.;
		real kdbasis = 0.;
		real rhov = 0.;
		real rhou = 0.;
		real rhow = 0.;
		real rhovdx = 0.; real rhoudx = 0.; real rhowdx = 0.;
		real rhovdy = 0.; real rhoudy = 0.; real rhowdy = 0.;
		real rhovdz = 0.; real rhoudz = 0.; real rhowdz = 0.;
		real tprime = 0.; real tdx = 0.; real tdy = 0.; real tdz = 0.;
		real rhoadx = 0.; real rhoady = 0.; real rhoadz = 0.;
		real qvdx = 0.; real qvdy = 0.; real qvdz = 0.;
		real pdx = 0.; real pdy = 0.; real pdz = 0.;
		real qvprime = 0.;
		real rhoprime = 0.;
		real qrprime = 0.;
		for (int var = 0; var < varDim; var++) {
		  for (int kkNode = (kk - 1); kkNode <= (kk + 2); ++kkNode) {
		    int kNode = kkNode;
		    if ((kNode < 0) or (kNode >= kDim)) continue;
		    for (int iiNode = (ii - 1); iiNode <= (ii + 2); ++iiNode) {
		      int iNode = iiNode;
		      if ((iNode < 0) or (iNode >= iDim)) continue;
		      for (int jjNode = (jj - 1); jjNode <= (jj + 2); ++jjNode) {
			int jNode = jjNode;
			if ((jNode < 0) or (jNode >= jDim)) continue;
			ibasis = Basis(iNode, i, iDim-1, iMin, DI, DIrecip, 0, iBCL[var], iBCR[var]);
			jbasis = Basis(jNode, j, jDim-1, jMin, DJ, DJrecip, 0, jBCL[var], jBCR[var]);
			kbasis = Basis(kNode, k, kDim-1, kMin, DK, DKrecip, 0, kBCL[var], kBCR[var]);
			idbasis = Basis(iNode, i, iDim-1, iMin, DI, DIrecip, 1, iBCL[var], iBCR[var]);
			jdbasis = Basis(jNode, j, jDim-1, jMin, DJ, DJrecip, 1, jBCL[var], jBCR[var]);
			kdbasis = Basis(kNode, k, kDim-1, kMin, DK, DKrecip, 1, kBCL[var], kBCR[var]);
			real basis3x = ibasis*jbasis*kbasis;
			int aIndex = varDim * iDim * jDim * kNode + varDim * iDim * jNode + varDim * iNode;
			switch (var) {
			case 0:
			  rhou   += Astate[aIndex] * basis3x;
			  rhoudx += Astate[aIndex] * idbasis * jbasis * kbasis;
			  rhoudy += Astate[aIndex] * ibasis * jdbasis * kbasis;
			  rhoudz += Astate[aIndex] * ibasis * jbasis * kdbasis;
			  break;
			case 1:
			  rhov   += Astate[aIndex + 1] * basis3x;
			  rhovdx += Astate[aIndex + 1] * idbasis * jbasis * kbasis;
			  rhovdy += Astate[aIndex + 1] * ibasis * jdbasis * kbasis;
			  rhovdz += Astate[aIndex + 1] * ibasis * jbasis * kdbasis;
			  break;
			case 2:
			  rhow   += Astate[aIndex + 2] * basis3x;
			  rhowdx += Astate[aIndex + 2] * idbasis * jbasis * kbasis;
			  rhowdy += Astate[aIndex + 2] * ibasis * jdbasis * kbasis;
			  rhowdz += Astate[aIndex + 2] * ibasis * jbasis * kdbasis;
			  break;
			case 3:
			  tprime += Astate[aIndex + 3] * basis3x;
			  tdx += Astate[aIndex + 3] * idbasis * jbasis * kbasis;
			  tdy += Astate[aIndex + 3] * ibasis * jdbasis * kbasis;
			  tdz += Astate[aIndex + 3] * ibasis * jbasis * kdbasis;
			  break;
			case 4:
			  qvprime += Astate[aIndex + 4] * basis3x;
			  qvdx += Astate[aIndex + 4] * idbasis * jbasis * kbasis;
			  qvdy += Astate[aIndex + 4] * ibasis * jdbasis * kbasis;
			  qvdz += Astate[aIndex + 4] * ibasis * jbasis * kdbasis;
			  break;
			case 5:
			  rhoprime += Astate[aIndex + 5] * basis3x;
			  rhoadx += Astate[aIndex + 5] * idbasis * jbasis * kbasis;
			  rhoady += Astate[aIndex + 5] * ibasis * jdbasis * kbasis;
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
		std::string gridref = (*configHash)["qr_variable"];
		if ((imu != 0) and (jmu != 0) and (kmu != 0)) {
		  int uJ = jIndex * 2 + (jmu+1)/2;
		  int uI = iIndex * 2 + (imu+1)/2;
		  int uK = kIndex * 2 + (kmu+1)/2;
		  int uIndex = varDim * (iDim - 1) * 2 * (jDim - 1) * 2 * uK
		    + varDim * (iDim - 1) * 2 * uJ + varDim * uI;

		  bgFields[uIndex] = rhou;
		  bgFields[uIndex + 1] = rhov;
		  bgFields[uIndex + 2] = rhow;
		  bgFields[uIndex + 3] = tprime;
		  bgFields[uIndex + 4] = qvprime;
		  bgFields[uIndex + 5] = rhoprime;
		  bgFields[uIndex + 6] = qrprime;
		}

		if (((*configHash)["output_mish"] == "false")
		    and (ihalf or jhalf or khalf)) continue;

		// Output it
		real rhoa = rhoBar + rhoprime / 100;
		real qv = refstate->bhypInvTransform(qBar + qvprime);
		real qbardz = 1000. * refstate->getReferenceVariable(ReferenceVariable::qvbhypref, heightm, 1);
		// qv derivatives multipled by 2 to account for hyperbolic transform
		qvdx = 2.0*qvdx;
		qvdy = 2.0*qvdy;
		qvdz = 2.0*(qbardz + qvdz);

		real qr;
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

		real h = 1005.7 * temp + 2.501e3 * qv + 9.81 * heightm;
		real airpress = temp * rhoa * 287. / 100.;
		//real tempc = temp - 273.15;
		//real satvp = 6.112 * exp((17.67 * tempc)/(243.5 + tempc));
		real satvp =  exp(-6096.9385 / temp + 16.635794 - 2.711193e-2 * temp
				  + 1.673952e-5 * temp * temp + 2.433502 * log(temp));
		real vp = temp * rhoq * 461. / 100.;
		//real vp = airpress * qv / (622 + qv);
		real press = airpress + vp;

		real pprime = press - refstate->getReferenceVariable(ReferenceVariable::pressref, heightm) / 100.;
		real hprime = h - refstate->getReferenceVariable(ReferenceVariable::href, heightm);

		real RoverCp = 0.2854 * (1 - 0.00028 * qv);
		real theta = temp * pow((1000 / press), RoverCp);
		real lcl = 2840/(3.5 * log(temp) - log(vp) - 4.805) + 55.0;
		real thetae = theta * exp(((3.376 / lcl) - 0.00254) * qv * (1 + 0.00081 * qv));
		real qvsat = 622 * satvp / airpress;
		real relhum = -999.;
		real thetaes = -999.;
		if (satvp != 0) {
		  relhum = 100 * vp / satvp;
		  lcl = 2840 / (3.5 * log(temp) - log(satvp) - 4.805) + 55.0;
		  thetaes = theta * exp(((3.376 / lcl) - 0.00254) * qvsat * (1 + 0.00081 * qvsat));
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
		  dewp = 237.3 * log(vp / 6.1078) / (17.2694 - log(vp / 6.1078)) + 273.15;
		}
		// Calculate the kinematic derivatives
		// rhoa derivatives divided by 100
		rhoadx /= 100.;
		rhoady /= 100.;
		rhoadz /= 100.;
		real rhodx = rhoadx * (1. + qv/1000.) + rhoa * qvdx/1000.;
		real rhody = rhoady * (1. + qv/1000.) + rhoa * qvdy/1000.;
		real rhodz = rhoadz * (1. + qv/1000.) + rhoa * qvdz/1000.;
		real rhobardz = 1000 * refstate->getReferenceVariable(ReferenceVariable::rhoref, heightm, 1);
		rhodz += rhobardz;
		real rhoabardz = 1000 * refstate->getReferenceVariable(ReferenceVariable::rhoaref, heightm, 1);
		rhoadz += rhoabardz;

		// Units 10-5
		real udx = 100. * (rhoudx - u * rhodx) / rho;
		real udy = 100. * (rhoudy - u * rhody) / rho;
		real udz = 100. * (rhoudz - u * rhodz) / rho;

		real vdx = 100. * (rhovdx - v * rhodx) / rho;
		real vdy = 100. * (rhovdy - v * rhody) / rho;
		real vdz = 100. * (rhovdz - v * rhodz) / rho;

		real wdx = 100. * (rhowdx - w * rhodx) / rho;
		real wdy = 100. * (rhowdy - w * rhody) / rho;
		real wdz = 100. * (rhowdz - w * rhodz) / rho;

		// Thermodynamic derivatives
		pdx = (tdx*rhoa + rhoadx*temp)*287./100. + (tdx*rhoq + (rhoadx*qv + qvdx*rhoa)*temp/1000.0)*461./100.;
		pdy = (tdy*rhoa + rhoady*temp)*287./100. + (tdy*rhoq + (rhoady*qv + qvdy*rhoa)*temp/1000.0)*461./100.;
		pdz = (tdz*rhoa + rhoadz*temp)*287./100. + (tdz*rhoq + (rhoadz*qv + qvdz*rhoa)*temp/1000.0)*461./100.;

		// Vorticity units are 10-5
		real vorticity = (vdx - udy);
		real divergence = (udx + vdy);
		real s1 = (udx - vdy);
		real s2 = (vdx + udy);
		real strain = sqrt(s1*s1 + s2*s2);
		real okuboweiss = vorticity*vorticity - s1*s1 -s2*s2;
		real mcresidual = rhoudx + rhovdy + rhowdz;

		// Add Coriolis parameter to relative vorticity
		real latReference = configHash->value("ref_lat").toFloat();
		real Coriolisf = 2 * 7.2921 * sin(latReference*acos(-1.) / 180); // Units 10^-5 s-1
		real absVorticity = vorticity + Coriolisf;

		std::string refmask = (*configHash)["mask_reflectivity"];
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
		    udx = -999.; udy = -999.; udz = -999.;
		    vdx = -999.; vdy = -999.; vdz = -999.;
		    wdx = -999.; wdy = -999.; wdz = -999.;
		    tdx = -999.; tdy = -999.; tdz = -999.;
		    qvdx = -999.; qvdy = -999.; qvdz = -999.;
		    pdx = -999.; pdy = -999.; pdz = -999.;
		    rhodx = -999.; rhody = -999.; rhodz = -999.;
		    dewp = -999.;
		    theta = -999.; thetae = -999.; thetaes = -999.;
		  }
		}

		if ((*configHash)["output_txt"] == "true") {
		  samuraistream << scientific << i << "\t" << j << "\t"  << k
				<< "\t" << u << "\t" << v << "\t" << w << "\t" << vorticity << "\t" << divergence
				<< "\t" << qv << "\t" << rho << "\t" << temp << "\t" << press
				<< "\t" << theta << "\t" << thetae << "\t" << thetaes << "\t"
				<< udx << "\t" << udy << "\t" << udz << "\t"
				<< vdx << "\t" << vdy << "\t" << vdz << "\t"
				<< wdx << "\t" << wdy << "\t" << wdz << "\t"
				<< rhowdz * 100. << "\t" << mcresidual << "\t" << qr << "\n";
		}

		// Sum up the TPW in the vertical, top level is tpw
		tpw += qv * rhoa * DK;

		// On the nodes
		if (!ihalf and !jhalf and !khalf) {
		  int fIndex = (iDim - 2) * (jDim - 2) * (sDim - 2);
		  int posIndex = (iDim - 2) * (jDim - 2) * (kIndex - 1) + (iDim - 2) * (jIndex - 1) + (iIndex - 1);

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
		  finalAnalysis[fIndex * 29 + posIndex] = udx;
		  finalAnalysis[fIndex * 30 + posIndex] = vdx;
		  finalAnalysis[fIndex * 31 + posIndex] = wdx;
		  finalAnalysis[fIndex * 32 + posIndex] = udy;
		  finalAnalysis[fIndex * 33 + posIndex] = vdy;
		  finalAnalysis[fIndex * 34 + posIndex] = wdy;
		  finalAnalysis[fIndex * 35 + posIndex] = udz;
		  finalAnalysis[fIndex * 36 + posIndex] = vdz;
		  finalAnalysis[fIndex * 37 + posIndex] = wdz;
		  finalAnalysis[fIndex * 38 + posIndex] = tdx;
		  finalAnalysis[fIndex * 39 + posIndex] = tdy;
		  finalAnalysis[fIndex * 40 + posIndex] = tdz;
		  finalAnalysis[fIndex * 41 + posIndex] = qvdx;
		  finalAnalysis[fIndex * 42 + posIndex] = qvdy;
		  finalAnalysis[fIndex * 43 + posIndex] = qvdz;
		  finalAnalysis[fIndex * 44 + posIndex] = pdx;
		  finalAnalysis[fIndex * 45 + posIndex] = pdy;
		  finalAnalysis[fIndex * 46 + posIndex] = pdz;
		  finalAnalysis[fIndex * 47 + posIndex] = rhodx;
		  finalAnalysis[fIndex * 48 + posIndex] = rhody;
		  finalAnalysis[fIndex * 49 + posIndex] = rhodz;
		  finalAnalysis[fIndex * 50 + posIndex] = mcresidual;
		}
	      }
	    }
	  }
	}
      }
    }
  }


  std::string fileName = "samurai_COAMPS_" + suffix;
  std::string outFileName = outputPath+ "/" + fileName;

  // Write the Obs to a summary text file
  if ((*configHash)["output_qc"] == "true") {
    std::string qcout = "samurai_QC_" + suffix + ".out";
    std::string qcFileName = outputPath + "/" + qcout;
    ofstream qcstream(qcFileName);
    ostream_iterator<string> os(qcstream, "\t ");
    *os++ = "Observation";
    *os++ = "Inverse Error";
    *os++ = "X";
    *os++ = "Y";
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
      real i = obsVector[mi + 2];
      real j = obsVector[mi + 3];
      real k = obsVector[mi + 4];
      real tempsum = 0;
      int ii = (int)((i - iMin) * DIrecip);
      int jj = (int)((j - jMin) * DJrecip);
      int kk = (int)((k - kMin) * DKrecip);
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
      datetime obtime;
      obtime.setTime_t(unixtime);
      //obtime.setTimeSpec(Qt::UTC); 
      //std::string timestring = obtime.toString("hh:mm:ss.zzz");
      qcstream << obTtime << "\t";

      // Multiply the weight by the ob -- Observations.in has individual weights already
      // Only non-derivative for now
      for (int t=7; t<14; t++) {
	*od++ = obsVector[mi+t] * obsVector[mi];
      }

      *od++ = tempsum;
      *od++ = obsVector[mi]-innovation[m];
      qcstream << endl;

    }
  }

  // Write out COAMPS flatfiles
  std::string flatFileName = outFileName + "_uuwind_sig_029385";
  writeFlatfile(flatFileName, 0);

  flatFileName = outFileName + "_vvwind_sig_029385";
  writeFlatfile(flatFileName, 1);

  adjustInternalDomain(-1);
  sDim -= 2;

  // Write out to a netCDF file
  if ((*configHash)["output_netcdf"] == "true") {
    std::string cdfFileName = outFileName + ".nc";
    if ( ! writeNetCDF(cdfFileName))
      std::cout << "Error writing netcdf file " << cdfFileName << std::endl;
  }
  // Write out to an asi file
  if ((*configHash)["output_asi"] == "true") {
    std::string asiFileName = outFileName + ".asi";
    if (!writeAsi(asiFileName))
      std::cout << "Error writing asi file " << asiFileName << std::endl;
  }

  // Set the domain back
  adjustInternalDomain(1);
  sDim += 2;

  // Free the memory for the analysis variables


  if (suffix != "analysis") 	// TODO: Need the keep analysis for return arrays
    delete[] finalAnalysis;
*/
  return true;

}

bool CostFunctionCOAMPS::writeFlatfile(const std::string& flatFileName, const int var)
{
/*fixme
	// Write out a binary file for use with COAMPS
	std::ofstream file(flatFileName);
	if (! file.is_open()) {
	  std::cerr << "CostFunctionCOAMPS::writeFlatfile: Failed to open " << flatFileName << std::endl;
	  return false;
	}

	for (int kIndex = 1; kIndex < sDim-1; kIndex++) {
		for (int jIndex = 1; jIndex < jDim-1; jIndex++) {
			for (int iIndex = 1; iIndex < iDim-1; iIndex++) {
				int fIndex = (iDim-2)*(jDim-2)*(sDim-2);
				int posIndex = (iDim-2)*(jDim-2)*(kIndex-1) + (iDim-2)*(jIndex-1) + (iIndex-1);
				float value = (float)finalAnalysis[fIndex * var + posIndex];
				out << value;
			}
		}
	}
	file.close();
*/
	return true;

}

bool CostFunctionCOAMPS::writeNetCDF(const std::string& netcdfFileName)
{
/* fixme
	Nc3Error err(Nc3Error::verbose_nonfatal);
	int NC_ERR = 0;

	// Create the file.
	std::cout << "Writing " << netcdfFileName << std::endl;
	Nc3File dataFile(netcdfFileName, Nc3File::Replace);

	// Check to see if the file was created.
	if(!dataFile.is_valid())
		return NC_ERR;

	// Define the dimensions. NetCDF will hand back an ncDim object for
	// each.
	Nc3Dim *lvlDim, *latDim, *lonDim, *timeDim;
	if (!(lonDim = dataFile.add_dim("longitude", iDim)))
		return NC_ERR;
	if (!(latDim = dataFile.add_dim("latitude", jDim)))
		return NC_ERR;
	if (!(lvlDim = dataFile.add_dim("altitude", sDim)))
		return NC_ERR;
	// Add an unlimited dimension...
	if (!(timeDim = dataFile.add_dim("time")))
		return NC_ERR;

	// Define the coordinate variables.
	Nc3Var *latVar, *lonVar, *lvlVar, *timeVar, *xVar, *yVar;
	if (!(lonVar = dataFile.add_var("longitude", nc3Float, lonDim)))
		return NC_ERR;
	if (!(latVar = dataFile.add_var("latitude", nc3Float, latDim)))
		return NC_ERR;
	if (!(xVar = dataFile.add_var("x", nc3Float, lonDim)))
		return NC_ERR;
	if (!(yVar = dataFile.add_var("y", nc3Float, latDim)))
		return NC_ERR;
	if (!(lvlVar = dataFile.add_var("altitude", nc3Float, lvlDim)))
		return NC_ERR;
	if (!(timeVar = dataFile.add_var("time", nc3Int, timeDim)))
		return NC_ERR;

	// Define units attributes for coordinate vars. This attaches a
	// text attribute to each of the coordinate variables, containing
	// the units.
	if (!latVar->add_att("units", "degrees_north"))
		return NC_ERR;
	if (!lonVar->add_att("units", "degrees_east"))
		return NC_ERR;
	if (!xVar->add_att("units", "km"))
		return NC_ERR;
	if (!yVar->add_att("units", "km"))
		return NC_ERR;
	if (!lvlVar->add_att("units", "km"))
		return NC_ERR;
	if (!timeVar->add_att("units", "seconds since 1970-01-01 00:00:00 +0000"))
		return NC_ERR;

	// Define the netCDF variables
	Nc3Var *u, *v, *w, *wspd, *relhum, *hprime, *qvprime, *rhoprime, *tprime, *pprime;
	Nc3Var *vorticity, *divergence, *okuboweiss, *strain, *tpw, *rhou, *rhov, *rhow;
	Nc3Var *rho, *press, *temp, *qv, *h, *qr, *absVorticity;
	Nc3Var *dudx, *dvdx, *dwdx, *dudy, *dvdy, *dwdy, *dudz, *dvdz, *dwdz;
	Nc3Var *dtdx, *dqdx, *dpdx, *dtdy, *dqdy, *dpdy, *dtdz, *dqdz, *dpdz;
	Nc3Var *drhodx, *drhody, *drhodz;
	Nc3Var *dewp, *theta, *thetae, *thetaes, *mcresidual;

	if (!(u = dataFile.add_var("U", nc3Float, timeDim,
                               lvlDim, latDim, lonDim)))
		return NC_ERR;
	if (!(v = dataFile.add_var("V", nc3Float, timeDim,
				   lvlDim, latDim, lonDim)))
		return NC_ERR;
	if (!(w = dataFile.add_var("W", nc3Float, timeDim,
				   lvlDim, latDim, lonDim)))
		return NC_ERR;
	if (!(wspd = dataFile.add_var("WSPD", nc3Float, timeDim,
				      lvlDim, latDim, lonDim)))
	  return NC_ERR;
	if (!(relhum = dataFile.add_var("RH", nc3Float, timeDim,
					lvlDim, latDim, lonDim)))
		return NC_ERR;
	if (!(hprime = dataFile.add_var("HP", nc3Float, timeDim,
                                    lvlDim, latDim, lonDim)))
		return NC_ERR;
	if (!(qvprime = dataFile.add_var("QVP", nc3Float, timeDim,
                                     lvlDim, latDim, lonDim)))
		return NC_ERR;
	if (!(rhoprime = dataFile.add_var("RHOAP", nc3Float, timeDim,
                                      lvlDim, latDim, lonDim)))
		return NC_ERR;
	if (!(tprime = dataFile.add_var("TP", nc3Float, timeDim,
									lvlDim, latDim, lonDim)))
		return NC_ERR;
	if (!(pprime = dataFile.add_var("PP", nc3Float, timeDim,
									lvlDim, latDim, lonDim)))
		return NC_ERR;
	if (!(vorticity = dataFile.add_var("VORT", nc3Float, timeDim,
                                       lvlDim, latDim, lonDim)))
		return NC_ERR;
	if (!(divergence = dataFile.add_var("DIV", nc3Float, timeDim,
                                        lvlDim, latDim, lonDim)))
		return NC_ERR;
	if (!(okuboweiss = dataFile.add_var("OW", nc3Float, timeDim,
                                        lvlDim, latDim, lonDim)))
		return NC_ERR;
	if (!(strain = dataFile.add_var("STRAIN", nc3Float, timeDim,
                                    lvlDim, latDim, lonDim)))
		return NC_ERR;
	if (!(tpw = dataFile.add_var("TPW", nc3Float, timeDim,
                                 lvlDim, latDim, lonDim)))
		return NC_ERR;
	if (!(rhou = dataFile.add_var("RHOU", nc3Float, timeDim,
                                  lvlDim, latDim, lonDim)))
		return NC_ERR;
	if (!(rhov = dataFile.add_var("RHOV", nc3Float, timeDim,
                                  lvlDim, latDim, lonDim)))
		return NC_ERR;
	if (!(rhow = dataFile.add_var("RHOW", nc3Float, timeDim,
                                  lvlDim, latDim, lonDim)))
		return NC_ERR;
	if (!(rho = dataFile.add_var("RHOA", nc3Float, timeDim,
                                 lvlDim, latDim, lonDim)))
		return NC_ERR;
	if (!(press = dataFile.add_var("P", nc3Float, timeDim,
                                   lvlDim, latDim, lonDim)))
		return NC_ERR;
	if (!(temp = dataFile.add_var("T", nc3Float, timeDim,
								  lvlDim, latDim, lonDim)))
		return NC_ERR;
	if (!(qv = dataFile.add_var("QV", nc3Float, timeDim,
								lvlDim, latDim, lonDim)))
		return NC_ERR;
	if (!(h = dataFile.add_var("H", nc3Float, timeDim,
                               lvlDim, latDim, lonDim)))
		return NC_ERR;
	if ((*configHash)["qr_variable"] == "dbz") {
	  if (!(qr = dataFile.add_var("DBZ", nc3Float, timeDim,
                                    lvlDim, latDim, lonDim)))
            return NC_ERR;
	} else {
        if (!(qr = dataFile.add_var("QR", nc3Float, timeDim,
                                    lvlDim, latDim, lonDim)))
            return NC_ERR;
	}
	if (!(absVorticity = dataFile.add_var("ABSVORT", nc3Float, timeDim,
                                          lvlDim, latDim, lonDim)))
		return NC_ERR;
	if (!(dewp = dataFile.add_var("DEWPOINT", nc3Float, timeDim,
								lvlDim, latDim, lonDim)))
		return NC_ERR;
	if (!(theta = dataFile.add_var("THETA", nc3Float, timeDim,
								lvlDim, latDim, lonDim)))
		return NC_ERR;
	if (!(thetae = dataFile.add_var("THETAE", nc3Float, timeDim,
								lvlDim, latDim, lonDim)))
		return NC_ERR;
	if (!(thetaes = dataFile.add_var("THETAES", nc3Float, timeDim,
								lvlDim, latDim, lonDim)))
		return NC_ERR;
	if (!(dudx = dataFile.add_var("DUDX", nc3Float, timeDim,
                                  lvlDim, latDim, lonDim)))
		return NC_ERR;
	if (!(dvdx = dataFile.add_var("DVDX", nc3Float, timeDim,
                                  lvlDim, latDim, lonDim)))
		return NC_ERR;
	if (!(dwdx = dataFile.add_var("DWDX", nc3Float, timeDim,
                                  lvlDim, latDim, lonDim)))
		return NC_ERR;
	if (!(dudy = dataFile.add_var("DUDY", nc3Float, timeDim,
                                  lvlDim, latDim, lonDim)))
		return NC_ERR;
	if (!(dvdy = dataFile.add_var("DVDY", nc3Float, timeDim,
                                  lvlDim, latDim, lonDim)))
		return NC_ERR;
	if (!(dwdy = dataFile.add_var("DWDY", nc3Float, timeDim,
                                  lvlDim, latDim, lonDim)))
		return NC_ERR;
	if (!(dudz = dataFile.add_var("DUDZ", nc3Float, timeDim,
                                  lvlDim, latDim, lonDim)))
		return NC_ERR;
	if (!(dvdz = dataFile.add_var("DVDZ", nc3Float, timeDim,
                                  lvlDim, latDim, lonDim)))
		return NC_ERR;
	if (!(dwdz = dataFile.add_var("DWDZ", nc3Float, timeDim,
                                  lvlDim, latDim, lonDim)))
		return NC_ERR;
	if (!(dtdx = dataFile.add_var("DTDX", nc3Float, timeDim,
                                  lvlDim, latDim, lonDim)))
		return NC_ERR;
	if (!(dqdx = dataFile.add_var("DQVDX", nc3Float, timeDim,
                                  lvlDim, latDim, lonDim)))
		return NC_ERR;
	if (!(dpdx = dataFile.add_var("DPDX", nc3Float, timeDim,
                                  lvlDim, latDim, lonDim)))
		return NC_ERR;
	if (!(dtdy = dataFile.add_var("DTDY", nc3Float, timeDim,
                                  lvlDim, latDim, lonDim)))
		return NC_ERR;
	if (!(dqdy = dataFile.add_var("DQVDY", nc3Float, timeDim,
                                  lvlDim, latDim, lonDim)))
		return NC_ERR;
	if (!(dpdy = dataFile.add_var("DPDY", nc3Float, timeDim,
                                  lvlDim, latDim, lonDim)))
		return NC_ERR;
	if (!(dtdz = dataFile.add_var("DTDZ", nc3Float, timeDim,
                                  lvlDim, latDim, lonDim)))
		return NC_ERR;
	if (!(dqdz = dataFile.add_var("DQVDZ", nc3Float, timeDim,
                                  lvlDim, latDim, lonDim)))
		return NC_ERR;
	if (!(dpdz = dataFile.add_var("DPDZ", nc3Float, timeDim,
                                  lvlDim, latDim, lonDim)))
		return NC_ERR;
	if (!(drhodx = dataFile.add_var("DRHODX", nc3Float, timeDim,
                                    lvlDim, latDim, lonDim)))
		return NC_ERR;
	if (!(drhody = dataFile.add_var("DRHODY", nc3Float, timeDim,
                                    lvlDim, latDim, lonDim)))
		return NC_ERR;
	if (!(drhodz = dataFile.add_var("DRHODZ", nc3Float, timeDim,
                                    lvlDim, latDim, lonDim)))
		return NC_ERR;
	if (!(mcresidual = dataFile.add_var("MCRESIDUAL", nc3Float, timeDim,
                                    lvlDim, latDim, lonDim)))
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
	if ((*configHash)["qr_variable"] == "dbz") {
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
	if (!dudx->add_att("units", "10-5s-1"))
		return NC_ERR;
	if (!dvdx->add_att("units", "10-5s-1"))
		return NC_ERR;
	if (!dwdx->add_att("units", "10-5s-1"))
		return NC_ERR;
	if (!dudy->add_att("units", "10-5s-1"))
		return NC_ERR;
	if (!dvdy->add_att("units", "10-5s-1"))
		return NC_ERR;
	if (!dwdy->add_att("units", "10-5s-1"))
		return NC_ERR;
	if (!dudz->add_att("units", "10-5s-1"))
		return NC_ERR;
	if (!dvdz->add_att("units", "10-5s-1"))
		return NC_ERR;
	if (!dwdz->add_att("units", "10-5s-1"))
		return NC_ERR;
    if (!dtdx->add_att("units", "K km-1"))
		return NC_ERR;
	if (!dqdx->add_att("units", "g kg-1 km-1"))
		return NC_ERR;
	if (!dpdx->add_att("units", "hPa km-1"))
		return NC_ERR;
	if (!dtdy->add_att("units", "K km-1"))
		return NC_ERR;
	if (!dqdy->add_att("units", "g kg-1 km-1"))
		return NC_ERR;
	if (!dpdy->add_att("units", "hPa km-1"))
		return NC_ERR;
	if (!dtdz->add_att("units", "K km-1"))
		return NC_ERR;
	if (!dqdz->add_att("units", "g kg-1 km-1"))
		return NC_ERR;
	if (!dpdz->add_att("units", "hPa km-1"))
		return NC_ERR;
	if (!drhodx->add_att("units", "kg m-3 km-1"))
		return NC_ERR;
	if (!drhody->add_att("units", "kg m-3 km-1"))
		return NC_ERR;
	if (!drhodz->add_att("units", "kg m-3 km-1"))
		return NC_ERR;
	if (!mcresidual->add_att("units", "kg m-3 km-1"))
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
	if (!dudx->add_att("long_name", "wind gradient"))
		return NC_ERR;
	if (!dvdx->add_att("long_name", "wind gradient"))
		return NC_ERR;
	if (!dwdx->add_att("long_name", "wind gradient"))
		return NC_ERR;
	if (!dudy->add_att("long_name", "wind gradient"))
		return NC_ERR;
	if (!dvdy->add_att("long_name", "wind gradient"))
		return NC_ERR;
	if (!dwdy->add_att("long_name", "wind gradient"))
		return NC_ERR;
	if (!dudz->add_att("long_name", "wind gradient"))
		return NC_ERR;
	if (!dvdz->add_att("long_name", "wind gradient"))
		return NC_ERR;
	if (!dwdz->add_att("long_name", "wind gradient"))
		return NC_ERR;
	if (!dtdx->add_att("long_name", "temperature gradient"))
		return NC_ERR;
	if (!dqdx->add_att("long_name", "moisture gradient"))
		return NC_ERR;
	if (!dpdx->add_att("long_name", "pressure gradient"))
		return NC_ERR;
	if (!dtdy->add_att("long_name", "temperature gradient"))
		return NC_ERR;
	if (!dqdy->add_att("long_name", "moisture gradient"))
		return NC_ERR;
	if (!dpdy->add_att("long_name", "pressure gradient"))
		return NC_ERR;
	if (!dtdz->add_att("long_name", "temperature gradient"))
		return NC_ERR;
	if (!dqdz->add_att("long_name", "moisture gradient"))
		return NC_ERR;
	if (!dpdz->add_att("long_name", "pressure gradient"))
		return NC_ERR;
	if (!drhodx->add_att("long_name", "density gradient"))
		return NC_ERR;
	if (!drhody->add_att("long_name", "density gradient"))
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
	if (!dudx->add_att("missing_value", -999.f))
		return NC_ERR;
	if (!dvdx->add_att("missing_value", -999.f))
		return NC_ERR;
	if (!dwdx->add_att("missing_value", -999.f))
		return NC_ERR;
	if (!dudy->add_att("missing_value", -999.f))
		return NC_ERR;
	if (!dvdy->add_att("missing_value", -999.f))
		return NC_ERR;
	if (!dwdy->add_att("missing_value", -999.f))
		return NC_ERR;
	if (!dudz->add_att("missing_value", -999.f))
		return NC_ERR;
	if (!dvdz->add_att("missing_value", -999.f))
		return NC_ERR;
	if (!dwdz->add_att("missing_value", -999.f))
		return NC_ERR;
	if (!dtdx->add_att("missing_value", -999.f))
		return NC_ERR;
	if (!dqdx->add_att("missing_value", -999.f))
		return NC_ERR;
	if (!dpdx->add_att("missing_value", -999.f))
		return NC_ERR;
	if (!dtdy->add_att("missing_value", -999.f))
		return NC_ERR;
	if (!dqdy->add_att("missing_value", -999.f))
		return NC_ERR;
	if (!dpdy->add_att("missing_value", -999.f))
		return NC_ERR;
	if (!dtdz->add_att("missing_value", -999.f))
		return NC_ERR;
	if (!dqdz->add_att("missing_value", -999.f))
		return NC_ERR;
	if (!dpdz->add_att("missing_value", -999.f))
		return NC_ERR;
	if (!drhodx->add_att("missing_value", -999.f))
		return NC_ERR;
	if (!drhody->add_att("missing_value", -999.f))
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
	if (!dudx->add_att("_FillValue", -999.f))
		return NC_ERR;
	if (!dvdx->add_att("_FillValue", -999.f))
		return NC_ERR;
	if (!dwdx->add_att("_FillValue", -999.f))
		return NC_ERR;
	if (!dudy->add_att("_FillValue", -999.f))
		return NC_ERR;
	if (!dvdy->add_att("_FillValue", -999.f))
		return NC_ERR;
	if (!dwdy->add_att("_FillValue", -999.f))
		return NC_ERR;
	if (!dudz->add_att("_FillValue", -999.f))
		return NC_ERR;
	if (!dvdz->add_att("_FillValue", -999.f))
		return NC_ERR;
	if (!dwdz->add_att("_FillValue", -999.f))
		return NC_ERR;
	if (!dtdx->add_att("_FillValue", -999.f))
		return NC_ERR;
	if (!dqdx->add_att("_FillValue", -999.f))
		return NC_ERR;
	if (!dpdx->add_att("_FillValue", -999.f))
		return NC_ERR;
	if (!dtdy->add_att("_FillValue", -999.f))
		return NC_ERR;
	if (!dqdy->add_att("_FillValue", -999.f))
		return NC_ERR;
	if (!dpdy->add_att("_FillValue", -999.f))
		return NC_ERR;
	if (!dtdz->add_att("_FillValue", -999.f))
		return NC_ERR;
	if (!dqdz->add_att("_FillValue", -999.f))
		return NC_ERR;
	if (!dpdz->add_att("_FillValue", -999.f))
		return NC_ERR;
	if (!drhodx->add_att("_FillValue", -999.f))
		return NC_ERR;
	if (!drhody->add_att("_FillValue", -999.f))
		return NC_ERR;
	if (!drhodz->add_att("_FillValue", -999.f))
		return NC_ERR;
	if (!mcresidual->add_att("_FillValue", -999.f))
		return NC_ERR;

	// Write the coordinate variable data to the file.
	real *lons = new real[iDim];
	real *lats = new real[jDim];
	real *levs = new real[sDim];
	real *x = new real[iDim];
	real *y = new real[jDim];
	int time[2];

	// Reference time and position from center file
	time[0] = std::stoi((*configHash)["ref_time"]);
	real latReference = std::stof((*configHash)["ref_lat"]);
	real lonReference = std::stof((*configHash)["ref_lon"]);
	real refX, refY;

	// GeographicLib::LambertConformalConic tm = GeographicLib::LambertConformalConic::Mercator();
	projection.Forward(lonReference, latReference, lonReference, refX, refY);
	for (int iIndex = 0; iIndex < iDim; iIndex++) {
		real i = (iMin + DI * iIndex)*1000;
		real j = (jMin + DJ * (jDim/2))*1000;
		real latnull = 0;
		projection.Reverse(lonReference,refX + i, refY + j, latnull, lons[iIndex]);
		x[iIndex] = i/1000;
	}

	for (int jIndex = 0; jIndex < jDim; jIndex++) {
		real i = (iMin + DI * (iDim/2))*1000;
		real j = (jMin + DJ * jIndex)*1000;
		real lonnull = 0;
		projection.Reverse(lonReference,refX + i, refY + j, lats[jIndex], lonnull);
		y[jIndex] = j/1000;
	}

	if (!lonVar->put(lons, iDim))
		return NC_ERR;

	if (!latVar->put(lats, jDim))
		return NC_ERR;

	if (!xVar->put(x, iDim))
		return NC_ERR;

	if (!yVar->put(y, jDim))
		return NC_ERR;

        const real sigma[42] = { 31885.0, 29385.0, 25085.0, 21985.0, 19685.0, 17885.0,
                16385.0, 15072.5, 13910.0, 12860.0, 11885.0, 10955.0, 10065.0,
                9215.0, 8405.0, 7635.0, 6905.0, 6215.0, 5565.0, 4955.0, 4385.0,
                3855.0, 3365.0, 2915.0, 2505.0, 2135.0, 1805.0, 1515.0, 1265.0,
                1050.0, 860.0, 690.0, 540.0, 410.0, 300.0, 210.0, 140.0, 90.0,
                55.0, 30.0, 10.0, 0.0 };

	for (int kIndex = 1; kIndex < sDim - 1; kIndex++) {
		levs[kIndex - 1] = (sigma[kIndex] * (sigma[0] - sigma[41])
				    / sigma[0] + sigma[41]) / 1000.0;
	}

	if (!lvlVar->put(levs, sDim))
		return NC_ERR;

	if (!timeVar->put(time, 1))
		return NC_ERR;

	// Write the data.
	for (int rec = 0; rec < 1; rec++)
	{
		if (!u->put_rec(&finalAnalysis[0], rec))
			return NC_ERR;
		if (!v->put_rec(&finalAnalysis[iDim*jDim*sDim*1], rec))
			return NC_ERR;
		if (!w->put_rec(&finalAnalysis[iDim*jDim*sDim*2], rec))
			return NC_ERR;
		if (!wspd->put_rec(&finalAnalysis[iDim*jDim*sDim*3], rec))
			return NC_ERR;
		if (!relhum->put_rec(&finalAnalysis[iDim*jDim*sDim*4], rec))
			return NC_ERR;
		if (!hprime->put_rec(&finalAnalysis[iDim*jDim*sDim*5], rec))
			return NC_ERR;
		if (!qvprime->put_rec(&finalAnalysis[iDim*jDim*sDim*6], rec))
			return NC_ERR;
		if (!rhoprime->put_rec(&finalAnalysis[iDim*jDim*sDim*7], rec))
			return NC_ERR;
		if (!tprime->put_rec(&finalAnalysis[iDim*jDim*sDim*8], rec))
			return NC_ERR;
		if (!pprime->put_rec(&finalAnalysis[iDim*jDim*sDim*9], rec))
			return NC_ERR;
		if (!vorticity->put_rec(&finalAnalysis[iDim*jDim*sDim*10], rec))
			return NC_ERR;
		if (!divergence->put_rec(&finalAnalysis[iDim*jDim*sDim*11], rec))
			return NC_ERR;
		if (!okuboweiss->put_rec(&finalAnalysis[iDim*jDim*sDim*12], rec))
			return NC_ERR;
		if (!strain->put_rec(&finalAnalysis[iDim*jDim*sDim*13], rec))
			return NC_ERR;
		if (!tpw->put_rec(&finalAnalysis[iDim*jDim*sDim*14], rec))
			return NC_ERR;
		if (!rhou->put_rec(&finalAnalysis[iDim*jDim*sDim*15], rec))
			return NC_ERR;
		if (!rhov->put_rec(&finalAnalysis[iDim*jDim*sDim*16], rec))
			return NC_ERR;
		if (!rhow->put_rec(&finalAnalysis[iDim*jDim*sDim*17], rec))
			return NC_ERR;
		if (!rho->put_rec(&finalAnalysis[iDim*jDim*sDim*18], rec))
			return NC_ERR;
		if (!press->put_rec(&finalAnalysis[iDim*jDim*sDim*19], rec))
			return NC_ERR;
		if (!temp->put_rec(&finalAnalysis[iDim*jDim*sDim*20], rec))
			return NC_ERR;
		if (!qv->put_rec(&finalAnalysis[iDim*jDim*sDim*21], rec))
			return NC_ERR;
		if (!h->put_rec(&finalAnalysis[iDim*jDim*sDim*22], rec))
			return NC_ERR;
		if (!qr->put_rec(&finalAnalysis[iDim*jDim*sDim*23], rec))
			return NC_ERR;
		if (!absVorticity->put_rec(&finalAnalysis[iDim*jDim*sDim*24], rec))
			return NC_ERR;
		if (!dewp->put_rec(&finalAnalysis[iDim*jDim*sDim*25], rec))
			return NC_ERR;
		if (!theta->put_rec(&finalAnalysis[iDim*jDim*sDim*26], rec))
			return NC_ERR;
		if (!thetae->put_rec(&finalAnalysis[iDim*jDim*sDim*27], rec))
			return NC_ERR;
		if (!thetaes->put_rec(&finalAnalysis[iDim*jDim*sDim*28], rec))
			return NC_ERR;
		if (!dudx->put_rec(&finalAnalysis[iDim*jDim*sDim*29], rec))
			return NC_ERR;
		if (!dvdx->put_rec(&finalAnalysis[iDim*jDim*sDim*30], rec))
			return NC_ERR;
		if (!dwdx->put_rec(&finalAnalysis[iDim*jDim*sDim*31], rec))
			return NC_ERR;
		if (!dudy->put_rec(&finalAnalysis[iDim*jDim*sDim*32], rec))
			return NC_ERR;
		if (!dvdy->put_rec(&finalAnalysis[iDim*jDim*sDim*33], rec))
			return NC_ERR;
		if (!dwdy->put_rec(&finalAnalysis[iDim*jDim*sDim*34], rec))
			return NC_ERR;
		if (!dudz->put_rec(&finalAnalysis[iDim*jDim*sDim*35], rec))
			return NC_ERR;
		if (!dvdz->put_rec(&finalAnalysis[iDim*jDim*sDim*36], rec))
			return NC_ERR;
		if (!dwdz->put_rec(&finalAnalysis[iDim*jDim*sDim*37], rec))
			return NC_ERR;
		if (!dtdx->put_rec(&finalAnalysis[iDim*jDim*sDim*38], rec))
			return NC_ERR;
		if (!dtdy->put_rec(&finalAnalysis[iDim*jDim*sDim*39], rec))
			return NC_ERR;
		if (!dtdz->put_rec(&finalAnalysis[iDim*jDim*sDim*40], rec))
			return NC_ERR;
		if (!dqdx->put_rec(&finalAnalysis[iDim*jDim*sDim*41], rec))
			return NC_ERR;
		if (!dqdy->put_rec(&finalAnalysis[iDim*jDim*sDim*42], rec))
			return NC_ERR;
		if (!dqdz->put_rec(&finalAnalysis[iDim*jDim*sDim*43], rec))
			return NC_ERR;
		if (!dpdx->put_rec(&finalAnalysis[iDim*jDim*sDim*44], rec))
			return NC_ERR;
		if (!dpdy->put_rec(&finalAnalysis[iDim*jDim*sDim*45], rec))
			return NC_ERR;
		if (!dpdz->put_rec(&finalAnalysis[iDim*jDim*sDim*46], rec))
			return NC_ERR;
		if (!drhodx->put_rec(&finalAnalysis[iDim*jDim*sDim*47], rec))
			return NC_ERR;
		if (!drhody->put_rec(&finalAnalysis[iDim*jDim*sDim*48], rec))
			return NC_ERR;
		if (!drhodz->put_rec(&finalAnalysis[iDim*jDim*sDim*49], rec))
			return NC_ERR;
		if (!mcresidual->put_rec(&finalAnalysis[iDim*jDim*sDim*50], rec))
			return NC_ERR;

	}

	// The file is automatically closed by the destructor. This frees
	// up any internal netCDF resources associated with the file, and
	// flushes any buffers.
	delete[] lats;
	delete[] lons;
	delete[] levs;
	delete[] x;
	delete[] y;
*/
	return true;

}

bool CostFunctionCOAMPS::writeAsi(const std::string& asiFileName)
{
	// Initialize header
	int id[511];
	for (int n = 1; n <= 510; n++) {
		id[n]=-999;
	}

	// Calculate headers
	/*fixme
	std::vector<std::string> fieldNames;
	fieldNames.push_back("U");
	fieldNames.push_back("V");
	fieldNames.push_back("W");
	fieldNames.push_back("WS");
	fieldNames.push_back("RH");
	fieldNames.push_back("HP");
	fieldNames.push_back("QP");
	fieldNames.push_back("RP");
	fieldNames.push_back("TP");
	fieldNames.push_back("PP");
	fieldNames.push_back("VO");
	fieldNames.push_back("DV");
	fieldNames.push_back("OW");
	fieldNames.push_back("S");
	fieldNames.push_back("PW");
	fieldNames.push_back("MU");
	fieldNames.push_back("MV");
	fieldNames.push_back("MW");
	fieldNames.push_back("RO");
	fieldNames.push_back("PS");
	fieldNames.push_back("TK");
	fieldNames.push_back("QV");
	fieldNames.push_back("HH");
	fieldNames.push_back("DZ");
	fieldNames.push_back("AV");
	fieldNames.push_back("DP");
	fieldNames.push_back("TH");
	fieldNames.push_back("TE");
	fieldNames.push_back("TS");
	id[175] = fieldNames.size();
	for(int n = 0; n < id[175]; n++) {
		std::string name_1 = fieldNames.at(n).substr(1);
		std::string name_2 = fieldNames.at(n).substr(fieldNames.at(n).length()/2,1);  // Possibly wrong; unsure how Qt's 'mid' corresponds to a substr in the middle (C++11)
		int int_1 = *name_1.toLatin1().data();
		int int_2 = *name_2.toLatin1().data();
		id[176 + (5 * n)] = (int_1 * 256) + int_2;
		id[177 + (5 * n)] = 8224;
		id[178 + (5 * n)] = 8224;
		id[179 + (5 * n)] = 8224;
		id[180 + (5 * n)] = 1;
	}*/

	// Cartesian file
	id[16] = 17217;
	id[17] = 21076;

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
	id[166] = (int)jMax*100;
	id[167] = (int)jDim;
	id[168] = (int)(DJ * 1000);
	id[169] = 2;

	// Z Header
	id[170] = (int)kMin*1000;
	id[171] = (int)kMax*1000;
	id[172] = (int)sDim;
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
/*fixme
	std::ofstream out(asiFileName);
	if(!out.is_open()) {
		std::cout << "Can't open CAPPI file for writing" << std::endl;
		return false;
	}

	// Write header
    int line = 0;
	for (int n = 1; n <= 510; n++) {
		line++;
		out.width(8) << id[n];
		if (line == 10) {
			out << endl;
            line = 0;
		}
	}

	// Write data
	for(int k = 0; k < kDim; k++) {
		out << reset << "level" << qSetFieldWidth(2) << k+1 << endl;
		for(int j = 0; j < jDim; j++) {
			out << reset << "azimuth" << qSetFieldWidth(3) << j+1 << endl;
			for(int n = 0; n < fieldNames.size(); n++) {
				out << reset << left << fieldNames.at(n) << endl;
				int line = 0;
				for (int i = 0; i < iDim;  i++){
					out << reset;
					out.precision(3) << scientific
					out.width(10) << qSetFieldWidth(10) <<
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
*/
	return true;
}
