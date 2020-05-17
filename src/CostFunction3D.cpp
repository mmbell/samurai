/*
 *  CostFunction3D.cpp
 *  samurai
 *
 *  Copyright 2008 Michael Bell. All rights reserved.
 *
 */

#include <cmath>
#include <GeographicLib/TransverseMercatorExact.hpp>

#include "CostFunction3D.h"
#include "MetObs.h"
#include "timing/gptl.h"

#define INDEX(i, j, k, idim, jdim, vdim, var) ((vdim) * ((idim) * ((jdim) * (k) + j) + i) + var)
#define KINDEX(i,dim,var) (dim * var + i)

CostFunction3D::CostFunction3D(const Projection& proj, const int& numObs, const int& stateSize)
  : CostFunction(proj, numObs, stateSize)
{

  // Set up the boundary condition hash
  bcHash["R0"] = R0;
  bcHash["R1T0"] = R1T0;
  bcHash["R1T1"] = R1T1;
  bcHash["R1T2"] = R1T2;
  bcHash["R1T10"] = R1T10;
  bcHash["R2T10"] = R2T10;
  bcHash["R2T20"] = R2T20;
  bcHash["R3"] = R3;
  bcHash["PERIODIC"] = PERIODIC;

  rankHash[R0] = 0;
  rankHash[R1T0] = 1;
  rankHash[R1T1] = 1;
  rankHash[R1T2] = 1;
  rankHash[R2T10] = 2;
  rankHash[R2T20] = 2;
  rankHash[PERIODIC] = 1;
  rankHash[R3] = 3;

  // Set the derivative array
  derivative[0][0] = 0;
  derivative[0][1] = 0;
  derivative[0][2] = 0;

  derivative[1][0] = 1;
  derivative[1][1] = 0;
  derivative[1][2] = 0;

  derivative[2][0] = 0;
  derivative[2][1] = 1;
  derivative[2][2] = 0;

  derivative[3][0] = 0;
  derivative[3][1] = 0;
  derivative[3][2] = 1;

  // Use the full basis unless otherwise specified
  basisappx = 0;

}

CostFunction3D::~CostFunction3D()
{
}

void CostFunction3D::finalize()
{

  delete iFilter;
  delete jFilter;
  delete kFilter;
  delete[] currState;
  delete[] currGradient;
  delete[] tempState;
  delete[] tempGradient;
  delete[] xt;
  delete[] df;
  delete[] CTHTd;
  delete[] stateU;
  delete[] obsVector;
  delete[] HCq;
  delete[] innovation;
  delete[] bgState;
  delete[] bgStdDev;
  delete[] stateA;
  delete[] stateB;
  delete[] stateC;
  delete[] H;
  delete[] JH;
  delete[] IH;

  if (basisappx > 0) {
    delete[] basis0;
    delete[] basis1;
  }
  for (int var = 0; var < varDim; ++var) {
    delete[] iGamma[var];
    delete[] jGamma[var];
    delete[] kGamma[var];
    delete[] iL[var];
    delete[] jL[var];
    delete[] kL[var];
  }


  delete[] kGammaL;
  delete[] kLL;
  fftw_destroy_plan(iForward);
  fftw_destroy_plan(iBackward);
  fftw_destroy_plan(jForward);
  fftw_destroy_plan(jBackward);
  fftw_destroy_plan(kForward);
  fftw_destroy_plan(kBackward);
  fftw_free(iFFTin);
  fftw_free(jFFTin);
  fftw_free(kFFTin);
  fftw_free(iFFTout);
  fftw_free(jFFTout);
  fftw_free(kFFTout);
}

void CostFunction3D::initialize(HashMap* config,
				real* bgU, real* obs, ReferenceState* ref)
{
  // Initialize number of variables
  varDim = 7;
  derivDim = 4;
  configHash = config;

  /* Set the output path */
  outputPath = (*configHash)["output_directory"];

  // Horizontal boundary conditions
  iBCL[0] = bcHash[(*configHash)["i_rhou_bcL"]];
  iBCR[0] = bcHash[(*configHash)["i_rhou_bcR"]];
  iBCL[1] = bcHash[(*configHash)["i_rhov_bcL"]];
  iBCR[1] = bcHash[(*configHash)["i_rhov_bcR"]];
  iBCL[2] = bcHash[(*configHash)["i_rhow_bcL"]];
  iBCR[2] = bcHash[(*configHash)["i_rhow_bcR"]];
  iBCL[3] = bcHash[(*configHash)["i_tempk_bcL"]];
  iBCR[3] = bcHash[(*configHash)["i_tempk_bcR"]];
  iBCL[4] = bcHash[(*configHash)["i_qv_bcL"]];
  iBCR[4] = bcHash[(*configHash)["i_qv_bcR"]];
  iBCL[5] = bcHash[(*configHash)["i_rhoa_bcL"]];
  iBCR[5] = bcHash[(*configHash)["i_rhoa_bcR"]];
  iBCL[6] = bcHash[(*configHash)["i_qr_bcL"]];
  iBCR[6] = bcHash[(*configHash)["i_qr_bcR"]];

  jBCL[0] = bcHash[(*configHash)["j_rhou_bcL"]];
  jBCR[0] = bcHash[(*configHash)["j_rhou_bcR"]];
  jBCL[1] = bcHash[(*configHash)["j_rhov_bcL"]];
  jBCR[1] = bcHash[(*configHash)["j_rhov_bcR"]];
  jBCL[2] = bcHash[(*configHash)["j_rhow_bcL"]];
  jBCR[2] = bcHash[(*configHash)["j_rhow_bcR"]];
  jBCL[3] = bcHash[(*configHash)["j_tempk_bcL"]];
  jBCR[3] = bcHash[(*configHash)["j_tempk_bcR"]];
  jBCL[4] = bcHash[(*configHash)["j_qv_bcL"]];
  jBCR[4] = bcHash[(*configHash)["j_qv_bcR"]];
  jBCL[5] = bcHash[(*configHash)["j_rhoa_bcL"]];
  jBCR[5] = bcHash[(*configHash)["j_rhoa_bcR"]];
  jBCL[6] = bcHash[(*configHash)["j_qr_bcL"]];
  jBCR[6] = bcHash[(*configHash)["j_qr_bcR"]];

  kBCL[0] = bcHash[(*configHash)["k_rhou_bcL"]];
  kBCR[0] = bcHash[(*configHash)["k_rhou_bcR"]];
  kBCL[1] = bcHash[(*configHash)["k_rhov_bcL"]];
  kBCR[1] = bcHash[(*configHash)["k_rhov_bcR"]];
  kBCL[2] = bcHash[(*configHash)["k_rhow_bcL"]];
  kBCR[2] = bcHash[(*configHash)["k_rhow_bcR"]];
  kBCL[3] = bcHash[(*configHash)["k_tempk_bcL"]];
  kBCR[3] = bcHash[(*configHash)["k_tempk_bcR"]];
  kBCL[4] = bcHash[(*configHash)["k_qv_bcL"]];
  kBCR[4] = bcHash[(*configHash)["k_qv_bcR"]];
  kBCL[5] = bcHash[(*configHash)["k_rhoa_bcL"]];
  kBCR[5] = bcHash[(*configHash)["k_rhoa_bcR"]];
  kBCL[6] = bcHash[(*configHash)["k_qr_bcL"]];
  kBCR[6] = bcHash[(*configHash)["k_qr_bcR"]];

  // Define the Reference state
  refstate = ref;

  // Assign local object pointers
  bgFields = bgU;
  rawObs = obs;
  iMin = std::stof((*configHash)["i_min"]);
  iMax = std::stof((*configHash)["i_max"]);
  DI = std::stof((*configHash)["i_incr"]);
  iDim = (int)((iMax - iMin)/DI) + 1;
  jMin = std::stof((*configHash)["j_min"]);
  jMax = std::stof((*configHash)["j_max"]);
  DJ = std::stof((*configHash)["j_incr"]);
  jDim = (int)((jMax - jMin)/DJ) + 1;
  kMin = std::stof((*configHash)["k_min"]);
  kMax = std::stof((*configHash)["k_max"]);
  DK = std::stof((*configHash)["k_incr"]);
  kDim = (int)((kMax - kMin)/DK) + 1;
	std::cout << "done with stof conversion..." << std::endl;

  DIrecip = 1. / DI;
  DJrecip = 1. / DJ;
  DKrecip = 1. / DK;

  // Adjust the internal, variable domain to include boundaries
  adjustInternalDomain(1);

  // Define nodes with internal domain
  int nodes = iDim * jDim * kDim;

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
  obsVector = new real[mObs*(7+varDim*derivDim)];
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
  if (kBCL[0] == PERIODIC) {
    kLDim = kDim-2;
  } else {
    kLDim = 4;
  }
  kRankMax=0;
  for (int var = 0; var < varDim; ++var) {
    iRank[var] = iDim - rankHash[iBCL[var]] - rankHash[iBCR[var]];
    if (iBCL[var] == PERIODIC) iRank[var]--;
    jRank[var] = jDim - rankHash[jBCL[var]] - rankHash[jBCR[var]];
    if (jBCL[var] == PERIODIC) jRank[var]--;
    kRank[var] = kDim - rankHash[kBCL[var]] - rankHash[kBCR[var]];
    if (kBCL[var] == PERIODIC) kRank[var]--;
    // Need to verify that Rank is sufficient (>2?)
    iL[var] = new real[iRank[var]*iLDim];
    jL[var] = new real[jRank[var]*jLDim];
    kL[var] = new real[kRank[var]*kLDim];
    iGamma[var] = new real[iRank[var] * iDim];
    jGamma[var] = new real[jRank[var] * jDim];
    kGamma[var] = new real[kRank[var] * kDim];
    kRankMax = max(kRank[var],kRankMax);
  }
  
  kGammaL = new real[kRankMax * varDim * kDim];
  kLL     = new real[kRankMax * varDim * kLDim];
  cout << "kRankMax: " << kRankMax << "\n"; 

  /* Precalculate the basis functions for lookup table option
     basisappx = configHash->value("spline_approximation").toInt();
     if (basisappx > 0) {
     basis0 = new real[2000000];
     basis1 = new real[2000000];
     fillBasisLookup();
     } */

  // Initialize the Fourier transforms
  iFFTin = (double*) fftw_malloc(sizeof(double) * iDim);
  iFFTout = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * iDim);
  iForward = fftw_plan_dft_r2c_1d(iDim, iFFTin, iFFTout, FFTW_MEASURE);
  iBackward = fftw_plan_dft_c2r_1d(iDim, iFFTout, iFFTin, FFTW_MEASURE);
  jFFTin = (double*) fftw_malloc(sizeof(double) * jDim);
  jFFTout = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * jDim);
  jForward = fftw_plan_dft_r2c_1d(jDim, jFFTin, jFFTout, FFTW_MEASURE);
  jBackward = fftw_plan_dft_c2r_1d(jDim, jFFTout, jFFTin, FFTW_MEASURE);
  kFFTin = (double*) fftw_malloc(sizeof(double) * kDim);
  kFFTout = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * kDim);
  kForward = fftw_plan_dft_r2c_1d(kDim, kFFTin, kFFTout, FFTW_MEASURE);
  kBackward = fftw_plan_dft_c2r_1d(kDim, kFFTout, kFFTin, FFTW_MEASURE);

}

void CostFunction3D::initState(const int iteration)
{
  GPTLstart("CostFunction3D::initState");
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

    
  // Set up the recursive filter
  iFilterScale = std::stof((*configHash)["i_filter_length"]);
  jFilterScale = std::stof((*configHash)["j_filter_length"]);
  kFilterScale = std::stof((*configHash)["k_filter_length"]);
  iFilter->setFilterLengthScale(iFilterScale);
  jFilter->setFilterLengthScale(jFilterScale);
  kFilter->setFilterLengthScale(kFilterScale);

  // Set up the Fourier filter
  for (int i = 0; i < 7; i++) {
    iMaxWavenumber[i] = -1.0;
    jMaxWavenumber[i] = -1.0;
    kMaxWavenumber[i] = -1.0;
  }
	if (configHash->exists("i_max_wavenumber")) {
    // Set all the variables to the same filter
    for (int i = 0; i < 7; i++) {
      iMaxWavenumber[i] = std::stof((*configHash)["i_max_wavenumber"]);
    }
  } else {
    // Set Fourier filter for individual variables
    iMaxWavenumber[0] = std::stof((*configHash)["i_max_wavenumber_rhou"]);
    iMaxWavenumber[1] = std::stof((*configHash)["i_max_wavenumber_rhov"]);
    iMaxWavenumber[2] = std::stof((*configHash)["i_max_wavenumber_rhow"]);
    iMaxWavenumber[3] = std::stof((*configHash)["i_max_wavenumber_tempk"]);
    iMaxWavenumber[4] = std::stof((*configHash)["i_max_wavenumber_qv"]);
    iMaxWavenumber[5] = std::stof((*configHash)["i_max_wavenumber_rhoa"]);
    iMaxWavenumber[6] = std::stof((*configHash)["i_max_wavenumber_qr"]);
  }

	if (configHash->exists("j_max_wavenumber")) {
    for (int i = 0; i < 7; i++) {
      jMaxWavenumber[i] = std::stof((*configHash)["j_max_wavenumber"]);
    }
  } else {
    jMaxWavenumber[0] = std::stof((*configHash)["j_max_wavenumber_rhou"]);
    jMaxWavenumber[1] = std::stof((*configHash)["j_max_wavenumber_rhov"]);
    jMaxWavenumber[2] = std::stof((*configHash)["j_max_wavenumber_rhow"]);
    jMaxWavenumber[3] = std::stof((*configHash)["j_max_wavenumber_tempk"]);
    jMaxWavenumber[4] = std::stof((*configHash)["j_max_wavenumber_qv"]);
    jMaxWavenumber[5] = std::stof((*configHash)["j_max_wavenumber_rhoa"]);
    jMaxWavenumber[6] = std::stof((*configHash)["j_max_wavenumber_qr"]);
  }

	if (configHash->exists("k_max_wavenumber")) {
    for (int i = 0; i < 7; i++) {
      jMaxWavenumber[i] = std::stof((*configHash)["k_max_wavenumber"]);
    }
  } else {
    kMaxWavenumber[0] = std::stof((*configHash)["k_max_wavenumber_rhou"]);
    kMaxWavenumber[1] = std::stof((*configHash)["k_max_wavenumber_rhov"]);
    kMaxWavenumber[2] = std::stof((*configHash)["k_max_wavenumber_rhow"]);
    kMaxWavenumber[3] = std::stof((*configHash)["k_max_wavenumber_tempk"]);
    kMaxWavenumber[4] = std::stof((*configHash)["k_max_wavenumber_qv"]);
    kMaxWavenumber[5] = std::stof((*configHash)["k_max_wavenumber_rhoa"]);
    kMaxWavenumber[6] = std::stof((*configHash)["k_max_wavenumber_qr"]);
  }

  // Set up the spline matrices
  setupSplines();

  // Flag whether or not to print the subgrid information
  if (((*configHash)["output_mish"] == "true") or
      ((*configHash)["save_mish"] == "true")) {
    mishFlag = 1;
  } else {
    mishFlag = 0;
  }

  // Mass continuity weight
  mcWeight = std::stof((*configHash)["mc_weight"]);
  cout << "Mass continuity weight set to " << mcWeight << endl;

  if (iteration == 1) {
    
    cout << "Initializing background..." << endl;
    // Set up the background state
    for (int n = 0; n < nState; n++) {
      bgState[n] = 0.0;
      // Initialize the std. dev. to 1 for the initial SC transform
      bgStdDev[n] = 1.0;
    }

    // If the user specified an error file, read that in instead of using
    // the same fixed value everywhere, and run the same transformations as with
    // the background

    double *mishError = variance.init((*configHash)["fractl_nc_file"], configHash, varDim);
    if (mishError != NULL) {
      double *SBError = new double[nState];
      double *meshError = new double[nState];
      if(meshError == NULL)
	std::cout << "Error: Couldn't allocate mesh error table" << std::endl;
      else {
	// variance.writeDebugNc("debug.out/std_errors.nc", true, mishError);
	SBtransform(mishError, SBError);
	// variance.writeDebugNc("debug.out/std_errors_SB.nc", false, SBError);
	SAtransform(SBError, meshError);
	variance.setMeshData(meshError);
	delete[] SBError;
	// variance.writeDebugNc("debug.out/std_errors_SA.nc", false, meshError);
      }
    }    	// end of background error transformations
 
    if ((*configHash)["load_bg_coefficients"] == "true") {
      for (int n = 0; n < nState; n++) {
	stateA[n] = bgFields[n];
      }
    } else {
      // SB Transform on the original bg fields	(Mish to Mesh)
      // variance.writeDebugNc("debug.out/bgU_orig.nc", true, bgFields);		// mish sized DEBUG
      SBtransform(bgFields, stateB);
      // variance.writeDebugNc("debug.out/StateB.nc", false, stateB);			// mesh sized DEBUG
      // SA transform = bg B's -> bg A's
      SAtransform(stateB, stateA);
      // variance.writeDebugNc("debug.out/StateA.nc", false, stateA);			// mesh sized DEBUG
    }

    // FF transform to match background and increment
    FFtransform(stateA, bgState); 
    // variance.writeDebugNc("debug.out/FF.nc", false, bgState);			// mesh sized DEBUG

  } // end of iteration == 1
  
  for (int var = 0; var < varDim; var++) {
    // Init node variance
    for (int iIndex = 0; iIndex < iDim; iIndex++) {
      for (int jIndex = 0; jIndex < jDim; jIndex++) {
	for (int kIndex = 0; kIndex < kDim; kIndex++) {
	  int bIndex = varDim * iDim * jDim*kIndex
	    + varDim * iDim * jIndex
	    + varDim * iIndex + var;
	  bgStdDev[bIndex] = variance.meshValueAt(var, iIndex, jIndex, kIndex);
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
	  int bIndex = varDim * iDim * jDim * kIndex
	    + varDim * iDim * jIndex
	    + varDim * iIndex + var;
	  varScale += bgState[bIndex] * bgState[bIndex];
	}
      }
    }
    varScale = sqrt(varScale / (iDim * jDim * kDim));

    if (varScale) {
      real errPct = 100 * bgError[var] / varScale;
      cout << "Variable " << var << " RMS = " << varScale << "\t BG Error = " << bgError[var]
	   << " ( " << errPct << " %)" << endl;
    } else {
      cout << "Variable " << var << " RMS = " << varScale << "\t BG Error = " << bgError[var]
	   << " ( Infinite! %)" << endl;
    }
  }

  // Load the obs locally and weight the nonlinear observation operators by interpolated bg fields
  obAdjustments();

  // Calculate the H matrix operator
  calcHmatrix();

  // d = y - HXb
  calcInnovation();

  // Output the original background field
  outputAnalysis("background", bgState);

  cout << "Beginning analysis...\n";

  // HTd
  calcHTranspose(innovation, stateC);

  FFtransform(stateC, stateA);
  SAtransform(stateA, stateB);
  SCtransform(stateB, CTHTd);

  //Htransform(stateB);
  GPTLstop("CostFunction3D::initState");
}

real CostFunction3D::funcValue(real* state)
{
  real qIP, obIP;
  GPTLstart("CostFunction3D::funcValue");
  qIP = 0.;
  obIP = 0.;

  GPTLstart("CostFunction3D::funcValue:updateHCq");
  updateHCq(state);
  GPTLstop("CostFunction3D::funcValue:updateHCq");

  GPTLstart("CostFunction3D::funcValue:other");
  // Compute inner product of state vector
  #pragma omp parallel for reduction(+:qIP)
  for (int n = 0; n < nState; n++) {
    qIP += state[n]*state[n];
  }

  // Subtract d from HCq to yield mObs length vector and compute inner product
  #pragma omp parallel for reduction(+:obIP)
  for (int m = 0; m < mObs; m++) {
    int obIndex = m*(7+varDim*derivDim) + 1;
    obIP += (HCq[m]-innovation[m])*(obsVector[obIndex])*(HCq[m]-innovation[m]);
  }
  GPTLstop("CostFunction3D::funcValue:other");

  real J = 0.5*(qIP + obIP);
  GPTLstop("CostFunction3D::funcValue");
  return J;

}

void CostFunction3D::funcGradient(real* state, real* gradient)
{
  GPTLstart("CostFunction3D::funcGradient");

  GPTLstart("CostFunction3D::funcGradient:updateHCq");
  updateHCq(state);
  GPTLstop("CostFunction3D::funcGradient:updateHCq");

  GPTLstart("CostFunction3D::funcGradient:calcHTranspose");
  // HTHCq
  calcHTranspose(HCq, stateC);
  GPTLstop("CostFunction3D::funcGradient:calcHTranspose");

  GPTLstart("CostFunction3D::funcGradient:FFtransform");
  FFtransform(stateC, stateA);
  GPTLstop("CostFunction3D::funcGradient:FFtransform");
  GPTLstart("CostFunction3D::funcGradient:SAtransform");
  SAtransform(stateA, stateB);
  GPTLstop("CostFunction3D::funcGradient:SAtransform");
  GPTLstart("CostFunction3D::funcGradient:SCtransform");
  SCtransform(stateB, stateC);
  GPTLstop("CostFunction3D::funcGradient:SCtransform");

  GPTLstart("CostFunction3D::funcGradient:gradient");
  for (int n = 0; n < nState; n++) {
    gradient[n] = state[n] + stateC[n] - CTHTd[n];
  }
  GPTLstop("CostFunction3D::funcGradient:gradient");

  GPTLstop("CostFunction3D::funcGradient");
}

void CostFunction3D::updateHCq(real* state)
{
  GPTLstart("CostFunction3D::updateHCq");
  SCtransform(state, stateB);
  SAtransform(stateB, stateA);
  FFtransform(stateA, stateC);
  Htransform(stateC, HCq);

  GPTLstop("CostFunction3D::updateHCq");
}

void CostFunction3D::updateBG()
{
  GPTLstart("CostFunction3D::updateBG");

  // S (SA transform) yield A's
  SCtransform(currState, stateB);
  SAtransform(stateB, stateA);
  FFtransform(stateA, stateC);
  outputAnalysis("increment", stateC);

  // In BG update we are directly summing C + A
  std::string cFilename = outputPath + "/samurai_Coefficients.out";
  ofstream cstream(cFilename);

  cstream << "Variable\tI\tJ\tK\tBackground\tAnalysis\tIncrement\n";
  for (int var = 0; var < varDim; var++) {
    for (int iIndex = 0; iIndex < iDim; iIndex++) {
      for (int jIndex = 0; jIndex < jDim; jIndex++) {
	for (int kIndex = 0; kIndex < kDim; kIndex++) {
	  cstream << var << "\t" << iIndex << "\t" << jIndex << "\t" << kIndex << "\t";
	  int bgIndex = varDim*iDim*jDim*kIndex + varDim*iDim*jIndex +varDim*iIndex + var;
	  cstream << bgState[bgIndex] << "\t";
	  bgState[bgIndex] += stateC[bgIndex];
	  cstream << bgState[bgIndex] << "\t";
	  cstream << stateC[bgIndex] << endl;
	}
      }
    }
  }

  outputAnalysis("analysis", bgState);
  GPTLstop("CostFunction3D::updateBG");
}

void CostFunction3D::calcInnovation()
{
  GPTLstart("CostFunction3D::calcInnovation");
  // Initialize and fill the innovation vector
  cout << "Initializing innovation vector..." << endl;

  // Use HCq to hold the transform, but C is not applied for the innovation
  Htransform(bgState, HCq);

  real innovationRMS = 0.;

#pragma omp parallel for reduction(+:innovationRMS)
  for (int m = 0; m < mObs; m++) {
    innovation[m] = obsVector[m*(7+varDim*derivDim)] - HCq[m];
    innovationRMS += (innovation[m]*innovation[m]);
    HCq[m] = 0.0;
  }

  if (mObs) innovationRMS /= mObs;
  innovationRMS = sqrt(innovationRMS);
  cout << "Innovation RMS : " << innovationRMS << endl;

  GPTLstop("CostFunction3D::calcInnovation");
}

void CostFunction3D::calcHTranspose(const real* yhat, real* Astate)
{
  GPTLstart("CostFunction3D::calcHTranspose");

  // Clear the Astate
  for (int n = 0; n < nState; n++) {
    Astate[n] = 0.0;
  }

  // Calculate H Transpose
  //#pragma omp parallel for  This didn't work
  for(int m=0; m<mObs; ++m) {
    int mi = m*(7+varDim*derivDim)+1;
    const int begin = IH[m];
    const int end = IH[m + 1];
    for(int j=begin; j<end; ++j) {
      //#pragma omp atomic  This didn't work
      Astate[JH[j]] += H[j] * yhat[m] * obsVector[mi];
    }
  }

  GPTLstop("CostFunction3D::calcHTranspose");
}

bool CostFunction3D::SAtransform(const real* Bstate, real* Astate)
{
  real kB[kDim], xk[kDim];
  real jB[jDim], xj[jDim];
  real iB[iDim], xi[iDim];
  int kRankVar;
  int iIndex,jIndex,kIndex;
  int i,j,k,l,m;
  real tmp;

  GPTLstart("CostFunction3D::SAtransform");
  //cout << "SAtransform:: point #1 \n";
  for (int var = 0; var < varDim; var++) {
    kRankVar = kRank[var];
    //GPTLstart("IJK Loop");
    //std::cout << "kDim: " << kDim << " kRankVar: " << kRankVar  << " kLDim: " << kLDim << std::endl;

    //cout << "SAtransform:: point #2 \n";
    #pragma omp parallel for private(tmp,kB,xk,k,l,m,iIndex,jIndex,kIndex) //[5.0.1]
    for (iIndex = 0; iIndex < iDim; iIndex++) {
      for (jIndex = 0; jIndex < jDim; jIndex++) {
	      for (k = 0; k < kDim; k++) {
                kIndex = INDEX(iIndex, jIndex, k, iDim, jDim, varDim, var);
	        kB[k] = Bstate[kIndex];
	      }

      	      // Multiply by gamma
	      for (m = 0; m < kRankVar; m++) {
	        //bk[m] = 0;
                tmp = 0;
	        for (k = 0; k < kDim; k++) {
	          //tmp += kGamma[var][kDim * m + k] * kB[k];
	          tmp += kGammaL[KINDEX(kDim * m + k,kRankMax*kDim,var)] * kB[k];
	        }
	        // std::cout << "m: " << m << " bk[" << m << "]: " << bk[m] << std
	        // Solve for A's using compact storage
                //cout << "-1:" << -(kLDim-1) << "\n"; 
	        for (l=-1;l>=-(kLDim-1);l--) {
	          if ((m + l >= 0) and ((m * kLDim - l) >= 0)) {
	            //tmp -= kL[var][m * kLDim - l] *xk[m + l];
	            tmp -= kLL[KINDEX(m * kLDim - l,kRankMax*kLDim,var)] *xk[m + l];
                  }
	        }
	        //xk[m] = tmp / kL[var][m * kLDim];
	        xk[m] = tmp / kLL[KINDEX(m * kLDim,kRankMax*kLDim,var)];
	      }
	      for (k = kRankVar - 1; k >= 0; k--) {
                tmp = xk[k];
	        for (l = 1; l <= (kLDim - 1); l++) {
	          if ((k + l < kRankVar) and (((k + l) * kLDim + l) < kRankVar * kLDim)) {
	            //tmp -= kL[var][(k + l) * kLDim +l] * xk[k + l];
	            tmp -= kLL[KINDEX((k + l) * kLDim +l,kRankMax*kLDim,var)] * xk[k + l];
                  }
	        }
	        //xk[k] = tmp / kL[var][k * kLDim];
	        xk[k] = tmp / kLL[KINDEX(k * kLDim,kRankMax*kLDim,var)];
	      }

	      for (k = 0; k < kDim; k++) {
	        // Multiply by gammaT
	        tmp = 0;
	        for (m = 0; m < kRankVar; m++) {
	          //tmp += kGamma[var][kDim * m + k] * xk[m];
	          tmp += kGammaL[KINDEX(kDim * m + k,kRankMax*kDim,var)] * xk[m];
	        }
	        // std::cout << "k: " << k << " ak[" << k << "]: " << ak[k] << "\n";
                kIndex = INDEX(iIndex, jIndex, k, iDim, jDim, varDim, var);
	        //Astate[INDEX(iIndex, jIndex, k, iDim, jDim, varDim, var)] = tmp; 
	        Astate[kIndex] = tmp; 
	      }
      }
    }
    //GPTLstop("IJK Loop");
    // cout << "SAtransform:: point #3 \n";

    //GPTLstart("IKJ Loop");
    #pragma omp parallel for private(tmp,jB,xj,j,l,m,iIndex,kIndex) //[5.0.2]
    for (iIndex = 0; iIndex < iDim; iIndex++) {
      for (kIndex = 0; kIndex < kDim; kIndex++) {
	      for (j = 0; j < jDim; j++) {
	        jB[j] = Astate[INDEX(iIndex, j, kIndex, iDim, jDim, varDim, var)];
	      }
	      for (m = 0; m < jRank[var]; m++) {
	        // Multiply by gamma
                tmp = 0;
	        for (j = 0; j < jDim; j++) {
	          tmp += jGamma[var][jDim*m + j]*jB[j];
	        }
	        // Solve for A's using compact storage
	        for (l=-1;l>=-(jLDim-1);l--) {
	          if ((m+l >= 0) and ((m*jLDim-l) >= 0)) {
	            tmp -= jL[var][m*jLDim-l]*xj[m+l];
                  }
	        }
	        xj[m] = tmp/jL[var][m*jLDim];
	      }

	      for (j=jRank[var]-1;j>=0;j--) {
                tmp=xj[j];
	        for (l=1;l<=(jLDim-1);l++) {
	          if ((j+l < jRank[var]) and (((j+l)*jLDim+l) < jRank[var]*jLDim)) {
	            tmp -= jL[var][(j+l)*jLDim+l]*xj[j+l];
                  }
	        }
	        xj[j] = tmp/jL[var][j*jLDim];
	      }

	      for (j = 0; j < jDim; j++) {
	        // Multiply by gammaT
                tmp = 0;
	        for (m = 0; m < jRank[var]; m++) {
	          tmp += jGamma[var][jDim*m + j]*xj[m];
	        }
	        // std::cout << "i: " << i << " aj[" << j << "]: " << tmp << "\n";	  
	        Astate[INDEX(iIndex, j, kIndex, iDim, jDim, varDim, var)] = tmp;
	      }
      }
    }
    //GPTLstop("IKJ Loop");
    //cout << "SAtransform:: point #4 \n";
    //GPTLstart("JKI Loop");
    #pragma omp parallel for private(tmp,iB,xi,i,l,m,jIndex,kIndex) //[5.0.3]
    for (jIndex = 0; jIndex < jDim; jIndex++) {
      for (kIndex = 0; kIndex < kDim; kIndex++) {
	      for (i = 0; i < iDim; i++) {
	        iB[i] = Astate[INDEX(i, jIndex, kIndex, iDim, jDim, varDim, var)];
	      }
	      // Multiply by gamma
	      for (m = 0; m < iRank[var]; m++) {
	        //bi[m] = 0;
	        tmp = 0;
	        for (i = 0; i < iDim; i++) {
	          tmp += iGamma[var][iDim*m + i]*iB[i];
	        }
	        // Solve for A's using compact storage
	        for (l=-1;l>=-(iLDim-1);l--) {
	          if ((m+l >= 0) and ((m*iLDim-l) >= 0)) {
	            tmp -= iL[var][m*iLDim-l]*xi[m+l];
                  }
	        }
	        xi[m] = tmp/iL[var][m*iLDim];
	      }
	      for (i=iRank[var]-1;i>=0;i--) {
                tmp=xi[i];
	        for (l=1;l<=(iLDim-1);l++) {
	          if ((i+l < iRank[var]) and (((i+l)*iLDim+l) < iRank[var]*iLDim)) {
	            tmp -= iL[var][(i+l)*iLDim+l]*xi[i+l];
                  }
	        }
	        xi[i] = tmp/iL[var][i*iLDim];
	      }

	      // Multiply by gammaT
	      for (i = 0; i < iDim; i++) {
	        //ai[i] = 0;
                tmp=0;
	        for (m = 0; m < iRank[var]; m++) {
	          tmp += iGamma[var][iDim*m + i]*xi[m];
	        }
	        // std::cout << "i: " << i << " ai[" << i << "]: " << tmp << "\n";	  
	        Astate[INDEX(i, jIndex, kIndex, iDim, jDim, varDim, var)] = tmp;
	      }

      }
    }
    //GPTLstop("JKI Loop");
    //cout << "SAtransform:: point #5 \n";
    //delete[] kGammaTmp;

  }
  GPTLstop("CostFunction3D::SAtransform");

  return true;
}

bool CostFunction3D::SAtranspose(const real* Astate, real* Bstate)
{
  GPTLstart("CostFunction3D::SAtranspose");

  for (int var = 0; var < varDim; var++) {
    int l;
    real* iB = new real[iDim];
    real* x = new real[iDim];
    real* b = new real[iRank[var]];
    real* a = new real[iDim];
    for (int jIndex = 0; jIndex < jDim; jIndex++) {
     for (int kIndex = 0; kIndex < kDim; kIndex++) {
			for (int iIndex = 0; iIndex < iDim; iIndex++) {
  			iB[iIndex] = Astate[INDEX(iIndex, jIndex, kIndex, iDim, jDim, varDim, var)];
			}
			// Multiply by gamma
			for (int m = 0; m < iRank[var]; m++) {
  			b[m] = 0;
  			for (int i = 0; i < iDim; i++) {
    			b[m] += iGamma[var][iDim*m + i]*iB[i];
  			}
			}

			// Solve for A's using compact storage
			real sum = 0;
			for (int i = 0; i < iRank[var]; i++) {
  			for (sum=b[i], l=-1;l>=-(iLDim-1);l--) {
    			if ((i+l >= 0) and ((i*iLDim-l) >= 0))
      			sum -= iL[var][i*iLDim-l]*x[i+l];
  				}
  				x[i] = sum/iL[var][i*iLDim];
				}
				for (int i=iRank[var]-1;i>=0;i--) {
  				for (sum=x[i], l=1;l<=(iLDim-1);l++) {
    				if ((i+l < iRank[var]) and (((i+l)*iLDim+l) < iRank[var]*iLDim))
      				sum -= iL[var][(i+l)*iLDim+l]*x[i+l];
  				}
  				x[i] = sum/iL[var][i*iLDim];
				}

				// Multiply by gammaT
				for (int i = 0; i < iDim; i++) {
  				a[i] = 0;
  				for (int m = 0; m < iRank[var]; m++) {
    				a[i] += iGamma[var][iDim*m + i]*x[m];
  				}
				}

				for (int iIndex = 0; iIndex < iDim; iIndex++) {
  				Bstate[INDEX(iIndex, jIndex, kIndex, iDim, jDim, varDim, var)] = a[iIndex];
				}

     	}
		}
   	delete[] iB;
   	delete[] x;
   	delete[] b;
   	delete[] a;

   	real* jB = new real[jDim];
   	x = new real[jDim];
   	b = new real[jRank[var]];
   	a = new real[jDim];
   	for (int kIndex = 0; kIndex < kDim; kIndex++) {
     	for (int iIndex = 0; iIndex < iDim; iIndex++) {
				for (int jIndex = 0; jIndex < jDim; jIndex++) {
  				jB[jIndex] = Bstate[INDEX(iIndex, jIndex, kIndex, iDim, jDim, varDim, var)];
				}
				// Multiply by gamma
				for (int m = 0; m < jRank[var]; m++) {
  				b[m] = 0;
  				for (int j = 0; j < jDim; j++) {
    				b[m] += jGamma[var][jDim*m + j]*jB[j];
  				}
				}

				// Solve for A's using compact storage
				real sum = 0;
				for (int j = 0; j < jRank[var]; j++) {
  				for (sum=b[j], l=-1;l>=-(jLDim-1);l--) {
    				if ((j+l >= 0) and ((j*jLDim-l) >= 0))
      				sum -= jL[var][j*jLDim-l]*x[j+l];
  				}
  				x[j] = sum/jL[var][j*jLDim];
				}
				for (int j=jRank[var]-1;j>=0;j--) {
  				for (sum=x[j], l=1;l<=(jLDim-1);l++) {
    				if ((j+l < jRank[var]) and (((j+l)*jLDim+l) < jRank[var]*jLDim))
      				sum -= jL[var][(j+l)*jLDim+l]*x[j+l];
  				}
  				x[j] = sum/jL[var][j*jLDim];
				}

				// Multiply by gammaT
				for (int j = 0; j < jDim; j++) {
  				a[j] = 0;
  				for (int m = 0; m < jRank[var]; m++) {
    				a[j] += jGamma[var][jDim*m + j]*x[m];
  				}
				}
				for (int jIndex = 0; jIndex < jDim; jIndex++) {
  				Bstate[INDEX(iIndex, jIndex, kIndex, iDim, jDim, varDim, var)] = a[jIndex];
				}
     	}
   	}
   	delete[] jB;
   	delete[] x;
   	delete[] b;
   	delete[] a;

   	real* kB = new real[kDim];
   	x = new real[kDim];
   	b = new real[kRank[var]];
   	a = new real[kDim];
   	for (int iIndex = 0; iIndex < iDim; iIndex++) {
     	for (int jIndex = 0; jIndex < jDim; jIndex++) {
				for (int kIndex = 0; kIndex < kDim; kIndex++) {
  				kB[kIndex] = Bstate[INDEX(iIndex, jIndex, kIndex, iDim, jDim, varDim, var)];
				}
				// Multiply by gamma
				for (int m = 0; m < kRank[var]; m++) {
  				b[m] = 0;
  				for (int k = 0; k < kDim; k++) {
    				b[m] += kGamma[var][kDim*m + k]*kB[k];
  				}
  				//std::cout << m << " " << b[m] << "\n";
				}

				// Solve for A's using compact storage
				real sum = 0;
				for (int k = 0; k < kRank[var]; k++) {
  				for (sum=b[k], l=-1;l>=-(kLDim-1);l--) {
    				if ((k+l >= 0) and ((k*kLDim-l) >= 0))
      				sum -= kL[var][k*kLDim-l]*x[k+l];
  				}
  				x[k] = sum/kL[var][k*kLDim];
				}
				for (int k=kRank[var]-1;k>=0;k--) {
  				for (sum=x[k], l=1;l<=(kLDim-1);l++) {
    				if ((k+l < kRank[var]) and (((k+l)*kLDim+l) < kRank[var]*kLDim))
      				sum -= kL[var][(k+l)*kLDim+l]*x[k+l];
  				}
  				x[k] = sum/kL[var][k*kLDim];
				}

				// Multiply by gammaT
				for (int k = 0; k < kDim; k++) {
  				a[k] = 0;
  				for (int m = 0; m < kRank[var]; m++) {
    				a[k] += kGamma[var][kDim*m + k]*x[m];
  				}
  				//std::cout << k << " " << a[k] << "\n";
				}

				for (int kIndex = 0; kIndex < kDim; kIndex++) {
  				Bstate[INDEX(iIndex, jIndex, kIndex, iDim, jDim, varDim, var)] = a[kIndex];
				}
     	}
   	}
   	delete[] kB;
   	delete[] b;
   	delete[] x;
   	delete[] a;

 	}

  GPTLstop("CostFunction3D::SAtranspose");
  return true;
}

void CostFunction3D::SBtransform(const real* Ustate, real* Bstate)
{
  int n,var;
  int is,ie,js,je,ks,ke;

  int ii,iIndex,imu,iNode,uI;
  int jj,jIndex,jmu,jNode,uJ;
  int kk,kIndex,kmu,kNode;
  int ui,bi;
  int iis,iie,jjs,jje,kks,kke;

  real i, ibasis;
  real j, jbasis;
  real k, kbasis;
  GPTLstart("CostFunction3D::SBtransform");
  // Clear the Bstate
  for (n = 0; n < nState; n++) {
    Bstate[n] = 0.;
  }
  real gausspoint = 0.5*sqrt(1./3.);

  for (var = 0; var < varDim; var++) {
    cout << " min(rankHash[iBCL[var]],1): " <<  min(rankHash[iBCL[var]],1) << " max(iDim-1-rankHash[iBCR[var]],iDim-2): " << max(iDim-1-rankHash[iBCR[var]],iDim-2) << "\n";
    is = min(rankHash[iBCL[var]],1);
    ie = max(iDim-1-rankHash[iBCR[var]],iDim-2);
    js = min(rankHash[jBCL[var]],1);
    je =  max(jDim-1-rankHash[jBCR[var]],jDim-2);
    ks = min(rankHash[kBCL[var]],1);
    ke = max(kDim-1-rankHash[kBCR[var]],kDim-2);
    //#pragma acc parallel loop private(ii,iIndex,imu,iNode,uI,jj,jIndex,jmu,jNode,uJ,kk,kIndex,kmu,kNode,ui,bi,iis,iie,jjs,jje,kks,kke,i, ibasis,j, jbasis,k, kbasis)
    for (iIndex = is; iIndex < ie; iIndex++) {
      for (imu = -1; imu <= 1; imu += 2) {
				i = iMin + DI * (iIndex + (gausspoint * imu + 0.5));
				ii = (int)((i - iMin)*DIrecip);
        iis=max(0,ii-1);iie=min(ii+2,iDim);
				for (iNode = iis; iNode <= iie; ++iNode) {
	  			ibasis = Basis(iNode, i, iDim-1, iMin, DI, DIrecip, 0, iBCL[var], iBCR[var]);
	  			uI = iIndex*2 + (imu+1)/2;

	  			for (jIndex = js; jIndex < je; jIndex++) {
	    			for (jmu = -1; jmu <= 1; jmu += 2) {
	      			j = jMin + DJ * (jIndex + (gausspoint * jmu + 0.5));
	      			jj = (int)((j - jMin)*DJrecip);
              jjs=max(0,jj-1);jje=min(jj+2,jDim);
	      			for (jNode = jjs; jNode <= jje; ++jNode) {
								jbasis = Basis(jNode, j, jDim-1, jMin, DJ, DJrecip, 0, jBCL[var], jBCR[var]);
								uJ = jIndex*2 + (jmu+1)/2;

								for (kIndex = ks; kIndex < ke; kIndex++) {
		  						for (kmu = -1; kmu <= 1; kmu += 2) {
		    						k = kMin + DK * (kIndex + (gausspoint * kmu + 0.5));
		    						kk = (int)((k - kMin)*DKrecip);
                    kks=max(0,kk-1);kke=min(kk+2,kDim);
		    						for (kNode = kks; kNode <= kke; ++kNode) {

                      ui = INDEX(uI, uJ, kIndex*2 + (kmu+1)/2, (iDim-1)*2, (jDim-1)*2, varDim, var);
		      						if (Ustate[ui] == 0) continue;
		     							kbasis = Basis(kNode, k, kDim-1, kMin, DK, DKrecip, 0, kBCL[var], kBCR[var]);
		      						bi = INDEX(iNode, jNode, kNode, iDim, jDim, varDim, var);
		      						Bstate[bi] += Ustate[ui] * 0.125 * ibasis * jbasis * kbasis;
		    						}
		  						}
								}
	      			}
	    			}
	  			}
				}
      }
    } // iIndex loop
  } // var loop
  GPTLstop("CostFunction3D::SBtransform");
}


void CostFunction3D::SBtranspose(const real* Bstate, real* Ustate)
{
  GPTLstart("CostFunction3D::SBtranspose");

  // Clear the Ustate
  for (int n = 0; n < nState; n++) {
    Ustate[n] = 0;
  }
  real gausspoint = 0.5*sqrt(1./3.);

  //#pragma omp parallel for
  for (int var = 0; var < varDim; var++) {
    for (int iIndex = min(rankHash[iBCL[var]],1); iIndex < max(iDim-1-rankHash[iBCR[var]],iDim-2); iIndex++) {
      for (int imu = -1; imu <= 1; imu += 2) {
				real i = iMin + DI * (iIndex + (gausspoint * imu + 0.5));
				int ii = (int)((i - iMin)*DIrecip);
				for (int iiNode = (ii-1); iiNode <= (ii+2); ++iiNode) {
	  			int iNode = iiNode;
	  			if ((iNode < 0) or (iNode >= iDim)) continue;
	  			real ibasis = Basis(iNode, i, iDim-1, iMin, DI, DIrecip, 0, iBCL[var], iBCR[var]);
	  			int uI = iIndex*2 + (imu+1)/2;
	  			for (int jIndex = min(rankHash[jBCL[var]],1); jIndex < max(jDim-1-rankHash[jBCR[var]],jDim-2); jIndex++) {
	    			for (int jmu = -1; jmu <= 1; jmu += 2) {
	      			real j = jMin + DJ * (jIndex + (gausspoint * jmu + 0.5));
	      			int jj = (int)((j - jMin)*DJrecip);
	      			for (int jjNode = (jj-1); jjNode <= (jj+2); ++jjNode) {
								int jNode = jjNode;
								if ((jNode < 0) or (jNode >= jDim)) continue;
								real jbasis = Basis(jNode, j, jDim-1, jMin, DJ, DJrecip, 0, jBCL[var], jBCR[var]);
								int uJ = jIndex*2 + (jmu+1)/2;
								real ijbasis = ibasis * jbasis;
								for (int kIndex = min(rankHash[kBCL[var]],1); kIndex < max(kDim-1-rankHash[kBCR[var]],kDim-2); kIndex++) {
		  						for (int kmu = -1; kmu <= 1; kmu += 2) {
		    						real k = kMin + DK * (kIndex + (gausspoint * kmu + 0.5));
		    						int kk = (int)((k - kMin)*DKrecip);
		    						for (int kNode = (kk-1); kNode <= (kk+2); ++kNode) {
		      						if ((kNode < 0) or (kNode >= kDim)) continue;
              				int ui = INDEX(uI, uJ, kIndex*2 + (kmu+1)/2, (iDim-1)*2, (jDim-1)*2, varDim, var);
		      						if (Ustate[ui] == 0) continue;
		      						real kbasis = Basis(kNode, k, kDim-1, kMin, DK, DKrecip, 0, kBCL[var], kBCR[var]);
		      						int bi = INDEX(iNode, jNode, kNode, iDim, jDim, varDim, var);
		      						Ustate[ui] += Bstate[bi] * 0.125 * ijbasis * kbasis;
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
  GPTLstop("CostFunction3D::SBtranspose");
}

void CostFunction3D::SCtransform(const real* Astate, real* Cstate)
{
  real iTemp[iDim], iq[iDim],is[iDim];
  real jTemp[jDim], jq[jDim],js[jDim];
  real kTemp[kDim], kq[kDim],ks[kDim];
  int n,iIndex,jIndex,kIndex,index;

  GPTLstart("CostFunction3D::SCtransform");
  // Disable recursive filter if less than 1
  if ((iFilterScale < 0) and (jFilterScale < 0) and (kFilterScale < 0)) {
    #pragma omp parallel for private(n) //[5.1]
  	for (n = 0; n < nState; n++) {
      Cstate[n]= Astate[n] * bgStdDev[n];
    }
  } else {
    // Isotropic Recursive filter, no anisotropic "triad" working yet
    for (int var = 0; var < varDim; var++) {

      #pragma omp parallel for private(iIndex,jIndex,kIndex,index,kTemp,kq,ks) //[5.2]
      for (iIndex = 0; iIndex < iDim; iIndex++) {
				for (jIndex = 0; jIndex < jDim; jIndex++) {
	  			for (kIndex = 0; kIndex < kDim; kIndex++) {
            index = INDEX(iIndex, jIndex, kIndex, iDim, jDim, varDim, var);
	    			kTemp[kIndex] = Astate[index];
	  			}
	  			if (kFilterScale > 0) kFilter->filterArray(kTemp, kq,ks,kDim);
	  			for (kIndex = 0; kIndex < kDim; kIndex++) {
            index = INDEX(iIndex, jIndex, kIndex, iDim, jDim, varDim, var);
	    			Cstate[index] = kTemp[kIndex];
	  			}
				}
      }

			#pragma omp parallel for private(iIndex,jIndex,kIndex,index,jTemp,jq,js) //[5.3]
      for (iIndex = 0; iIndex < iDim; iIndex++) {
        for (kIndex = 0; kIndex < kDim; kIndex++) {
	  			for (jIndex = 0; jIndex < jDim; jIndex++) {
            index = INDEX(iIndex, jIndex, kIndex, iDim, jDim, varDim, var);
	    			jTemp[jIndex] = Cstate[index];
	  			}
	  			if (jFilterScale > 0) jFilter->filterArray(jTemp, jq, js, jDim);
	  			for (jIndex = 0; jIndex < jDim; jIndex++) {
            index = INDEX(iIndex, jIndex, kIndex, iDim, jDim, varDim, var);
	    			Cstate[index] = jTemp[jIndex];
	  			}
				}
      }

			#pragma omp parallel for private(iIndex,jIndex,kIndex,index,iTemp,iq,is) //[5.4]
      for (jIndex = 0; jIndex < jDim; jIndex++) {
				for (kIndex = 0; kIndex < kDim; kIndex++) {
	  			for (iIndex = 0; iIndex < iDim; iIndex++) {
            index = INDEX(iIndex, jIndex, kIndex, iDim, jDim, varDim, var);
	    			iTemp[iIndex] = Cstate[index];
	  			}
	  			if (iFilterScale > 0) iFilter->filterArray(iTemp, iq, is, iDim);
	  			for (iIndex = 0; iIndex < iDim; iIndex++) {
	    			// D
	    			index = INDEX(iIndex, jIndex, kIndex, iDim, jDim, varDim, var);
	    			Cstate[index] = iTemp[iIndex] * bgStdDev[index];
	  			}
				}
      }
    }  // end var for loop
  } // end else
  GPTLstop("CostFunction3D::SCtransform");
}

void CostFunction3D::SCtranspose(const real* Cstate, real* Astate)
{
  GPTLstart("CostFunction3D::SCtranspose");
  if ((iFilterScale < 0) and (jFilterScale < 0) and (kFilterScale < 0)) {
		//#pragma omp parallel for //[6]
		//#pragma acc parallel loop //[6]
    for (int n = 0; n < nState; n++) {
      Astate[n]= Cstate[n] * bgStdDev[n];
    }
  } else {
    // Isotropic Recursive filter, no anisotropic "triad" working yet
		//#pragma acc parallel loop //[7]
    for (int n = 0; n < nState; n++) {
      Astate[n]= Cstate[n];
    }
    //_#pragma omp parallel for
    for (int var = 0; var < varDim; var++) {

      // These are local for parallelization
      real* iTemp = new real[iDim];
      real* jTemp = new real[jDim];
      real* kTemp = new real[kDim];
      real* iPad = NULL;
      real* jPad = NULL;
      real* kPad = NULL;
      real* iq = NULL;
      real* is = NULL;
      real* jq = NULL;
      real* js = NULL;
      real* kq = NULL;
      real* ks = NULL;
      //if (iBCL[var] == PERIODIC) {
      iPad = new real[iDim*3];
      iq   = new real[iDim*3];
      is   = new real[iDim*3];
      //}
      //if (jBCL[var] == PERIODIC) {
      jPad = new real[jDim*3];
      jq   = new real[jDim*3];
      js   = new real[jDim*3];
      //}
      //if (kBCL[var] == PERIODIC) {
      kPad = new real[kDim*3];
      kq   = new real[kDim*3];
      ks   = new real[kDim*3];
      //}

      //FI & D
      for (int jIndex = 0; jIndex < jDim; jIndex++) {
				for (int kIndex = 0; kIndex < kDim; kIndex++) {
	  			for (int iIndex = 0; iIndex < iDim; iIndex++) {
	  				int cIndex = INDEX(iIndex, jIndex, kIndex, iDim, jDim, varDim, var);
	    			iTemp[iIndex] = Astate[cIndex] * bgStdDev[cIndex];
	  			}
	  			if (iFilterScale > 0) {
	    			if (iBCL[var] == PERIODIC) {
	      			// Pad the array to account for periodicity
	      			for (int iIndex = 0; iIndex < iDim; iIndex++) {
								iPad[iIndex] = iTemp[iIndex];
								iPad[iIndex+iDim] = iTemp[iIndex];
								iPad[iIndex+iDim*2] = iTemp[iIndex];
	      			}
	      			iFilter->filterArray(iPad, iq,is,iDim*3);
	      			for (int iIndex = 0; iIndex < iDim; iIndex++) {
								iTemp[iIndex] = iPad[iIndex+iDim];
	      			}
	    			} else {
	      			for (int iIndex = 0; iIndex < iDim; iIndex++) {
								iPad[iIndex] = 0.0;
								iPad[iIndex+iDim] = iTemp[iIndex];
								iPad[iIndex+iDim*2] = 0.0;
	      			}
	      			iFilter->filterArray(iPad, iq,is, iDim*3);
	      			for (int iIndex = 0; iIndex < iDim; iIndex++) {
									iTemp[iIndex] = iPad[iIndex+iDim];
	      			}
	    			}
	  			}
	  			for (int iIndex = 0; iIndex < iDim; iIndex++) {
	    			Astate[INDEX(iIndex, jIndex, kIndex, iDim, jDim, varDim, var)] = iTemp[iIndex];
	  			}
				}
      }
      //FJ
      for (int kIndex = 0; kIndex < kDim; kIndex++) {
				for (int iIndex = 0; iIndex < iDim; iIndex++) {
	  			for (int jIndex = 0; jIndex < jDim; jIndex++) {
	    			jTemp[jIndex] = Astate[INDEX(iIndex, jIndex, kIndex, iDim, jDim, varDim, var)];
	  			}
	  			if (jFilterScale > 0) {
	    			if (jBCL[var] == PERIODIC) {
	      			// Pad the array to account for periodicity
	      			for (int jIndex = 0; jIndex < jDim; jIndex++) {
								jPad[jIndex] = jTemp[jIndex];
								jPad[jIndex+jDim] = jTemp[jIndex];
								jPad[jIndex+jDim*2] = jTemp[jIndex];
	      			}
	      			jFilter->filterArray(jPad, jq,js,jDim*3);
	      			for (int jIndex = 0; jIndex < jDim; jIndex++) {
								jTemp[jIndex] = jPad[jIndex+jDim];
	      			}
	    			} else {
	      			for (int jIndex = 0; jIndex < jDim; jIndex++) {
								jPad[jIndex] = 0.0;
								jPad[jIndex+jDim] = jTemp[jIndex];
								jPad[jIndex+jDim*2] = 0.0;
	      			}
	      			jFilter->filterArray(jPad, jq,js, jDim*3);
	      			for (int jIndex = 0; jIndex < jDim; jIndex++) {
								jTemp[jIndex] = jPad[jIndex+jDim];
	      			}
	    			}
	  			}
	  			for (int jIndex = 0; jIndex < jDim; jIndex++) {
	    			Astate[INDEX(iIndex, jIndex, kIndex, iDim, jDim, varDim, var)] = jTemp[jIndex];
	  			}
				}
      }
      //FK
      for (int iIndex = 0; iIndex < iDim; iIndex++) {
				for (int jIndex = 0; jIndex < jDim; jIndex++) {
	  			for (int kIndex = 0; kIndex < kDim; kIndex++) {
	   				kTemp[kIndex] = Astate[INDEX(iIndex, jIndex, kIndex, iDim, jDim, varDim, var)];
	  			}
          if (kFilterScale > 0) {
	    			if (kBCL[var] == PERIODIC) {
	      			// Pad the array to account for periodicity
	      			for (int kIndex = 0; kIndex < kDim; kIndex++) {
								kPad[kIndex] = kTemp[kIndex];
								kPad[kIndex+kDim] = kTemp[kIndex];
								kPad[kIndex+kDim*2] = kTemp[kIndex];
	      			}
	      			kFilter->filterArray(kPad, kq,ks,kDim*3);
	      			for (int kIndex = 0; kIndex < kDim; kIndex++) {
								kTemp[kIndex] = kPad[kIndex+kDim];
	      			}
	    			} else {
	      			for (int kIndex = 0; kIndex < kDim; kIndex++) {
								kPad[kIndex] = 0.0;
								kPad[kIndex+kDim] = kTemp[kIndex];
								kPad[kIndex+kDim*2] = 0.0;
	      			}
	      			kFilter->filterArray(kPad, kq,ks, kDim*3);
	      			for (int kIndex = 0; kIndex < kDim; kIndex++) {
								kTemp[kIndex] = kPad[kIndex+kDim];
	      			}
	    			}
          }
          for (int kIndex = 0; kIndex < kDim; kIndex++) {
	    			Astate[INDEX(iIndex, jIndex, kIndex, iDim, jDim, varDim, var)] = kTemp[kIndex];
	  			}
				}
      }
      delete[] iTemp;
      delete[] jTemp;
      delete[] kTemp;
      //if (iBCL[var] == PERIODIC) {
      delete[] iPad;
      //}
      //if (jBCL[var] == PERIODIC) {
      delete[] jPad;
      //}
      //if (kBCL[var] == PERIODIC) {
      delete[] kPad;
      //}

    }
  }
  GPTLstop("CostFunction3D::SCtranspose");
}

bool CostFunction3D::setupSplines()
{
  int var,k;
  int index;
  GPTLstart("CostFunction3D::setupSplines");

  // Do the spline via a Cholesky decomposition
  // and manipulate the DC Filter
  real Pi = acos(-1.);

  real cutoff_wl = std::stof((*configHash)["i_spline_cutoff"]);
  cout << "i Spline cutoff set to " << cutoff_wl << endl;
  real eq = pow( (cutoff_wl/(2*Pi)) , 6);
  calcSplineCoefficients(iDim, eq, iBCL, iBCR, iMin, DI, DIrecip, iLDim, iL, iGamma);

  cutoff_wl = std::stof((*configHash)["j_spline_cutoff"]);
  cout << "j Spline cutoff set to " << cutoff_wl << endl;
  eq = pow( (cutoff_wl/(2*Pi)) , 6);
  calcSplineCoefficients(jDim, eq, jBCL, jBCR, jMin, DJ, DJrecip, jLDim, jL, jGamma);

  cutoff_wl = std::stof((*configHash)["k_spline_cutoff"]);
  cout << "k Spline cutoff set to " << cutoff_wl << endl;
  eq = pow( (cutoff_wl/(2*Pi)) , 6);
  calcSplineCoefficients(kDim, eq, kBCL, kBCR, kMin, DK, DKrecip, kLDim, kL, kGamma);

  cout << "Before linearlization of kGammaL array\n";
  for (var=0;var<varDim;var++) {
    // Linearlize the kGamma array 
    for (k=0;k<kRankMax*kDim;k++) {
     index = KINDEX(k,kRankMax*kDim,var);
     kGammaL[index]=kGamma[var][k];
    }
    // Linearlize the kL array 
    for (k=0;k<kRankMax*kLDim;k++) {
     index = KINDEX(k,kRankMax*kLDim,var);
     kLL[index]=kL[var][k];
    }
  }
  GPTLstop("CostFunction3D::setupSplines");
  return true;

}


void CostFunction3D::obAdjustments() {
  GPTLstart("CostFunction3D::obAdjustments");

  // Load the obs locally and weight the nonlinear observation operators by interpolated bg fields
  for (int m = 0; m < mObs; m++) {
    int mi = m*(7+varDim*derivDim);
    for (int ob = 0; ob < (7+varDim*derivDim); ob++) {
      obsVector[mi+ob] = rawObs[mi+ob];
    }
    real type = obsVector[mi+5];
    if (type <= 1) continue;

    real i = obsVector[mi+2];
    real j = obsVector[mi+3];
    real k = obsVector[mi+4];

    // Double check to make sure obs are in the domain
    if ((i < iMin) or (i > iMax)
	or (j < jMin) or (j > jMax)
	or (k < kMin) or (k > kMax)) {
      cout << "Error! Observations are found outside the domain where the spline is undefined.\n";
      cout << "This can only happen if you bypassed preprocessing -- check your samurai_Observations.in and re-run.\n";
    }
    real rhoprime = 0.;
    real qvprime = 0.;

    int ii = (int)((i - iMin)*DIrecip);
    int jj = (int)((j - jMin)*DJrecip);
    int kk = (int)((k - kMin)*DKrecip);
    real ibasis = 0.;
    real jbasis = 0.;
    real kbasis = 0.;
    for (int iiNode = (ii-1); iiNode <= (ii+2); ++iiNode) {
      int iNode = iiNode;
      if ((iNode < 0) or (iNode >= iDim)) continue;
      for (int jjNode = (jj-1); jjNode <= (jj+2); ++jjNode) {
	int jNode = jjNode;
	if ((jNode < 0) or (jNode >= jDim)) continue;
	for (int kkNode = (kk-1); kkNode <= (kk+2); ++kkNode) {
	  int kNode = kkNode;
	  if ((kNode < 0) or (kNode >= kDim)) continue;
	  ibasis = Basis(iNode, i, iDim-1, iMin, DI, DIrecip, 0, iBCL[4], iBCR[4]);
	  jbasis = Basis(jNode, j, jDim-1, jMin, DJ, DJrecip, 0, jBCL[4], jBCR[4]);
	  kbasis = Basis(kNode, k, kDim-1, kMin, DK, DKrecip, 0, kBCL[4], kBCR[4]);
	  qvprime += bgState[varDim*iDim*jDim*kNode + varDim*iDim*jNode +varDim*iNode + 4] * ibasis * jbasis * kbasis;
	  ibasis = Basis(iNode, i, iDim-1, iMin, DI, DIrecip, 0, iBCL[5], iBCR[5]);
	  jbasis = Basis(jNode, j, jDim-1, jMin, DJ, DJrecip, 0, jBCL[5], jBCR[5]);
	  kbasis = Basis(kNode, k, kDim-1, kMin, DK, DKrecip, 0, kBCL[5], kBCR[5]);
	  rhoprime += bgState[varDim*iDim*jDim*kNode + varDim*iDim*jNode +varDim*iNode + 5] * ibasis * jbasis * kbasis;
	}
      }
    }
    real heightm = 1000*k;
    real rhoBar = refstate->getReferenceVariable(ReferenceVariable::rhoaref, heightm);
    real qBar = refstate->getReferenceVariable(ReferenceVariable::qvbhypref, heightm);
    real qv = refstate->bhypInvTransform(qBar + qvprime);
    real rhoa = rhoBar + rhoprime / 100;
    real rhoBG = rhoa + rhoa*qv/1000.;

    // Adjust the relevant fields
    if (type == MetObs::sfmr) {
      obsVector[mi] *= rhoBG;
      // Still need to think carefully about nonlinear SFMR operator in Cartesian space
      //real wsBG = vBG*vBG; // + uBG*uBG;
      //obsVector[mi+2] = 2.*vBG/wsBG;
      //obsVector[mi+3] = 0.; //2.*uBG/wsBG;
    }
    if ((type == MetObs::radar) or (type == MetObs::qscat)
	or (type == MetObs::ascat) or (type == MetObs::AMV)
	or (type == MetObs::lidar)) {
      obsVector[mi] *= rhoBG;
    }
  }
  GPTLstop("CostFunction3D::obAdjustments");
}

void CostFunction3D::fillBasisLookup()
{

  real ONESIXTH = 1./6.;
  for (int i=0; i < 2000000; i++) {
    real z = 2.0 - real(i)/1000000.;
    real b = (z*z*z) * ONESIXTH;
    z -= 1.0;
    if (z > 0)
      b -= (z*z*z) * 4 * ONESIXTH;
    basis0[i] = b;

    z = 2.0 - real(i)/1000000.;
    b = (z*z) * ONESIXTH;
    z -= 1.0;
    if (z > 0)
      b -= (z*z) * 4 * ONESIXTH;
    basis1[i] = b;
  }

}

real CostFunction3D::Basis(const int& m, const real& x, const int& M,const real& xmin,
			   const real& DX, const real& DXrecip, const int& derivative,
			   const int& BL, const int& BR, const real& lambda)
{
  real b = 0;
  real xm = xmin + (m * DX);
  real delta = (x - xm) * DXrecip;
  real z = fabs(delta);
  real ONESIXTH = 1./6.; real FOURSIXTH = 4./6.;
  if (z < 2.0) {
    switch (derivative) {
    case 0:
      if (basisappx == NONE) {
	// Unapproximated
	z = 2.0 - z;
	b = (z*z*z) * ONESIXTH;
	z -= 1.0;
	if (z > 0)
	  b -= (z*z*z) * FOURSIXTH;
      } else if (basisappx == PARTIAL) {
	// Cheaper approximation
	real zi = z*1000000.;
	int z1 = int(zi);
	b = basis0[z1] + (basis0[z1+1]-basis0[z1])*(zi - z1);
      } else {
	// Cheapest approximation
	real zi = z*1000000.;
	int z1 = int(zi);
	b = basis0[z1];
      }
      break;
    case 1:
      if (basisappx == NONE) {
	z = 2.0 - z;
	b = (z*z) * ONESIXTH;
	z -= 1.0;
	if (z > 0)
	  b -= (z*z) * FOURSIXTH;
      } else if (basisappx == PARTIAL) {
	real zi = z*1000000.;
	int z1 = int(zi);
	b = basis1[z1] + (basis1[z1+1]-basis1[z1])*(zi - z1);

      } else {
	real zi = z*1000000.;
	int z1 = int(zi);
	b = basis1[z1];
      }
      b *= ((delta > 0) ? -1.0 : 1.0) * 3.0 * DXrecip;
      break;
    case 2:
      z = 2.0 - z;
      b = z;
      z -= 1.0;
      if (z > 0)
	b -= z * 4;
      b *= DXrecip * DXrecip;
      break;
    case 3:
      if (z > 1.0) {
	b = 1;
      } else if (z < 1.0) {
	b = -3.;
      }
      b *= ((delta > 0) ? -1.0 : 1.0) * DXrecip * DXrecip * DXrecip;
      break;
    }
  }
  /* if ((m > 2) and (m < M-2)) return b;
  // Add the boundary conditions if we get this far
  real bc = BasisBC(b, m, x, M, xmin, DX, DXrecip, derivative, BL, BR, lambda);
  return bc; */
  return b;
}

real CostFunction3D::BasisBC(real b, const int& m, const real& x, const int& M, const real& xmin,
                             const real& DX, const real& DXrecip, const int& derivative,
                             const int& BL, const int& BR, const real& lambda)
{
  real bmod = 0;
  int node = -2;
  real coeffmod = 0.;
  real ONESIXTH = 1./6.; real FOURSIXTH = 4./6.;
  if (m == 0) {
    // Left BC
    switch (BL) {
    case RX:
      // Absolutely no boundary condition
      return b;
    case R0:
      // No boundary condition, but buffered so use R1T2 on outer node
      node = -1;
      coeffmod = 2.;
      break;
    }
  } else if (m == 1) {
    switch (BL) {
    case RX:
      // Absolutely no boundary condition
      return b;
    case R0:
      // No boundary condition, but buffered so use R1T2 on outer node
      node = -1;
      coeffmod = -1.;
      break;
    case R1T0:
      node = 0;
      coeffmod = -4.;
      break;
    case R1T1:
      node = 0;
      coeffmod = 0.;
      break;
    case R1T2:
      node = 0;
      coeffmod = 2.;
      break;
    case R1T10:
      node = 0;
      coeffmod = -4./(3.*lambda + 1.);
      break;
    case R2T10:
      return 0;
    case R2T20:
      return 0;
    case R3:
      return 0;
    case PERIODIC:
      node = M-1;
      coeffmod = 1.;
      break;
    }
  } else if (m == 2) {
    // Left BC
    switch (BL) {
    case RX:
      return b;
    case R0:
      return b;
    case R1T0:
      node = 0;
      coeffmod = -1.;
      break;
    case R1T1:
      node = 0;
      coeffmod = 1.;
      break;
    case R1T2:
      node = 0;
      coeffmod = -1.;
      break;
    case R1T10:
      node = 0;
      coeffmod = (3.*lambda - 1.)/(3.*lambda + 1.);
      break;
    case R2T10:
      node = 0;
      coeffmod = 1.0;
    case R2T20:
      node = 0;
      coeffmod = -1.0;
    case R3:
      return 0;
    case PERIODIC:
      node = M;
      coeffmod = 1.;
      break;
    }
  } else if (m == M) {
    // Right BC
    switch (BR) {
    case RX:
      return b;
    case R0:
      node = M+1;
      coeffmod = 2.;
      break;
    }
  } else if (m == (M-1)) {
    switch (BR) {
    case RX:
      return b;
    case R0:
      node = M+1;
      coeffmod = -1.;
      break;
    case R1T0:
      node = M;
      coeffmod = -4.;
      break;
    case R1T1:
      node = M;
      coeffmod = 0.;
      break;
    case R1T2:
      node = M;
      coeffmod = 2.;
      break;
    case R1T10:
      node = M;
      coeffmod = -4./(3.*lambda + 1.);
      break;
    case R2T10:
      return 0;
    case R2T20:
      return 0;
    case R3:
      return 0;
    case PERIODIC:
      return 0;
    }
  } else if (m == (M-2)) {
    // Right BC
    switch (BR) {
    case RX:
      return b;
    case R0:
      return b;
    case R1T0:
      node = M;
      coeffmod = -1.;
      break;
    case R1T1:
      node = M;
      coeffmod = 1.;
      break;
    case R1T2:
      node = M;
      coeffmod = -1.;
      break;
    case R1T10:
      node = M;
      coeffmod = (3.*lambda - 1.)/(3.*lambda + 1.);
      break;
    case R2T10:
      node = M;
      coeffmod = 1.;
      break;
    case R2T20:
      node = M;
      coeffmod = -1.;
      break;
    case R3:
      return 0;
    case PERIODIC:
      node = 0;
      coeffmod = 1.;
      break;
    }
  }

  real xm = xmin + (node * DX);
  real delta = (x - xm) * DXrecip;
  real z = fabs(delta);
  if (z < 2.0) {
    switch (derivative) {
    case 0:
      if (basisappx == NONE) {
	// Unapproximated
	z = 2.0 - z;
	bmod = (z*z*z) * ONESIXTH;
	z -= 1.0;
	if (z > 0)
	  bmod -= (z*z*z) * FOURSIXTH;
      } else if (basisappx == PARTIAL) {
	// Cheaper approximation
	real zi = z*1000000.;
	int z1 = int(zi);
	bmod = basis0[z1] + (basis0[z1+1]-basis0[z1])*(zi - z1);
      } else {
	// Cheapest approximation
	real zi = z*1000000.;
	int z1 = int(zi);
	bmod = basis0[z1];
      }
      break;
    case 1:
      if (basisappx == NONE) {
	z = 2.0 - z;
	bmod = (z*z) * ONESIXTH;
	z -= 1.0;
	if (z > 0)
	  bmod -= (z*z) * FOURSIXTH;
      } else if (basisappx == PARTIAL) {
	real zi = z*1000000.;
	int z1 = int(zi);
	bmod = basis1[z1] + (basis1[z1+1]-basis1[z1])*(zi - z1);

      } else {
	real zi = z*1000000.;
	int z1 = int(zi);
	bmod = basis1[z1];
      }
      bmod *= ((delta > 0) ? -1.0 : 1.0) * 3.0 * DXrecip;
      break;
    case 2:
      z = 2.0 - z;
      bmod = z;
      z -= 1.0;
      if (z > 0)
	bmod -= z * 4;
      bmod *= DXrecip * DXrecip;
      break;
    case 3:
      if (z > 1.0) {
	bmod = 1;
      } else if (z < 1.0) {
	bmod = -3.;
      }
      bmod *= ((delta > 0) ? -1.0 : 1.0) * DXrecip * DXrecip * DXrecip;
      break;
    }
  }
  b += coeffmod * bmod;

  // R2 needs one more addition
  if ((BL == R2T10) and (m == 2)) {
    node = 1;
    coeffmod = -0.5;
  } else if ((BR == R2T10) and (m == M-2)) {
    node = M-1;
    coeffmod = -0.5;
  } else {
    return b;
  }

  xm = xmin + (node * DX);
  delta = (x - xm) * DXrecip;
  z = fabs(delta);
  if (z < 2.0) {
    switch (derivative) {
    case 0:
      if (basisappx == NONE) {
	// Unapproximated
	z = 2.0 - z;
	bmod = (z*z*z) * ONESIXTH;
	z -= 1.0;
	if (z > 0)
	  bmod -= (z*z*z) * FOURSIXTH;
      } else if (basisappx == PARTIAL) {
	// Cheaper approximation
	real zi = z*1000000.;
	int z1 = int(zi);
	bmod = basis0[z1] + (basis0[z1+1]-basis0[z1])*(zi - z1);
      } else {
	// Cheapest approximation
	real zi = z*1000000.;
	int z1 = int(zi);
	bmod = basis0[z1];
      }
      break;
    case 1:
      if (basisappx == NONE) {
	z = 2.0 - z;
	bmod = (z*z) * ONESIXTH;
	z -= 1.0;
	if (z > 0)
	  bmod -= (z*z) * FOURSIXTH;
      } else if (basisappx == PARTIAL) {
	real zi = z*1000000.;
	int z1 = int(zi);
	bmod = basis1[z1] + (basis1[z1+1]-basis1[z1])*(zi - z1);

      } else {
	real zi = z*1000000.;
	int z1 = int(zi);
	bmod = basis1[z1];
      }
      bmod *= ((delta > 0) ? -1.0 : 1.0) * 3.0 * DXrecip;
      break;
    case 2:
      z = 2.0 - z;
      bmod = z;
      z -= 1.0;
      if (z > 0)
	bmod -= z * 4;
      bmod *= DXrecip * DXrecip;
      break;
    case 3:
      if (z > 1.0) {
	bmod = 1;
      } else if (z < 1.0) {
	bmod = -3.;
      }
      bmod *= ((delta > 0) ? -1.0 : 1.0) * DXrecip * DXrecip * DXrecip;
      break;
    }
  }
  b += coeffmod * bmod;

  return b;

}

void CostFunction3D::adjustInternalDomain(int increment)
{

  iMin -= DI * increment;
  iMax += DI * increment;
  iDim += 2 * increment;
  jMin -= DJ * increment;
  jMax += DJ * increment;
  jDim += 2 * increment;;
  kMin -= DK * increment;
  kMax += DK * increment;
  kDim += 2 * increment;

}

void CostFunction3D::calcSplineCoefficients(const int& Dim, const real& eq, const int* BCL, const int* BCR,
                                            const real& xMin, const real& DX, const real& DXrecip, const int& LDim,
                                            real* L[7], real* gamma[7])
{

  for (int var = 0; var < varDim; var++) {
    int pDim = Dim;
    int mDim = Dim - rankHash[BCL[var]] - rankHash[BCR[var]];
    real pMin = xMin;

    // Subtract one for the periodic case rank mismatch
    if (BCL[var] == PERIODIC) mDim--;
    for (int i = 0; i < mDim*LDim; i++) {
      L[var][i] = 0;
    }

    // Allocate memory and initialize matrices
    real** P = new real*[mDim];
    real* p = new real[mDim];
    for (int i = 0; i < mDim; i++) {
      P[i] = new real[mDim];
      p[i] = 0.;
    }

    real** PP = new real*[pDim];
    real** tmp = new real*[pDim];
    real** G = new real*[mDim];
    real** GT = new real*[pDim];
    for (int i = 0; i < pDim; i++) {
      PP[i] = new real[pDim];
      tmp[i] = new real[mDim];
      GT[i] = new real[mDim];
      for (int j = 0; j < mDim; j++) {
	tmp[i][j] = 0;
	GT[i][j] = 0;
      }
      for (int j = 0; j < pDim; j++) {
	PP[i][j] = 0;
      }
    }

    for (int i = 0; i < mDim; i++) {
      G[i] = new real[pDim];
      for (int j = 0; j < mDim; j++) {
	P[i][j] = 0;
      }
      for (int j = 0; j < pDim; j++) {
	G[i][j] = 0;
      }
    }

    // Set boundary conditions
    switch (BCL[var]) {
      /* case R0: no BC enforced */
    case R1T0:
      G[0][0] = -4.0;
      G[1][0] = -1.0;
      break;
    case R1T1:
      G[1][0] =  1.0;
      break;
    case R1T2:
      G[0][0] =  2.0;
      G[1][0] = -1.0;
      break;
    case R2T10:
      G[0][0] =  1.0;
      G[0][1] = -0.5;
      break;
    case R2T20:
      G[0][0] = -1.0;
      break;
    case PERIODIC:
      G[mDim-1][0] = 1;
      break;
    default:
      break;
    }
    switch (BCR[var]) {
      /* case R0: no BC enforced */
    case R1T0:
      G[mDim-1][pDim-1] = -4.0;
      G[mDim-2][pDim-1] = -1.0;
      break;
    case R1T1:
      G[mDim-2][pDim-1] =  1.0;
      break;
    case R1T2:
      G[mDim-1][pDim-1] =  2.0;
      G[mDim-2][pDim-1] = -1.0;
      break;
    case R2T10:
      G[mDim-1][pDim-1] =  1.0;
      G[mDim-1][pDim-2] = -0.5;
      break;
    case R2T20:
      G[mDim-1][pDim-1] = -1.0;
      break;
    case PERIODIC:
      G[0][pDim-2] = 1;
      G[1][pDim-1] = 1;
      break;
    default:
      break;
    }

    for (int i = 0; i < mDim; i++) {
      G[i][i+rankHash[BCL[var]]] = 1;
    }
    for (int i = 0; i < mDim; i++) {
      for (int j = 0; j < pDim; j++) {
	GT[j][i] = G[i][j];
	gamma[var][Dim*i + j] = G[i][j];
	//std::cout << G[i][j] << " ";
      } //std::cout << "\n";
    } //std::cout << "\n";

    /* for (int i = 0; i < mDim; i++) {
       for (int j = 0; j < Dim; j++) {
       std::cout << gamma[var][Dim*i + j] << " ";
       } std::cout << "\n";
       } std::cout << "\n"; */

    /* for (int i = 0; i < pDim; i++) {
       for (int j = 0; j < mDim; j++) {
       //std::cout << GT[i][j] << " ";
       //std::cout << gamma[var][pDim*j + i] << " ";
       } //std::cout << "\n";
       } //std::cout << "\n"; */

    for (int Index = min(rankHash[BCL[var]],1); Index < max(pDim-1-rankHash[BCR[var]],pDim-2); Index++) {
      //for (int Index = 1; Index < pDim-2; Index++) {
      for (int mu = -1; mu <= 1; mu += 2) {
	real i = pMin + DX * (Index + (0.5*sqrt(1./3.) * mu + 0.5));
	int ii = (int)((i - pMin)*DXrecip);
	for (int Node = max(ii-1,0); Node <= min(ii+2,pDim-1); ++Node) {
	  real pm = Basis(Node, i, mDim-1, pMin, DX, DXrecip, 0, RX, RX);
	  real qm = Basis(Node, i, mDim-1, pMin, DX, DXrecip, 3, RX, RX);
	  real pn, qn;
	  PP[Node][Node] += 0.5 * ((pm * pm) + eq * (qm * qm));
	  if ((Node+1) < pDim) {
	    pn = Basis(Node+1, i, mDim-1, pMin, DX, DXrecip, 0, RX, RX);
	    qn = Basis(Node+1, i, mDim-1, pMin, DX, DXrecip, 3, RX, RX);
	    PP[Node][Node+1] += 0.5 * ((pm * pn) + eq * (qm * qn));
	    PP[Node+1][Node] += 0.5 * ((pm * pn) + eq * (qm * qn));
	  }
	  if ((Node+2) < pDim) {
	    pn = Basis(Node+2, i, mDim-1, pMin, DX, DXrecip, 0, RX, RX);
	    qn = Basis(Node+2, i, mDim-1, pMin, DX, DXrecip, 3, RX, RX);
	    PP[Node][Node+2] += 0.5 * ((pm * pn) + eq * (qm * qn));
	    PP[Node+2][Node] += 0.5 * ((pm * pn) + eq * (qm * qn));
	  }
	  if ((Node+3) < pDim) {
	    pn = Basis(Node+3, i, mDim-1, pMin, DX, DXrecip, 0, RX, RX);
	    qn = Basis(Node+3, i, mDim-1, pMin, DX, DXrecip, 3, RX, RX);
	    PP[Node][Node+3] += 0.5 * ((pm * pn) + eq * (qm * qn));
	    PP[Node+3][Node] += 0.5 * ((pm * pn) + eq * (qm * qn));
	  }
	}
      }
    }

    /* for (int i = 0; i < pDim; i++) {
       for (int j = 0; j < pDim; j++) {
       std::cout << PP[i][j] << " ";
       } std::cout << std::endl;
       }
       std::cout << std::endl; */

    for (int i = 0; i < pDim; i++) {
      for (int j = 0; j < mDim; j++) {
	//std::cout << PP[i][j] << " ";
	for (int k = 0; k < pDim; k++) {
	  tmp[i][j] += PP[i][k]*GT[k][j]; //PP[i][k]*gamma[var][pDim*j + k]
	}
      } //std::cout << std::endl;
    }
    //std::cout << std::endl;
    for (int i = 0; i < mDim; i++) {
      for (int j = 0; j < mDim; j++) {
	//std::cout << G[i][j] << " ";
	//P[i][j] = 0;
	for (int k = 0; k < pDim; k++) {
	  P[i][j] += G[i][k]*tmp[k][j]; //gamma[var][pDim*i + k]*tmp[k][j]
	}
	//std::cout << P[i][j] << " ";
      } //std::cout << std::endl;
    }
    //std::cout << std::endl;


    /* for (int i = 0; i < mDim; i++) {
       for (int j = 0; j < mDim; j++) {
       std::cout << P[i][j] << " ";
       } std::cout << std::endl;
       }
       std::cout << std::endl; */

    // Cholesky decomp of P+Q
    for (int i=0;i<mDim;i++) {
      for (int j=i;j<mDim;j++) {
	real sum=P[i][j];
	for (int k=i-1;k>=0;k--) {
	  sum -= P[i][k]*P[j][k];
	}
	if (i == j) {
	  if (sum <= 0.0) {
	    std::cout << "cholesky failed at i,j sum\n";
	    exit(1);
	    break;
	  } else {
	    p[i] = sqrt(sum);
	  }
	} else {
	  P[j][i]=sum/p[i];
	  if (p[i] == 0.) {
	    std::cout << "Problem! " << i << "\t" << j << "\n";
	    exit(1);
	  }
	}
      }
    }

    // Reduced representation of decomposed P+Q
    for (int i = 0; i < mDim; i++) {
      L[var][i*LDim] = p[i];
      //std::cout << L[var][i*LDim] << " ";
      for (int n=1;n<LDim;n++) {
	if ((i-n) >= 0) {
	  L[var][i*LDim+n] = P[i][i-n];
	}
	//std::cout << L[var][i*LDim + n] << " ";
      } //std::cout << endl;
    } //std::cout << endl;

    // Free memory
    for (int i = 0; i < pDim; i++) {
      delete[] PP[i];
      delete[] tmp[i];
      delete[] GT[i];
    }
    for (int i = 0; i < mDim; i++) {
      delete[] G[i];
    }
    delete[] PP;
    delete[] tmp;
    delete[] G;
    delete[] GT;

    for (int i = 0; i < mDim; i++) {
      delete[] P[i];
    }
    delete[] P;
    delete[] p;
  }
}

void CostFunction3D::FFtransform(const real* Astate, real* Cstate)
{
  GPTLstart("CostFunction3D::FFtransform");
  // Copy to the new state in case no FFT is enforced
  #pragma omp parallel for //[8]
  for (int n = 0; n < nState; n++) {
    Cstate[n]= Astate[n];
  }

  for (int var = 0; var < varDim; var++) {
    // Enforce max wavenumber
    if ((kBCL[var] == PERIODIC) and (kMaxWavenumber[var] >= 0)) {
      for (int iIndex = 0; iIndex < iDim; iIndex++) {
        for (int jIndex = 0; jIndex < jDim; jIndex++) {
          for (int kIndex = 0; kIndex < kDim; kIndex++) {
            kFFTin[kIndex] = Cstate[varDim*iDim*jDim*kIndex + varDim*iDim*jIndex + varDim*iIndex + var];
          }
          fftw_execute(kForward);
          for (int kIndex = kMaxWavenumber[var]+1; kIndex < (kDim/2)+1; kIndex++) {
            kFFTout[kIndex][0] = 0.0;
            kFFTout[kIndex][1] = 0.0;
          }
          fftw_execute(kBackward);
          for (int kIndex = 0; kIndex < kDim; kIndex++) {
            Cstate[varDim*iDim*jDim*kIndex + varDim*iDim*jIndex +varDim*iIndex + var] = kFFTin[kIndex]/kDim;
          }
        }
      }
    }

    // Enforce max wavenumber
    if ((jBCL[var] == PERIODIC) and (jMaxWavenumber[var] >= 0)) {
      for (int kIndex = 0; kIndex < kDim; kIndex++) {
	for (int iIndex = 0; iIndex < iDim; iIndex++) {
	  for (int jIndex = 0; jIndex < jDim; jIndex++) {
	    jFFTin[jIndex] = Cstate[varDim*iDim*jDim*kIndex + varDim*iDim*jIndex + varDim*iIndex + var];
	  }
	  fftw_execute(jForward);
	  for (int jIndex = jMaxWavenumber[var]+1; jIndex < (jDim/2)+1; jIndex++) {
	    jFFTout[jIndex][0] = 0.0;
	    jFFTout[jIndex][1] = 0.0;
	  }
	  fftw_execute(jBackward);
	  for (int jIndex = 0; jIndex < jDim; jIndex++) {
	    Cstate[varDim*iDim*jDim*kIndex + varDim*iDim*jIndex +varDim*iIndex + var] = jFFTin[jIndex]/jDim;
	  }
	}
      }
    }
    if ((iBCL[var] == PERIODIC) and (iMaxWavenumber[var] >= 0)) {
      for (int kIndex = 0; kIndex < kDim; kIndex++) {
	for (int jIndex = 0; jIndex < jDim; jIndex++) {
	  for (int iIndex = 0; iIndex < iDim; iIndex++) {
	    iFFTin[iIndex] = Cstate[varDim*iDim*jDim*kIndex + varDim*iDim*jIndex + varDim*iIndex + var];
	  }
	  fftw_execute(iForward);
	  for (int iIndex = iMaxWavenumber[var]+1; iIndex < (iDim/2)+1; iIndex++) {
	    iFFTout[iIndex][0] = 0.0;
	    iFFTout[iIndex][1] = 0.0;
	  }
	  fftw_execute(iBackward);
	  for (int iIndex = 0; iIndex < iDim; iIndex++) {
	    Cstate[varDim*iDim*jDim*kIndex + varDim*iDim*jIndex +varDim*iIndex + var] = iFFTin[iIndex]/iDim;
	  }
	}
      }
    }
  }
  GPTLstop("CostFunction3D::FFtransform");
}

void CostFunction3D::calcHmatrix()
{
  GPTLstart("CostFunction3D::calcHmatrix");

  std::cout << "Build H transform matrix...\n";

  // Allocate memory for the sparse matrix
  // 193 elements for columns assumes a maximum of 3 variables
  // will be used in an observation.
  real **Hbuild = new real*[mObs];
  real **JHbuild = new real*[mObs];
  for (int m = 0; m < mObs; m++) {
    Hbuild[m] = new real[193];
    JHbuild[m] = new real[193];
  }
  IH = new int[mObs+1];

  // H
#pragma omp parallel for
  for (int m = 0; m < mObs; m++) {
    //real* Hrow = new real[nState];
    //for (int i = 0; i < nState; i++) { Hrow[i] = 0.0; }
    int mi = m*(7+varDim*derivDim);
    real i = obsVector[mi+2];
    real j = obsVector[mi+3];
    real k = obsVector[mi+4];
    int ii = (int)((i - iMin)*DIrecip);
    int jj = (int)((j - jMin)*DJrecip);
    int kk = (int)((k - kMin)*DKrecip);
    real ibasis = 0;
    real jbasis = 0;
    real kbasis = 0;
    int hi = 1;
    for (int var = 0; var < varDim; var++) {
      for (int d = 0; d < derivDim; d++) {
        int wgt_index = mi + (7*(d+1)) + var;
        if (!obsVector[wgt_index]) continue;
        for (int iiNode = (ii-1); iiNode <= (ii+2); ++iiNode) {
          int iNode = iiNode;
          if ((iNode < 0) or (iNode >= iDim)) continue;
          ibasis = Basis(iNode, i, iDim-1, iMin, DI, DIrecip, derivative[d][0], iBCL[var], iBCR[var]);
          if (!ibasis) continue;
          for (int jjNode = (jj-1); jjNode <= (jj+2); ++jjNode) {
            int jNode = jjNode;
            if ((jNode < 0) or (jNode >= jDim)) continue;
            jbasis = Basis(jNode, j, jDim-1, jMin, DJ, DJrecip, derivative[d][1], jBCL[var], jBCR[var]);
            if (!jbasis) continue;
            for (int kkNode = (kk-1); kkNode <= (kk+2); ++kkNode) {
              int kNode = kkNode;
              if ((kNode < 0) or (kNode >= kDim)) continue;
              kbasis = Basis(kNode, k, kDim-1, kMin, DK, DKrecip, derivative[d][2], kBCL[var], kBCR[var]);
              if (!kbasis) continue;
              int cIndex = varDim*iDim*jDim*kNode + varDim*iDim*jNode +varDim*iNode + var;
              real weight = ibasis * jbasis * kbasis * obsVector[wgt_index];
              Hbuild[m][hi] = weight;
              JHbuild[m][hi] = cIndex;
              hi++;
            }
          }
        }
      }
    }
    Hbuild[m][0] = hi;
    if (hi > 193) {
      cout << "Overflow in H matrix calculation!" << hi << "\n";
    }
  }

  int nonzeros = 0;
  for (int m = 0; m < mObs; m++) {
    nonzeros += Hbuild[m][0] - 1;
  }

  IH[mObs] = nonzeros;
  //std::cout << "Non-zero entries in sparse H matrix: " << nonzeros << " = " << 100.0*float(nonzeros)/(float(mObs)*float(nState)) << " %\n";

  H = new real[nonzeros];
  JH = new int[nonzeros];

  int hi = 0;
  for (int m = 0; m < mObs; m++) {
    IH[m] = hi;
    for (int n = 1; n < Hbuild[m][0]; n++) {
      H[hi] = Hbuild[m][n];
      JH[hi] = JHbuild[m][n];
      hi++;
    }
  }

  for (int m = 0; m < mObs; m++) {
    delete[] Hbuild[m];
    delete[] JHbuild[m];
  }
  delete[] Hbuild;
  delete[] JHbuild;

  GPTLstop("CostFunction3D::calcHmatrix");
}

void CostFunction3D::Htransform(const real* Cstate, real* Hstate)
{
  GPTLstart("CostFunction3D::Htransform");
  // Multiply the state by the observation matrix
#pragma omp parallel for
  for(int i=0; i<mObs; ++i) {
    Hstate[i] = 0.0;
    const int begin = IH[i];
    const int end = IH[i + 1];
    for(int j=begin; j<end; ++j) {
      Hstate[i] += H[j] * Cstate[JH[j]];
    }
  }

  GPTLstop("CostFunction3D::Htransform");
}

// Copy the final results into the given arrays
// Source is row major order (C)
// Dest is column major order (Fortran)

// bool CostFunction3D::copy3DArrayCtoF(real *src, float *dest, int iDim, int jDim, int kDim)
bool CostFunction3D::copy3DArray(real *src, float *dest, int iDim, int jDim, int kDim)
{
  for(int i = 0; i < iDim; i++)
    for(int j = 0; j < jDim; j++)
      for(int k = 0; k < kDim; k++) {
	float f = *(src + k + kDim * (j + jDim * i));
	*(dest + i + iDim * (j + jDim * k)) = f;
      }
  return true;
}

bool CostFunction3D::copyResults(int iDim, int jDim, int kDim,
				 float *u, float *v, float *w, float *t, float *p)
{
  // TODO do some error checking on the dimentions
  // if bad, return false
  
  // question: Do the result sizes match the nx, ny, nsigmas passed in the VarDriver3d::run function?

  // See CostFunctionCOAMPS.cpp around line  355 for the meaning of the indices used here (0, 1, 2, 19, 26)
  copy3DArray(&finalAnalysis[0], u, iDim, jDim, kDim);
  copy3DArray(&finalAnalysis[iDim * jDim * kDim * 1],  v, iDim, jDim, kDim);
  copy3DArray(&finalAnalysis[iDim * jDim * kDim * 2],  w, iDim, jDim, kDim);
  copy3DArray(&finalAnalysis[iDim * jDim * kDim * 19], p, iDim, jDim, kDim);
  copy3DArray(&finalAnalysis[iDim * jDim * kDim * 26], t, iDim, jDim, kDim);
	  
  return true;
}

