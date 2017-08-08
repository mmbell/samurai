/*
 *  ReferenceState.cpp
 *  samurai
 *
 *  Copyright 2012 Michael Bell. All rights reserved.
 *
 */

#include "ReferenceState.h"
#include <QString>
#include <cmath>
#include <iostream>
#include <fstream>
#include <QVector>
using namespace ReferenceVariable;

ReferenceState::ReferenceState(const QString& config)
{
	QVector<real> altitude, theta, qv, pi;
	real sfcpress;
	std::ifstream refstream(config.toLatin1().data());
	//QFile refFile(config);
	//if(!refFile.open(QIODevice::ReadOnly)) {
	if (!refstream.good()) {
	  std::cout << "Can't open Reference State file for reading, "
		    << " '" << config.toLatin1().data()
		    << "' using default..." << std::endl;
		sfcpress = 1014.80;
		altitude.push_back(log(10.0)); theta.push_back(298.6949); qv.push_back(bhypTransform(18.63960));
		altitude.push_back(log(124.0)); theta.push_back(299.6500); qv.push_back(bhypTransform(18.58188));
		altitude.push_back(log(810.0)); theta.push_back(301.688); qv.push_back(bhypTransform(15.30626));
		altitude.push_back(log(1541.0)); theta.push_back(304.5541); qv.push_back(bhypTransform(11.98349));
		altitude.push_back(log(3178.0)); theta.push_back(312.2750); qv.push_back(bhypTransform(6.76311));
		altitude.push_back(log(4437.0)); theta.push_back(317.8749); qv.push_back(bhypTransform(4.15019));
		altitude.push_back(log(5887.0)); theta.push_back(324.8602); qv.push_back(bhypTransform(2.42535));
		altitude.push_back(log(7596.0)); theta.push_back(332.5846); qv.push_back(bhypTransform(1.11535));
		altitude.push_back(log(9690.0)); theta.push_back(339.6121); qv.push_back(bhypTransform(0.32924));
		altitude.push_back(log(10949.0)); theta.push_back(342.8986); qv.push_back(bhypTransform(0.13712));
		altitude.push_back(log(12418.0)); theta.push_back(346.4510); qv.push_back(bhypTransform(0.04282));
		altitude.push_back(log(14203.0)); theta.push_back(353.9290); qv.push_back(bhypTransform(0.03));
		altitude.push_back(log(16590.0)); theta.push_back(383.2672); qv.push_back(bhypTransform(0.03));
		altitude.push_back(log(20726.0)); theta.push_back(494.1519); qv.push_back(bhypTransform(0.010));
		altitude.push_back(log(40000.0)); theta.push_back(1010.8810); qv.push_back(bhypTransform(0.001));
	} else {
		altitude.push_back(log(10.0));
		real altin, thetain, qvin, uin, vin;
		refstream >> sfcpress >> thetain >> qvin;
		theta.push_back(thetain);
		qv.push_back(qvin);
		while (refstream >> altin >> thetain >> qvin >> uin >> vin) {
			altitude.push_back(log(altin));
			theta.push_back(thetain);
			qv.push_back(bhypTransform(qvin));
		}
	}

	// Create the splines
	if (altitude.size() == 1) {
		std::cout << "Only one level found in reference spline setup. Please check reference file and re-run.\n";
	}
	thetaSpline = new SplineD(&altitude.front(), altitude.size(), theta.data(), 0, SplineBase::BC_ZERO_FIRST);
	qvSpline = new SplineD(&altitude.front(), altitude.size(), qv.data(), 0, SplineBase::BC_ZERO_FIRST);

	// Integrate the hydrostatic equation
	real finalalt = exp(altitude.last());
	altitude.clear();

	real gamma = 287.04/1005.7;
	real qvsfc = bhypInvTransform(qv.first())/1000.0;
	real pressa = sfcpress - sfcpress * (qvsfc / (qvsfc + 0.622));
	real pisfc = pow((pressa/1000.0),gamma);
	pi.push_back(pisfc);
	altitude.push_back(log(10.0));
	for (float i = 11.0; i<= finalalt; i++) {
		real height_u = log(i);
		real height_d = log(i-1);
		real theta_u = thetaSpline->evaluate(height_u);
		real theta_d = thetaSpline->evaluate(height_d);
		real newpi = pi.last() - 9.81/(1005.7 * 0.5*(theta_u + theta_d));
		pi.push_back(newpi);
		altitude.push_back(height_u);
	}
	// Unclear why this needs a cutoff wavelength, guessing it has to do with very fine height increments
	piSpline = new SplineD(&altitude.front(), altitude.size(), pi.data(), 2, SplineBase::BC_ZERO_SECOND);

}

ReferenceState::~ReferenceState()
{

	delete thetaSpline;
	delete qvSpline;
	delete piSpline;

}

/* Cubic B-spline from file defines the background reference state
	which is in hydrostatic balance. Dunion (2010) moist tropical sounding is default
	qvbhypref,
	rhoaref,
	rhoref,
	href,
	tempref,
	pressref */
real ReferenceState::getReferenceVariable(const int& refVariable, const real& heightm, const int& dz)
{
	real logheight = 0.0;
	real invheight = 0.0;
	if (heightm < 10.0) {
		logheight = log(10.0);
		invheight = 0.10;
	} else {
		logheight = log(heightm);
		invheight = 1.0/heightm;
	}
	if (refVariable == qvbhypref) {
		real qvbhyp = 0.;
		if (dz == 0) {
			qvbhyp = qvSpline->evaluate(logheight);
		} else {
			qvbhyp = qvSpline->slope(logheight)*invheight;
		}
		return qvbhyp;
	} else {
		real theta = 0.;
		real pi = 0.;
		real qv = 0.;
		pi = piSpline->evaluate(logheight);
		theta = thetaSpline->evaluate(logheight);
		qv = qvSpline->evaluate(logheight);
		qv = bhypInvTransform(qv)/1000.0;
		real temp = pi * theta;
		real pressa = 100000.0 * pow(pi, (1005.7/287.04));
		real rhoa = pressa/(temp*287.04);
		real rho = rhoa + qv * rhoa;
		real press = pressa * (1.0 + qv/0.622);
		real h = 1005.7*temp + 9.81*heightm + 2.5e6*qv;
		if (dz == 0) {
			switch (refVariable) {
					case rhoaref:
						return rhoa;
					case rhoref:
						return rho;
					case href:
						return h;
					case tempref:
						return temp;
					case pressref:
						return press;
					default:
						break;
			}
		} else {
			real dthetadz = thetaSpline->slope(logheight)*invheight;
			real qvdz = 0.002*qvSpline->slope(logheight)*invheight;
			real dpidz = piSpline->slope(logheight)*invheight;
			real dpdz = -rho * 9.81;
			real dtdz = pi*dthetadz +theta*dpidz;
			real dpadz = (dpdz - 0.622*pressa*qvdz)/(1.0 + 0.622*qv);
			real drhoadz = (1.0/287.04)*(dpadz/temp - pressa*dtdz/(temp*temp));
			real drhodz = (1.0/287.04)*(dpdz/temp - press*dtdz/(temp*temp));
			real dhdz = 1005.7*dtdz + 9.81 + 2.5e3*qvdz;
			switch (refVariable) {
					case rhoaref:
						return drhoadz;
					case rhoref:
						return drhodz;
					case href:
						return dhdz;
					case tempref:
						return dtdz;
					case pressref:
						return dpadz;
					default:
						break;
			}
		}
	}
	return 0;
}

/* Biased Hyperbolic transform for positive definite quanitity
 See Ooyama (2001) Journal of Atmospheric Sciences */

real ReferenceState::bhypTransform(const real& qv)
{

	real qvbhyp = 0.5*((qv + 1.e-7) - 1.e-14/(qv + 1.e-7));
	return qvbhyp;

}

/* Quasi-Inverse of Biased Hyperbolic transform for positive definite quanitity */

real ReferenceState::bhypInvTransform(const real& qvbhyp)
{
	real qv = 1.0e-6;
	if (qvbhyp > 0) {
		qv = sqrt(qvbhyp*qvbhyp + 1.e-14) + qvbhyp - 1.e-7;
	}
	return qv;
}
