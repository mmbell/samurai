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
	std::ifstream refstream(config.toAscii().data());
	//QFile refFile(config);
	//if(!refFile.open(QIODevice::ReadOnly) {
	if (!refstream.good()) {
		std::cout << "Can't open Reference State file for reading, using default..." << std::endl;
		sfcpress = 1014.80;
		altitude.push_back(0.0); theta.push_back(298.6949); qv.push_back(bhypTransform(18.63960));
		altitude.push_back(0.124); theta.push_back(299.6500); qv.push_back(bhypTransform(18.58188));
		altitude.push_back(0.810); theta.push_back(301.688); qv.push_back(bhypTransform(15.30626));
		altitude.push_back(1.541); theta.push_back(304.5541); qv.push_back(bhypTransform(11.98349));
		altitude.push_back(3.178); theta.push_back(312.2750); qv.push_back(bhypTransform(6.76311));
		altitude.push_back(4.437); theta.push_back(317.8749); qv.push_back(bhypTransform(4.15019));
		altitude.push_back(5.887); theta.push_back(324.8602); qv.push_back(bhypTransform(2.42535));
		altitude.push_back(7.596); theta.push_back(332.5846); qv.push_back(bhypTransform(1.11535));
		altitude.push_back(9.690); theta.push_back(339.6121); qv.push_back(bhypTransform(0.32924));
		altitude.push_back(10.949); theta.push_back(342.8986); qv.push_back(bhypTransform(0.13712));
		altitude.push_back(12.418); theta.push_back(346.4510); qv.push_back(bhypTransform(0.04282));
		altitude.push_back(14.203); theta.push_back(353.9290); qv.push_back(bhypTransform(0.03));
		altitude.push_back(16.590); theta.push_back(383.2672); qv.push_back(bhypTransform(0.03));
		altitude.push_back(20.726); theta.push_back(494.1519); qv.push_back(bhypTransform(0.010));
		//altitude.push_back(40.000); theta.push_back(1010.8810); qv.push_back(bhypTransform(0.001));
	} else {
		altitude.push_back(0.0);
		real altin, thetain, qvin, uin, vin;
		refstream >> sfcpress >> thetain >> qvin;
		theta.push_back(thetain);
		qv.push_back(bhypTransform(qvin));
		while (refstream >> altin >> thetain >> qvin >> uin >> vin) {
			altitude.push_back(altin/1000.0);
			theta.push_back(thetain);
			qv.push_back(bhypTransform(qvin));
		}
	}

	// Create the splines
	if (altitude.size() == 1) {
		std::cout << "Only one level found in reference spline setup. Please check reference file and re-run.\n";
	}
	thetaSpline = new SplineD(&altitude.front(), altitude.size(), theta.data(), 2, SplineBase::BC_ZERO_FIRST);
	qvSpline = new SplineD(&altitude.front(), altitude.size(), qv.data(), 2, SplineBase::BC_ZERO_FIRST);

	// Integrate the hydrostatic equation
	real gamma = 287.04/1005.7;
	real reps = 461.5/287.04;
	real pisfc = pow((sfcpress/1000.0),gamma);
	pi.push_back(pisfc);
	for (float i = 1; i< altitude.size(); i++) {
		real height_u = altitude[i];
		real height_d = altitude[i-1];
		real theta_u = thetaSpline->evaluate(height_u);
		real theta_d = thetaSpline->evaluate(height_d);
		real qv_u = qvSpline->evaluate(height_u);
		qv_u = bhypInvTransform(qv_u)/1000.0;
		real qv_d = qvSpline->evaluate(height_d);
		qv_d = bhypInvTransform(qv_d)/1000.0;
		real thetav_u = theta_u * (1.0 + qv_u*reps)/(1.0+qv_u);
		real thetav_d = theta_d * (1.0 + qv_d*reps)/(1.0+qv_d);
		real newpi = pi.last() - 9.81*1000.0*(height_u - height_d)/(1005.7 * 0.5*(thetav_u + thetav_d));
		pi.push_back(newpi);
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
	real height = 0.0;
	real invheightm = 0.0;
	if (heightm < 10.0) {
		height = 0.01;
		invheightm = 0.1;
	} else {
		height = heightm/1000.0;
		invheightm = 1.0/heightm;
	}
	if (refVariable == qvbhypref) {
		real qvbhyp = 0.;
		if (dz == 0) {
			qvbhyp = qvSpline->evaluate(height);
		} else {
			qvbhyp = qvSpline->slope(height)*invheightm/1000.0;
		}
		return qvbhyp;
	} else {
		real theta = 0.;
		real pi = 0.;
		real qv = 0.;
		pi = piSpline->evaluate(height);
		theta = thetaSpline->evaluate(height);
		qv = qvSpline->evaluate(height);
		qv = bhypInvTransform(qv)/1000.0;
		real temp = pi * theta;
		real press = 100000.0 * pow(pi, (1005.7/287.04));
		real pressa = press - press * (qv / (qv + 0.622));
		real rhoa = pressa/(temp*287.04);
		real rho = rhoa + qv * rhoa;
		//real press = pressa * (1.0 + qv/0.622);
		real h = 1005.7*temp + 9.81*height/1000.0 + 2.5e6*qv;
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
			real dthetadz = thetaSpline->slope(height)/1000.0;
			real qvdz = 0.002*qvSpline->slope(height)/1000.0;
			real dpidz = piSpline->slope(height)/1000.0;
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
						return dpdz;
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
