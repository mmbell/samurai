/*
 *  ReferenceState.cpp
 *  samurai
 *
 *  Copyright 2012 Michael Bell. All rights reserved.
 *
 */

#include "ReferenceState.h"
#include <QString>
#include <QFile>
#include <QTextStream>
#include <cmath>
#include <iostream>
using namespace ReferenceVariable;

ReferenceState::ReferenceState(const QString& config)
{
	QFile refFile(config);
	if(!refFile.open(QIODevice::ReadOnly)) {
		std::cout << "Can't open Reference State file for reading, using default..." << std::endl;
		qvbhypcoeff[0] = 9.4826;
		qvbhypcoeff[1] = -0.0026721;
		qvbhypcoeff[2] = 2.8312e-07;
		qvbhypcoeff[3] = -1.3217e-11;
		qvbhypcoeff[4] = 2.2749e-16;
		
		rhoacoeff[0] = 1.1439;
		rhoacoeff[1] = -0.00010117;
		rhoacoeff[2] = 3.2486e-09;
		rhoacoeff[3] = -3.4898e-14;
		rhoacoeff[4] = -2.6925e-19;
		
		dpdzcoeff[0] = -11.432;
		dpdzcoeff[1] = 0.0010633;
		dpdzcoeff[2] = -4.0545e-08;
		dpdzcoeff[3] =  7.9634e-13;
		dpdzcoeff[4] = -5.8778e-18;

		sfcpress = 1015.10;
	} else {
		QTextStream in(&refFile);
		QString coeff;
		in >> coeff >> sfcpress;
		in >> coeff >> qvbhypcoeff[0] >> qvbhypcoeff[1] >> qvbhypcoeff[2] >> qvbhypcoeff[3] >> qvbhypcoeff[4];
		in >> coeff >> rhoacoeff[0] >> rhoacoeff[1] >> rhoacoeff[2] >> rhoacoeff[3] >> rhoacoeff[4];
		in >> coeff >> dpdzcoeff[0] >> dpdzcoeff[1] >> dpdzcoeff[2] >> dpdzcoeff[3] >> dpdzcoeff[4];
	}
}

ReferenceState::~ReferenceState()
{
}

/* Fourth order polynomial coefficients define the background reference state
	which is in hydrostatic balance. Dunion (2010) moist tropical sounding is only
	current implementation */

real ReferenceState::getReferenceVariable(const int& refVariable, const real& heightm, const int& dz)
{
	
	if (refVariable == qvbhypref) {
		real qvbhyp = 0.;
		for (int i = dz; i < 5; i++) {
			real power = pow(heightm, i-dz);
			if (dz) {
				qvbhyp += qvbhypcoeff[i] * power * i;
			} else {
				qvbhyp += qvbhypcoeff[i] * power;
			}
		}
		if (!dz) {
			if (qvbhyp < 0.) qvbhyp = 0.;
		}
		return qvbhyp;
	} else if (refVariable == rhoaref) {
		real rhoa = 0.;
		for (int i = dz; i < 5; i++) {
			real power = pow(heightm, i-dz); 
			if (dz) {
				rhoa += rhoacoeff[i] * power * i;
			} else {
				rhoa += rhoacoeff[i] * power;
			}
		}
		return rhoa;
	} else if (refVariable == rhoref) {
		real rho = 0.;
		real qvbhyp = 0.;
		real rhoa = 0.;
		for (int i = 0; i < 5; i++) {
			real power = pow(heightm, i); 
			rhoa += rhoacoeff[i] * power;
			qvbhyp += qvbhypcoeff[i] * power;
		}
		if (qvbhyp < 0.) qvbhyp = 0.;
		real qv = bhypInvTransform(qvbhyp);
		rho = rhoa*qv/1000. + rhoa;
		if (dz) {
			real rhoadz = 0.;
			real qvdz = 0.;
			for (int i = dz; i < 5; i++) {
				real power = pow(heightm, i-dz); 
				rhoadz += rhoacoeff[i] * power * i;
				qvdz += qvbhypcoeff[i] * power * i;
			}
			rho = rhoadz * (1 + qv/1000.) + rhoa * (qvdz/500.);
			return rho;
		}
		return rho;
	} else if ((refVariable == href) or (refVariable == tempref) or (refVariable == pressref)) {
		// Integrate hydrostatic equation to get pressure and/or solve for T or h
		real press = 0.;
		real temp = 0.;
		real qvbhyp = 0.;
		real rhoa = 0.;
		real dpdz = 0.;
		for (int i = 0; i < 5; i++) {
			real power = pow(heightm, i);
			real power1 = pow(heightm, i+1);
			press += dpdzcoeff[i] * power1 / (i+1);
			dpdz += dpdzcoeff[i] * power;
			rhoa += rhoacoeff[i] * power;
			qvbhyp += qvbhypcoeff[i] * power;
		}
		if (qvbhyp < 0.) qvbhyp = 0.;
		real qv = bhypInvTransform(qvbhyp);
		press += sfcpress*100.0;
		if (dz) {
			real rhoadz = 0.;
			real qvdz = 0.;
			for (int i = dz; i < 5; i++) {
				real power = pow(heightm, i-dz); 
				rhoadz += rhoacoeff[i] * power * i;
				qvdz += qvbhypcoeff[i] * power * i;
			}
			real alphadz = -(rhoadz*(1 + qv/1000.) + rhoa*qvdz/500.)/(rhoa*(1+qv/1000.)*rhoa*(1+qv/1000.));
			real dtdz = (press*alphadz + dpdz/(rhoa*(1+qv/1000.)))/286.9;
			//real alphadz = 1/(rhoadz * (286.9 + 461.5*qv/1000.) + 461.5 * rhoa * (qvdz/500.));
			//real dtdz = press*alphadz + dpdz/(286.9*rhoa + 461.5*rhoa*qv/1000.);
			real dhdz = 1005.7*dtdz + 9.81 + 2.5e3*qvdz;
			switch (refVariable) {
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
		temp = press/(286.9*rhoa + 461.5*rhoa*qv/1000.);
		real h = 1005.7*temp + 9.81*heightm + 2.5e3*qv;
		switch (refVariable) {
			case href:
				return h;
			case tempref:
				return temp;
			case pressref:
				return press;
			default:
				break;
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
	real qv = 0.;
	if (qvbhyp > 0) {
		qv = sqrt(qvbhyp*qvbhyp + 1.e-14) + qvbhyp - 1.e-7;
	}
	return qv;
}

