/*
 *  ReferenceState.h
 *  samurai
 *
 *  Copyright 2012 Michael Bell. All rights reserved.
 *
 */

#ifndef REFSTATE_H
#define REFSTATE_H

#include "precision.h"
#include <QString>

class ReferenceState
{
    
public:
    ReferenceState(const QString& config);
    ~ReferenceState();
    real getReferenceVariable(const int& refVariable, const real& heightm, const int& dz = 0);
    real bhypTransform(const real& qv);
    real bhypInvTransform(const real& qvbhyp);

private:
    real qvbhypcoeff[5];
	real rhoacoeff[5];
	real dpdzcoeff[5];
};

namespace ReferenceVariable {
    
    enum variables {
        qvbhypref,
        rhoaref,
        rhoref,
        href,
        tempref,
        pressref
    };	
};

#endif
