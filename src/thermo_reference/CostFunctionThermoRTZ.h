#ifndef COSTFUNCTIONTHERMORTZ_H
#define COSTFUNCTIONTHERMORTZ_H

#include "CostFunctionThermo.h"
#include <QString>

class CostFunctionThermoRTZ: public CostFunctionThermo
{
public:
    CostFunctionThermoRTZ(const int& numObs = 0, const int& stateSize = 0);
    ~CostFunctionThermoRTZ();

private:
    bool outputAnalysis(const QString& suffix, real* Astate);
    bool writeAsi(const QString& asiFileName);
    bool writeNetCDF(const QString& netcdfFileName);
};

#endif // COSTFUNCTIONTHERMORTZ_H
