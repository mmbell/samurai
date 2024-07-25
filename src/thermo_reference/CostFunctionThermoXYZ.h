#ifndef COSTFUNCTIONTHERMOXYZ_H
#define COSTFUNCTIONTHERMOXYZ_H

#include "CostFunctionThermo.h"
#include <QString>

class CostFunctionThermoXYZ: public CostFunctionThermo
{
public:
    CostFunctionThermoXYZ(const int& numObs = 0, const int& stateSize = 0);
    ~CostFunctionThermoXYZ();

private:
    bool outputAnalysis(const QString& suffix, real* Astate);
    bool writeAsi(const QString& asiFileName);
    bool writeNetCDF(const QString& netcdfFileName);
};

#endif // COSTFUNCTIONTHERMOXYZ_H





