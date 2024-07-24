#include "CostFunctionThermoXYZ.h"
#include <cmath>
#include <QTextStream>
#include <QDir>
#include <QDateTime>
#include <netcdfcpp.h>
#include <GeographicLib/TransverseMercatorExact.hpp>

CostFunctionThermoXYZ::CostFunctionThermoXYZ(const int& numObs, const int& stateSize)
    : CostFunctionThermo(numObs, stateSize)
{
}

CostFunctionThermoXYZ::~CostFunctionThermoXYZ()
{
}

bool CostFunctionThermoXYZ::outputAnalysis(const QString& suffix, real* Astate)
{

    cout << "Outputting " << suffix.toStdString() << "...\n";
    // H --> to Mish for output
    QString samuraiout = "samurai_XYZ_" + suffix + ".out";
    ofstream samuraistream;
    if (configHash->value("output_txt") == "true") {
        samuraistream.open(outputPath.absoluteFilePath(samuraiout).toAscii().data());
        samuraistream << "X\tY\tZ\tpip\tthetarhop\tftheta\n";
        samuraistream.precision(10);
    }

    int analysisDim = 12;
    int analysisSize = (iDim-2)*(jDim-2)*(kDim-2);
    finalAnalysis = new real[analysisSize*analysisDim];
    real gausspoint = 0.5*sqrt(1./3.);

    for (int iIndex = 1; iIndex < iDim-1; iIndex++) {
        for (int ihalf = 0; ihalf <= mishFlag; ihalf++) {
            for (int imu = -ihalf; imu <= ihalf; imu++) {
                real i = iMin + DI * (iIndex + (gausspoint * imu + 0.5*ihalf));
                if (i > ((iDim-1)*DI + iMin)) continue;

                for (int jIndex = 1; jIndex < jDim-1; jIndex++) {
                    for (int jhalf =0; jhalf <=mishFlag; jhalf++) {
                        for (int jmu = -jhalf; jmu <= jhalf; jmu++) {
                            real j = jMin + DJ * (jIndex + (gausspoint * jmu + 0.5*jhalf));
                            if (j > ((jDim-1)*DJ + jMin)) continue;

                            for (int kIndex = 1; kIndex < kDim-1; kIndex++) {
                                for (int khalf =0; khalf <=mishFlag; khalf++) {
                                    for (int kmu = -khalf; kmu <= khalf; kmu++) {
                                        real k = kMin + DK * (kIndex + (gausspoint * kmu + 0.5*khalf));
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
                                        real pipdx = 0.; real thetarhopdx = 0.; real fthetadx = 0.;
                                        real pipdy = 0.; real thetarhopdy = 0.; real fthetady = 0.;
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
                                                                pipdx += Astate[aIndex] * idbasis * jbasis * kbasis;
                                                                pipdy += Astate[aIndex] * ibasis * jdbasis * kbasis;
                                                                pipdz += Astate[aIndex] * ibasis * jbasis * kdbasis;
                                                                break;
                                                            case 1:
                                                                thetarhop +=  Astate[aIndex + 1] * basis3x;
                                                                thetarhopdx += Astate[aIndex + 1] * idbasis * jbasis * kbasis;
                                                                thetarhopdy += Astate[aIndex + 1] * ibasis * jdbasis * kbasis;
                                                                thetarhopdz += Astate[aIndex + 1] * ibasis * jbasis * kdbasis;
                                                                break;
                                                            case 2:
                                                                ftheta += Astate[aIndex + 2] * basis3x;
                                                                fthetadx += Astate[aIndex + 2] * idbasis * jbasis * kbasis;
                                                                fthetady += Astate[aIndex + 2] * ibasis * jdbasis * kbasis;
                                                                fthetadz += Astate[aIndex + 2] * ibasis * jbasis * kdbasis;
                                                                break;
                                                        }
                                                    }
                                                }
                                            }
                                        }


                                        if (configHash->value("output_txt") == "true") {
                                            samuraistream << scientific << i << "\t" << j << "\t"  << k
                                            << "\t" << pip << "\t" << thetarhop << "\t" << ftheta << "\t" << "\n";
                                        }

                                        // On the nodes
                                        if (!ihalf and !jhalf and !khalf){
                                            int fIndex = (iDim-2)*(jDim-2)*(kDim-2);
                                            int posIndex = (iDim-2)*(jDim-2)*(kIndex-1) + (iDim-2)*(jIndex-1) + (iIndex-1);
                                            finalAnalysis[fIndex * 0 + posIndex] = pip;
                                            finalAnalysis[fIndex * 1 + posIndex] = thetarhop;
                                            finalAnalysis[fIndex * 2 + posIndex] = ftheta;
                                            finalAnalysis[fIndex * 3 + posIndex] = pipdx;
                                            finalAnalysis[fIndex * 4 + posIndex] = thetarhopdx;
                                            finalAnalysis[fIndex * 5 + posIndex] = fthetadx;
                                            finalAnalysis[fIndex * 6 + posIndex] = pipdy;
                                            finalAnalysis[fIndex * 7 + posIndex] = thetarhopdy;
                                            finalAnalysis[fIndex * 8 + posIndex] = fthetady;
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

    QString fileName = "samurai_XYZ_" + suffix;
    QString outFileName = outputPath.absoluteFilePath(fileName);

    // Write the Obs to a summary text file
    if (configHash->value("output_qc") == "true") {
        QString qcout = "samurai_QC_" + suffix + ".out";
        QString qcFileName = outputPath.absoluteFilePath(qcout);
        ofstream qcstream(qcFileName.toAscii().data());
        ostream_iterator<string> os(qcstream, "\t ");
        *os++ = "Observation";
        *os++ = "Inverse Error";
        *os++ = "X";
        *os++ = "Y";
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

            // Multiply the weight by the ob -- Observations.in has individual weights already
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

bool CostFunctionThermoXYZ::writeNetCDF(const QString& netcdfFileName)
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
    NcDim *lvlDim, *latDim, *lonDim, *timeDim;
    if (!(lonDim = dataFile.add_dim("longitude", iDim)))
        return NC_ERR;
    if (!(latDim = dataFile.add_dim("latitude", jDim)))
        return NC_ERR;
    if (!(lvlDim = dataFile.add_dim("altitude", kDim)))
        return NC_ERR;
    // Add an unlimited dimension...
    if (!(timeDim = dataFile.add_dim("time")))
        return NC_ERR;

    // Define the coordinate variables.
    NcVar *latVar, *lonVar, *lvlVar, *timeVar, *xVar, *yVar;
    if (!(lonVar = dataFile.add_var("longitude", ncFloat, lonDim)))
        return NC_ERR;
    if (!(latVar = dataFile.add_var("latitude", ncFloat, latDim)))
        return NC_ERR;
    if (!(xVar = dataFile.add_var("x", ncFloat, lonDim)))
        return NC_ERR;
    if (!(yVar = dataFile.add_var("y", ncFloat, latDim)))
        return NC_ERR;
    if (!(lvlVar = dataFile.add_var("altitude", ncFloat, lvlDim)))
        return NC_ERR;
    if (!(timeVar = dataFile.add_var("time", ncInt, timeDim)))
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
    NcVar *pip, *thetarhop, *ftheta;
    NcVar *dpipdx, *dthetarhopdx, *dfthetadx, *dpipdy, *dthetarhopdy, *dfthetady, *dpipdz, *dthetarhopdz, *dfthetadz;


    if (!(pip = dataFile.add_var("PIP", ncFloat, timeDim,
                               lvlDim, latDim, lonDim)))
        return NC_ERR;
    if (!(thetarhop = dataFile.add_var("THETARHOP", ncFloat, timeDim,
                               lvlDim, latDim, lonDim)))
        return NC_ERR;
    if (!(ftheta = dataFile.add_var("FTHETA", ncFloat, timeDim,
                               lvlDim, latDim, lonDim)))
        return NC_ERR;
    if (!(dpipdx = dataFile.add_var("DPIPDX", ncFloat, timeDim,
                                  lvlDim, latDim, lonDim)))
        return NC_ERR;
    if (!(dthetarhopdx = dataFile.add_var("DTHETARHOPDX", ncFloat, timeDim,
                                  lvlDim, latDim, lonDim)))
        return NC_ERR;
    if (!(dfthetadx = dataFile.add_var("DFTHETADX", ncFloat, timeDim,
                                  lvlDim, latDim, lonDim)))
        return NC_ERR;
    if (!(dpipdy = dataFile.add_var("DPIPDY", ncFloat, timeDim,
                                  lvlDim, latDim, lonDim)))
        return NC_ERR;
    if (!(dthetarhopdy = dataFile.add_var("DTHETARHOPDY", ncFloat, timeDim,
                                  lvlDim, latDim, lonDim)))
        return NC_ERR;
    if (!(dfthetady = dataFile.add_var("DFTHETADY", ncFloat, timeDim,
                                  lvlDim, latDim, lonDim)))
        return NC_ERR;
    if (!(dpipdz = dataFile.add_var("DPIPDZ", ncFloat, timeDim,
                                  lvlDim, latDim, lonDim)))
        return NC_ERR;
    if (!(dthetarhopdz = dataFile.add_var("DTHETARHOPDZ", ncFloat, timeDim,
                                  lvlDim, latDim, lonDim)))
        return NC_ERR;
    if (!(dfthetadz = dataFile.add_var("DFTHETADZ", ncFloat, timeDim,
                                  lvlDim, latDim, lonDim)))
        return NC_ERR;

    // Define units attributes for data variables.
    if (!pip->add_att("units", "???"))
        return NC_ERR;
    if (!thetarhop->add_att("units", "???"))
        return NC_ERR;
    if (!ftheta->add_att("units", "???"))
        return NC_ERR;
    if (!dpipdx->add_att("units", "10-5s-1"))
        return NC_ERR;
    if (!dthetarhopdx->add_att("units", "10-5s-1"))
        return NC_ERR;
    if (!dfthetadx->add_att("units", "10-5s-1"))
        return NC_ERR;
    if (!dpipdy->add_att("units", "10-5s-1"))
        return NC_ERR;
    if (!dthetarhopdy->add_att("units", "10-5s-1"))
        return NC_ERR;
    if (!dfthetady->add_att("units", "10-5s-1"))
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
    if (!dpipdx->add_att("long_name", "pi prime gradient"))
        return NC_ERR;
    if (!dthetarhopdx->add_att("long_name", "theta rho prime gradient"))
        return NC_ERR;
    if (!dfthetadx->add_att("long_name", "f theta gradient"))
        return NC_ERR;
    if (!dpipdy->add_att("long_name", "pi prime gradient"))
        return NC_ERR;
    if (!dthetarhopdy->add_att("long_name", "theta rho prime gradient"))
        return NC_ERR;
    if (!dfthetady->add_att("long_name", "f theta gradient"))
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
    if (!dpipdx->add_att("missing_value", -999.f))
        return NC_ERR;
    if (!dthetarhopdx->add_att("missing_value", -999.f))
        return NC_ERR;
    if (!dfthetadx->add_att("missing_value", -999.f))
        return NC_ERR;
    if (!dpipdy->add_att("missing_value", -999.f))
        return NC_ERR;
    if (!dthetarhopdy->add_att("missing_value", -999.f))
        return NC_ERR;
    if (!dfthetady->add_att("missing_value", -999.f))
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
    if (!dpipdx->add_att("_FillValue", -999.f))
        return NC_ERR;
    if (!dthetarhopdx->add_att("_FillValue", -999.f))
        return NC_ERR;
    if (!dfthetadx->add_att("_FillValue", -999.f))
        return NC_ERR;
    if (!dpipdy->add_att("_FillValue", -999.f))
        return NC_ERR;
    if (!dthetarhopdy->add_att("_FillValue", -999.f))
        return NC_ERR;
    if (!dfthetady->add_att("_FillValue", -999.f))
        return NC_ERR;
    if (!dpipdz->add_att("_FillValue", -999.f))
        return NC_ERR;
    if (!dthetarhopdz->add_att("_FillValue", -999.f))
        return NC_ERR;
    if (!dfthetadz->add_att("_FillValue", -999.f))
        return NC_ERR;

    // Write the coordinate variable data to the file.
    real *lons = new real[iDim];
    real *lats = new real[jDim];
    real *levs = new real[kDim];
    real *x = new real[iDim];
    real *y = new real[jDim];
    int time[2];

    // Reference time and position from center file
    time[0] = configHash->value("ref_time").toInt();
    real latReference = configHash->value("ref_lat").toFloat();
    real lonReference = configHash->value("ref_lon").toFloat();
    real refX, refY;

    GeographicLib::TransverseMercatorExact tm = GeographicLib::TransverseMercatorExact::UTM();
    tm.Forward(lonReference, latReference, lonReference, refX, refY);
    for (int iIndex = 0; iIndex < iDim; iIndex++) {
        real i = (iMin + DI * iIndex)*1000;
        real j = (jMin + DJ * (jDim/2))*1000;
        real latnull = 0;
        tm.Reverse(lonReference,refX + i, refY + j, latnull, lons[iIndex]);
        x[iIndex] = i/1000;
    }

    for (int jIndex = 0; jIndex < jDim; jIndex++) {
        real i = (iMin + DI * (iDim/2))*1000;
        real j = (jMin + DJ * jIndex)*1000;
        real lonnull = 0;
        tm.Reverse(lonReference,refX + i, refY + j, lats[jIndex], lonnull);
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
        if (!dpipdx->put_rec(&finalAnalysis[iDim*jDim*kDim*3], rec))
            return NC_ERR;
        if (!dthetarhopdx->put_rec(&finalAnalysis[iDim*jDim*kDim*4], rec))
            return NC_ERR;
        if (!dfthetadx->put_rec(&finalAnalysis[iDim*jDim*kDim*5], rec))
            return NC_ERR;
        if (!dpipdy->put_rec(&finalAnalysis[iDim*jDim*kDim*6], rec))
            return NC_ERR;
        if (!dthetarhopdy->put_rec(&finalAnalysis[iDim*jDim*kDim*7], rec))
            return NC_ERR;
        if (!dfthetady->put_rec(&finalAnalysis[iDim*jDim*kDim*8], rec))
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
    delete[] lats;
    delete[] lons;
    delete[] levs;
    delete[] x;
    delete[] y;
    return true;

}

bool CostFunctionThermoXYZ::writeAsi(const QString& asiFileName)
{
    cout << "Writing Asi not implemented yet" << endl;
    return true;
}
