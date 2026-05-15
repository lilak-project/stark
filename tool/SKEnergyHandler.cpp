#include "SKEnergyHandler.h"
#include "LKLogger.h"
#include <string>
#include <fstream>

using namespace std;

ClassImp(SKEnergyHandler)

SKEnergyHandler::SKEnergyHandler(TString energyCalibrationName, TString positionCalibrationName)
{
    int det, side, strip;
    string dummy;

    ifstream energyCalibrationFile(energyCalibrationName);
    if (energyCalibrationFile.is_open()==false) {
        e_error << "[SKEnergyHandler] Cannot open " << energyCalibrationName << endl;
        return;
    }
    else {
        e_info << "[SKEnergyHandler] Reading energy calibration file " << energyCalibrationName << endl;
    }

    double c0_itcpt, c0_slope, c1_itcptL, c1_slopeL, c1_itcptR, c1_slopeR, c2_par0, c2_par1, c2_par2, c3_itcpt, c3_slope;
    getline(energyCalibrationFile,dummy);
    while (energyCalibrationFile>>det>>side>>strip>>c0_itcpt>>c0_slope>>c1_itcptL>>c1_slopeL>>c1_itcptR>>c1_slopeR>>c2_par0>>c2_par1>>c2_par2>>c3_itcpt>>c3_slope)
    {
        fC0Parameters[det][side][strip][0] = c0_itcpt;
        fC0Parameters[det][side][strip][1] = c0_slope;
        fC1Parameters[det][side][strip][0][0] = c1_itcptL;
        fC1Parameters[det][side][strip][0][1] = c1_slopeL;
        fC1Parameters[det][side][strip][1][0] = c1_itcptR;
        fC1Parameters[det][side][strip][1][1] = c1_slopeR;
        fC2Parameters[det][side][strip][0] = c2_par0;
        fC2Parameters[det][side][strip][1] = c2_par1;
        fC2Parameters[det][side][strip][2] = c2_par2;
        fC3Parameters[det][side][strip][0] = c3_itcpt;
        fC3Parameters[det][side][strip][1] = c3_slope;
    }

    if (positionCalibrationName.IsNull())
        return;

    /////////////////////////////////////////////////////////////////////

    ifstream positionCalibrationFile(positionCalibrationName);
    if (positionCalibrationFile.is_open()==false) {
        e_error << "[SKEnergyHandler] Cannot open " << positionCalibrationName << endl;
        return;
    }
    else {
        e_info << "[SKEnergyHandler] Reading position calibration file " << positionCalibrationName << endl;
    }

    double x00 ,x01 ,x02 ,x03 ,x04 ,x10 ,x11 ,x12 ,x13 ,x14;
    getline(positionCalibrationFile,dummy);
    while (positionCalibrationFile>>det>>side>>strip >> x00 >> x01 >> x02 >> x03 >> x04 >> x10 >> x11 >> x12 >> x13 >> x14)
    {
        fCPParameters[det][side][strip][0][0] = x00;
        fCPParameters[det][side][strip][0][1] = x01;
        fCPParameters[det][side][strip][0][2] = x02;
        fCPParameters[det][side][strip][0][3] = x03;
        fCPParameters[det][side][strip][0][4] = x04;
        fCPParameters[det][side][strip][1][0] = x10;
        fCPParameters[det][side][strip][1][1] = x11;
        fCPParameters[det][side][strip][1][2] = x12;
        fCPParameters[det][side][strip][1][3] = x13;
        fCPParameters[det][side][strip][1][4] = x14;
    }
}

double SKEnergyHandler::RestoreEnergy(int det, int side, int strip, double &energy)
{
    energy = fC0Parameters[det][side][strip][1] * energy + fC0Parameters[det][side][strip][0];
    return energy;
}

void SKEnergyHandler::RestoreEnergyPosition(int det, int side, int strip, double energyL, double energyR, double &position, double &energySum, bool skipPositionCalibration)
{
    // slope correction
    energyL = fC1Parameters[det][side][strip][0][1] * energyL + fC1Parameters[det][side][strip][0][0];
    energyR = fC1Parameters[det][side][strip][1][1] * energyR + fC1Parameters[det][side][strip][1][0];
    // ballistic correction
    energySum = energyL + energyR;
    position = (energyR - energyL) / energySum;
    energySum = energySum / (fC2Parameters[det][side][strip][0] + fC2Parameters[det][side][strip][1]*position + fC2Parameters[det][side][strip][2]*position*position) * fGate1Energy;
    // energy correction
    energySum = fC3Parameters[det][side][strip][1] * energySum + fC3Parameters[det][side][strip][0];
    // position correction
    if (!skipPositionCalibration) RestorePosition(det, side, strip, position, energySum);
}

double SKEnergyHandler::RestorePosition(int det, int side, int strip, double &position, double energySum)
{
    double x01 = fCPParameters[det][side][strip][0][1];
    double x02 = fCPParameters[det][side][strip][0][2];
    double x03 = fCPParameters[det][side][strip][0][3];
    double x11 = fCPParameters[det][side][strip][1][1];
    double x12 = fCPParameters[det][side][strip][1][2];
    double x13 = fCPParameters[det][side][strip][1][3];
    double x1 = (x11-x01)/(fGate1Energy-fGate0Energy) * (energySum-fGate0Energy) + x01;
    double x2 = (x12-x02)/(fGate1Energy-fGate0Energy) * (energySum-fGate0Energy) + x02;
    double x3 = (x13-x03)/(fGate1Energy-fGate0Energy) * (energySum-fGate0Energy) + x03;
    if (position>x2) position =  0.5/(x3-x2)*(position-x2);
    else             position = -0.5/(x1-x2)*(position-x2);
    return position;
}
