#include "SKEnergyHandler.h"
#include "LKLogger.h"
#include <string>
#include <fstream>

using namespace std;

ClassImp(SKEnergyHandler)

void SKEnergyHandler::InitializeParameters()
{
    for (auto det = 0; det < kMaxDetectors; ++det) {
        for (auto side = 0; side < kMaxSides; ++side) {
            for (auto strip = 0; strip < kMaxStrips; ++strip) {
                fC0Parameters[det][side][strip][0] = 0;
                fC0Parameters[det][side][strip][1] = 1;
                for (auto lr = 0; lr < 2; ++lr) {
                    fC1Parameters[det][side][strip][lr][0] = 0;
                    fC1Parameters[det][side][strip][lr][1] = 1;
                }
                for (auto cpSide = 0; cpSide < 2; ++cpSide) {
                    for (auto i = 0; i < 5; ++i)
                        fCPParameters[det][side][strip][cpSide][i] = 0;
                }
                fC2Parameters[det][side][strip][0] = fGate1Energy;
                fC2Parameters[det][side][strip][1] = 0;
                fC2Parameters[det][side][strip][2] = 0;
                fC3Parameters[det][side][strip][0] = 0;
                fC3Parameters[det][side][strip][1] = 1;
            }
        }
    }
}

bool SKEnergyHandler::IsValidChannel(int det, int side, int strip) const
{
    return det >= 0 && det < kMaxDetectors
        && side >= 0 && side < kMaxSides
        && strip >= 0 && strip < kMaxStrips;
}

SKEnergyHandler::SKEnergyHandler(TString energyCalibrationName, TString positionCalibrationName)
{
    InitializeParameters();

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
        if (!IsValidChannel(det, side, strip)) {
            e_warning << "[SKEnergyHandler] Skip invalid energy calibration channel "
                      << det << " " << side << " " << strip << endl;
            continue;
        }
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
        if (!IsValidChannel(det, side, strip)) {
            e_warning << "[SKEnergyHandler] Skip invalid position calibration channel "
                      << det << " " << side << " " << strip << endl;
            continue;
        }
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
        fHasPositionCalibration = true;
    }
}

double SKEnergyHandler::RestoreEnergy(int det, int side, int strip, double &energy)
{
    if (!IsValidChannel(det, side, strip))
        return energy;
    energy = fC0Parameters[det][side][strip][1] * energy + fC0Parameters[det][side][strip][0];
    return energy;
}

void SKEnergyHandler::RestoreEnergyPosition(int det, int side, int strip, double energyL, double energyR, double &position, double &energySum, bool skipPositionCalibration)
{
    position = 0;
    energySum = 0;
    if (!IsValidChannel(det, side, strip))
        return;

    // slope correction
    energyL = fC1Parameters[det][side][strip][0][1] * energyL + fC1Parameters[det][side][strip][0][0];
    energyR = fC1Parameters[det][side][strip][1][1] * energyR + fC1Parameters[det][side][strip][1][0];
    // ballistic correction
    energySum = energyL + energyR;
    if (energySum <= 0)
        return;
    position = (energyR - energyL) / energySum;
    auto scale = fC2Parameters[det][side][strip][0] + fC2Parameters[det][side][strip][1]*position + fC2Parameters[det][side][strip][2]*position*position;
    if (scale == 0) {
        energySum = 0;
        return;
    }
    energySum = energySum / scale * fGate1Energy;
    // energy correction
    energySum = fC3Parameters[det][side][strip][1] * energySum + fC3Parameters[det][side][strip][0];
    // position correction
    if (!skipPositionCalibration && fHasPositionCalibration)
        RestorePosition(det, side, strip, position, energySum);
}

double SKEnergyHandler::RestorePosition(int det, int side, int strip, double &position, double energySum)
{
    if (!IsValidChannel(det, side, strip))
        return position;

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
