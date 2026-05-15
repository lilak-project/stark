#ifndef SKENERGYHANDLER_HH
#define SKENERGYHANDLER_HH

#include "TObject.h"
#include "TString.h"

class SKEnergyHandler : public TObject
{
    public:
        SKEnergyHandler(TString energyCalibrationName, TString positionCalibrationName="");
        virtual ~SKEnergyHandler() {}

        double RestoreEnergy(int det, int side, int strip, double &energy);
        void RestoreEnergyPosition(int det, int side, int strip, double energyL, double energyR, double &position, double &energySum, bool skipPositionCalibration=false); ///< for Resistive Strip
        double RestorePosition(int det, int side, int strip, double &position, double energySum);

    private:
        double fC0Parameters[40][2][8][2];
        double fC1Parameters[40][2][8][2][2];
        double fCPParameters[40][2][8][2][5];
        double fC2Parameters[40][2][8][3];
        double fC3Parameters[40][2][8][2];
        const double fGate1Energy = 5.486;
        const double fGate0Energy = 3.1822;

    ClassDef(SKEnergyHandler,1);
};

#endif
