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
        enum { kMaxDetectors = 128, kMaxSides = 2, kMaxStrips = 8 };
        bool IsValidChannel(int det, int side, int strip) const;
        void InitializeParameters();

        double fC0Parameters[kMaxDetectors][kMaxSides][kMaxStrips][2];
        double fC1Parameters[kMaxDetectors][kMaxSides][kMaxStrips][2][2];
        double fCPParameters[kMaxDetectors][kMaxSides][kMaxStrips][2][5];
        double fC2Parameters[kMaxDetectors][kMaxSides][kMaxStrips][3];
        double fC3Parameters[kMaxDetectors][kMaxSides][kMaxStrips][2];
        bool fHasPositionCalibration = false;
        const double fGate1Energy = 5.486;
        const double fGate0Energy = 3.1822;

    ClassDef(SKEnergyHandler,1);
};

#endif
