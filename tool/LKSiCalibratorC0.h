#ifndef LKSICALIBRATORC0_HH
#define LKSICALIBRATORC0_HH

#include "TNamed.h"

#include <vector>

class TH1;

class LKSiCalibratorC0 : public TNamed
{
    public:
        struct Result
        {
            bool success = false;
            double entries = 0;
            double intercept = 0;
            double slope = 0;
            std::vector<double> peakMeans;
            std::vector<double> peakSigmas;
        };

    public:
        LKSiCalibratorC0();
        virtual ~LKSiCalibratorC0() {}

        Result Fit(TH1 *hist, const std::vector<double> &gateEnergies, double expectedResolution=0.01, double entriesCut=100.0);

    ClassDef(LKSiCalibratorC0, 1);
};

#endif
