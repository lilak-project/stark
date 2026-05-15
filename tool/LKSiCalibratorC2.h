#ifndef LKSICALIBRATORC2_HH
#define LKSICALIBRATORC2_HH

#include "TNamed.h"

#include <vector>

class TH2;
class TF1;

class LKSiCalibratorC2 : public TNamed
{
    public:
        struct GateResult
        {
            bool success = false;
            double entries = 0;
            double b0 = 0;
            double b1 = 0;
            double b2 = 0;
        };

        struct Result
        {
            bool success = false;
            int selectedGate = -1;
            double entries = 0;
            double b0 = 0;
            double b1 = 0;
            double b2 = 0;
            std::vector<GateResult> gates;
        };

    public:
        LKSiCalibratorC2();
        virtual ~LKSiCalibratorC2() {}

        void Clear(Option_t *option="");
        GateResult FitSingle(TH2 *hist);
        Result Fit(const std::vector<TH2 *> &positionEnergyHists, int chooseGate=1, double entriesCut=100.0);
        TF1 *GetLastFit() const { return fLastFit; }

    private:
        TF1 *fLastFit = nullptr;

    ClassDef(LKSiCalibratorC2, 1);
};

#endif
