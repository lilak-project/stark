#ifndef LKSICALIBRATORC1_HH
#define LKSICALIBRATORC1_HH

#include "TNamed.h"

#include <vector>

class TH1;
class TH2;
class TF1;
class TGraphErrors;

class LKSiCalibratorC1 : public TNamed
{
    public:
        struct GateResult
        {
            bool success = false;
            double entries = 0;
            double meanLeft = 0;
            double sigmaLeft = 0;
            double meanRight = 0;
            double sigmaRight = 0;
        };

        struct Result
        {
            bool success = false;
            int nPoints = 0;
            double interceptLeft = 0;
            double slopeLeft = 0;
            double interceptRight = 0;
            double slopeRight = 0;
            std::vector<GateResult> gates;
        };

    public:
        LKSiCalibratorC1();
        virtual ~LKSiCalibratorC1() {}

        void Clear(Option_t *option="");
        Result Fit(const std::vector<TH2 *> &leftRightHists, const std::vector<double> &gateEnergies, double entriesCut=100.0);
        Result Fit(const std::vector<TH1 *> &leftHists, const std::vector<TH1 *> &rightHists, const std::vector<double> &gateEnergies, double entriesCut=100.0);

        TGraphErrors *GetLeftGraph() const { return fLeftGraph; }
        TGraphErrors *GetRightGraph() const { return fRightGraph; }
        TF1 *GetLeftFit() const { return fLeftFit; }
        TF1 *GetRightFit() const { return fRightFit; }

    private:
        GateResult FitProjectionPair(TH1 *leftHist, TH1 *rightHist) const;

    private:
        TGraphErrors *fLeftGraph = nullptr;
        TGraphErrors *fRightGraph = nullptr;
        TF1 *fLeftFit = nullptr;
        TF1 *fRightFit = nullptr;

    ClassDef(LKSiCalibratorC1, 1);
};

#endif
