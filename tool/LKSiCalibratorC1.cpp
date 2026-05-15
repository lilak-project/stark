#include "LKSiCalibratorC1.h"

#include "TF1.h"
#include "TGraphErrors.h"
#include "TH1.h"
#include "TH2.h"

ClassImp(LKSiCalibratorC1)

LKSiCalibratorC1::LKSiCalibratorC1()
    : TNamed("LKSiCalibratorC1", "")
{
}

void LKSiCalibratorC1::Clear(Option_t *)
{
    delete fLeftGraph;
    delete fRightGraph;
    delete fLeftFit;
    delete fRightFit;
    fLeftGraph = nullptr;
    fRightGraph = nullptr;
    fLeftFit = nullptr;
    fRightFit = nullptr;
}

LKSiCalibratorC1::GateResult LKSiCalibratorC1::FitProjectionPair(TH1 *leftHist, TH1 *rightHist) const
{
    GateResult result;
    if (leftHist == nullptr || rightHist == nullptr)
        return result;

    result.entries = leftHist->GetEntries() + rightHist->GetEntries();
    TF1 fitLeft("fitLeftPeak", "gaus(0)", leftHist->GetXaxis()->GetXmin(), leftHist->GetXaxis()->GetXmax());
    TF1 fitRight("fitRightPeak", "gaus(0)", rightHist->GetXaxis()->GetXmin(), rightHist->GetXaxis()->GetXmax());

    auto maxBinLeft = leftHist->GetMaximumBin();
    auto maxBinRight = rightHist->GetMaximumBin();
    auto meanLeft = leftHist->GetBinCenter(maxBinLeft);
    auto meanRight = rightHist->GetBinCenter(maxBinRight);
    auto sigmaLeft = 0.05 * (leftHist->GetXaxis()->GetXmax() - leftHist->GetXaxis()->GetXmin());
    auto sigmaRight = 0.05 * (rightHist->GetXaxis()->GetXmax() - rightHist->GetXaxis()->GetXmin());

    fitLeft.SetRange(meanLeft - sigmaLeft, meanLeft + sigmaLeft);
    fitLeft.SetParameters(leftHist->GetBinContent(maxBinLeft), meanLeft, sigmaLeft / 3.);
    leftHist->Fit(&fitLeft, "Q0NR");

    fitRight.SetRange(meanRight - sigmaRight, meanRight + sigmaRight);
    fitRight.SetParameters(rightHist->GetBinContent(maxBinRight), meanRight, sigmaRight / 3.);
    rightHist->Fit(&fitRight, "Q0NR");

    result.meanLeft = fitLeft.GetParameter(1);
    result.sigmaLeft = fitLeft.GetParameter(2);
    result.meanRight = fitRight.GetParameter(1);
    result.sigmaRight = fitRight.GetParameter(2);
    result.success = true;
    return result;
}

LKSiCalibratorC1::Result LKSiCalibratorC1::Fit(const std::vector<TH2 *> &leftRightHists, const std::vector<double> &gateEnergies, double entriesCut)
{
    std::vector<TH1 *> leftHists;
    std::vector<TH1 *> rightHists;
    leftHists.reserve(leftRightHists.size());
    rightHists.reserve(leftRightHists.size());
    for (size_t iGate = 0; iGate < leftRightHists.size(); ++iGate) {
        auto hist = leftRightHists[iGate];
        if (hist == nullptr) {
            leftHists.push_back(nullptr);
            rightHists.push_back(nullptr);
            continue;
        }
        leftHists.push_back(hist->ProjectionX(Form("%s_px_%zu", hist->GetName(), iGate)));
        rightHists.push_back(hist->ProjectionY(Form("%s_py_%zu", hist->GetName(), iGate)));
    }

    auto result = Fit(leftHists, rightHists, gateEnergies, entriesCut);
    for (auto hist : leftHists)
        delete hist;
    for (auto hist : rightHists)
        delete hist;
    return result;
}

LKSiCalibratorC1::Result LKSiCalibratorC1::Fit(const std::vector<TH1 *> &leftHists, const std::vector<TH1 *> &rightHists, const std::vector<double> &gateEnergies, double entriesCut)
{
    Clear();
    Result result;
    fLeftGraph = new TGraphErrors();
    fRightGraph = new TGraphErrors();

    auto numGates = std::min(gateEnergies.size(), std::min(leftHists.size(), rightHists.size()));
    for (size_t iGate = 0; iGate < numGates; ++iGate) {
        auto gateResult = FitProjectionPair(leftHists[iGate], rightHists[iGate]);
        if (gateResult.entries < entriesCut || !gateResult.success) {
            result.gates.push_back(gateResult);
            continue;
        }
        auto idxL = fLeftGraph->GetN();
        fLeftGraph->SetPoint(idxL, gateResult.meanLeft, gateEnergies[iGate]);
        fLeftGraph->SetPointError(idxL, gateResult.sigmaLeft, gateEnergies[iGate] * 0.01);
        auto idxR = fRightGraph->GetN();
        fRightGraph->SetPoint(idxR, gateResult.meanRight, gateEnergies[iGate]);
        fRightGraph->SetPointError(idxR, gateResult.sigmaRight, gateEnergies[iGate] * 0.01);
        result.gates.push_back(gateResult);
    }

    if (fLeftGraph->GetN() < 2 || fRightGraph->GetN() < 2)
        return result;

    fLeftFit = new TF1("fitC1Left", "pol1", 0, 100000);
    fRightFit = new TF1("fitC1Right", "pol1", 0, 100000);
    fLeftGraph->Fit(fLeftFit, "Q0R");
    fRightGraph->Fit(fRightFit, "Q0R");

    result.interceptLeft = fLeftFit->GetParameter(0);
    result.slopeLeft = fLeftFit->GetParameter(1);
    result.interceptRight = fRightFit->GetParameter(0);
    result.slopeRight = fRightFit->GetParameter(1);
    result.nPoints = std::min(fLeftGraph->GetN(), fRightGraph->GetN());
    result.success = true;
    return result;
}
