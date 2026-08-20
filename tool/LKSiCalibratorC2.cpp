#include "LKSiCalibratorC2.h"

#include "TF1.h"
#include "TGraphErrors.h"
#include "TH1.h"
#include "TH2.h"

#include <algorithm>
#include <cmath>

ClassImp(LKSiCalibratorC2)

LKSiCalibratorC2::LKSiCalibratorC2()
    : TNamed("LKSiCalibratorC2", "")
{
}

void LKSiCalibratorC2::Clear(Option_t *)
{
    delete fLastFit;
    fLastFit = nullptr;
}

LKSiCalibratorC2::GateResult LKSiCalibratorC2::FitSingle(TH2 *hist, double sourceEnergy, double expectedResolution, double rangeScale)
{
    GateResult result;
    if (hist == nullptr)
        return result;

    result.entries = hist->GetEntries();
    if (result.entries <= 0)
        return result;

    Clear();
    auto xAxis = hist -> GetXaxis();
    auto yAxis = hist -> GetYaxis();
    auto yMin = yAxis -> GetXmin();
    auto yMax = yAxis -> GetXmax();
    if (sourceEnergy > 0) {
        auto halfRange = std::max(10.0 * expectedResolution * sourceEnergy,
                                  rangeScale * 5.0 * expectedResolution * sourceEnergy);
        yMin = std::max(yMin, sourceEnergy - halfRange);
        yMax = std::min(yMax, sourceEnergy + halfRange);
    }

    TGraphErrors graph;
    for (auto ix = 1; ix <= hist->GetNbinsX(); ++ix) {
        auto projection = hist -> ProjectionY(Form("%s_py_%d", hist->GetName(), ix), ix, ix);
        projection -> SetDirectory(nullptr);

        auto entries = projection -> Integral(projection->FindBin(yMin), projection->FindBin(yMax));
        if (entries < 10) {
            delete projection;
            continue;
        }

        auto yBin1 = projection -> FindBin(yMin);
        auto yBin2 = projection -> FindBin(yMax);
        auto maxBin = yBin1;
        auto maxContent = projection -> GetBinContent(maxBin);
        for (auto iy = yBin1 + 1; iy <= yBin2; ++iy) {
            auto content = projection -> GetBinContent(iy);
            if (content > maxContent) {
                maxContent = content;
                maxBin = iy;
            }
        }

        auto peak = projection -> GetBinCenter(maxBin);
        auto binWidth = projection -> GetXaxis() -> GetBinWidth(1);
        auto sigma0 = std::max(peak * expectedResolution, binWidth);
        TF1 fitPeak(Form("%s_fit_peak_%d", hist->GetName(), ix), "gaus(0)",
                    std::max(yMin, peak - 5.0 * sigma0),
                    std::min(yMax, peak + 5.0 * sigma0));
        fitPeak.SetParameters(maxContent, peak, sigma0);
        projection -> Fit(&fitPeak, "Q0NR");

        auto mean = fitPeak.GetParameter(1);
        auto sigma = std::abs(fitPeak.GetParameter(2));
        if (sigma <= 0 || mean < yMin || mean > yMax) {
            mean = peak;
            sigma = binWidth;
        }

        auto iPoint = graph.GetN();
        graph.SetPoint(iPoint, xAxis->GetBinCenter(ix), mean);
        graph.SetPointError(iPoint, 0.5 * xAxis->GetBinWidth(ix), sigma / std::sqrt(std::max(1.0, entries)));
        delete projection;
    }

    result.points = graph.GetN();
    if (result.points < 3)
        return result;

    fLastFit = new TF1(Form("fit_%s", hist->GetName()), "pol2", xAxis->GetXmin(), xAxis->GetXmax());
    fLastFit->SetParameters(sourceEnergy > 0 ? sourceEnergy : graph.GetY()[0], 0.0, 0.2);
    graph.Fit(fLastFit, "Q0NR");

    result.b0 = fLastFit->GetParameter(0);
    result.b1 = fLastFit->GetParameter(1);
    result.b2 = fLastFit->GetParameter(2);
    result.success = true;
    return result;
}

LKSiCalibratorC2::Result LKSiCalibratorC2::Fit(const std::vector<TH2 *> &positionEnergyHists, int chooseGate, double entriesCut)
{
    Result result;
    result.selectedGate = chooseGate;
    for (auto hist : positionEnergyHists)
        result.gates.push_back(FitSingle(hist));

    if (chooseGate < 0 || chooseGate >= (int) result.gates.size())
        return result;
    auto gate = result.gates[chooseGate];
    if (!gate.success || gate.entries < entriesCut)
        return result;

    result.entries = gate.entries;
    result.b0 = gate.b0;
    result.b1 = gate.b1;
    result.b2 = gate.b2;
    result.success = true;
    return result;
}
