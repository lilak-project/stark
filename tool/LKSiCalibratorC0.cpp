#include "LKSiCalibratorC0.h"

#include "TF1.h"
#include "TGraphErrors.h"
#include "TH1.h"
#include "TSpectrum.h"

#include <algorithm>

ClassImp(LKSiCalibratorC0)

LKSiCalibratorC0::LKSiCalibratorC0()
    : TNamed("LKSiCalibratorC0", "")
{
}

LKSiCalibratorC0::Result LKSiCalibratorC0::Fit(TH1 *hist, const std::vector<double> &gateEnergies, double expectedResolution, double entriesCut)
{
    Result result;
    if (hist == nullptr)
        return result;

    result.entries = hist->GetEntries();
    if (result.entries < entriesCut)
        return result;

    TSpectrum spectrum((int) gateEnergies.size() + 2);
    auto numPeaks = spectrum.Search(hist, 5, "goff nodraw");
    if (numPeaks < (int) gateEnergies.size())
        return result;

    std::vector<double> peaks(spectrum.GetPositionX(), spectrum.GetPositionX() + numPeaks);
    std::sort(peaks.begin(), peaks.end());

    TGraphErrors graph;
    TF1 fitPeak("fitPeak", "gaus(0)", hist->GetXaxis()->GetXmin(), hist->GetXaxis()->GetXmax());
    for (size_t iGate = 0; iGate < gateEnergies.size(); ++iGate) {
        auto peak = peaks[iGate];
        fitPeak.SetRange(peak - 5.0 * peak * expectedResolution, peak + 5.0 * peak * expectedResolution);
        fitPeak.SetParameters(hist->GetBinContent(hist->FindBin(peak)), peak, peak * expectedResolution);
        hist->Fit(&fitPeak, "Q0NR");

        auto mean = fitPeak.GetParameter(1);
        auto sigma = fitPeak.GetParameter(2);
        fitPeak.SetRange(mean - 1.0 * sigma, mean + 2.5 * sigma);
        hist->Fit(&fitPeak, "Q0NR");

        mean = fitPeak.GetParameter(1);
        sigma = fitPeak.GetParameter(2);
        result.peakMeans.push_back(mean);
        result.peakSigmas.push_back(sigma);
        graph.SetPoint((int) iGate, mean, gateEnergies[iGate]);
        graph.SetPointError((int) iGate, sigma, gateEnergies[iGate] * 0.01);
    }

    TF1 fitLine("fitLine", "pol1", hist->GetXaxis()->GetXmin(), hist->GetXaxis()->GetXmax());
    graph.Fit(&fitLine, "Q0R");
    result.intercept = fitLine.GetParameter(0);
    result.slope = fitLine.GetParameter(1);
    result.success = (graph.GetN() >= 2);
    return result;
}
