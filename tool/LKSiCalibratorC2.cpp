#include "LKSiCalibratorC2.h"

#include "TF1.h"
#include "TH2.h"

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

LKSiCalibratorC2::GateResult LKSiCalibratorC2::FitSingle(TH2 *hist)
{
    GateResult result;
    if (hist == nullptr)
        return result;

    result.entries = hist->GetEntries();
    if (result.entries <= 0)
        return result;

    Clear();
    fLastFit = new TF1(Form("fit_%s", hist->GetName()), "pol2", -1, 1);
    fLastFit->SetParameters(5.48556, 0.0, 0.2);
    hist->Fit(fLastFit, "Q0R");

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
