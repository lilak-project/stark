#include "SKSiCalibrationTask.h"

#include "LKLogger.h"
#include "LKDrawing.h"
#include "LKParameterContainer.h"
#include "LKRun.h"
#include "LKDrawingGroup.h"
#include "LKSiCalibratorC0.h"
#include "LKSiCalibratorC1.h"
#include "LKSiCalibratorC2.h"
#include "LKSiChannel.h"
#include "LKSiliconArray.h"

#include "TFile.h"
#include "TF1.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TList.h"
#include "TLine.h"
#include "TPaveText.h"
#include "TProfile.h"
#include "TSystem.h"

#include <algorithm>
#include <fstream>
#include <iomanip>
#include <map>
#include <sstream>
#include <string>
#include <vector>

namespace {
// A TF1 registers itself in gROOT->GetListOfFunctions(), which ROOT deletes at the end of the
// process. The histograms that are drawn in the saved canvases keep the same pointers in their
// own list of functions, so that deletion leaves dangling entries behind and the clean-up sweep
// of ROOT crashes. Hand the ownership to the histogram, which deletes its functions itself.
void AttachFunction(TH1 *hist, TF1 *function)
{
    function -> AddToGlobalList(false);
    hist -> GetListOfFunctions() -> Add(function);
}

template <typename HistMap>
void AddHistByDetector(LKDrawingGroup *group, HistMap &histMap, const char *prefix, Option_t *option = "")
{
    std::map<int, LKDrawingGroup *> groupMap;
    for (auto &entry : histMap) {
        auto det = std::get<0>(entry.first);
        auto detGroup = groupMap[det];
        if (detGroup == nullptr) {
            detGroup = group -> CreateGroup(Form("%s%d", prefix, det));
            groupMap[det] = detGroup;
        }

        std::vector<TObject *> paveTextArray;
        auto functions = entry.second -> GetListOfFunctions();
        for (auto object : *functions) {
            if (object -> InheritsFrom(TPaveText::Class()))
                paveTextArray.push_back(object);
        }

        if (paveTextArray.empty()) {
            detGroup -> AddHist(entry.second, option);
            continue;
        }

        auto drawing = new LKDrawing();
        drawing -> SetPaveCorner(1);
        drawing -> Add(entry.second, option);
        for (auto object : paveTextArray) {
            functions -> Remove(object);
            drawing -> Add(object, "same");
        }
        detGroup -> AddDrawing(drawing);
    }
}
}

ClassImp(SKSiCalibrationTask)

SKSiCalibrationTask::SKSiCalibrationTask()
    : LKTask("SKSiCalibrationTask", "SKSiCalibrationTask")
{
}

SKSiCalibrationTask::~SKSiCalibrationTask()
{
    for (auto &entry : fEnergyHistMap) delete entry.second;
    for (auto &entry : fEnergySumHistMap) delete entry.second;
    for (auto &entry : fLeftRightHistMap) delete entry.second;
    for (auto &entry : fPositionEnergyHistMap) delete entry.second;
    for (auto &entry : fFitHistArray) delete entry.second;
}

bool SKSiCalibrationTask::Init()
{
    fSiChannelArray = fRun -> GetBranchA("SiChannel", "LKSiChannel");
    if (fSiChannelArray == nullptr) {
        lk_error << "Branch SiChannel does not exist. Run LKSetSiChannelTask first or use a reco file with SiChannel." << endl;
        return false;
    }

    fSiArray = (LKSiliconArray *) fRun -> FindDetectorPlane("LKSiliconArray");

    fPar -> UpdatePar(fStage, "SKSiCalibrationTask/stage");
    fPar -> UpdatePar(fSourceEnergy, "SKSiCalibrationTask/sourceEnergy");
    fPar -> UpdatePar(fSourceEnergyText, "SKSiCalibrationTask/sourceEnergies");
    fPar -> UpdatePar(fC2ReferenceEnergy, "SKSiCalibrationTask/c2ReferenceEnergy");
    fPar -> UpdatePar(fExpectedResolution, "SKSiCalibrationTask/expectedResolution");
    fPar -> UpdatePar(fGateDisplayRangeScale, "SKSiCalibrationTask/gateDisplayRangeScale");
    fPar -> UpdatePar(fEntriesCut, "SKSiCalibrationTask/entriesCut");
    fPar -> UpdatePar(fOutputPath, "SKSiCalibrationTask/outputPath");
    fPar -> UpdatePar(fOutputPrefix, "SKSiCalibrationTask/outputPrefix");
    fPar -> UpdatePar(fInputC0FileName, "SKSiCalibrationTask/inputC0File");
    fPar -> UpdatePar(fInputC1FileName, "SKSiCalibrationTask/inputC1File");
    fPar -> UpdatePar(fInputC2FileName, "SKSiCalibrationTask/inputC2File");
    fPar -> UpdatePar(fInputC3FileName, "SKSiCalibrationTask/inputC3File");
    fPar -> UpdatePar(fListFileName, "SKSiCalibrationTask/listFile");

    if (!fSourceEnergyText.IsNull() && !ParseSourceEnergies(fSourceEnergyText))
        return false;

    if (fSourceEnergies.empty() && fSourceEnergy > 0)
        fSourceEnergies.push_back(fSourceEnergy);

    if (fSourceEnergies.empty() && !fListFileName.IsNull()) {
        std::ifstream listFile(fListFileName.Data());
        int runID = -1;
        TString detName;
        TString mappingName;
        TString voltage;
        TString sourceName;
        double sourceEnergy = 0;
        while (listFile >> runID >> detName >> mappingName >> voltage >> sourceName >> sourceEnergy) {
            if (runID == fRun->GetRunID()) {
                fSourceEnergy = sourceEnergy;
                fSourceEnergies.push_back(sourceEnergy);
                break;
            }
        }
    }

    if (fSourceEnergies.empty()) {
        lk_error << "Set SKSiCalibrationTask/sourceEnergies, SKSiCalibrationTask/sourceEnergy, or SKSiCalibrationTask/listFile." << endl;
        return false;
    }
    std::sort(fSourceEnergies.begin(), fSourceEnergies.end());
    fSourceEnergy = fSourceEnergies.back();

    if (fOutputPrefix.IsNull())
        fOutputPrefix = Form("pre_li7_%04d", fRun->GetRunID());

    gSystem -> mkdir(fOutputPath, true);
    if (!LoadParametersForStage())
        return false;
    InitializeExpectedChannels();

    lk_info << "SKSiCalibrationTask stage=" << fStage << ", sourceEnergy=" << fSourceEnergy
            << ", sourceEnergies=" << GetSourceEnergyListText()
            << ", outputPrefix=" << fOutputPrefix << endl;
    return true;
}

void SKSiCalibrationTask::Exec(Option_t *)
{
    auto numChannels = fSiChannelArray -> GetEntriesFast();
    for (auto iChannel = 0; iChannel < numChannels; ++iChannel) {
        auto channel = (LKSiChannel *) fSiChannelArray -> At(iChannel);
        if (channel == nullptr)
            continue;

        auto det = channel -> GetDetNum();
        auto side = channel -> GetSide();
        auto strip = channel -> GetStrip();

        if (channel -> IsStandaloneChannel()) {
            if (!(fStage == 0 || fStage == 1 || fStage >= 6))
                continue;
            auto energy = channel -> GetEnergy();
            if (fStage >= 6) {
                auto it = fC0Parameters.find(StripKey(det, side, strip));
                if (it != fC0Parameters.end())
                    energy = ApplyLinear(it->second, energy);
            }
            GetEnergyHist(det, side, strip) -> Fill(energy);
            continue;
        }

        if (!(fStage == 0 || fStage >= 2))
            continue;

        if (channel -> GetDirection() != 0 || channel -> GetEnergy2() < 0)
            continue;

        auto energyR = channel -> GetEnergy();
        auto energyL = channel -> GetEnergy2();
        if (fStage >= 3) {
            auto itL = fC1Parameters.find(LRKey(det, side, strip, 0));
            auto itR = fC1Parameters.find(LRKey(det, side, strip, 1));
            if (itL != fC1Parameters.end()) energyL = ApplyLinear(itL->second, energyL);
            if (itR != fC1Parameters.end()) energyR = ApplyLinear(itR->second, energyR);
        }

        auto energySum = energyL + energyR;
        if (energySum <= 0)
            continue;
        auto position = (energyR - energyL) / energySum;
        if (fStage >= 4)
            ApplyC2(det, side, strip, position, energySum);
        if (fStage >= 6)
            ApplyC3(det, side, strip, energySum);

        GetLeftRightHist(det, side, strip) -> Fill(energyR, energyL);
        if (fStage != 2) {
            GetEnergySumHist(det, side, strip) -> Fill(energySum);
            GetPositionEnergyHist(det, side, strip) -> Fill(position, energySum);
        }
    }
}

bool SKSiCalibrationTask::EndOfRun()
{
    WriteChannelSummary();

    if (fStage == 1) {
        FitAndWriteC0();
    }
    else if (fStage == 2) {
        FitAndWriteC1();
    }
    else if (fStage == 3) {
        FitAndWriteC2();
        ApplyStageEnergyZoom();
    }
    else if (fStage == 4) {
        ApplyStageEnergyZoom();
    }
    else if (fStage == 5) {
        FitAndWriteC3();
        ApplyStageEnergyZoom();
    }
    else if (fStage >= 6) {
        WriteCombinedCalibrationFile();
        ApplyStageEnergyZoom();
    }

    WriteHistograms();

    return true;
}

TH1D *SKSiCalibrationTask::GetEnergyHist(int det, int side, int strip)
{
    StripKey key(det, side, strip);
    auto it = fEnergyHistMap.find(key);
    if (it != fEnergyHistMap.end())
        return it->second;

    auto name = MakeHistName("energy", det, side, strip);
    auto eMax = fStage >= 6 ? GetCalibratedEnergyAxisMax() : 6000.;
    auto hist = new TH1D(name, Form("%s;energy;counts", name.Data()), 500, 0, eMax);
    hist -> SetDirectory(nullptr);
    fEnergyHistMap[key] = hist;
    return hist;
}

TH1D *SKSiCalibrationTask::GetEnergySumHist(int det, int side, int strip)
{
    StripKey key(det, side, strip);
    auto it = fEnergySumHistMap.find(key);
    if (it != fEnergySumHistMap.end())
        return it->second;

    auto name = MakeHistName("esum", det, side, strip);
    auto eMax = fStage >= 3 ? GetCalibratedEnergyAxisMax() : 6000.;
    auto hist = new TH1D(name, Form("%s;energy sum;counts", name.Data()), 500, 0, eMax);
    hist -> SetDirectory(nullptr);
    fEnergySumHistMap[key] = hist;
    return hist;
}

TH2D *SKSiCalibrationTask::GetLeftRightHist(int det, int side, int strip)
{
    StripKey key(det, side, strip);
    auto it = fLeftRightHistMap.find(key);
    if (it != fLeftRightHistMap.end())
        return it->second;

    auto name = MakeHistName("left_right", det, side, strip);
    auto eMax = fStage >= 3 ? GetCalibratedEnergyAxisMax() : 6000.;
    auto hist = new TH2D(name, Form("%s;right;left", name.Data()), 500, 0, eMax, 500, 0, eMax);
    hist -> SetDirectory(nullptr);
    fLeftRightHistMap[key] = hist;
    return hist;
}

TH2D *SKSiCalibrationTask::GetPositionEnergyHist(int det, int side, int strip)
{
    StripKey key(det, side, strip);
    auto it = fPositionEnergyHistMap.find(key);
    if (it != fPositionEnergyHistMap.end())
        return it->second;

    auto name = MakeHistName("rpos_esum", det, side, strip);
    auto eMax = fStage >= 3 ? GetCalibratedEnergyAxisMax() : 6000.;
    auto hist = new TH2D(name, Form("%s;relative position;energy sum", name.Data()), 240, -1.2, 1.2, 500, 0, eMax);
    hist -> SetDirectory(nullptr);
    fPositionEnergyHistMap[key] = hist;
    return hist;
}

bool SKSiCalibrationTask::LoadC0Parameters(TString fileName)
{
    std::ifstream input(fileName.Data());
    if (!input.is_open())
        return false;
    std::string line;
    std::getline(input, line);
    int det = -1, side = -1, strip = -1;
    double entries = 0, intercept = 0, slope = 0;
    while (input >> det >> side >> strip >> entries >> intercept >> slope)
        fC0Parameters[StripKey(det, side, strip)] = {intercept, slope};
    return true;
}

bool SKSiCalibrationTask::LoadC1Parameters(TString fileName)
{
    std::ifstream input(fileName.Data());
    if (!input.is_open())
        return false;
    std::string line;
    std::getline(input, line);
    int det = -1, side = -1, strip = -1, nPoints = 0;
    double interceptL = 0, slopeL = 0, interceptR = 0, slopeR = 0;
    while (input >> det >> side >> strip >> nPoints >> interceptL >> slopeL >> interceptR >> slopeR) {
        fC1Parameters[LRKey(det, side, strip, 0)] = {interceptL, slopeL};
        fC1Parameters[LRKey(det, side, strip, 1)] = {interceptR, slopeR};
    }
    return true;
}

bool SKSiCalibrationTask::LoadC2Parameters(TString fileName)
{
    std::ifstream input(fileName.Data());
    if (!input.is_open())
        return false;
    std::string line;
    std::getline(input, line);
    int det = -1, side = -1, strip = -1;
    double entries = 0, b0 = 0, b1 = 0, b2 = 0;
    while (input >> det >> side >> strip >> entries >> b0 >> b1 >> b2)
        fC2Parameters[StripKey(det, side, strip)] = {b0, b1, b2};
    return true;
}

bool SKSiCalibrationTask::LoadC3Parameters(TString fileName)
{
    std::ifstream input(fileName.Data());
    if (!input.is_open())
        return false;
    std::string line;
    std::getline(input, line);
    int det = -1, side = -1, strip = -1;
    double entries = 0, intercept = 0, slope = 0;
    while (input >> det >> side >> strip >> entries >> intercept >> slope)
        fC3Parameters[StripKey(det, side, strip)] = {intercept, slope};
    return true;
}

bool SKSiCalibrationTask::LoadParametersForStage()
{
    if (fInputC0FileName.IsNull()) fInputC0FileName = MakeFileName("c0.par");
    if (fInputC1FileName.IsNull()) fInputC1FileName = MakeFileName("c1.par");
    if (fInputC2FileName.IsNull()) fInputC2FileName = MakeFileName("c2.par");
    if (fInputC3FileName.IsNull()) fInputC3FileName = MakeFileName("c3.par");

    if (fStage >= 2) {
        if (!LoadC0Parameters(fInputC0FileName)) {
            lk_error << "Cannot load C0 parameter file " << fInputC0FileName
                     << ". Run stage 1 first or set SKSiCalibrationTask/inputC0File." << endl;
            return false;
        }
        lk_info << "Loaded C0 parameter file " << fInputC0FileName << endl;
    }
    if (fStage >= 3) {
        if (!LoadC1Parameters(fInputC1FileName)) {
            lk_error << "Cannot load C1 parameter file " << fInputC1FileName
                     << ". Run stage 2 first or set SKSiCalibrationTask/inputC1File." << endl;
            return false;
        }
        lk_info << "Loaded C1 parameter file " << fInputC1FileName << endl;
    }
    if (fStage >= 4) {
        if (!LoadC2Parameters(fInputC2FileName)) {
            lk_error << "Cannot load C2 parameter file " << fInputC2FileName
                     << ". Run stage 3 first or set SKSiCalibrationTask/inputC2File." << endl;
            return false;
        }
        lk_info << "Loaded C2 parameter file " << fInputC2FileName << endl;
    }
    if (fStage >= 6) {
        if (!LoadC3Parameters(fInputC3FileName)) {
            lk_error << "Cannot load C3 parameter file " << fInputC3FileName
                     << ". Run stage 5 first or set SKSiCalibrationTask/inputC3File." << endl;
            return false;
        }
        lk_info << "Loaded C3 parameter file " << fInputC3FileName << endl;
    }

    return true;
}

double SKSiCalibrationTask::ApplyLinear(const std::array<double,2> &par, double value) const
{
    return par[0] + par[1] * value;
}

bool SKSiCalibrationTask::ApplyC2(int det, int side, int strip, double position, double &energy) const
{
    auto it = fC2Parameters.find(StripKey(det, side, strip));
    if (it == fC2Parameters.end())
        return false;
    auto scale = it->second[0] + it->second[1] * position + it->second[2] * position * position;
    if (scale == 0)
        return false;
    energy = energy / scale * fC2ReferenceEnergy;
    return true;
}

bool SKSiCalibrationTask::ApplyC3(int det, int side, int strip, double &energy) const
{
    auto it = fC3Parameters.find(StripKey(det, side, strip));
    if (it == fC3Parameters.end())
        return false;
    energy = ApplyLinear(it->second, energy);
    return true;
}

bool SKSiCalibrationTask::ParseSourceEnergies(TString sourceEnergies)
{
    sourceEnergies.ReplaceAll(",", " ");
    sourceEnergies.ReplaceAll(":", " ");
    sourceEnergies.ReplaceAll(";", " ");

    std::stringstream stream(sourceEnergies.Data());
    double energy = 0;
    while (stream >> energy) {
        if (energy > 0)
            fSourceEnergies.push_back(energy);
    }

    if (fSourceEnergies.empty()) {
        lk_error << "Cannot parse SKSiCalibrationTask/sourceEnergies = " << sourceEnergies << endl;
        return false;
    }
    return true;
}

double SKSiCalibrationTask::FitPeak(TH1D *hist, TString fitName) const
{
    if (hist == nullptr || hist->GetEntries() < fEntriesCut)
        return 0;

    hist -> SetStats(0);
    auto peak = hist -> GetBinCenter(hist -> GetMaximumBin());
    if (peak <= 0)
        return 0;

    auto sigma0 = std::max(peak * fExpectedResolution, hist->GetXaxis()->GetBinWidth(1));
    TF1 fit(fitName, "gaus", peak - 5 * sigma0, peak + 5 * sigma0);
    fit.SetParameters(hist->GetMaximum(), peak, sigma0);
    hist -> Fit(&fit, "Q0NR");

    auto mean = fit.GetParameter(1);
    auto sigma = fit.GetParameter(2);
    if (sigma > 0) {
        fit.SetRange(mean - sigma, mean + 2.5 * sigma);
        hist -> Fit(&fit, "Q0NR");
        mean = fit.GetParameter(1);
        sigma = fit.GetParameter(2);
    }

    auto fitDraw = (TF1 *) fit.Clone(fitName);
    fitDraw -> SetLineColor(kRed+1);
    fitDraw -> SetLineWidth(2);
    fitDraw -> SetNpx(500);
    AttachFunction(hist, fitDraw);

    auto text = new TPaveText(0.55, 0.66, 0.88, 0.88, "NDC");
    text -> SetName(Form("fit_info_%s", hist->GetName()));
    text -> SetTextFont(132);
    text -> SetTextAlign(12);
    text -> AddText(Form("entries = %.0f", hist->GetEntries()));
    text -> AddText(Form("amp. = %.2f", fit.GetParameter(0)));
    text -> AddText(Form("mean = %.2f", mean));
    text -> AddText(Form("sigma = %.4f", sigma));
    if (mean != 0)
        text -> AddText(Form("res. = %.4f", sigma / mean));
    ApplyPaveTextStyle(text);
    hist -> GetListOfFunctions() -> Add(text);

    return mean;
}

void SKSiCalibrationTask::ApplyPaveTextStyle(TPaveText *text) const
{
    if (text == nullptr)
        return;

    text -> SetFillStyle(0);
    text -> SetFillColorAlpha(kWhite, 0);
    text -> SetBorderSize(1);
}

void SKSiCalibrationTask::AddC0FitAnnotations(TH1D *hist, const std::vector<double> &amps, const std::vector<double> &means, const std::vector<double> &sigmas) const
{
    if (hist == nullptr)
        return;

    hist -> SetStats(0);
    auto numPeaks = std::min(amps.size(), std::min(means.size(), sigmas.size()));
    if (numPeaks == 0)
        return;

    auto text = new TPaveText(0.43, 0.58, 0.88, 0.88, "NDC");
    text -> SetName(Form("fit_info_%s", hist->GetName()));
    text -> SetTextFont(132);
    text -> SetTextAlign(12);
    text -> AddText(Form("entries = %.0f", hist->GetEntries()));

    for (size_t iPeak = 0; iPeak < numPeaks; ++iPeak) {
        auto fit = new TF1(Form("fit_c0_%s_peak%zu", hist->GetName(), iPeak), "gaus(0)",
                           means[iPeak] - sigmas[iPeak], means[iPeak] + 2.5 * sigmas[iPeak]);
        fit -> SetParameters(amps[iPeak], means[iPeak], sigmas[iPeak]);
        fit -> SetLineColor(kRed + (int) iPeak);
        fit -> SetLineWidth(2);
        fit -> SetNpx(500);
        AttachFunction(hist, fit);

        text -> AddText(Form("p%zu: m=%.2f s=%.4f r=%.4f",
                             iPeak, means[iPeak], sigmas[iPeak],
                             means[iPeak] != 0 ? sigmas[iPeak] / means[iPeak] : 0));
    }

    ApplyPaveTextStyle(text);
    hist -> GetListOfFunctions() -> Add(text);
}

void SKSiCalibrationTask::AddC2FitAnnotation(TH2D *hist, const std::array<double,3> &fitPar, const std::array<double,3> &c2Par, double entries) const
{
    if (hist == nullptr)
        return;

    hist -> SetStats(0);
    auto fit = new TF1(Form("fit_%s", hist->GetName()), "pol2", -1, 1);
    fit -> SetParameters(fitPar[0], fitPar[1], fitPar[2]);
    fit -> SetLineColor(kRed+1);
    fit -> SetLineWidth(2);
    fit -> SetNpx(500);
    AttachFunction(hist, fit);

    auto text = new TPaveText(0.44, 0.52, 0.88, 0.88, "NDC");
    text -> SetName(Form("fit_info_%s", hist->GetName()));
    text -> SetTextFont(132);
    text -> SetTextAlign(12);
    text -> AddText(Form("entries = %.0f", entries));
    text -> AddText(Form("fit: %.5f %+.5fx %+.5fx^{2}", fitPar[0], fitPar[1], fitPar[2]));
    text -> AddText(Form("C2 b0 = %.5f", c2Par[0]));
    text -> AddText(Form("C2 b1 = %.5f", c2Par[1]));
    text -> AddText(Form("C2 b2 = %.5f", c2Par[2]));
    ApplyPaveTextStyle(text);
    hist -> GetListOfFunctions() -> Add(text);
}

bool SKSiCalibrationTask::GetGateDisplayRange(double axisMin, double axisMax, double &rangeMin, double &rangeMax) const
{
    if (fSourceEnergies.empty())
        return false;

    auto minEnergy = *std::min_element(fSourceEnergies.begin(), fSourceEnergies.end());
    auto maxEnergy = *std::max_element(fSourceEnergies.begin(), fSourceEnergies.end());
    auto gateWidth = maxEnergy * fExpectedResolution * 5.;
    auto span = maxEnergy - minEnergy;
    auto margin = gateWidth;
    if (span <= 0)
        margin = std::max(1., fGateDisplayRangeScale) * gateWidth;
    else
        margin = 0.5 * std::max(1., fGateDisplayRangeScale - 1.) * span;

    rangeMin = std::max(axisMin, minEnergy - margin);
    rangeMax = std::min(axisMax, maxEnergy + margin);
    return rangeMax > rangeMin;
}

void SKSiCalibrationTask::SetEnergyGateDisplayRange(TH1D *hist) const
{
    if (hist == nullptr)
        return;

    double rangeMin = 0;
    double rangeMax = 0;
    if (GetGateDisplayRange(hist->GetXaxis()->GetXmin(), hist->GetXaxis()->GetXmax(), rangeMin, rangeMax))
        hist -> GetXaxis() -> SetRangeUser(rangeMin, rangeMax);
}

void SKSiCalibrationTask::SetEnergyGateDisplayRange(TH2D *hist) const
{
    if (hist == nullptr)
        return;

    double rangeMin = 0;
    double rangeMax = 0;
    if (GetGateDisplayRange(hist->GetYaxis()->GetXmin(), hist->GetYaxis()->GetXmax(), rangeMin, rangeMax))
        hist -> GetYaxis() -> SetRangeUser(rangeMin, rangeMax);
}

void SKSiCalibrationTask::SetEnergyGateDisplayRangeXY(TH2D *hist) const
{
    if (hist == nullptr)
        return;

    double rangeMin = 0;
    double rangeMax = 0;
    if (GetGateDisplayRange(hist->GetXaxis()->GetXmin(), hist->GetXaxis()->GetXmax(), rangeMin, rangeMax))
        hist -> GetXaxis() -> SetRangeUser(rangeMin, rangeMax);
    if (GetGateDisplayRange(hist->GetYaxis()->GetXmin(), hist->GetYaxis()->GetXmax(), rangeMin, rangeMax))
        hist -> GetYaxis() -> SetRangeUser(rangeMin, rangeMax);
}

void SKSiCalibrationTask::AddEnergyPositionGateLines(TH2D *hist) const
{
    if (hist == nullptr)
        return;

    SetEnergyGateDisplayRange(hist);

    auto x1 = hist -> GetXaxis() -> GetXmin();
    auto x2 = hist -> GetXaxis() -> GetXmax();
    for (auto iEnergy = 0u; iEnergy < fSourceEnergies.size(); ++iEnergy) {
        auto energy = fSourceEnergies[iEnergy];
        auto halfWidth = std::max(energy * fExpectedResolution * 5., hist->GetYaxis()->GetBinWidth(1));
        auto line = new TLine(x1, energy, x2, energy);
        line -> SetLineColor(kRed + (int) iEnergy);
        line -> SetLineWidth(2);
        hist -> GetListOfFunctions() -> Add(line);

        auto lineLow = new TLine(x1, energy - halfWidth, x2, energy - halfWidth);
        lineLow -> SetLineColor(kRed + (int) iEnergy);
        lineLow -> SetLineStyle(2);
        hist -> GetListOfFunctions() -> Add(lineLow);

        auto lineHigh = new TLine(x1, energy + halfWidth, x2, energy + halfWidth);
        lineHigh -> SetLineColor(kRed + (int) iEnergy);
        lineHigh -> SetLineStyle(2);
        hist -> GetListOfFunctions() -> Add(lineHigh);
    }

    std::vector<TObject *> textObjects;
    auto functions = hist -> GetListOfFunctions();
    for (auto object : *functions) {
        if (object -> InheritsFrom(TPaveText::Class()))
            textObjects.push_back(object);
    }
    for (auto object : textObjects) {
        functions -> Remove(object);
        functions -> Add(object);
    }
}

void SKSiCalibrationTask::AddEnergyPositionGateLines()
{
    for (auto const &entry : fPositionEnergyHistMap)
        AddEnergyPositionGateLines(entry.second);
}

void SKSiCalibrationTask::ApplyStageEnergyZoom()
{
    for (auto const &entry : fEnergyHistMap)
        SetEnergyGateDisplayRange(entry.second);
    for (auto const &entry : fEnergySumHistMap)
        SetEnergyGateDisplayRange(entry.second);
    for (auto const &entry : fPositionEnergyHistMap)
        AddEnergyPositionGateLines(entry.second);
}

void SKSiCalibrationTask::FitAndWriteC0()
{
    std::ofstream output(MakeFileName("c0.par").Data());
    output << "# det side strip entries intercept slope\n";
    output << std::setprecision(12);

    LKSiCalibratorC0 calibrator;
    for (auto const &entry : fEnergyHistMap) {
        auto det = std::get<0>(entry.first);
        auto side = std::get<1>(entry.first);
        auto strip = std::get<2>(entry.first);
        auto result = calibrator.Fit(entry.second, fSourceEnergies, fExpectedResolution, fEntriesCut);
        if (!result.success) {
            auto peak = FitPeak(entry.second, Form("fit_c0_det%03d_side%d_strip%02d", det, side, strip));
            result.entries = entry.second -> GetEntries();
            result.intercept = 0;
            result.slope = peak > 0 ? fSourceEnergy / peak : 0;
        }
        else
            AddC0FitAnnotations(entry.second, result.peakAmps, result.peakMeans, result.peakSigmas);
        fC0Parameters[entry.first] = {result.intercept, result.slope};
        output << det << " " << side << " " << strip << " "
               << result.entries << " " << result.intercept << " " << result.slope << "\n";
    }
    WriteC0ParameterContainer();
}

void SKSiCalibrationTask::FitAndWriteC1()
{
    std::ofstream output(MakeFileName("c1.par").Data());
    output << "# det side strip n_points intercept_left slope_left intercept_right slope_right\n";
    output << std::setprecision(12);

    for (auto const &entry : fLeftRightHistMap) {
        auto det = std::get<0>(entry.first);
        auto side = std::get<1>(entry.first);
        auto strip = std::get<2>(entry.first);
        auto hist = entry.second;
        auto entries = hist -> GetEntries();

        double intercept = 0;
        double slope = 0;
        double xIntercept = 0;
        double yIntercept = 0;
        if (entries >= fEntriesCut) {
            auto profile = hist -> ProfileX(Form("%s_profile", hist->GetName()));
            profile -> SetDirectory(nullptr);
            TF1 fit(Form("fit_c1_slope_det%03d_side%d_strip%02d", det, side, strip), "pol1",
                    hist->GetXaxis()->GetXmin(), hist->GetXaxis()->GetXmax());
            profile -> Fit(&fit, "Q0NR");
            intercept = fit.GetParameter(0);
            slope = fit.GetParameter(1);
            if (slope != 0) {
                yIntercept = intercept;
                xIntercept = -intercept / slope;
            }

            auto fitDraw = (TF1 *) fit.Clone(Form("fit_c1_slope_det%03d_side%d_strip%02d", det, side, strip));
            fitDraw -> SetLineColor(kRed+1);
            fitDraw -> SetLineWidth(2);
            fitDraw -> SetNpx(500);
            hist -> SetStats(0);
            AttachFunction(hist, fitDraw);

            auto text = new TPaveText(0.48, 0.66, 0.88, 0.88, "NDC");
            text -> SetName(Form("fit_info_%s", hist->GetName()));
            text -> SetTextFont(132);
            text -> SetTextAlign(12);
            text -> AddText(Form("entries = %.0f", entries));
            text -> AddText(Form("left = %.2f + %.4f right", intercept, slope));
            text -> AddText(Form("right intercept = %.2f", xIntercept));
            text -> AddText(Form("left intercept = %.2f", yIntercept));
            ApplyPaveTextStyle(text);
            hist -> GetListOfFunctions() -> Add(text);
            delete profile;
        }

        auto slopeL = yIntercept > 0 ? fSourceEnergy / yIntercept : 0;
        auto slopeR = xIntercept > 0 ? fSourceEnergy / xIntercept : 0;
        auto nPoints = (slopeL > 0 && slopeR > 0) ? 1 : 0;
        fC1Parameters[LRKey(det, side, strip, 0)] = {0, slopeL};
        fC1Parameters[LRKey(det, side, strip, 1)] = {0, slopeR};
        output << det << " " << side << " " << strip << " " << nPoints << " "
               << 0 << " " << slopeL << " " << 0 << " " << slopeR << "\n";
    }
    WriteC1ParameterContainer();
}

void SKSiCalibrationTask::FitAndWriteC2()
{
    std::ofstream output(MakeFileName("c2.par").Data());
    output << "# det side strip entries b0 b1 b2\n";
    output << std::setprecision(12);

    LKSiCalibratorC2 calibrator;
    for (auto const &entry : fPositionEnergyHistMap) {
        auto det = std::get<0>(entry.first);
        auto side = std::get<1>(entry.first);
        auto strip = std::get<2>(entry.first);
        auto result = calibrator.FitSingle(entry.second, fSourceEnergy, fExpectedResolution, fGateDisplayRangeScale);
        auto fitSucceeded = result.success && result.entries >= fEntriesCut;
        if (!fitSucceeded) {
            result.entries = entry.second -> GetEntries();
            result.b0 = fSourceEnergy > 0 ? fSourceEnergy : fC2ReferenceEnergy;
            result.b1 = 0;
            result.b2 = 0;
        }

        std::array<double,3> fitPar = {result.b0, result.b1, result.b2};
        if (fitSucceeded && fSourceEnergy > 0) {
            auto factor = fC2ReferenceEnergy / fSourceEnergy;
            result.b0 *= factor;
            result.b1 *= factor;
            result.b2 *= factor;
        }
        else if (!fitSucceeded) {
            result.b0 = fC2ReferenceEnergy;
            result.b1 = 0;
            result.b2 = 0;
        }
        std::array<double,3> c2Par = {result.b0, result.b1, result.b2};
        AddC2FitAnnotation(entry.second, fitPar, c2Par, result.entries);
        fC2Parameters[entry.first] = c2Par;
        output << det << " " << side << " " << strip << " "
               << result.entries << " " << result.b0 << " " << result.b1 << " " << result.b2 << "\n";
    }
    WriteC2ParameterContainer();
}

void SKSiCalibrationTask::FitAndWriteC3()
{
    std::ofstream output(MakeFileName("c3.par").Data());
    output << "# det side strip entries intercept slope\n";
    output << std::setprecision(12);

    LKSiCalibratorC0 calibrator;
    for (auto const &entry : fEnergySumHistMap) {
        auto det = std::get<0>(entry.first);
        auto side = std::get<1>(entry.first);
        auto strip = std::get<2>(entry.first);
        auto hist = entry.second;

        auto result = calibrator.Fit(hist, fSourceEnergies, fExpectedResolution, fEntriesCut);
        if (!result.success) {
            auto peak = FitPeak(hist, Form("fit_c3_det%03d_side%d_strip%02d", det, side, strip));
            result.entries = hist -> GetEntries();
            result.intercept = 0;
            result.slope = peak > 0 ? fSourceEnergy / peak : 1;
        }
        else
            AddC0FitAnnotations(hist, result.peakAmps, result.peakMeans, result.peakSigmas);

        fC3Parameters[entry.first] = {result.intercept, result.slope};
        output << det << " " << side << " " << strip << " "
               << result.entries << " " << result.intercept << " " << result.slope << "\n";
    }

    WriteC3ParameterContainer();
}

void SKSiCalibrationTask::WriteHistograms()
{
    auto taskGroup = fRun -> GetTopDrawingGroup();
    taskGroup -> SetName("top");

    if (!fEnergyHistMap.empty()) {
        auto energyGroup = taskGroup -> CreateGroup(Form("stage%d_Energy", fStage));
        AddHistByDetector(energyGroup, fEnergyHistMap, "E");
    }

    if (!fEnergySumHistMap.empty()) {
        auto esumGroup = taskGroup -> CreateGroup(Form("stage%d_EnergySum", fStage));
        AddHistByDetector(esumGroup, fEnergySumHistMap, "ES");
    }

    if (!fLeftRightHistMap.empty()) {
        auto leftRightGroup = taskGroup -> CreateGroup(Form("stage%d_LeftRight", fStage));
        AddHistByDetector(leftRightGroup, fLeftRightHistMap, "LR", "colz");
    }

    if (!fPositionEnergyHistMap.empty()) {
        auto rposEsumGroup = taskGroup -> CreateGroup(Form("stage%d_EnergyPosition", fStage));
        AddHistByDetector(rposEsumGroup, fPositionEnergyHistMap, "EP", "colz");
    }

    if (!fFitHistArray.empty()) {
        auto fitGroup = taskGroup -> CreateGroup(Form("stage%d_Fit", fStage));
        std::map<int, LKDrawingGroup *> groupMap;
        for (auto &entry : fFitHistArray) {
            auto det = std::get<0>(entry.first);
            auto detGroup = groupMap[det];
            if (detGroup == nullptr) {
                detGroup = fitGroup -> CreateGroup(Form("F%d", det));
                groupMap[det] = detGroup;
            }
            detGroup -> AddHist(entry.second);
        }
    }
}

void SKSiCalibrationTask::WriteC0ParameterContainer() const
{
    LKParameterContainer par;
    par.SetName("SKSiCalibrationTask_C0");
    par.AddPar("SKSiCalibrationTask/stage", 1);
    par.AddPar("SKSiCalibrationTask/sourceEnergy", fSourceEnergy);
    par.AddPar("SKSiCalibrationTask/sourceEnergies", GetSourceEnergyListText());
    for (auto const &entry : fC0Parameters) {
        auto det = std::get<0>(entry.first);
        auto side = std::get<1>(entry.first);
        auto strip = std::get<2>(entry.first);
        par.AddPar(
            Form("SKSiCalibrationTask/C0/det%03d_side%d_strip%02d", det, side, strip),
            Form("%0.12g,%0.12g", entry.second[0], entry.second[1]),
            "intercept,slope"
        );
    }
    par.SaveAs(MakeFileName("c0.lkpar"));
    fRun -> RegisterObject("SKSiCalibrationTask_C0", par.CloneParameterContainer("SKSiCalibrationTask_C0"));
}

void SKSiCalibrationTask::WriteC1ParameterContainer() const
{
    LKParameterContainer par;
    par.SetName("SKSiCalibrationTask_C1");
    par.AddPar("SKSiCalibrationTask/stage", 2);
    par.AddPar("SKSiCalibrationTask/sourceEnergy", fSourceEnergy);
    par.AddPar("SKSiCalibrationTask/sourceEnergies", GetSourceEnergyListText());
    for (auto const &entry : fC1Parameters) {
        auto det = std::get<0>(entry.first);
        auto side = std::get<1>(entry.first);
        auto strip = std::get<2>(entry.first);
        auto lr = std::get<3>(entry.first);
        par.AddPar(
            Form("SKSiCalibrationTask/C1/det%03d_side%d_strip%02d_lr%d", det, side, strip, lr),
            Form("%0.12g,%0.12g", entry.second[0], entry.second[1]),
            "intercept,slope"
        );
    }
    par.SaveAs(MakeFileName("c1.lkpar"));
    fRun -> RegisterObject("SKSiCalibrationTask_C1", par.CloneParameterContainer("SKSiCalibrationTask_C1"));
}

void SKSiCalibrationTask::WriteC2ParameterContainer() const
{
    LKParameterContainer par;
    par.SetName("SKSiCalibrationTask_C2");
    par.AddPar("SKSiCalibrationTask/stage", 3);
    par.AddPar("SKSiCalibrationTask/sourceEnergy", fSourceEnergy);
    par.AddPar("SKSiCalibrationTask/sourceEnergies", GetSourceEnergyListText());
    par.AddPar("SKSiCalibrationTask/c2ReferenceEnergy", fC2ReferenceEnergy);
    for (auto const &entry : fC2Parameters) {
        auto det = std::get<0>(entry.first);
        auto side = std::get<1>(entry.first);
        auto strip = std::get<2>(entry.first);
        par.AddPar(
            Form("SKSiCalibrationTask/C2/det%03d_side%d_strip%02d", det, side, strip),
            Form("%0.12g,%0.12g,%0.12g", entry.second[0], entry.second[1], entry.second[2]),
            "b0,b1,b2"
        );
    }
    par.SaveAs(MakeFileName("c2.lkpar"));
    fRun -> RegisterObject("SKSiCalibrationTask_C2", par.CloneParameterContainer("SKSiCalibrationTask_C2"));
}

void SKSiCalibrationTask::WriteC3ParameterContainer() const
{
    LKParameterContainer par;
    par.SetName("SKSiCalibrationTask_C3");
    par.AddPar("SKSiCalibrationTask/stage", 5);
    par.AddPar("SKSiCalibrationTask/sourceEnergy", fSourceEnergy);
    par.AddPar("SKSiCalibrationTask/sourceEnergies", GetSourceEnergyListText());
    for (auto const &entry : fC3Parameters) {
        auto det = std::get<0>(entry.first);
        auto side = std::get<1>(entry.first);
        auto strip = std::get<2>(entry.first);
        par.AddPar(
            Form("SKSiCalibrationTask/C3/det%03d_side%d_strip%02d", det, side, strip),
            Form("%0.12g,%0.12g", entry.second[0], entry.second[1]),
            "intercept,slope"
        );
    }
    par.SaveAs(MakeFileName("c3.lkpar"));
    fRun -> RegisterObject("SKSiCalibrationTask_C3", par.CloneParameterContainer("SKSiCalibrationTask_C3"));
}

void SKSiCalibrationTask::WriteCombinedCalibrationFile()
{
    std::ofstream output(MakeFileName("energy_calibration.dat").Data());
    output << "# det side strip c0_itcpt c0_slope c1_itcptL c1_slopeL c1_itcptR c1_slopeR c2_b0 c2_b1 c2_b2 c3_itcpt c3_slope\n";
    output << std::setprecision(12);

    std::map<StripKey, bool> keys;
    for (auto const &entry : fC0Parameters) keys[entry.first] = true;
    for (auto const &entry : fC2Parameters) keys[entry.first] = true;
    for (auto const &entry : fC3Parameters) keys[entry.first] = true;
    for (auto const &entry : fEnergyHistMap) keys[entry.first] = true;
    for (auto const &entry : fEnergySumHistMap) keys[entry.first] = true;

    for (auto const &entry : keys) {
        auto det = std::get<0>(entry.first);
        auto side = std::get<1>(entry.first);
        auto strip = std::get<2>(entry.first);
        auto c0 = fC0Parameters.count(entry.first) ? fC0Parameters[entry.first] : std::array<double,2>{0, 1};
        auto c1L = fC1Parameters.count(LRKey(det, side, strip, 0)) ? fC1Parameters[LRKey(det, side, strip, 0)] : std::array<double,2>{0, 1};
        auto c1R = fC1Parameters.count(LRKey(det, side, strip, 1)) ? fC1Parameters[LRKey(det, side, strip, 1)] : std::array<double,2>{0, 1};
        auto c2 = fC2Parameters.count(entry.first) ? fC2Parameters[entry.first] : std::array<double,3>{fC2ReferenceEnergy, 0, 0};
        auto c3 = fC3Parameters.count(entry.first) ? fC3Parameters[entry.first] : std::array<double,2>{0, 1};
        output << det << " " << side << " " << strip << " "
               << c0[0] << " " << c0[1] << " "
               << c1L[0] << " " << c1L[1] << " "
               << c1R[0] << " " << c1R[1] << " "
               << c2[0] << " " << c2[1] << " " << c2[2] << " "
               << c3[0] << " " << c3[1] << "\n";
    }
    WriteCombinedParameterContainer();
}

void SKSiCalibrationTask::WriteCombinedParameterContainer() const
{
    LKParameterContainer par;
    par.SetName("SKSiCalibrationTask_EnergyCalibration");
    par.AddPar("SKSiCalibrationTask/stage", 6);
    par.AddPar("SKSiCalibrationTask/sourceEnergy", fSourceEnergy);
    par.AddPar("SKSiCalibrationTask/sourceEnergies", GetSourceEnergyListText());
    par.AddPar("SKSiCalibrationTask/c2ReferenceEnergy", fC2ReferenceEnergy);

    std::map<StripKey, bool> keys;
    for (auto const &entry : fC0Parameters) keys[entry.first] = true;
    for (auto const &entry : fC2Parameters) keys[entry.first] = true;
    for (auto const &entry : fC3Parameters) keys[entry.first] = true;
    for (auto const &entry : fEnergyHistMap) keys[entry.first] = true;
    for (auto const &entry : fEnergySumHistMap) keys[entry.first] = true;

    for (auto const &entry : keys) {
        auto det = std::get<0>(entry.first);
        auto side = std::get<1>(entry.first);
        auto strip = std::get<2>(entry.first);
        auto c0 = fC0Parameters.count(entry.first) ? fC0Parameters.at(entry.first) : std::array<double,2>{0, 1};
        auto c1L = fC1Parameters.count(LRKey(det, side, strip, 0)) ? fC1Parameters.at(LRKey(det, side, strip, 0)) : std::array<double,2>{0, 1};
        auto c1R = fC1Parameters.count(LRKey(det, side, strip, 1)) ? fC1Parameters.at(LRKey(det, side, strip, 1)) : std::array<double,2>{0, 1};
        auto c2 = fC2Parameters.count(entry.first) ? fC2Parameters.at(entry.first) : std::array<double,3>{fC2ReferenceEnergy, 0, 0};
        auto c3 = fC3Parameters.count(entry.first) ? fC3Parameters.at(entry.first) : std::array<double,2>{0, 1};
        par.AddPar(
            Form("SKSiCalibrationTask/EnergyCalibration/det%03d_side%d_strip%02d", det, side, strip),
            Form("%0.12g,%0.12g,%0.12g,%0.12g,%0.12g,%0.12g,%0.12g,%0.12g,%0.12g,%0.12g,%0.12g",
                c0[0], c0[1],
                c1L[0], c1L[1],
                c1R[0], c1R[1],
                c2[0], c2[1], c2[2],
                c3[0], c3[1]),
            "c0_itcpt,c0_slope,c1_itcptL,c1_slopeL,c1_itcptR,c1_slopeR,c2_b0,c2_b1,c2_b2,c3_itcpt,c3_slope"
        );
    }

    par.SaveAs(MakeFileName("energy_calibration.lkpar"));
    fRun -> RegisterObject("SKSiCalibrationTask_EnergyCalibration", par.CloneParameterContainer("SKSiCalibrationTask_EnergyCalibration"));
}

void SKSiCalibrationTask::InitializeExpectedChannels()
{
    if (fSiArray == nullptr) {
        lk_warning << "LKSiliconArray is not added. SKSiCalibrationTask will summarize only channels seen in SiChannel." << endl;
        return;
    }

    for (auto iDetector = 0; iDetector < fSiArray->GetNumSiDetectors(); ++iDetector) {
        auto detector = fSiArray -> GetSiDetector(iDetector);
        if (detector == nullptr)
            continue;

        auto numRegisteredChannels = detector -> GetNumRegisteredChannels();
        for (auto iChannel = 0; iChannel < numRegisteredChannels; ++iChannel) {
            auto channel = detector -> GetRegisteredChannel(iChannel);
            if (channel == nullptr)
                continue;

            auto det = channel -> GetDetNum();
            auto side = channel -> GetSide();
            auto strip = channel -> GetStrip();
            StripKey key(det, side, strip);
            if (channel -> IsPairedChannel()) {
                fExpectedPairedMap[key] = true;
                if (fStage == 0 || fStage >= 3) {
                    GetEnergySumHist(det, side, strip);
                    GetPositionEnergyHist(det, side, strip);
                }
                if (fStage == 0 || fStage >= 2)
                    GetLeftRightHist(det, side, strip);
            }
            else {
                fExpectedStandaloneMap[key] = true;
                if (fStage == 0 || fStage == 1 || fStage >= 6)
                    GetEnergyHist(det, side, strip);
            }
        }
    }

    lk_info << "Expected standalone strips = " << fExpectedStandaloneMap.size()
            << ", expected paired strips = " << fExpectedPairedMap.size() << endl;
}

void SKSiCalibrationTask::WriteChannelSummary() const
{
    std::ofstream output(MakeFileName(Form("stage%d.summary.dat", fStage)).Data());
    output << "# det side strip channel_type expected entries status\n";

    LKParameterContainer par;
    par.SetName(Form("SKSiCalibrationTask_Stage%dSummary", fStage));
    par.AddPar("SKSiCalibrationTask/stage", fStage);
    par.AddPar("SKSiCalibrationTask/sourceEnergies", GetSourceEnergyListText());
    par.AddPar("SKSiCalibrationTask/expectedStandalone", (int) fExpectedStandaloneMap.size());
    par.AddPar("SKSiCalibrationTask/expectedPaired", (int) fExpectedPairedMap.size());

    int foundStandalone = 0;
    int missingStandalone = 0;
    int foundPaired = 0;
    int missingPaired = 0;

    auto writeSummaryLine = [&](StripKey key, TString type, bool expected, double entries) {
        auto det = std::get<0>(key);
        auto side = std::get<1>(key);
        auto strip = std::get<2>(key);
        TString status = entries > 0 ? "found" : "missing";
        output << det << " " << side << " " << strip << " "
               << type << " " << (expected ? 1 : 0) << " "
               << entries << " " << status << "\n";
        par.AddPar(
            Form("SKSiCalibrationTask/Summary/stage%d/det%03d_side%d_strip%02d_%s", fStage, det, side, strip, type.Data()),
            Form("%d,%0.12g,%s", expected ? 1 : 0, entries, status.Data()),
            "expected,entries,status"
        );
    };

    for (auto const &entry : fExpectedStandaloneMap) {
        auto it = fEnergyHistMap.find(entry.first);
        auto entries = it == fEnergyHistMap.end() ? 0 : it->second->GetEntries();
        if (entries > 0) ++foundStandalone;
        else ++missingStandalone;
        writeSummaryLine(entry.first, "single", true, entries);
    }

    for (auto const &entry : fExpectedPairedMap) {
        auto it = fEnergySumHistMap.find(entry.first);
        auto entries = it == fEnergySumHistMap.end() ? 0 : it->second->GetEntries();
        if (entries > 0) ++foundPaired;
        else ++missingPaired;
        writeSummaryLine(entry.first, "paired", true, entries);
    }

    if (fExpectedStandaloneMap.empty() && fExpectedPairedMap.empty()) {
        for (auto const &entry : fEnergyHistMap)
            writeSummaryLine(entry.first, "single", false, entry.second->GetEntries());
        for (auto const &entry : fEnergySumHistMap)
            writeSummaryLine(entry.first, "paired", false, entry.second->GetEntries());
    }

    par.AddPar("SKSiCalibrationTask/foundStandalone", foundStandalone);
    par.AddPar("SKSiCalibrationTask/missingStandalone", missingStandalone);
    par.AddPar("SKSiCalibrationTask/foundPaired", foundPaired);
    par.AddPar("SKSiCalibrationTask/missingPaired", missingPaired);
    par.SaveAs(MakeFileName(Form("stage%d.summary.lkpar", fStage)));
    fRun -> RegisterObject(Form("SKSiCalibrationTask_Stage%dSummary", fStage), par.CloneParameterContainer(Form("SKSiCalibrationTask_Stage%dSummary", fStage)));
}

double SKSiCalibrationTask::GetCalibratedEnergyAxisMax() const
{
    auto maxEnergy = fSourceEnergy;
    for (auto energy : fSourceEnergies)
        if (energy > maxEnergy)
            maxEnergy = energy;

    if (maxEnergy <= 0)
        return 20.;
    return std::max(20., 1.6 * maxEnergy);
}

TString SKSiCalibrationTask::GetSourceEnergyListText() const
{
    TString text;
    for (auto iEnergy = 0u; iEnergy < fSourceEnergies.size(); ++iEnergy) {
        if (iEnergy > 0)
            text += ",";
        text += Form("%0.12g", fSourceEnergies[iEnergy]);
    }
    return text;
}

TString SKSiCalibrationTask::MakeFileName(TString suffix) const
{
    return Form("%s/%s.%s", fOutputPath.Data(), fOutputPrefix.Data(), suffix.Data());
}

TString SKSiCalibrationTask::MakeHistName(TString quantity, int det, int side, int strip) const
{
    return Form("stage%d_%s_det%03d_side%d_strip%02d", fStage, quantity.Data(), det, side, strip);
}
