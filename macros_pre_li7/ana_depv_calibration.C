#include "LKCompiled.h"
#include "LKDrawingGroup.h"
#include "LKBinning.h"
#include "LKParameterContainer.h"
#include "LKRun.h"
#include "LKSiliconArray.h"
#include "LKSiChannel.h"
#include "TCanvas.h"
#include "TClonesArray.h"
#include "TF1.h"
#include "TH2D.h"
#include "TH1D.h"
#include "TLine.h"
#include "TPaveText.h"
#include "TString.h"
#include "TSystem.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <map>
#include <tuple>
#include <vector>

struct FitResult {
    double amp = 0;
    double energy = 0;
    double fwhm = 0;
};

void Fill(LKRun *run, TClonesArray *siChannelArray, LKSiliconArray *siArrayPlane, std::map<std::pair<int, int>, TH2D *> &histMap, const std::map<std::tuple<int, int, int>, double> *correctionMap = nullptr);
void AddCorrectionFactors(TH2D *hist, int detID, int side, std::map<std::tuple<int, int, int>, double> &correctionMap, std::map<std::pair<int, int>, std::vector<TH1D *>> *fitDrawingMap = nullptr, double targetEnergy = 3.182, double fitStdDevRange = 5.0, std::map<std::tuple<int, int, int>, FitResult> *fitResultMap = nullptr);
TH2D *CreateCorrectedHistogram(TH2D *hist, double targetEnergy = 3.182);
void WriteCorrectedFitData(TString detName, TString voltage, const std::map<std::tuple<int, int, int>, FitResult> &correctedFitResultMap, const std::map<std::tuple<int, int, int>, FitResult> &rawFitResultMap);
TString SafeFileToken(TString token);

LKBinning bnnADC(500,0,5000);
LKBinning bnnEnergy(500,0,7);
double targetEnergy = 1;//sourceEnergy.Atof();

void ana_depv_calibration(
    int runID = -1,
    TString parameterFileName = "ana_depv_calibration.conf",
    bool saveImages = false,
    bool draw_fittings = false,
    bool draw_corrected_fit = false,
    bool drawTop = false
    )
{
    auto par = LKParameterContainer(parameterFileName);
    if (runID<0)
        runID = par.GetParInt("LKRun/RunID");
    TString inputFileName = Form("/Users/jungwoo/data/stark/reco/pre_li7_%04d.reco.root",runID);
    double fitStdDevRange = 5.0;
    TString detName, mappingDir, voltage, sourceName, sourceEnergy;
    par.FindAndRetrieveColumnValue("list.txt", 0, runID, 1, detName);
    par.FindAndRetrieveColumnValue("list.txt", 0, runID, 2, mappingDir);
    par.FindAndRetrieveColumnValue("list.txt", 0, runID, 3, voltage);
    par.FindAndRetrieveColumnValue("list.txt", 0, runID, 4, sourceName);
    par.FindAndRetrieveColumnValue("list.txt", 0, runID, 5, sourceEnergy);
    e_info << inputFileName << endl;
    e_info << detName <<endl;
    e_info << mappingDir <<endl;
    e_info << voltage <<endl;
    e_info << sourceName <<endl;
    e_info << sourceEnergy <<endl;
    par.UpdatePar(bnnADC,"bnnADC");
    par.UpdatePar(bnnEnergy,"bnnEnergy");
    par.UpdatePar(fitStdDevRange,"fitStdDevRange");
    targetEnergy = sourceEnergy.Atof();

    auto run = new LKRun();
    run->SetInputFile(inputFileName, "event");
    if (!run->Init()) {
        Error("draw_detector_channels", "Failed to initialize LKRun with input %s", inputFileName.Data());
        return;
    }

    auto siPar = new LKParameterContainer();
    siPar->AddLine(Form("si_array/MappingPath %s # mapping directory for LKSiliconArray", mappingDir.Data()));
    auto siArrayPlane = new LKSiliconArray();
    siArrayPlane->SetPar(siPar);
    if (!siArrayPlane->Init()) {
        Error("draw_detector_channels", "Failed to initialize LKSiliconArray with mapping %s", mappingDir.Data());
        return;
    }

    auto siChannelArray = run->GetBranchA("SiChannel", "LKSiChannel");
    if (siChannelArray == nullptr) {
        Error("draw_detector_channels", "Branch %s does not exist or is not an LKSiChannel TClonesArray.","SiChannel");
        return;
    }

    std::map<std::pair<int, int>, int> maxStrip;
    std::map<int, TString> detTitle;
    for (int iChannel = 0; iChannel < siArrayPlane->GetNumSiDetectors() * 256; ++iChannel) {
        auto dummy = siArrayPlane->GetSiChannel(iChannel);
        if (dummy == nullptr)
            continue;
        const int detIndex = dummy->GetDetID();
        const int detID = dummy->GetDetNum();
        const int side = dummy->GetSide();
        const int strip = dummy->GetStrip();
        maxStrip[{detID, side}] = std::max(maxStrip[{detID, side}], strip);
        auto detector = siArrayPlane->GetSiDetector(detIndex);
        detTitle[detID] = detector ? detector->GetName() : Form("det%d", detID);
    }

    // histogram
    std::map<std::pair<int, int>, TH2D *> histMap;
    for (auto const &entry : maxStrip) {
        const int detID = entry.first.first;
        const int side = entry.first.second;
        const int nStrips = std::max(1, entry.second + 1);
        TString sideName = side == 0 ? "junction" : "ohmic";
        TString histName = Form("hist_det%d_side%d", detID, side);
        TString histTitle = Form("det=%d %s %s;strip;fEnergy", detID, detTitle[detID].Data(), sideName.Data());
        LKBinning bnnStrips(nStrips, -0.5, nStrips - 0.5); 
        histMap[entry.first] = (bnnStrips*bnnADC).NewH2(histName, histTitle);
    }

    Fill(run, siChannelArray, siArrayPlane, histMap);

    std::map<std::tuple<int, int, int>, double> correctionMap;
    std::map<std::pair<int, int>, TH2D *> correctedHistMap;
    std::map<std::pair<int, int>, std::vector<TH1D *>> fitDrawingMap;
    std::map<std::pair<int, int>, std::vector<TH1D *>> correctedFitDrawingMap;
    std::map<std::tuple<int, int, int>, FitResult> rawFitResultMap;
    std::map<std::tuple<int, int, int>, FitResult> correctedFitResultMap;
    for (auto const &entry : histMap) {
        const int detID = entry.first.first;
        const int side = entry.first.second;
        AddCorrectionFactors(
            entry.second,
            detID,
            side,
            correctionMap,
            draw_fittings ? &fitDrawingMap : nullptr,
            targetEnergy,
            fitStdDevRange,
            &rawFitResultMap
            );
        correctedHistMap[entry.first] = CreateCorrectedHistogram(entry.second, targetEnergy);
    }
    Fill(run, siChannelArray, siArrayPlane, correctedHistMap, &correctionMap);
    if (draw_corrected_fit) {
        std::map<std::tuple<int, int, int>, double> correctedFitCorrectionMap;
        for (auto const &entry : correctedHistMap) {
            const int detID = entry.first.first;
            const int side = entry.first.second;
            AddCorrectionFactors(
                entry.second,
                detID,
                side,
                correctedFitCorrectionMap,
                &correctedFitDrawingMap,
                targetEnergy,
                fitStdDevRange,
                &correctedFitResultMap
                );
        }
    }
    else {
        std::map<std::tuple<int, int, int>, double> correctedFitCorrectionMap;
        for (auto const &entry : correctedHistMap) {
            const int detID = entry.first.first;
            const int side = entry.first.second;
            AddCorrectionFactors(
                entry.second,
                detID,
                side,
                correctedFitCorrectionMap,
                nullptr,
                targetEnergy,
                fitStdDevRange,
                &correctedFitResultMap
                );
        }
    }

    WriteCorrectedFitData(detName, voltage, correctedFitResultMap, rawFitResultMap);

    auto top = new LKDrawingGroup();
    for (int iDetector = 0; iDetector < siArrayPlane->GetNumSiDetectors(); ++iDetector)
    {
        auto detector = siArrayPlane->GetSiDetector(iDetector);
        const int detNumber = detector != nullptr ? detector->GetDetNum() : iDetector;
        auto group = top -> CreateGroup(Form("det%d",detNumber));
        auto histJ = histMap[{detNumber, 0}];
        auto histO = histMap[{detNumber, 1}];
        auto histJCorrected = correctedHistMap[{detNumber, 0}];
        auto histOCorrected = correctedHistMap[{detNumber, 1}];
        auto histGroup = group;
        if (draw_fittings)
            histGroup = group -> CreateGroup(Form("det%d_histograms", detNumber));
        if (histJ != nullptr)
            histGroup -> Add(histJ);
        if (histO != nullptr)
            histGroup -> Add(histO);
        if (histJCorrected != nullptr)
            histGroup -> Add(histJCorrected);
        if (histOCorrected != nullptr)
            histGroup -> Add(histOCorrected);

        int increase = 0;
        LKDrawingGroup *fitGroup = group -> CreateGroup(Form("det%d_fittings_%d", detNumber, ++increase));
        if (draw_fittings || draw_corrected_fit) {
            if (draw_fittings) {
                for (auto side : {0, 1}) {
                    for (auto projection : fitDrawingMap[{detNumber, side}]) {
                        if (fitGroup -> GetNumDrawings()>=12)
                            fitGroup = group -> CreateGroup(Form("det%d_fittings_%d", detNumber, ++increase));
                        fitGroup -> Add(projection);
                    }
                }
            }
            fitGroup = group -> CreateGroup(Form("det%d_fit_corrected", detNumber));
            if (draw_corrected_fit) {
                for (auto side : {0, 1}) {
                    for (auto projection : correctedFitDrawingMap[{detNumber, side}]) {
                        if (fitGroup -> GetNumDrawings()>=12)
                            fitGroup = group -> CreateGroup(Form("det%d_fittings_%d", detNumber, ++increase));
                        fitGroup -> Add(projection);
                    }
                }
            }
        }
    }

    if (drawTop)
        top -> Draw("viewer");

    Info("draw_detector_channels", "Drew detector channels using LKSiliconArray mapping from %s (%lld entries).",
         mappingDir.Data(), run->GetEntries());
}

void Fill(
    LKRun *run,
    TClonesArray *siChannelArray,
    LKSiliconArray *siArrayPlane,
    std::map<std::pair<int, int>, TH2D *> &histMap,
    const std::map<std::tuple<int, int, int>, double> *correctionMap)
{
    const auto nEntries = run->GetEntries();
    for (Long64_t iEntry = 0; iEntry < nEntries; ++iEntry)
    {
        run->GetEntry(iEntry);
        if (siChannelArray == nullptr)
            continue;

        const int nChannels = siChannelArray->GetEntriesFast();
        for (int iChannel = 0; iChannel < nChannels; ++iChannel) {
            auto channel = (LKSiChannel *) siChannelArray->At(iChannel);
            if (channel == nullptr)
                continue;

            auto dummy = siArrayPlane->GetSiChannel(channel->GetCobo(), channel->GetAsad(), channel->GetAget(), channel->GetChan());
            if (dummy == nullptr)
                continue;

            const int detID = dummy->GetDetNum();
            const int side = dummy->GetSide();
            const int strip = dummy->GetStrip();
            auto hist = histMap[{detID, side}];
            if (hist == nullptr)
                continue;

            double energy = channel->GetEnergy();
            if (correctionMap != nullptr) {
                auto correction = correctionMap->find(std::make_tuple(detID, side, strip));
                if (correction != correctionMap->end())
                    energy *= correction->second;
            }
            hist->Fill(strip, energy);
        }
    }
}

void AddCorrectionFactors(
    TH2D *hist,
    int detID,
    int side,
    std::map<std::tuple<int, int, int>, double> &correctionMap,
    std::map<std::pair<int, int>, std::vector<TH1D *>> *fitDrawingMap,
    double targetEnergy,
    double fitStdDevRange,
    std::map<std::tuple<int, int, int>, FitResult> *fitResultMap
    )
{
    if (hist == nullptr)
        return;

    const int nXBins = hist->GetNbinsX();
    for (int xBin = 1; xBin <= nXBins; ++xBin) {
        auto projection = hist->ProjectionY(Form("%s_py_%d", hist->GetName(), xBin), xBin, xBin);
        if (projection == nullptr || projection->GetEntries() < 10) {
            delete projection;
            continue;
        }

        const int strip = lround(hist->GetXaxis()->GetBinCenter(xBin));
        const int maxBin = projection->GetMaximumBin();
        const double peak = projection->GetBinCenter(maxBin);
        const double amp = projection->GetBinContent(maxBin);
        const double rms = 0.01*peak;//projection->GetRMS();
        const double fitHalfWidth = std::max(0.1, fitStdDevRange * rms);
        const double fitMin = std::max(projection->GetXaxis()->GetXmin(), peak - fitHalfWidth);
        const double fitMax = std::min(projection->GetXaxis()->GetXmax(), peak + fitHalfWidth);
        projection->SetName(Form("%s_fit_det%d_side%d_strip%d", hist->GetName(), detID, side, strip));
        projection->SetTitle(Form("det=%d side=%d strip=%d;fEnergy;counts", detID, side, strip));
        TF1 fit(Form("gaus_det%d_side%d_strip%d", detID, side, strip), "gaus", fitMin, fitMax);
        fit.SetParameters(amp, peak, rms);
        const char *fitOption = fitDrawingMap == nullptr ? "QNR" : "QR";
        if (projection->Fit(&fit, fitOption) == 0) {
            const double mean = fit.GetParameter(1);
            const double sigma = std::abs(fit.GetParameter(2));
            if (mean > 0) {
                correctionMap[std::make_tuple(detID, side, strip)] = targetEnergy / mean;
                if (fitResultMap != nullptr) {
                    (*fitResultMap)[std::make_tuple(detID, side, strip)] = {
                        fit.GetParameter(0),
                        mean,
                        2.354820045 * sigma
                    };
                }
            }
        }
        //projection -> GetXaxis() -> SetRangeUser(fit.GetParameter(1)-5*fit.GetParameter(2),fit.GetParameter(1)+5*fit.GetParameter(2));
        projection -> GetXaxis() -> SetRangeUser(fitMin,fitMax);
        if (fitDrawingMap != nullptr) {
            auto initialFit = new TF1(Form("initial_gaus_det%d_side%d_strip%d", detID, side, strip), "gaus", fitMin, fitMax);
            initialFit->SetParameters(amp, peak, rms);
            initialFit->SetLineColor(kBlue + 1);
            initialFit->SetLineStyle(2);
            initialFit->SetLineWidth(2);

            const double yMaxLine = std::max(projection->GetMaximum() * 1.08, 1.0);
            auto lineMin = new TLine(fitMin, 0, fitMin, yMaxLine);
            auto lineMax = new TLine(fitMax, 0, fitMax, yMaxLine);
            lineMin->SetLineColor(kGreen + 2);
            lineMax->SetLineColor(kGreen + 2);
            lineMin->SetLineStyle(7);
            lineMax->SetLineStyle(7);

            auto text = new TPaveText(0.12, 0.62, 0.48, 0.88, "NDC");
            text->SetFillColor(0);
            text->SetFillStyle(0);
            text->SetBorderSize(1);
            text->SetTextAlign(12);
            text->AddText(Form("init A = %.4g", amp));
            text->AddText(Form("init #mu = %.4g", peak));
            text->AddText(Form("init #sigma = %.4g", rms));
            text->AddText(Form("fit range = %.4g - %.4g", fitMin, fitMax));

            projection->GetListOfFunctions()->Add(initialFit);
            projection->GetListOfFunctions()->Add(lineMin);
            projection->GetListOfFunctions()->Add(lineMax);
            projection->GetListOfFunctions()->Add(text);
            (*fitDrawingMap)[{detID, side}].push_back(projection);
        }
        else
            delete projection;
    }
}

TH2D *CreateCorrectedHistogram(TH2D *hist, double targetEnergy)
{
    if (hist == nullptr)
        return nullptr;

    auto bnnStrips = LKBinning(hist);//.GetBinningX();
    auto corrected = (bnnStrips*bnnEnergy).NewH2(Form("%s_corrected", hist->GetName()), Form("%s corrected to %.3f MeV", hist->GetTitle(), targetEnergy));
    return corrected;
}

TString SafeFileToken(TString token)
{
    token = token.Strip(TString::kBoth);
    if (token.IsNull())
        token = "unknown";
    token.ReplaceAll(" ", "_");
    token.ReplaceAll("/", "_");
    token.ReplaceAll("\\", "_");
    return token;
}

void WriteCorrectedFitData(
    TString detName,
    TString voltage,
    const std::map<std::tuple<int, int, int>, FitResult> &correctedFitResultMap,
    const std::map<std::tuple<int, int, int>, FitResult> &rawFitResultMap)
{
    gSystem->mkdir("data_depv", true);
    TString fileName = Form("data_depv/%s_%s_er.dat", SafeFileToken(detName).Data(), SafeFileToken(voltage).Data());
    std::ofstream output(fileName.Data());
    if (!output.is_open()) {
        Error("ana_depv_calibration", "Cannot open output file %s", fileName.Data());
        return;
    }

    output << "# det_number side channel corrected_fit_amp corrected_fit_energy corrected_fit_fwhm raw_fit_mean\n";
    output << std::setprecision(10);
    for (auto const &entry : correctedFitResultMap) {
        const int detNumber = std::get<0>(entry.first);
        const int side = std::get<1>(entry.first);
        const int channel = std::get<2>(entry.first);
        const auto &fit = entry.second;
        const auto rawFitIt = rawFitResultMap.find(entry.first);
        const double rawFitMean = rawFitIt != rawFitResultMap.end() ? rawFitIt->second.energy : 0;
        output << detNumber << " "
               << side << " "
               << channel << " "
               << fit.amp << " "
               << fit.energy << " "
               << fit.fwhm << " "
               << rawFitMean << "\n";
    }

    Info("ana_depv_calibration", "Wrote corrected fit data to %s (%zu channels).", fileName.Data(), correctedFitResultMap.size());
}
