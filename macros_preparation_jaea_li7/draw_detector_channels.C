#include "TCanvas.h"
#include "TClonesArray.h"
#include "TF1.h"
#include "TH2D.h"
#include "TH1D.h"
#include "TLine.h"
#include "TPaveText.h"
#include "TString.h"

#include <algorithm>
#include <cmath>
#include <map>
#include <tuple>
#include <vector>

void Fill(
    LKRun *run,
    TClonesArray *siChannelArray,
    LKSiliconArray *siArrayPlane,
    std::map<std::pair<int, int>, TH2D *> &histMap,
    const std::map<std::tuple<int, int, int>, double> *correctionMap = nullptr)
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

            const int detID = dummy->GetDetID();
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
    std::map<std::pair<int, int>, std::vector<TH1D *>> *fitDrawingMap = nullptr,
    double targetEnergy = 3.182,
    double fitStdDevRange = 5.0,
    int fitProjectionRebin = 1)
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
        if (fitProjectionRebin > 1)
            projection->Rebin(fitProjectionRebin);

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
            if (mean > 0) {
                correctionMap[std::make_tuple(detID, side, strip)] = targetEnergy / mean;
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

TH2D *CreateCorrectedHistogram(TH2D *hist, double targetEnergy = 3.182, int correctedEnergyBins = 1000)
{
    if (hist == nullptr)
        return nullptr;

    const int nXBins = hist->GetNbinsX();
    const double xMin = hist->GetXaxis()->GetXmin();
    const double xMax = hist->GetXaxis()->GetXmax();
    const double yMax = targetEnergy*2.0;//std::max(5.0, targetEnergy * 2.5);
    auto corrected = new TH2D(
        Form("%s_corrected", hist->GetName()),
        Form("%s corrected to %.3f MeV", hist->GetTitle(), targetEnergy),
        nXBins,
        xMin,
        xMax,
        correctedEnergyBins,
        0,
        yMax);
    return corrected;
}

void draw_detector_channels(
    TString inputFileName = "/Users/jungwoo/data/stark/reco/pre_li7_0111.reco.root",
    TString mappingDir = "/Users/jungwoo/Research/lilak/stark/macros_preparation_jaea_li7/mapping_QB_DepV_Measure",
    TString treeName = "event",
    TString branchName = "SiChannel",
    bool saveImages = false,
    bool draw_fittings = true,
    bool draw_corrected_fit = true,
    int energyBins = 500,
    int correctedEnergyBins = 500,
    double fitStdDevRange = 5.0,
    int fitProjectionRebin = 1,
    bool drawTop = true)
{
    auto run = new LKRun();
    run->SetInputFile(inputFileName, treeName);

    auto par = new LKParameterContainer();
    par->AddPar("stark/ForceMapping", mappingDir, "mapping directory for LKSiliconArray");
    run->AddPar(par);

    auto siArrayPlane = new LKSiliconArray();
    run->Add(siArrayPlane);
    if (!run->Init()) {
        Error("draw_detector_channels", "Failed to initialize LKRun with input %s and mapping %s",
              inputFileName.Data(), mappingDir.Data());
        return;
    }

    auto siChannelArray = run->GetBranchA(branchName, "LKSiChannel");
    if (siChannelArray == nullptr) {
        Error("draw_detector_channels", "Branch %s does not exist or is not an LKSiChannel TClonesArray.",
              branchName.Data());
        return;
    }

    std::map<std::pair<int, int>, int> maxStrip;
    std::map<int, TString> detTitle;
    for (int iChannel = 0; iChannel < siArrayPlane->GetNumSiDetectors() * 256; ++iChannel) {
        auto dummy = siArrayPlane->GetSiChannel(iChannel);
        if (dummy == nullptr)
            continue;
        const int detID = dummy->GetDetID();
        const int side = dummy->GetSide();
        const int strip = dummy->GetStrip();
        maxStrip[{detID, side}] = std::max(maxStrip[{detID, side}], strip);
        auto detector = siArrayPlane->GetSiDetector(detID);
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
        histMap[entry.first] = new TH2D(histName, histTitle, nStrips, -0.5, nStrips - 0.5, energyBins, 0, 2000);
    }

    Fill(run, siChannelArray, siArrayPlane, histMap);

    std::map<std::tuple<int, int, int>, double> correctionMap;
    std::map<std::pair<int, int>, TH2D *> correctedHistMap;
    std::map<std::pair<int, int>, std::vector<TH1D *>> fitDrawingMap;
    std::map<std::pair<int, int>, std::vector<TH1D *>> correctedFitDrawingMap;
    for (auto const &entry : histMap) {
        const int detID = entry.first.first;
        const int side = entry.first.second;
        AddCorrectionFactors(
            entry.second,
            detID,
            side,
            correctionMap,
            draw_fittings ? &fitDrawingMap : nullptr,
            3.182,
            fitStdDevRange,
            fitProjectionRebin);
        correctedHistMap[entry.first] = CreateCorrectedHistogram(entry.second, 3.182, correctedEnergyBins);
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
                3.182,
                fitStdDevRange,
                fitProjectionRebin);
        }
    }

    auto top = new LKDrawingGroup();
    for (int iDetector = 0; iDetector < siArrayPlane->GetNumSiDetectors(); ++iDetector)
    {
        auto group = top -> CreateGroup(Form("det%d",iDetector));
        auto histJ = histMap[{iDetector, 0}];
        auto histO = histMap[{iDetector, 1}];
        auto histJCorrected = correctedHistMap[{iDetector, 0}];
        auto histOCorrected = correctedHistMap[{iDetector, 1}];
        auto histGroup = group;
        if (draw_fittings)
            histGroup = group -> CreateGroup(Form("det%d_histograms", iDetector));
        if (histJ != nullptr)
            histGroup -> Add(histJ);
        if (histO != nullptr)
            histGroup -> Add(histO);
        if (histJCorrected != nullptr)
            histGroup -> Add(histJCorrected);
        if (histOCorrected != nullptr)
            histGroup -> Add(histOCorrected);

        if (draw_fittings || draw_corrected_fit) {
            auto fitGroup = group -> CreateGroup(Form("det%d_fittings", iDetector));
            if (draw_fittings) {
                //auto rawFitGroup = fitGroup -> CreateGroup(Form("det%d_raw_fittings", iDetector));
                for (auto side : {0, 1}) {
                    for (auto projection : fitDrawingMap[{iDetector, side}])
                        fitGroup -> Add(projection);
                }
            }
            fitGroup = group -> CreateGroup(Form("det%d_fit_corrected", iDetector));
            if (draw_corrected_fit) {
                //auto correctedFitGroup = fitGroup -> CreateGroup(Form("det%d_corrected_fittings", iDetector));
                for (auto side : {0, 1}) {
                    for (auto projection : correctedFitDrawingMap[{iDetector, side}])
                        fitGroup -> Add(projection);
                }
            }
        }
    }

    if (drawTop)
        top -> Draw("viewer");

    Info("draw_detector_channels", "Drew detector channels using LKSiliconArray mapping from %s (%lld entries).",
         mappingDir.Data(), run->GetEntries());
}
