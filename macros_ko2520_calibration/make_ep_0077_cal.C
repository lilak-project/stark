#include "TCanvas.h"
#include "TFile.h"
#include "TH2D.h"
#include "TKey.h"
#include "TTree.h"

#include <array>
#include <fstream>
#include <map>
#include <tuple>
#include <vector>

void make_ep_0077_cal(const char *inputPath="data_reco2", const char *outputName="ep_0077_cal.root")
{
    TFile reference(Form("%s/ep_0077.root", inputPath), "READ");
    TFile calibratedData("ko2520_0077.sihit_cal.root", "READ");
    TFile recoData(Form("%s/ko2520_0077.reco.root", inputPath), "READ");
    auto tree = dynamic_cast<TTree *>(calibratedData.Get("event"));
    auto recoTree = dynamic_cast<TTree *>(recoData.Get("event"));
    if (reference.IsZombie() || calibratedData.IsZombie() || recoData.IsZombie()
        || tree == nullptr || recoTree == nullptr) {
        Error("make_ep_0077_cal", "Input file is missing");
        return;
    }

    // QQQ5 junction strips are standalone channels, just like X6 ohmic
    // channels. Apply the same C0 calibration: E = intercept + slope * ADC.
    using ChannelKey = std::tuple<int, int, int>;
    std::map<ChannelKey, std::array<double, 2>> c0;
    std::ifstream calibration("data_calibration/run_0077/ko2520_0077.energy_calibration.dat");
    std::string header;
    std::getline(calibration, header);
    int calDet, calSide, calStrip;
    double p[11];
    while (calibration >> calDet >> calSide >> calStrip
        >> p[0] >> p[1] >> p[2] >> p[3] >> p[4] >> p[5]
        >> p[6] >> p[7] >> p[8] >> p[9] >> p[10])
        c0[{calDet, calSide, calStrip}] = {p[0], p[1]};

    struct Plot {
        TString canvasName;
        TString histName;
        TString title;
        int det = -1;
        int detNumber = -1;
        bool qqq5 = false;
        TH2D *hist = nullptr;
    };
    std::vector<Plot> plots;

    TIter nextKey(reference.GetListOfKeys());
    while (auto key = dynamic_cast<TKey *>(nextKey())) {
        auto canvas = dynamic_cast<TCanvas *>(key->ReadObj());
        if (canvas == nullptr)
            continue;
        TIter nextPad(canvas->GetListOfPrimitives());
        while (auto object = nextPad()) {
            auto pad = dynamic_cast<TPad *>(object);
            if (pad == nullptr)
                continue;
            TIter nextObject(pad->GetListOfPrimitives());
            while (auto primitive = nextObject()) {
                auto source = dynamic_cast<TH2 *>(primitive);
                if (source == nullptr || !TString(source->GetName()).BeginsWith("hist_e_vs_position_det_"))
                    continue;
                int det = -1;
                sscanf(source->GetName(), "hist_e_vs_position_det_%d", &det);
                Plot plot;
                plot.canvasName = canvas->GetName();
                plot.canvasName.ReplaceAll("ep_0077", "ep_0077_cal");
                plot.histName = source->GetName();
                plot.histName.ReplaceAll("hist_e_", "hist_e_cal_");
                plot.title = source->GetTitle();
                plot.title += "; calibrated";
                plot.det = det;
                plot.qqq5 = TString(source->GetTitle()).BeginsWith("QQQ5");
                if (plot.qqq5)
                    sscanf(source->GetTitle(), "QQQ5 det#%d", &plot.detNumber);
                if (plot.qqq5)
                    plot.hist = new TH2D(plot.histName, plot.title, 32, 0, 32, 1000, 0, 20);
                else
                    plot.hist = new TH2D(plot.histName, plot.title, 200, -1, 1, 1000, 0, 20);
                plot.hist->SetDirectory(nullptr);
                plots.push_back(plot);
            }
        }
    }

    gROOT->cd();
    for (auto &plot : plots) {
        plot.hist->SetDirectory(gDirectory);
        if (plot.qqq5) {
            // Build one strip-dependent expression and scan the tree once per
            // detector.  The former loop scanned the 1.9 GB tree 32 times.
            TString energy = "0";
            for (auto strip = 31; strip >= 0; --strip) {
                auto it = c0.find({plot.detNumber, 0, strip});
                if (it == c0.end() || it->second[1] <= 0)
                    continue;
                auto intercept = it->second[0];
                auto slope = it->second[1];
                energy = Form("(SiChannel.fStrip==%d?(%0.17g+%0.17g*SiChannel.fEnergy):(%s))",
                    strip, intercept, slope, energy.Data());
            }
            TString cut = Form("SiChannel.fDetNum==%d&&SiChannel.fSide==0"
                "&&SiChannel.fEnergy>0&&%s>0", plot.detNumber, energy.Data());
            recoTree->Draw(Form("%s:SiChannel.fStrip>>%s", energy.Data(), plot.histName.Data()),
                cut, "goff");
            Info("make_ep_0077_cal", "QQQ5 det#%d (idx %d): %.0f calibrated entries",
                plot.detNumber, plot.det, plot.hist->GetEntries());
            continue;
        }
        auto xE = plot.qqq5 ? "SiHit.fJunctionStrip" : "SiHit.fRelativeZ";
        auto xdE = plot.qqq5 ? "SiHit.fJunctionStrip" : "SiHit.fRelativeZdE";
        tree->Draw(Form("SiHit.fEnergy:%s>>%s", xE, plot.histName.Data()),
            Form("SiHit.fDetID==%d&&SiHit.fEnergy>0", plot.det), "goff");
        tree->Draw(Form("SiHit.fdE:%s>>+%s", xdE, plot.histName.Data()),
            Form("SiHit.fdEDetID==%d&&SiHit.fdE>0", plot.det), "goff");
        Info("make_ep_0077_cal", "det %d: %.0f entries", plot.det, plot.hist->GetEntries());
    }

    std::map<TString, std::vector<TH2D *>> groups;
    for (auto &plot : plots)
        groups[plot.canvasName].push_back(plot.hist);

    TFile output(outputName, "RECREATE");
    for (auto &group : groups) {
        TString referenceName = group.first;
        referenceName.ReplaceAll("ep_0077_cal", "ep_0077");
        auto referenceCanvas = dynamic_cast<TCanvas *>(reference.Get(referenceName));
        if (referenceCanvas == nullptr)
            continue;
        auto canvas = dynamic_cast<TCanvas *>(referenceCanvas->Clone(group.first));
        canvas->SetName(group.first);
        canvas->SetTitle(group.first);
        for (auto i = 0u; i < group.second.size(); ++i) {
            auto pad = dynamic_cast<TPad *>(canvas->GetListOfPrimitives()->At(i));
            if (pad == nullptr)
                continue;
            pad->cd();
            pad->Clear();
            group.second[i]->Draw("colz");
            group.second[i]->Write();
        }
        canvas->Modified();
        canvas->Update();
        canvas->Write();
    }
    output.Close();
    Info("make_ep_0077_cal", "Created %s with %zu detector histograms", outputName, plots.size());
}
