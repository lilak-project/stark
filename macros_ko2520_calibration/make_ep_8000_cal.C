#include "TCanvas.h"
#include "TChain.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2D.h"
#include "TKey.h"
#include "TTreeReader.h"
#include "TTreeReaderArray.h"

#include <array>
#include <fstream>
#include <map>
#include <string>
#include <tuple>
#include <vector>

namespace {
struct Cal {
    double c0i=0, c0s=0, c1iL=0, c1sL=1, c1iR=0, c1sR=1;
    double c2b0=5.486, c2b1=0, c2b2=0, c3i=0, c3s=1;
};
using Key = std::tuple<int,int,int>;

std::map<Key,Cal> load_calibration(const char *name)
{
    std::map<Key,Cal> result;
    std::ifstream input(name);
    std::string header;
    std::getline(input,header);
    int det,side,strip;
    Cal p;
    while (input>>det>>side>>strip>>p.c0i>>p.c0s
        >>p.c1iL>>p.c1sL>>p.c1iR>>p.c1sR
        >>p.c2b0>>p.c2b1>>p.c2b2>>p.c3i>>p.c3s)
        result[{det,side,strip}]=p;
    return result;
}

TH2 *find_hist(TPad *pad)
{
    if (!pad) return nullptr;
    TIter next(pad->GetListOfPrimitives());
    while (auto object=next())
        if (auto hist=dynamic_cast<TH2*>(object)) return hist;
    return nullptr;
}
}

void make_ep_8000_cal(const char *inputPath="data_reco",
    const char *outputName="ep_8000_cal.root")
{
    TString referenceName=Form("%s/ep_8000.root",inputPath);
    const char *calibrationName="data_calibration/run_0077/ko2520_0077.energy_calibration.dat";

    TFile reference(referenceName,"READ");
    auto calibration=load_calibration(calibrationName);
    if (reference.IsZombie() || calibration.empty()) {
        Error("make_ep_8000_cal","Cannot read reference or calibration file");
        return;
    }

    TChain chain("event");
    for (int run=83;run<=89;++run)
        chain.Add(Form("%s/ko2520_%04d.reco.root",inputPath,run));

    struct CanvasSet { TCanvas *canvas=nullptr; std::vector<TH2D*> histograms; };
    std::vector<CanvasSet> canvases;
    std::map<int,TH2D*> junction;
    std::map<int,TH2D*> ohmic;

    TIter nextKey(reference.GetListOfKeys());
    while (auto key=dynamic_cast<TKey*>(nextKey())) {
        auto original=dynamic_cast<TCanvas*>(key->ReadObj());
        if (!original) continue;
        TString newName=original->GetName();
        newName.ReplaceAll("8000","8000_cal");
        auto canvas=dynamic_cast<TCanvas*>(original->Clone(newName));
        canvas->SetName(newName);
        canvas->SetTitle(newName);
        CanvasSet set;
        set.canvas=canvas;

        bool isMultiplicity=TString(original->GetName()).Contains("multiplicity");
        if (!isMultiplicity) {
            TIter nextPad(canvas->GetListOfPrimitives());
            while (auto object=nextPad()) {
                auto pad=dynamic_cast<TPad*>(object);
                auto source=find_hist(pad);
                if (!source) continue;
                int detIndex=-1;
                if (sscanf(source->GetName(),"hist_e_vs_position_det_%d",&detIndex)!=1)
                    sscanf(source->GetName(),"hist_ohmic_energy_det_%d",&detIndex);
                TString histName=source->GetName();
                histName.ReplaceAll("hist_e_","hist_e_cal_");
                histName.ReplaceAll("hist_ohmic_energy_","hist_ohmic_energy_cal_");
                TString title=source->GetTitle();
                auto hist=new TH2D(histName,title,
                    source->GetNbinsX(),source->GetXaxis()->GetXmin(),source->GetXaxis()->GetXmax(),
                    1000,0,20);
                hist->SetDirectory(nullptr);
                hist->GetXaxis()->SetTitle(source->GetXaxis()->GetTitle());
                hist->GetYaxis()->SetTitle("Calibrated energy (MeV)");
                pad->cd();
                pad->Clear();
                hist->Draw("colz");
                set.histograms.push_back(hist);
                if (TString(original->GetName()).BeginsWith("cep_")) junction[detIndex]=hist;
                else ohmic[detIndex]=hist;
            }
        }
        canvases.push_back(set);
    }

    TTreeReader reader(&chain);
    TTreeReaderArray<Int_t> detIndex(reader,"SiChannel.fPadID");
    TTreeReaderArray<Int_t> detNumber(reader,"SiChannel.fDetNum");
    TTreeReaderArray<Bool_t> side(reader,"SiChannel.fSide");
    TTreeReaderArray<Int_t> strip(reader,"SiChannel.fStrip");
    TTreeReaderArray<Bool_t> direction(reader,"SiChannel.fDirection");
    TTreeReaderArray<Double_t> energyRRaw(reader,"SiChannel.fEnergy");
    TTreeReaderArray<Double_t> energyLRaw(reader,"SiChannel.fEnergy2");

    Long64_t events=0, junctionFills=0, ohmicFills=0;
    while (reader.Next()) {
        ++events;
        for (auto i=0u;i<detIndex.GetSize();++i) {
            auto key=Key{detNumber[i],(int)side[i],strip[i]};
            auto found=calibration.find(key);
            if (found==calibration.end()) continue;
            const auto &p=found->second;

            if (side[i]==1) {
                auto h=ohmic.find(detIndex[i]);
                auto energy=p.c0i+p.c0s*energyRRaw[i];
                if (h!=ohmic.end() && p.c0s>0 && energy>0) {
                    h->second->Fill(strip[i],energy);
                    ++ohmicFills;
                }
                continue;
            }

            auto h=junction.find(detIndex[i]);
            if (h==junction.end()) continue;
            if (detNumber[i]>=1 && detNumber[i]<=4) {
                auto energy=p.c0i+p.c0s*energyRRaw[i];
                if (p.c0s>0 && energy>0) {
                    h->second->Fill(strip[i],energy);
                    ++junctionFills;
                }
            }
            else if (direction[i]==0 && energyRRaw[i]>0 && energyLRaw[i]>0) {
                auto energyL=p.c1iL+p.c1sL*energyLRaw[i];
                auto energyR=p.c1iR+p.c1sR*energyRRaw[i];
                auto sum=energyL+energyR;
                if (sum<=0) continue;
                auto position=(energyR-energyL)/sum;
                auto scale=p.c2b0+p.c2b1*position+p.c2b2*position*position;
                if (scale==0) continue;
                auto energy=p.c3i+p.c3s*(sum/scale*5.486);
                if (energy>0) {
                    h->second->Fill(position,energy);
                    ++junctionFills;
                }
            }
        }
        if (events%100000==0)
            Info("make_ep_8000_cal","event %lld / %lld",events,chain.GetEntries());
    }

    TFile output(outputName,"RECREATE");
    for (auto &set:canvases) {
        for (auto hist:set.histograms) hist->Write();
        set.canvas->Modified();
        set.canvas->Update();
        set.canvas->Write();
    }
    TNamed sourceInfo("calibration_source",calibrationName);
    sourceInfo.Write();
    output.Close();
    Info("make_ep_8000_cal","Created %s: events=%lld, junction=%lld, ohmic=%lld",
        outputName,events,junctionFills,ohmicFills);
}
