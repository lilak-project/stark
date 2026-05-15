#include "LKLogger.h"

TString fAnaName, fMainName;
TTree* fTree = nullptr;
ELARK* fELARK = nullptr;
LKDrawingGroup* fTop = nullptr;

void SetExperimentEnergyIndex();

void draw_phi_distribution(int energyIndex=8, int optionIndex=3)
{
    TString path = "data_ana";
    TString cut_pid_name = "all";
    if (optionIndex==0) cut_pid_name = "cutX_ProtonAndStop";
    else if (optionIndex==1) cut_pid_name = "cutX_ProtonAndERange";
    else if (optionIndex==2) cut_pid_name = "cutX_Stopped";
    else if (optionIndex==3) cut_pid_name = "cutX_ProtonInERange";

    fAnaName = Form("KO2421_e%d_%s",energyIndex,cut_pid_name.Data());
    if (energyIndex==1) fAnaName = "KO2421_blank";
    SetExperimentEnergyIndex();

    TFile* file;
    if (energyIndex==1) file = new TFile("/home/ejungwoo/lilak/ko2421/macros/data_ana/KO2421_blank_cutX_pidAll_NT360.ana_eve.root");
    else file = new TFile(Form("%s/%s_NT180.ana_eve.root",path.Data(),fAnaName.Data()));
    if (file->IsOpen()==false) { e_error << file->GetName() << endl; return; }
    auto runHeader = (LKParameterContainer*) file -> Get("RunHeader");
    fMainName = runHeader -> GetParString("MainName");
    if (energyIndex==1) fMainName = fMainName(0,13);
    else fMainName = fMainName(0,9);
    fMainName.ReplaceAll("_e","-E");
    fTree = (TTree*) file -> Get("event");

    auto hist1 = new TH1D("hist1",fMainName+" 12-Ring",12,0,360);
    for (auto detID=0; detID<12; ++detID) {
        auto detector = fELARK -> GetSiDetector(detID);
        auto detID2 = fELARK -> FindEPairDetectorID(detID);
        auto detector2 = fELARK -> GetSiDetector(detID2);
        TString content = Form("%d-%s/%s",detID,detector->GetDetTypeName().Data(),detector2->GetDetTypeName().Data());;
        auto phi = detector -> GetPhi0();
        auto bin = hist1 -> FindBin(phi);
        hist1 -> GetXaxis() -> SetBinLabel(bin,content);
    }

    auto hist2 = new TH1D("hist2",fMainName+" 16-Ring",16,0,360);
    for (auto detID=12; detID<12+16; ++detID) {
        auto detector = fELARK -> GetSiDetector(detID);
        TString content = Form("%d-%s",detID,detector->GetDetTypeName().Data());
        auto phi = detector -> GetPhi0();
        auto bin = hist2 -> FindBin(phi);
        hist2 -> GetXaxis() -> SetBinLabel(bin,content);
    }

    TString selection = "fKeyEnergy>0&&fKeyEnergy<30&&fTheta>10&&fTheta<90";
    fTree -> Draw("fPhi>>hist1",selection+"&&"+"fIsEPairDetector","goff");
    fTree -> Draw("fPhi>>hist2",selection+"&&"+"!fIsEPairDetector","goff");

    auto graphAtt = new TGraph();
    graphAtt -> SetLineWidth(2);
    graphAtt -> SetLineColor(kBlue);
    graphAtt -> SetFillColor(kAzure-4);

    hist1 -> Draw();
    auto draw1 = LKSAM::GetSAM()->MakeTH1Polar(hist1,6,10,0,12,true,true,0,0.1,0.1,graphAtt);
    auto draw2 = LKSAM::GetSAM()->MakeTH1Polar(hist2,6,10,0,16,true,true,0,0.1,0.1,graphAtt);

    fTop -> Add(draw1);
    fTop -> Add(draw2);
    //fTop -> SetCanvasSize(1200,600);
    fTop -> SetCanvasSize(1.5*1200,1.5*600);
    fTop -> Draw();
    fTop -> Save();
}

void SetExperimentEnergyIndex()
{
    fELARK = new ELARK();
    fELARK -> AddPar("config_stark.mac");
    fELARK -> Init();

    TString topName = fAnaName;
    topName.ReplaceAll("KO2421","PhiDist");
    fTop = new LKDrawingGroup(topName);
}
