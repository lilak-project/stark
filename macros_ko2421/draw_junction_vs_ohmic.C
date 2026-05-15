#include "LKLogger.h"

void draw_junction_vs_ohmic()
{
    gStyle -> SetOptFit(1);
    auto elark = new ELARK();
    elark -> AddPar("config_stark.mac");
    elark -> Init();

    auto top = new LKDrawingGroup("si_junction_vs_ohmic");
    top -> AddOption("wide_canvas");

    //TString tag = "reco5";
    TString tag = "reco_xmulth";

    auto tree = new TChain("event");
    for (auto i=20;  i<=38;  ++i) { TString name = Form("/home/ejungwoo/lilak/ko2421/macros/data_reco/stark_%04d.%s.root",i,tag.Data()); cout << name << endl; tree -> AddFile(name); }
    for (auto i=39;  i<=57;  ++i) { TString name = Form("/home/ejungwoo/lilak/ko2421/macros/data_reco/stark_%04d.%s.root",i,tag.Data()); cout << name << endl; tree -> AddFile(name); }
    for (auto i=59;  i<=72;  ++i) { TString name = Form("/home/ejungwoo/lilak/ko2421/macros/data_reco/stark_%04d.%s.root",i,tag.Data()); cout << name << endl; tree -> AddFile(name); }
    for (auto i=81;  i<=83;  ++i) { TString name = Form("/home/ejungwoo/lilak/ko2421/macros/data_reco/stark_%04d.%s.root",i,tag.Data()); cout << name << endl; tree -> AddFile(name); }
    for (auto i=92;  i<=92;  ++i) { TString name = Form("/home/ejungwoo/lilak/ko2421/macros/data_reco/stark_%04d.%s.root",i,tag.Data()); cout << name << endl; tree -> AddFile(name); }
    for (auto i=94;  i<=97;  ++i) { TString name = Form("/home/ejungwoo/lilak/ko2421/macros/data_reco/stark_%04d.%s.root",i,tag.Data()); cout << name << endl; tree -> AddFile(name); }
    for (auto i=99;  i<=107; ++i) { TString name = Form("/home/ejungwoo/lilak/ko2421/macros/data_reco/stark_%04d.%s.root",i,tag.Data()); cout << name << endl; tree -> AddFile(name); }
    for (auto i=110; i<=111; ++i) { TString name = Form("/home/ejungwoo/lilak/ko2421/macros/data_reco/stark_%04d.%s.root",i,tag.Data()); cout << name << endl; tree -> AddFile(name); }

    TClonesArray* array1 = nullptr;
    TClonesArray* array2 = nullptr;
    TClonesArray* array3 = nullptr;
    tree -> SetBranchAddress("SiHit_12RdE",&array1);
    tree -> SetBranchAddress("SiHit_12RE",&array2);
    tree -> SetBranchAddress("SiHit_16R",&array3);

    TH2D* histJvsO[40][12];
    for (auto detID=0; detID<40; ++detID)
    {
        auto detector = elark -> GetSiDetector(detID);
        TString name = detector -> GetDetTypeName();
        if (detID<12) name = name + " 12-Ring E";
        else if (detID>=28) name = name + " 12-Ring dE";
        else name = name + " 16-Ring";
        for (auto strip=0; strip<12; ++strip)
        {
            if (strip<8)
                histJvsO[detID][strip] = new TH2D(Form("hist_%d_jstrip_%d",detID,strip),Form("Det-%d (%s) Junction Strip-%d;Ohmic Energy;Junction Energy",detID,name.Data(),strip),200,0,50,200,0,50);
            else
                histJvsO[detID][strip] = new TH2D(Form("hist_%d_ostrip_%d",detID,strip-8),Form("Det-%d (%s) Ohmic Strip-%d;Ohmic Energy;Junction Energy",detID,name.Data(),strip-8),200,0,50,200,0,50);
        }
    }

    auto group = top -> CreateGroup("junction_vs_ohmic");
    //for (auto detID : {15})
    for (auto detID=0; detID<40; ++detID)
    {
        auto sub = group -> CreateGroup(Form("det-%d",detID));
        for (auto strip=0; strip<12; ++strip)
        {
            auto f1 = new TF1(Form("f1%d%d",detID,strip),"x",0,50);
            f1 -> SetLineWidth(1);
            //f1 -> SetLineStyle(2);
            f1 -> SetLineColor(kBlack);
            //auto fit = new TF1(Form("fit%d%d",detID,strip),"[0]*x",0,50);
            //fit -> SetParameter(0,1);
            //fit -> SetLineWidth(1);
            //fit -> SetLineColor(kRed);
            //histJvsO[detID][strip] -> Fit(fit,"R0Q");
            auto drawing = sub -> CreateDrawing();
            drawing -> Add(histJvsO[detID][strip]);
            drawing -> Add(f1);
            //drawing -> Add(fit);
            drawing -> AddOption("stats_tl_corner");
            drawing -> SetOptStat(10);
            drawing -> SetOptFit(1);
            drawing -> SetLogz();
        }
    }

    auto numEvents = tree -> GetEntries();
    for (auto iEvent=0; iEvent<numEvents; ++iEvent)
    {
        if (iEvent%100000==0) e_info << iEvent << endl;
        tree -> GetEntry(iEvent);
        for (auto array : {array1,array2,array3})
        {
            auto n = array -> GetEntries();
            for (auto iHit=0; iHit<n; ++iHit)
            {
                auto hit = (EKSiHit*) array -> At(iHit);
                auto detID = hit -> GetDetID();
                auto z = hit -> GetRelativeZ();
                auto ej = hit -> GetEnergy();
                auto eo = hit -> GetEnergyOhmic();
                auto sj = hit -> GetJunctionStrip();
                auto so = hit -> GetOhmicStrip();
                if (detID<0) 
                    continue;
                histJvsO[detID][sj] -> Fill(eo,ej);
                histJvsO[detID][so+8] -> Fill(eo,ej);
            }
        }
    }

    //top -> Draw("v:save_all");
    top -> Draw("v");
}
