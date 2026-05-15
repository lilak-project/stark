#include "LKLogger.h"

void draw_relative_z()
{
    bool findPoint = false;
    auto elark = new ELARK();
    elark -> AddPar("config_stark.mac");
    elark -> Init();

    auto top = new LKDrawingGroup("relative_position");
    //if (findPoint==false) top -> AddOption("vertical_canvas");

    auto tree = new TChain("event");

    // 4
    for (auto i=20; i<=21; ++i) { TString name = Form("/Users/ejungwoo/lilak/ko2421/macros/data_reco/stark_%04d.reco8.root",i); cout << name << endl; tree -> AddFile(name); }
    for (auto i=20; i<=38; ++i) { TString name = Form("/Users/ejungwoo/lilak/ko2421/macros/data_reco/stark_%04d.reco8.root",i); cout << name << endl; tree -> AddFile(name); }

    // 5
    //for (auto i=39; i<=57; ++i) { TString name = Form("/home/ejungwoo/lilak/ko2421/macros/data_reco/stark_%04d.reco8.root",i); cout << name << endl; tree -> AddFile(name); }
    //for (auto i=59; i<=72; ++i) { TString name = Form("/home/ejungwoo/lilak/ko2421/macros/data_reco/stark_%04d.reco8.root",i); cout << name << endl; tree -> AddFile(name); }

    // 8
    //for (auto i=81; i<=83; ++i) { TString name = Form("/home/ejungwoo/lilak/ko2421/macros/data_reco/stark_%04d.reco8.root",i); cout << name << endl; tree -> AddFile(name); }
    //for (auto i=92; i<=92; ++i) { TString name = Form("/home/ejungwoo/lilak/ko2421/macros/data_reco/stark_%04d.reco8.root",i); cout << name << endl; tree -> AddFile(name); }
    //for (auto i=94; i<=97; ++i) { TString name = Form("/home/ejungwoo/lilak/ko2421/macros/data_reco/stark_%04d.reco8.root",i); cout << name << endl; tree -> AddFile(name); }
    //for (auto i=99; i<=107; ++i) { TString name = Form("/home/ejungwoo/lilak/ko2421/macros/data_reco/stark_%04d.reco8.root",i); cout << name << endl; tree -> AddFile(name); }
    //for (auto i=110; i<=111; ++i) { TString name = Form("/home/ejungwoo/lilak/ko2421/macros/data_reco/stark_%04d.reco8.root",i); cout << name << endl; tree -> AddFile(name); }

    TClonesArray* channelArray = nullptr;
    //TClonesArray* array1 = nullptr;
    //TClonesArray* array2 = nullptr;
    //TClonesArray* array3 = nullptr;
    tree -> SetBranchAddress("SiChannel",&channelArray);
    //tree -> SetBranchAddress("SiHit_12RdE",&array1);
    //tree -> SetBranchAddress("SiHit_12RE",&array2);
    //tree -> SetBranchAddress("SiHit_16R",&array3);

    auto tree2 = new TChain("event");
    tree2 -> AddFile("/Users/ejungwoo/lilak/ko2421/macros/data_ana/KO2421_e4_cutX_ProtonCut_NT180.ana_eve.root");
    tree2 -> AddFile("/Users/ejungwoo/lilak/ko2421/macros/data_ana/KO2421_e5_cutX_ProtonCut_NT180.ana_eve.root");
    tree2 -> AddFile("/Users/ejungwoo/lilak/ko2421/macros/data_ana/KO2421_e8_cutX_ProtonCut_NT180.ana_eve.root");
    TClonesArray* recoArray = nullptr;
    tree2 -> SetBranchAddress("RecoHit",&recoArray);

    int countChannels = 0;
    TH2D* histStrip[40];
    for (auto detID=0; detID<40; ++detID)
    {
        auto detector = elark -> GetSiDetector(detID);
        TString name = detector -> GetDetTypeName();
        if (detID<12) name = name + " 12-Ring E";
        else if (detID>=28) name = name + " 12-Ring dE";
        else name = name + " 16-Ring";
        histStrip[detID] = new TH2D(Form("hist_%d_junction",detID),Form("Det-%d (%s);Relative-position;E_{Ohmic}",detID,name.Data()),200,-1.5,1.5,200,0,30);
    }

    auto group1 = top -> CreateGroup("relative_position",!findPoint);
    auto group2 = top -> CreateGroup("channels",findPoint);
    auto sub1 = group1 -> CreateGroup("ring_12_1");
    auto sub2 = group1 -> CreateGroup("ring_12_2");
    auto sub3 = group1 -> CreateGroup("ring_16_1");
    auto sub4 = group1 -> CreateGroup("ring_16_2");
    for (auto detID=0; detID<6; ++detID)
    {
        auto drawing = sub1 -> CreateDrawing();
        drawing -> Add(histStrip[detID]);
        drawing -> SetGridx();
        drawing -> SetGridy();
        drawing -> SetOptStat(10);
    }
    for (auto detID=6; detID<12; ++detID)
    {
        auto drawing = sub2 -> CreateDrawing();
        drawing -> Add(histStrip[detID]);
        drawing -> SetGridx();
        drawing -> SetGridy();
        drawing -> SetOptStat(10);
    }
    for (auto detID=12; detID<20; ++detID)
    {
        auto drawing = sub3 -> CreateDrawing();
        drawing -> Add(histStrip[detID]);
        drawing -> SetGridx();
        drawing -> SetGridy();
        drawing -> SetOptStat(10);
    }
    for (auto detID=20; detID<28; ++detID)
    {
        auto drawing = sub4 -> CreateDrawing();
        drawing -> Add(histStrip[detID]);
        drawing -> SetGridx();
        drawing -> SetGridy();
        drawing -> SetOptStat(10);
    }

    //auto numEvents = tree -> GetEntries();
    auto numEvents = tree2 -> GetEntries();
    for (auto iEvent=0; iEvent<numEvents; ++iEvent)
    {
        if (iEvent%100000==0) e_info << iEvent << endl;
        //tree -> GetEntry(iEvent);
        tree2 -> GetEntry(iEvent);
        TString note0;
        TString note1;
        auto numHits = recoArray -> GetEntries();
        bool foundPoint = false;
        for (auto iHit=0; iHit<numHits; ++iHit)
        {
            auto hit = (EKRecoHit*) recoArray -> At(iHit);
            auto z = hit -> GetRelativeZ();
            auto eo = hit -> GetEOhmic();
            auto detID = hit -> GetDetID();
            histStrip[detID] -> Fill(z,eo);
            if (findPoint&&eo>25.) {
                foundPoint = true;
                note0 = Form("hit_detID=%d",detID);
                note1 = Form("hit_energy=%.2f",eo);
            }
        }
        if (foundPoint)
        {
            auto sub = group2 -> CreateGroup(Form("event-%d",iEvent));
            auto nHits = channelArray -> GetEntries();
            for (auto iHit=0; iHit<nHits; ++iHit)
            {
                auto channel = (LKSiChannel*) channelArray -> At(iHit);
                auto drawing = sub -> CreateDrawing();
                drawing -> Add(channel->GetHist(Form("channel%d%d",iEvent,iHit)));
                auto pv = new TPaveText();
                pv -> AddText(note0);
                pv -> AddText(note1);
                pv -> AddText(Form("ch_detID=%d",channel->GetDetID()));
                pv -> AddText(Form("ch_energy=%.2f",channel->GetEnergy()));
                drawing -> Add(pv);
                drawing -> SetPaveCorner(2);
            }
            if (countChannels++>1)
                break;
        }
    }

    //top -> Draw("v");
    top -> Draw();
    top -> Save();
}
