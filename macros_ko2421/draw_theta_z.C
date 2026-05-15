#include "LKLogger.h"

void draw_theta_z()
{
    bool findPoint = true;

    auto elark = new ELARK();
    elark -> AddPar("config_stark.mac");
    elark -> Init();

    auto top = new LKDrawingGroup("strip_position");
    if (findPoint==false)top -> AddOption("vertical_canvas");

    auto tree = new TChain("event");
    for (auto i=20;  i<=21;  ++i) { TString name = Form("/home/ejungwoo/lilak/ko2421/macros/data_reco/stark_%04d.reco5.root",i); cout << name << endl; tree -> AddFile(name); }
    for (auto i=20;  i<=38;  ++i) { TString name = Form("/home/ejungwoo/lilak/ko2421/macros/data_reco/stark_%04d.reco5.root",i); cout << name << endl; tree -> AddFile(name); }
    for (auto i=39;  i<=57;  ++i) { TString name = Form("/home/ejungwoo/lilak/ko2421/macros/data_reco/stark_%04d.reco5.root",i); cout << name << endl; tree -> AddFile(name); }
    for (auto i=59;  i<=72;  ++i) { TString name = Form("/home/ejungwoo/lilak/ko2421/macros/data_reco/stark_%04d.reco5.root",i); cout << name << endl; tree -> AddFile(name); }
    for (auto i=81;  i<=83;  ++i) { TString name = Form("/home/ejungwoo/lilak/ko2421/macros/data_reco/stark_%04d.reco5.root",i); cout << name << endl; tree -> AddFile(name); }
    //for (auto i=92;  i<=92;  ++i) { TString name = Form("/home/ejungwoo/lilak/ko2421/macros/data_reco/stark_%04d.reco5.root",i); cout << name << endl; tree -> AddFile(name); }
    //for (auto i=94;  i<=97;  ++i) { TString name = Form("/home/ejungwoo/lilak/ko2421/macros/data_reco/stark_%04d.reco5.root",i); cout << name << endl; tree -> AddFile(name); }
    //for (auto i=99;  i<=107; ++i) { TString name = Form("/home/ejungwoo/lilak/ko2421/macros/data_reco/stark_%04d.reco5.root",i); cout << name << endl; tree -> AddFile(name); }
    //for (auto i=110; i<=111; ++i) { TString name = Form("/home/ejungwoo/lilak/ko2421/macros/data_reco/stark_%04d.reco5.root",i); cout << name << endl; tree -> AddFile(name); }

    TClonesArray* array0 = nullptr;
    TClonesArray* array1 = nullptr;
    TClonesArray* array2 = nullptr;
    TClonesArray* array3 = nullptr;
    tree -> SetBranchAddress("SiChannel",&array0);
    tree -> SetBranchAddress("SiHit_12RdE",&array1);
    tree -> SetBranchAddress("SiHit_12RE",&array2);
    tree -> SetBranchAddress("SiHit_16R",&array3);

    int countChannels = 0;

    TH2D* histStrip[40][12];
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
                histStrip[detID][strip] = new TH2D(Form("hist_%d_junction_%d",detID,strip),Form("Det-%d (%s) Junction Strip-%d",detID,name.Data(),strip),200,-1.5,1.5,200,0,30);
            else
                histStrip[detID][strip] = new TH2D(Form("hist_%d_ohmic_%d",detID,strip-8),Form("Det-%d (%s) Ohmic Strip-%d",detID,name.Data(),strip-8),200,-1.5,1.5,200,0,30);
        }
    }

    //for (auto detID : {15})
    auto group = top -> CreateGroup("e_vs_relativeZ",!findPoint);
    auto group2 = top -> CreateGroup("channels",findPoint);
    for (auto detID=0; detID<40; ++detID)
    {
        auto sub = group -> CreateGroup(Form("det-%d",detID));
        //for (auto strip=0; strip<12; ++strip) {
        for (auto strip=0; strip<8; ++strip) {
            auto drawing = sub -> CreateDrawing();
            drawing -> Add(histStrip[detID][strip]);
            drawing -> SetGridx();
            drawing -> SetGridy();
            drawing -> SetOptStat(10);
        }
    }

    auto numEvents = tree -> GetEntries();
    for (auto iEvent=0; iEvent<numEvents; ++iEvent)
    {
        if (iEvent%100000==0) e_info << iEvent << endl;
        tree -> GetEntry(iEvent);
        //for (auto array : {array1,array2,array3})
        bool foundPoint = false;
        TString note0;
        TString note1;
        //for (auto array : {array1,array2,array3})
        if (array1->GetEntries()+array2->GetEntries()+array3->GetEntries()!=1)
            continue;
        for (auto array : {array1,array3})
        {
            auto numHits = array -> GetEntries();
            if (numHits!=1)
                continue;
            for (auto iHit=0; iHit<numHits; ++iHit)
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
                //histStrip[detID][sj] -> Fill(z,ej);
                histStrip[detID][sj] -> Fill(z,eo);
                //histStrip[detID][so+8] -> Fill(z,eo);
                if (findPoint&&eo>15&&eo<16) {
                    foundPoint=true;
                    note0 = Form("hit_detID=%d",detID);
                    note1 = Form("hit_energy=%.2f",ej);
                }
            }
        }
        if (foundPoint)
        {
            auto sub = group2 -> CreateGroup(Form("event-%d",iEvent));

            auto nHits = array0 -> GetEntries();
            for (auto iHit=0; iHit<nHits; ++iHit)
            {
                auto channel = (LKSiChannel*) array0 -> At(iHit);
                auto drawing = sub -> CreateDrawing();
                auto hist = channel->GetHist(Form("channel%d_%d_%d",iEvent,iHit,countChannels));
                drawing -> Add(hist);
                auto pv = new TPaveText();
                pv -> AddText(note0);
                pv -> AddText(note1);
                //pv -> AddText(Form("ch_detID=%d",channel->GetDetID()));
                pv -> AddText(Form("ch_energy=%.2f",channel->GetEnergy()));
                drawing -> Add(pv);
                drawing -> SetPaveCorner(2);
            }
            if (countChannels++>10)
                break;
        }
    }

    //top -> Draw("v:save_all");
    top -> Draw("v");
    //top -> Draw();
    //top -> Save();
}
