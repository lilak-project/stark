void make_all_histograms(int energyIndex=5)
{
    int histType = 2;
    TString tag = "hist8";
    TString tagTitle = "Multiplicity cut";
    TCut selection = "(RecoHit.fPairID<12) ? ((RecoHit.fPairID>=4&&RecoHit.fPairID<12)?RecoHit.fdEJunction:RecoHit.fdEOhmic)>1 : 1";
    TCut cutb4 = "RecoHit.fPairID<4";
    TCut cuta4 = "RecoHit.fPairID>=4";
    TString topName = "";

    if      (histType==0) { topName = "ETdEEXM";  tag = "reco7"; tagTitle = "No multiplicity cut"; }
    else if (histType==1) { topName = "ETdEEYM";  tag = "hist8"; tagTitle = "Multiplicity cut"; }
    else if (histType==2) { topName = "ETdEEYMG"; tag = "hist8"; tagTitle = "Multiplicity cut, Proton"; selection = "RecoHit.fIsInGate"; }
    else if (histType==3) { topName = "ETdEEXMG"; tag = "reco7"; tagTitle = "Multiplicity cut, Proton"; selection = "RecoHit.fIsInGate"; }
    else
        return;

    auto top = new LKDrawingGroup(topName);
    int energyIndices[] = {4,5,8};
    TString energies[] = {"4.41","5.838","8.13"};

    LKBinning bnn_tta("tta","#theta_{Lab}",400,25,75);
    LKBinning bnn_e1("e1","Energy",300,0,30);
    LKBinning bnn_de("de","dE",300,0,30);
    LKBinning bnn_e2("e2","Energy",500,0,50);

    for (auto i : {0,1,2})
    {
        TString title1 = Form("[%s MeV] (%s)",energies[i].Data(),tagTitle.Data());

        int energyIndex = energyIndices[i];
        TString name = Form("e%d",energyIndex);
        auto runList = LKParameterContainer("config_short_cut.mac").GetParIntRange(Form("e%d/RunIDRange",energyIndex));
        auto tree = new TChain("event");
        for (auto id : runList)
            tree -> Add(Form("data_reco/stark_%04d.%s.root",id,tag.Data()));
        cout << tree -> GetEntries() << endl;

        TString name1 = name+"_hist_b4";
        TString name2 = name+"_hist_a4";
        TString name3 = name+"_hdee_b4";
        TString name4 = name+"_hdee_a4";
        auto hist1 = (bnn_tta*bnn_e2).NewH2(name1,title1 + " X6 pair");
        auto hist2 = (bnn_tta*bnn_e2).NewH2(name2,title1 + " CSD pair");
        auto hist3 = (bnn_e1*bnn_de) .NewH2(name3,title1 + " X6 pair");
        auto hist4 = (bnn_e1*bnn_de) .NewH2(name4,title1 + " CSD pair");
        tree -> Draw(Form("RecoHit.fKeyEnergy:RecoHit.fTheta>>%s",name1.Data()),selection&&cutb4,"colz goff");
        tree -> Draw(Form("RecoHit.fKeyEnergy:RecoHit.fTheta>>%s",name2.Data()),selection&&cuta4,"colz goff");
        tree -> Draw(Form("((fPairID>=4&&fPairID<12)?RecoHit.fdEJunction:RecoHit.fdEOhmic):RecoHit.fKeyEnergy>>%s",name3.Data()),selection&&cutb4&&TCut("RecoHit.fIsEPairDetector"),"colz goff");
        tree -> Draw(Form("((fPairID>=4&&fPairID<12)?RecoHit.fdEJunction:RecoHit.fdEOhmic):RecoHit.fKeyEnergy>>%s",name4.Data()),selection&&cuta4&&TCut("RecoHit.fIsEPairDetector"),"colz goff");

        auto draw0 = top -> CreateDrawing(); draw0 -> Add(hist3); draw0 -> SetStatCorner(1); if (histType!=2) draw0 -> SetLogz();
        auto draw1 = top -> CreateDrawing(); draw1 -> Add(hist1); draw1 -> SetStatCorner(1);
        auto draw2 = top -> CreateDrawing(); draw2 -> Add(hist4); draw2 -> SetStatCorner(1); if (histType!=2) draw2 -> SetLogz();
        auto draw3 = top -> CreateDrawing(); draw3 -> Add(hist2); draw3 -> SetStatCorner(1); draw3 -> SetLogz();
    }

    top -> Draw("v:wx=1200:wy=750");
    top -> Save();
}
