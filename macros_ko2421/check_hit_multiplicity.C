void check_hit_multiplicity()
{
    auto tree = new TChain("event");
    for (auto i=20;  i<=21;  ++i) { TString name = Form("/home/ejungwoo/lilak/ko2421/macros/data_reco/stark_%04d.reco5.root",i); cout << name << endl; tree -> AddFile(name); }
    //for (auto i=20;  i<=38;  ++i) { TString name = Form("/home/ejungwoo/lilak/ko2421/macros/data_reco/stark_%04d.reco5.root",i); cout << name << endl; tree -> AddFile(name); }
    //for (auto i=39;  i<=57;  ++i) { TString name = Form("/home/ejungwoo/lilak/ko2421/macros/data_reco/stark_%04d.reco5.root",i); cout << name << endl; tree -> AddFile(name); }
    //for (auto i=59;  i<=72;  ++i) { TString name = Form("/home/ejungwoo/lilak/ko2421/macros/data_reco/stark_%04d.reco5.root",i); cout << name << endl; tree -> AddFile(name); }
    //for (auto i=81;  i<=83;  ++i) { TString name = Form("/home/ejungwoo/lilak/ko2421/macros/data_reco/stark_%04d.reco5.root",i); cout << name << endl; tree -> AddFile(name); }
    //for (auto i=92;  i<=92;  ++i) { TString name = Form("/home/ejungwoo/lilak/ko2421/macros/data_reco/stark_%04d.reco5.root",i); cout << name << endl; tree -> AddFile(name); }
    //for (auto i=94;  i<=97;  ++i) { TString name = Form("/home/ejungwoo/lilak/ko2421/macros/data_reco/stark_%04d.reco5.root",i); cout << name << endl; tree -> AddFile(name); }
    //for (auto i=99;  i<=107; ++i) { TString name = Form("/home/ejungwoo/lilak/ko2421/macros/data_reco/stark_%04d.reco5.root",i); cout << name << endl; tree -> AddFile(name); }
    //for (auto i=110; i<=111; ++i) { TString name = Form("/home/ejungwoo/lilak/ko2421/macros/data_reco/stark_%04d.reco5.root",i); cout << name << endl; tree -> AddFile(name); }

    TString tag[3] = { "12RdE", "12RE", "16R" };
    TClonesArray* array[3];
    TH1D* hist[3];
    TH1D* histDouble = new TH1D("histDouble","",4,0,4);
    auto top = new LKDrawingGroup();

    for (auto i : {0,1,2}) {
        array[i] = nullptr;
        tree -> SetBranchAddress(Form("SiHit_%s",tag[i].Data()),&array[i]);
        hist[i] = new TH1D(Form("hist%d",i),"",4,0,4);
        top -> AddHist(hist[0]);
    }
    top -> AddHist(histDouble);

    int numEvents = tree -> GetEntries();
    for (auto iEvent=0; iEvent<numEvents; ++iEvent)
    {
        tree -> GetEntry(iEvent);
        int numHits[3];
        for (auto i : {0,1,2}) {
            numHits[i] = array[i] -> GetEntries();
            hist[i] -> Fill(numHits[i]);
        }
        int numHits01 = ((numHits[0]>0||numHits[1]>0)?(numHits[0]+numHits[1]):0);
        if (numHits01>0 && numHits[2]>0)
            histDouble -> Fill(0);
        else if (numHits01>0)
            histDouble -> Fill(1);
        else if (numHits[2]>0)
            histDouble -> Fill(2);
        else
            histDouble -> Fill(4);
    }

    top -> Draw();
}
