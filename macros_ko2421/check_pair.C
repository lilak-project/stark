void check_pair()
{
    auto run = new LKRun();
    run -> AddDetector(new STARK());
    run -> AddPar("config_check.mac");
    run -> Init();
    run -> SetAutoTermination(false);
    auto runID = run -> GetRunID();
    auto tree = run -> GetInputTree();
    TClonesArray* array = nullptr;
    tree -> SetBranchAddress("SiChannel",&array);
    auto numEvents = tree -> GetEntries();
    //cout << "Number of events: " << numEvents << endl;

    auto cvs2 = new TCanvas("cvs2","",1600,900);
    cvs2 -> Divide(3,2);
    tree -> Draw("SiChannel.fEnergy:SiChannel.GetStrip()>>(8,0,8,200,0,8000)","SiChannel.GetDetID()==36","colz");//
    return;
    int icvs = 1;
    //for (auto sel : {3})
    //for (auto sel : {4,5,6})
    //for (auto sel : {7,8,9,10,11})
    for (auto sel : {36})
    {
        //auto hist = new TH2D(Form("hist%d",sel),Form("%d",sel),40,0,40,200,0,5000);
        auto hist2 = new TH2D(Form("hist2_%d",sel),Form("%d",sel),8,0,8,200,0,5000);
        for (auto iEvent=0; iEvent<numEvents; ++iEvent)
        {
            tree -> GetEntry(iEvent);
            auto numHits = array -> GetEntries();
            if (iEvent%20000==0) cout << "Event-" << iEvent << ", Number of hits: " << numHits << endl;
            bool pass = false;
            for (auto iHit=0; iHit<numHits; ++iHit)
            {
                auto channel = (LKSiChannel *) array -> At(iHit);
                //if (channel -> GetDetID()==sel) {
                if (channel -> GetDetID()==sel) {
                    hist2 -> Fill(channel->GetStrip(),channel->GetEnergy());
                    pass = true;
                    //break;
                }
            }
            //if (!pass)
            //    continue;

            //for (auto iHit=0; iHit<numHits; ++iHit)
            //{
            //    auto channel = (LKSiChannel *) array -> At(iHit);
            //    if (channel->GetDetID()==sel)
            //        continue;
            //    //if (channel->GetEnergy()>5)
            //    hist -> Fill(channel->GetDetID(),channel->GetEnergy());
            //}
        }

        cvs2 -> cd(icvs++);
        //hist -> Draw("colz");
        hist2 -> Draw("colz");
    }
}
