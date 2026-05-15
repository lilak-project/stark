void check_event()
{
    auto run = new LKRun();
    run -> AddDetector(new STARK());
    run -> AddPar("config_check.mac");
    run -> Init();
    run -> SetAutoTermination(false);
    auto runID = run -> GetRunID();
    auto tree = run -> GetInputTree();
    TClonesArray* array = nullptr;

    tree -> SetBranchAddress("SiHit",&array);
    auto numEvents = tree -> GetEntries();
    //cout << "Number of events: " << numEvents << endl;

    for (auto iEvent=0; iEvent<numEvents; ++iEvent)
    {
        tree -> GetEntry(iEvent);
        auto numHits = array -> GetEntries();
        if (iEvent%20000==0) cout << "Event-" << iEvent << ", Number of hits: " << numHits << endl;
        bool pass = false;
        for (auto iHit=0; iHit<numHits; ++iHit)
        {
            auto hit = (SKSiHit *) array -> At(iHit);
            if (hit -> GetdEDetID()==37 && hit -> GetdE()>10) {
                lk_debug << iEvent << " " << hit -> GetdE() << endl;
                //hist2 -> Fill(hit->GetStrip(),hit->GetEnergy());
                //pass = true;
                //break;
            }
        }
    }
}
