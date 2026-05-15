void create_flat_data(int runNo=95)
{
    vector<int> runNo4 = {20,21,22,23,24,25,29,30,32,33,34,35,36,38};
    vector<int> runNo5 = {39,40,41,42,44,45,46,47,48,49,50,51,52,53,54,56,57,58,59,61,62,63,65,66,67,68,69,72};
    vector<int> runNo6 = {92,94,95,96,97,99,100,101,102,103,104,105,106,107,110,111};
    vector<int> runNoB = {73,76};

    for (auto list : {runNo4, runNo5, runNo6, runNoB})
    {
        for (auto runNo : list)
        {
            auto file = new TFile(Form("data_reco/stark_%04d.reco8.root",runNo));
            auto tree = (TTree*) file -> Get("event");
            if (tree==nullptr) {
                cout << "ERROR" << endl;
                continue;
            }
            TClonesArray *RecoHit = nullptr;
            tree -> SetBranchAddress("RecoHit", &RecoHit);
            tree -> SetBranchStatus("SiChannel", false);
            tree -> SetBranchStatus("SiHit_12RdE", false);
            tree -> SetBranchStatus("SiHit_12RE", false);
            tree -> SetBranchStatus("SiHit_16R", false);
            tree -> SetBranchStatus("RecoHit", true);

            auto file2 = new TFile(Form("data_flat/run_%04d.root",runNo),"recreate");
            auto tree2 = new TTree("event","triggered events");

            int    dEIndex, eIndex, detID, pairID, isInGate, junctionStrip, ohmicStrip, numHits;
            bool   isEPairDetector;
            double keyEnergy, qValue, eJunction, eOhmic, dEJunction, dEOhmic, phi, theta, solidAngle, relativeZ, realZ;

            //tree2 -> Branch("n"              , &numHits);
            //tree2 -> Branch("dEIndex"        , &dEIndex);
            //tree2 -> Branch("eIndex"         , &eIndex);
            //tree2 -> Branch("detID"          , &detID);
            tree2 -> Branch("pairID"         , &pairID);
            //tree2 -> Branch("isInGate"       , &isInGate);
            //tree2 -> Branch("isEPairDetector", &isEPairDetector);
            tree2 -> Branch("keyEnergy"      , &keyEnergy);
            //tree2 -> Branch("qValue"         , &qValue);
            tree2 -> Branch("eJunction"      , &eJunction);
            tree2 -> Branch("eOhmic"         , &eOhmic);
            tree2 -> Branch("dEJunction"     , &dEJunction);
            tree2 -> Branch("dEOhmic"        , &dEOhmic);
            tree2 -> Branch("phi"            , &phi);
            tree2 -> Branch("theta"          , &theta);
            //tree2 -> Branch("solidAngle"     , &solidAngle);
            tree2 -> Branch("relativeZ"      , &relativeZ);
            tree2 -> Branch("realZ"          , &realZ);
            tree2 -> Branch("junctionStrip"  , &junctionStrip);
            tree2 -> Branch("ohmicStrip"     , &ohmicStrip);

            auto numEvents = tree -> GetEntries();
            for (auto iEvent=0; iEvent<numEvents; ++iEvent)
            {
                if (iEvent%5000==0) cout << runNo << ": " << iEvent << " / " << numEvents << " " << numHits << " (" << 100*double(iEvent)/numEvents << " %)" << endl;
                tree -> GetEntry(iEvent);

                numHits = RecoHit -> GetEntries();
                if (numHits!=1)
                    continue;

                auto hit = (EKRecoHit*) RecoHit -> At(0); 
                dEIndex         = hit->GetdEIndex();
                eIndex          = hit->GetEIndex();
                detID           = hit->GetDetID();
                pairID          = hit->GetPairID();
                isInGate        = hit->GetIsInGate();
                isEPairDetector = hit->GetIsEPairDetector();
                keyEnergy       = hit->GetKeyEnergy();
                qValue          = hit->GetQValue();
                eJunction       = hit->GetEJunction();
                eOhmic          = hit->GetEOhmic();
                dEJunction      = hit->GetdEJunction();
                dEOhmic         = hit->GetdEOhmic();
                phi             = hit->GetPhi();
                theta           = hit->GetTheta();
                solidAngle      = hit->GetSolidAngle();
                relativeZ       = hit->GetRelativeZ();
                realZ           = hit->GetRealZ();
                junctionStrip   = hit->GetJunctionStrip();
                ohmicStrip      = hit->GetOhmicStrip();
                tree2 -> Fill();
            }
            file2 -> cd();
            tree2 -> Write();
            cout << file2 -> GetName() << " " << tree2 -> GetEntries() << endl;

            file -> Close();
            file2 -> Close();
        }
    }
}
