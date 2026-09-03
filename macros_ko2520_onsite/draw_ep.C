void draw_ep()
{
    constexpr int maxDetectorIndex = 256;

    double energyMax=4200;
    Long64_t maxEvents=-1;

    auto par = new LKParameterContainer("get_common.mac");
    auto runID = par -> GetParInt("LKRun/RunID");
    cout << runID << endl;

    TChain *tree = nullptr;

    if (runID==8000) {
        tree = new TChain("event");
        tree -> Add(Form("data_reco/ko2520_%04d.reco.root",83));
        tree -> Add(Form("data_reco/ko2520_%04d.reco.root",84));
        tree -> Add(Form("data_reco/ko2520_%04d.reco.root",85));
        tree -> Add(Form("data_reco/ko2520_%04d.reco.root",86));
        tree -> Add(Form("data_reco/ko2520_%04d.reco.root",87));
        tree -> Add(Form("data_reco/ko2520_%04d.reco.root",88));
        tree -> Add(Form("data_reco/ko2520_%04d.reco.root",89));
        tree -> Add(Form("data_reco/ko2520_%04d.reco.root",90));
    }
    else {
        TString fileName = Form("data_reco/ko2520_%04d.reco.root",runID);
        cout << fileName << endl;
        auto file = new TFile(fileName);
        tree = (TChain*) file -> Get("event");
    }

    //TString inputFile = Form("data_reco/ko2520_%04d.reco.root",runID);

    LKSiliconMapping mapping;
    if (!mapping.Load("mapping_ko2520")) {
        cerr << "Failed to load mapping_ko2520" << endl;
        return;
    }

    //auto tree = new TChain("event");
    //if (tree->Add(inputFile) == 0) {
    //    cerr << "No input file matched: " << inputFile << endl;
    //    return;
    //}

    TString nameTop = Form("ep_%04d",runID);
    auto top1         = new LKDrawingGroup(nameTop);
    auto drawings16   = top1->CreateGroup("X6_16_ring");
    auto drawings12DE = top1->CreateGroup("X6_12dE_ring");
    auto drawings12E  = top1->CreateGroup("X6_12E_ring");
    auto drawingsQQQ5 = top1->CreateGroup("QQQ5");

    TString nameTop2 = Form("ohmic_energy_%04d",runID);
    auto top2              = new LKDrawingGroup(nameTop2);
    auto drawingsOhmic16   = top2->CreateGroup("X6_16_ring");
    auto drawingsOhmic12DE = top2->CreateGroup("X6_12dE_ring");
    auto drawingsOhmic12E  = top2->CreateGroup("X6_12E_ring");
    auto drawingsOhmicQQQ5 = top2->CreateGroup("QQQ5");

    TString nameTop3 = Form("event_multiplicity_%04d",runID);
    auto top3 = new LKDrawingGroup(nameTop3);

    enum { kAll=0, kRing16=1, kRing12DE=2, kRing12E=3, kQQQ5=4, kNumMult=5 };
    TString multName[kNumMult] = {"all", "X6_16", "X6_12dE", "X6_12E", "QQQ5"};
    TH1D *histMultiplicity[kNumMult] = {nullptr};
    for (auto i=0; i<kNumMult; ++i) {
        histMultiplicity[i] = new TH1D(
            Form("hist_multiplicity_%s",multName[i].Data()),
            Form("%s event multiplicity;number of fired channels;events",multName[i].Data()),
            100,0,100);
        auto drawing = top3->CreateDrawing(multName[i]);
        drawing->Add(histMultiplicity[i],"hist");
    }

    TH2D *histPosition[maxDetectorIndex] = {nullptr};
    TH2D *histOhmic[maxDetectorIndex] = {nullptr};
    for (auto iDet=0; iDet<mapping.GetNumDetectors(); ++iDet) {
        auto det = mapping.GetDetectorByVectorIndex(iDet);
        if (det == nullptr || det->detIndex < 0 || det->detIndex >= maxDetectorIndex)
            continue;
        if (det->detType == "X6") {
            histPosition[det->detIndex] = new TH2D(
                Form("hist_e_vs_position_det_%03d",det->detIndex),
                Form("X6 det#%d (idx=%d, %s);(E_{L}-E_{R})/(E_{L}+E_{R});E_{L}+E_{R}",
                     det->detNumber,det->detIndex,det->ringType.Data()),
                200,-1,1,300,0,energyMax);
        }
        else if (det->detType == "QQQ5") {
            histPosition[det->detIndex] = new TH2D(
                Form("hist_e_vs_position_det_%03d",det->detIndex),
                Form("QQQ5 det#%d (idx=%d);junction strip;Energy",
                     det->detNumber,det->detIndex),
                32,0,32,300,0,energyMax);
        }
        LKDrawingGroup *ringDrawings = nullptr;
        if      (det->ringType == "16")   ringDrawings = drawings16;
        else if (det->ringType == "12dE") ringDrawings = drawings12DE;
        else if (det->ringType == "12E")  ringDrawings = drawings12E;
        else if (det->ringType == "Q")    ringDrawings = drawingsQQQ5;
        if (ringDrawings != nullptr && histPosition[det->detIndex] != nullptr) {
            auto drawing = ringDrawings->CreateDrawing(Form("idx%d_det%d",det->detIndex,det->detNumber));
            drawing->Add(histPosition[det->detIndex],"colz");
            //drawing->SetLogz();
        }

        histOhmic[det->detIndex] = new TH2D(
            Form("hist_ohmic_energy_det_%03d",det->detIndex),
            Form("%s det#%d (idx=%d, %s);ohmic strip;Energy",
                 det->detType.Data(),det->detNumber,det->detIndex,det->ringType.Data()),
            4,0,4,300,0,energyMax);
        LKDrawingGroup *ohmicRingDrawings = nullptr;
        if      (det->ringType == "16")   ohmicRingDrawings = drawingsOhmic16;
        else if (det->ringType == "12dE") ohmicRingDrawings = drawingsOhmic12DE;
        else if (det->ringType == "12E")  ohmicRingDrawings = drawingsOhmic12E;
        else if (det->ringType == "Q")    ohmicRingDrawings = drawingsOhmicQQQ5;
        if (ohmicRingDrawings != nullptr) {
            auto drawing = ohmicRingDrawings->CreateDrawing(Form("idx%d_det%d",det->detIndex,det->detNumber));
            drawing->Add(histOhmic[det->detIndex],"colz");
        }
    }

    TClonesArray *siChannelArray = nullptr;
    tree->SetBranchAddress("SiChannel",&siChannelArray);
    Long64_t numEvents = tree->GetEntries();
    if (maxEvents >= 0 && maxEvents < numEvents)
        numEvents = maxEvents;

    Long64_t numFilled[maxDetectorIndex] = {0};
    for (Long64_t event=0; event<numEvents; ++event) {
        tree->GetEntry(event);
        int multiplicity[kNumMult] = {0};
        for (auto iChannel=0; iChannel<siChannelArray->GetEntriesFast(); ++iChannel)
        {
            auto channel = (LKSiChannel*) siChannelArray->At(iChannel);
            int detIndex = channel->GetDetID();
            if (detIndex < 0 || detIndex >= maxDetectorIndex)
                continue;
            auto det = mapping.FindDetectorByIndex(detIndex);
            if (det == nullptr)
                continue;

            double energy1 = channel->GetEnergy1();
            if (energy1 > 0) {
                ++multiplicity[kAll];
                if      (det->ringType == "16")   ++multiplicity[kRing16];
                else if (det->ringType == "12dE") ++multiplicity[kRing12DE];
                else if (det->ringType == "12E")  ++multiplicity[kRing12E];
                else if (det->ringType == "Q")    ++multiplicity[kQQQ5];
            }

            if (channel->GetSide() == 1 && energy1 > 0 && histOhmic[detIndex] != nullptr)
                histOhmic[detIndex]->Fill(channel->GetStrip(),energy1);

            if (channel->GetSide() != 0 || histPosition[detIndex] == nullptr)
                continue;
            if (det->detType == "X6") {
                double energySum = channel->GetEnergySum();
                if (energySum > 0) {
                    histPosition[detIndex]->Fill(channel->GetEnergyPos(),energySum);
                    ++numFilled[detIndex];
                }
            }
            else if (det->detType == "QQQ5" && channel->GetEnergy1() > 0) {
                histPosition[detIndex]->Fill(channel->GetStrip(),channel->GetEnergy1());
                ++numFilled[detIndex];
            }
        }
        for (auto i=0; i<kNumMult; ++i)
            histMultiplicity[i]->Fill(multiplicity[i]);
        if (event%10000 == 0)
            cout << "event " << event << " / " << numEvents << endl;
    }

    TString outName = TString("data_reco/")+nameTop+".root";
    cout << outName << endl;
    auto out = new TFile(outName,"recreate");

    top1->Draw();
    drawings16   -> GetCanvas() -> Write();
    drawings12DE -> GetCanvas() -> Write();
    drawings12E  -> GetCanvas() -> Write();
    drawingsQQQ5 -> GetCanvas() -> Write();

    top2->Draw();
    drawingsOhmic16   -> GetCanvas() -> Write();
    drawingsOhmic12DE -> GetCanvas() -> Write();
    drawingsOhmic12E  -> GetCanvas() -> Write();
    drawingsOhmicQQQ5 -> GetCanvas() -> Write();

    top3->Draw();
    top3->GetCanvas()->Write();
    out->Close();

    cout << endl << "Entries by detector" << endl;
    for (auto iDet=0; iDet<mapping.GetNumDetectors(); ++iDet) {
        auto det = mapping.GetDetectorByVectorIndex(iDet);
        if (det != nullptr)
            printf("idx=%3d det#=%4d %-5s %-4s entries=%lld\n",
                   det->detIndex,det->detNumber,det->detType.Data(),
                   det->ringType.Data(),numFilled[det->detIndex]);
    }
}
