void draw_sihit(TString inputName="", TString outputName="")
{
    constexpr int maxDetectorIndex = 256;

    double energyMax = 30;
    Long64_t maxEvents = -1;

    auto par = new LKParameterContainer("get_common.mac");
    auto runID = par->GetParInt("LKRun/RunID");
    cout << "run " << runID << endl;

    auto tree = new TChain("event");
    if (!inputName.IsNull()) {
        cout << inputName << endl;
        tree->Add(inputName);
    }
    else if (runID == 8000) {
        for (int inputRun=83; inputRun<=90; ++inputRun)
            tree->Add(Form("data_reco/ko2520_%04d.sihit.root",inputRun));
    }
    else {
        inputName = Form("data_reco/ko2520_%04d.sihit.root",runID);
        cout << inputName << endl;
        tree->Add(inputName);
    }
    if (tree->GetNtrees() == 0) {
        cerr << "No sihit input file was found for run " << runID << endl;
        return;
    }
    if (tree->GetBranch("SiHit") == nullptr) {
        cerr << "SiHit branch is not present in the input file" << endl;
        return;
    }

    LKSiliconMapping mapping;
    if (!mapping.Load("mapping_ko2520")) {
        cerr << "Failed to load mapping_ko2520" << endl;
        return;
    }

    TString nameTop = Form("sihit_energy_position_%04d",runID);
    auto topPosition = new LKDrawingGroup(nameTop);
    auto position16   = topPosition->CreateGroup("X6_16_ring");
    auto position12DE = topPosition->CreateGroup("X6_12dE_ring");
    auto position12E  = topPosition->CreateGroup("X6_12E_ring");
    auto positionQQQ5 = topPosition->CreateGroup("QQQ5");
    auto position12Sum = topPosition->CreateGroup("X6_12_dE_plus_E");
    auto pid12TotalVsDE = topPosition->CreateGroup("X6_12_dE_vs_E_plus_dE");
    auto pid12DEVsE = topPosition->CreateGroup("X6_12_dE_vs_E");

    TString nameTop2 = Form("sihit_junction_ohmic_%04d",runID);
    auto topOhmic = new LKDrawingGroup(nameTop2);
    auto ohmic16   = topOhmic->CreateGroup("X6_16_ring");
    auto ohmic12DE = topOhmic->CreateGroup("X6_12dE_ring");
    auto ohmic12E  = topOhmic->CreateGroup("X6_12E_ring");
    auto ohmicQQQ5 = topOhmic->CreateGroup("QQQ5");

    TString nameTopOhmicPosition = Form("sihit_ohmic_energy_position_%04d",runID);
    auto topOhmicPosition = new LKDrawingGroup(nameTopOhmicPosition);
    auto ohmicPosition16   = topOhmicPosition->CreateGroup("X6_16_ring");
    auto ohmicPosition12DE = topOhmicPosition->CreateGroup("X6_12dE_ring");
    auto ohmicPosition12E  = topOhmicPosition->CreateGroup("X6_12E_ring");
    auto ohmicPositionQQQ5 = topOhmicPosition->CreateGroup("QQQ5");
    auto pid12OhmicDEVsE = topOhmicPosition->CreateGroup("X6_12_ohmic_dE_vs_E");
    auto ohmic12TotalVsPosition = topOhmicPosition->CreateGroup(
        "X6_12_ohmic_dE_plus_E_vs_position");

    TString nameTopRingSum = Form("sihit_ring_sum_energy_position_%04d",runID);
    auto topRingSum = new LKDrawingGroup(nameTopRingSum);

    TString nameTopOhmicRingSum = Form("sihit_ring_sum_ohmic_energy_position_%04d",runID);
    auto topOhmicRingSum = new LKDrawingGroup(nameTopOhmicRingSum);

    TString nameTop3 = Form("sihit_multiplicity_%04d",runID);
    auto topMultiplicity = new LKDrawingGroup(nameTop3);

    enum { kAll=0, kRing16=1, kRing12DE=2, kRing12E=3, kQQQ5=4, kNumMult=5 };
    enum { kSum16=0, kSum12DE=1, kSum12E=2, kSumQQQ5=3, kNumRingSums=4 };
    TString ringSumName[kNumRingSums] = {"X6_16", "X6_12dE", "X6_12E", "QQQ5"};
    TH2D *histRingSum[kNumRingSums] = {nullptr};
    TH2D *histOhmicRingSum[kNumRingSums] = {nullptr};
    for (int i=0; i<kNumRingSums; ++i) {
        bool isQQQ5 = (i == kSumQQQ5);
        histRingSum[i] = new TH2D(
            Form("hist_sihit_ring_sum_energy_position_%s",ringSumName[i].Data()),
            Form("%s detector-summed junction energy vs position;%s;junction energy",
                 ringSumName[i].Data(),isQQQ5 ? "junction strip" : "(E_{R}-E_{L})/(E_{R}+E_{L})"),
            isQQQ5 ? 32 : 200,isQQQ5 ? 0 : -1.2,isQQQ5 ? 32 : 1.2,
            300,0,energyMax);
        auto drawing = topRingSum->CreateDrawing(ringSumName[i]);
        drawing->Add(histRingSum[i],"colz");

        histOhmicRingSum[i] = new TH2D(
            Form("hist_sihit_ring_sum_ohmic_energy_position_%s",ringSumName[i].Data()),
            Form("%s detector-summed ohmic energy vs junction position;%s;ohmic energy",
                 ringSumName[i].Data(),isQQQ5 ? "junction strip" : "(E_{R}-E_{L})/(E_{R}+E_{L})"),
            isQQQ5 ? 32 : 200,isQQQ5 ? 0 : -1.2,isQQQ5 ? 32 : 1.2,
            300,0,energyMax);
        auto ohmicDrawing = topOhmicRingSum->CreateDrawing(ringSumName[i]);
        ohmicDrawing->Add(histOhmicRingSum[i],"colz");
    }
    auto histRingSumTotalEnergy = new TH2D(
        "hist_sihit_ring_sum_total_energy_position_X6_12",
        "12-ring detector-summed coincidence;E-detector (E_{R}-E_{L})/(E_{R}+E_{L});dE+E",
        200,-1.2,1.2,600,0,2*energyMax);
    auto ringSumTotalDrawing = topRingSum->CreateDrawing("X6_12_dE_plus_E");
    ringSumTotalDrawing->Add(histRingSumTotalEnergy,"colz");
    auto histRingSumDEVsE = new TH2D(
        "hist_sihit_ring_sum_de_vs_e_X6_12",
        "12-ring detector-summed junction coincidence;E;dE",
        300,0,energyMax,300,0,energyMax);
    auto ringSumDEVsEDrawing = topRingSum->CreateDrawing("X6_12_dE_vs_E");
    ringSumDEVsEDrawing->Add(histRingSumDEVsE,"colz");
    auto histRingSumDEVsTotal = new TH2D(
        "hist_sihit_ring_sum_de_vs_total_X6_12",
        "12-ring detector-summed junction coincidence;E+dE;dE",
        600,0,2*energyMax,300,0,energyMax);
    auto ringSumDEVsTotalDrawing = topRingSum->CreateDrawing("X6_12_dE_vs_E_plus_dE");
    ringSumDEVsTotalDrawing->Add(histRingSumDEVsTotal,"colz");
    auto histRingSumOhmicDEVsE = new TH2D(
        "hist_sihit_ring_sum_ohmic_de_vs_e_X6_12",
        "12-ring detector-summed ohmic coincidence;E_{ohmic};dE_{ohmic}",
        300,0,energyMax,300,0,energyMax);
    auto ringSumOhmicDEVsEDrawing = topOhmicRingSum->CreateDrawing("X6_12_ohmic_dE_vs_E");
    ringSumOhmicDEVsEDrawing->Add(histRingSumOhmicDEVsE,"colz");
    auto histRingSumOhmicTotalVsPosition = new TH2D(
        "hist_sihit_ring_sum_ohmic_total_vs_position_X6_12",
        "12-ring detector-summed ohmic coincidence;E-detector (E_{R}-E_{L})/(E_{R}+E_{L});dE_{ohmic}+E_{ohmic}",
        200,-1.2,1.2,600,0,2*energyMax);
    auto ringSumOhmicTotalVsPositionDrawing = topOhmicRingSum->CreateDrawing(
        "X6_12_ohmic_dE_plus_E_vs_position");
    ringSumOhmicTotalVsPositionDrawing->Add(histRingSumOhmicTotalVsPosition,"colz");

    TString multName[kNumMult] = {"all", "X6_16", "X6_12dE", "X6_12E", "QQQ5"};
    TH1D *histMultiplicity[kNumMult] = {nullptr};
    for (int i=0; i<kNumMult; ++i) {
        histMultiplicity[i] = new TH1D(
            Form("hist_sihit_multiplicity_%s",multName[i].Data()),
            Form("%s SiHit multiplicity;number of detector hits;events",multName[i].Data()),
            50,0,50);
        auto drawing = topMultiplicity->CreateDrawing(multName[i]);
        drawing->Add(histMultiplicity[i],"hist");
    }

    TH2D *histPosition[maxDetectorIndex] = {nullptr};
    TH2D *histOhmic[maxDetectorIndex] = {nullptr};
    TH2D *histOhmicPosition[maxDetectorIndex] = {nullptr};
    TH2D *histTotalEnergy[maxDetectorIndex] = {nullptr};
    TH2D *histTotalVsDE[maxDetectorIndex] = {nullptr};
    TH2D *histDEVsE[maxDetectorIndex] = {nullptr};
    TH2D *histOhmicDEVsE[maxDetectorIndex] = {nullptr};
    TH2D *histOhmicTotalVsPosition[maxDetectorIndex] = {nullptr};
    int dEDetectorByE[maxDetectorIndex];
    Long64_t numFilled[maxDetectorIndex] = {0};
    Long64_t numCoincidence[maxDetectorIndex] = {0};
    Long64_t numOhmicCoincidence[maxDetectorIndex] = {0};
    for (int id=0; id<maxDetectorIndex; ++id)
        dEDetectorByE[id] = -1;

    for (int iDet=0; iDet<mapping.GetNumDetectors(); ++iDet) {
        auto det = mapping.GetDetectorByVectorIndex(iDet);
        if (det == nullptr || det->detIndex < 0 || det->detIndex >= maxDetectorIndex)
            continue;

        int id = det->detIndex;
        if (det->detType == "QQQ5") {
            histPosition[id] = new TH2D(
                Form("hist_sihit_energy_position_det_%03d",id),
                Form("QQQ5 det#%d (idx=%d);junction strip;restored energy",
                     det->detNumber,id),
                32,0,32,300,0,energyMax);
        }
        else {
            histPosition[id] = new TH2D(
                Form("hist_sihit_energy_position_det_%03d",id),
                Form("X6 det#%d (idx=%d, %s);(E_{R}-E_{L})/(E_{R}+E_{L});E_{R}+E_{L}",
                     det->detNumber,id,det->ringType.Data()),
                200,-1.2,1.2,300,0,energyMax);
        }
        histOhmic[id] = new TH2D(
            Form("hist_sihit_junction_ohmic_det_%03d",id),
            Form("%s det#%d (idx=%d, %s);junction energy;ohmic energy",
                 det->detType.Data(),det->detNumber,id,det->ringType.Data()),
            300,0,energyMax,300,0,energyMax);
        bool isQQQ5 = (det->detType == "QQQ5");
        histOhmicPosition[id] = new TH2D(
            Form("hist_sihit_ohmic_energy_position_det_%03d",id),
            Form("%s det#%d (idx=%d, %s);%s;ohmic energy",
                 det->detType.Data(),det->detNumber,id,det->ringType.Data(),
                 isQQQ5 ? "junction strip" : "(E_{R}-E_{L})/(E_{R}+E_{L})"),
            isQQQ5 ? 32 : 200,isQQQ5 ? 0 : -1.2,isQQQ5 ? 32 : 1.2,
            300,0,energyMax);

        LKDrawingGroup *positionGroup = nullptr;
        LKDrawingGroup *ohmicGroup = nullptr;
        LKDrawingGroup *ohmicPositionGroup = nullptr;
        if (det->ringType == "16") {
            positionGroup = position16;
            ohmicGroup = ohmic16;
            ohmicPositionGroup = ohmicPosition16;
        }
        else if (det->ringType == "12dE") {
            positionGroup = position12DE;
            ohmicGroup = ohmic12DE;
            ohmicPositionGroup = ohmicPosition12DE;
        }
        else if (det->ringType == "12E") {
            positionGroup = position12E;
            ohmicGroup = ohmic12E;
            ohmicPositionGroup = ohmicPosition12E;
        }
        else if (det->ringType == "Q") {
            positionGroup = positionQQQ5;
            ohmicGroup = ohmicQQQ5;
            ohmicPositionGroup = ohmicPositionQQQ5;
        }

        TString drawingName = Form("idx%d_det%d",id,det->detNumber);
        if (positionGroup != nullptr) {
            auto drawing = positionGroup->CreateDrawing(drawingName);
            drawing->Add(histPosition[id],"colz");
        }
        if (ohmicGroup != nullptr) {
            auto drawing = ohmicGroup->CreateDrawing(drawingName);
            drawing->Add(histOhmic[id],"colz");
        }
        if (ohmicPositionGroup != nullptr) {
            auto drawing = ohmicPositionGroup->CreateDrawing(drawingName);
            drawing->Add(histOhmicPosition[id],"colz");
        }

        if (det->ringType == "12E") {
            const LKSiliconMapping::DetectorInfo *dEDet = nullptr;
            for (int jDet=0; jDet<mapping.GetNumDetectors(); ++jDet) {
                auto candidate = mapping.GetDetectorByVectorIndex(jDet);
                if (candidate != nullptr && candidate->ringType == "12dE" &&
                    candidate->phiNumber == det->phiNumber) {
                    dEDet = candidate;
                    break;
                }
            }
            if (dEDet != nullptr) {
                dEDetectorByE[id] = dEDet->detIndex;
                histTotalEnergy[id] = new TH2D(
                    Form("hist_sihit_total_energy_pair_%03d_%03d",dEDet->detIndex,id),
                    Form("12-ring coincidence dE det#%d (idx=%d) + E det#%d (idx=%d);E relative z;dE+E",
                         dEDet->detNumber,dEDet->detIndex,det->detNumber,id),
                    200,-1.2,1.2,600,0,2*energyMax);
                auto drawing = position12Sum->CreateDrawing(
                    Form("phi%d_dEidx%d_Eidx%d",det->phiNumber,dEDet->detIndex,id));
                drawing->Add(histTotalEnergy[id],"colz");

                histTotalVsDE[id] = new TH2D(
                    Form("hist_sihit_total_vs_de_pair_%03d_%03d",dEDet->detIndex,id),
                    Form("12-ring coincidence dE det#%d (idx=%d) + E det#%d (idx=%d);E+dE;dE",
                         dEDet->detNumber,dEDet->detIndex,det->detNumber,id),
                    600,0,2*energyMax,300,0,energyMax);
                auto pidDrawing = pid12TotalVsDE->CreateDrawing(
                    Form("phi%d_dEidx%d_Eidx%d",det->phiNumber,dEDet->detIndex,id));
                pidDrawing->Add(histTotalVsDE[id],"colz");

                histDEVsE[id] = new TH2D(
                    Form("hist_sihit_de_vs_e_pair_%03d_%03d",dEDet->detIndex,id),
                    Form("12-ring junction coincidence dE det#%d + E det#%d;E;dE",
                         dEDet->detNumber,det->detNumber),
                    300,0,energyMax,300,0,energyMax);
                auto deVsEDrawing = pid12DEVsE->CreateDrawing(
                    Form("phi%d_dEidx%d_Eidx%d",det->phiNumber,dEDet->detIndex,id));
                deVsEDrawing->Add(histDEVsE[id],"colz");

                histOhmicDEVsE[id] = new TH2D(
                    Form("hist_sihit_ohmic_de_vs_e_pair_%03d_%03d",dEDet->detIndex,id),
                    Form("12-ring ohmic coincidence dE det#%d + E det#%d;E_{ohmic};dE_{ohmic}",
                         dEDet->detNumber,det->detNumber),
                    300,0,energyMax,300,0,energyMax);
                auto ohmicDEVsEDrawing = pid12OhmicDEVsE->CreateDrawing(
                    Form("phi%d_dEidx%d_Eidx%d",det->phiNumber,dEDet->detIndex,id));
                ohmicDEVsEDrawing->Add(histOhmicDEVsE[id],"colz");

                histOhmicTotalVsPosition[id] = new TH2D(
                    Form("hist_sihit_ohmic_total_vs_position_pair_%03d_%03d",dEDet->detIndex,id),
                    Form("12-ring ohmic coincidence dE det#%d + E det#%d;E-detector (E_{R}-E_{L})/(E_{R}+E_{L});dE_{ohmic}+E_{ohmic}",
                         dEDet->detNumber,det->detNumber),
                    200,-1.2,1.2,600,0,2*energyMax);
                auto ohmicTotalVsPositionDrawing = ohmic12TotalVsPosition->CreateDrawing(
                    Form("phi%d_dEidx%d_Eidx%d",det->phiNumber,dEDet->detIndex,id));
                ohmicTotalVsPositionDrawing->Add(histOhmicTotalVsPosition[id],"colz");
            }
        }
    }

    TClonesArray *siHitArray = nullptr;
    tree->SetBranchAddress("SiHit",&siHitArray);
    Long64_t numEvents = tree->GetEntries();
    if (maxEvents >= 0 && maxEvents < numEvents)
        numEvents = maxEvents;

    // A paired SiHit contains both detectors: fDetID is E and fdEDetID is dE.
    // Fill both detector histograms so that the 12dE ring is not omitted.
    for (Long64_t event=0; event<numEvents; ++event) {
        tree->GetEntry(event);
        int multiplicity[kNumMult] = {0};
        double eventEnergy[maxDetectorIndex] = {0};
        double eventPosition[maxDetectorIndex];
        double eventOhmicEnergy[maxDetectorIndex] = {0};
        for (int id=0; id<maxDetectorIndex; ++id)
            eventPosition[id] = -999;

        for (int iHit=0; iHit<siHitArray->GetEntriesFast(); ++iHit) {
            auto hit = (SKSiHit*) siHitArray->At(iHit);

            int detIDs[2] = {hit->GetDetID(), hit->GetdEDetID()};
            double energies[2] = {hit->GetE(), hit->GetdE()};
            double ohmicEnergies[2] = {hit->GetEnergyOhmic(), hit->GetdEOhmic()};
            double relativeZ[2] = {hit->GetRelativeZ(), hit->GetRelativeZdE()};

            for (int component=0; component<2; ++component) {
                int id = detIDs[component];
                double energy = energies[component];
                if (id < 0 || id >= maxDetectorIndex)
                    continue;

                auto det = mapping.FindDetectorByIndex(id);
                if (det == nullptr || histPosition[id] == nullptr)
                    continue;

                double ohmicEnergy = ohmicEnergies[component];
                if (ohmicEnergy > eventOhmicEnergy[id])
                    eventOhmicEnergy[id] = ohmicEnergy;

                // X6 position and energy must come from a reconstructed
                // left-right junction pair. Standalone/ohmic hits are excluded.
                if (det->detType == "X6") {
                    if (component != 0 ||
                        hit->GetEnergyLeft() <= 0 || hit->GetEnergyRight() <= 0)
                        continue;
                    double energyLeft = hit->GetEnergyLeft();
                    double energyRight = hit->GetEnergyRight();
                    energy = energyLeft+energyRight;
                    relativeZ[component] =
                        (energyRight-energyLeft)/(energyRight+energyLeft);
                }
                // In current sihit files an ohmic QQQ5 channel is a separate
                // SiHit. Do not draw it as a junction-strip hit.
                else if (ohmicEnergy > 0) {
                    continue;
                }

                if (energy <= 0)
                    continue;

                ++multiplicity[kAll];
                if      (det->ringType == "16")   ++multiplicity[kRing16];
                else if (det->ringType == "12dE") ++multiplicity[kRing12DE];
                else if (det->ringType == "12E")  ++multiplicity[kRing12E];
                else if (det->ringType == "Q")    ++multiplicity[kQQQ5];

                double position = (det->detType == "QQQ5")
                                ? hit->GetJunctionStrip() : relativeZ[component];
                histPosition[id]->Fill(position,energy);

                int ringSumIndex = -1;
                if      (det->ringType == "16")   ringSumIndex = kSum16;
                else if (det->ringType == "12dE") ringSumIndex = kSum12DE;
                else if (det->ringType == "12E")  ringSumIndex = kSum12E;
                else if (det->ringType == "Q")    ringSumIndex = kSumQQQ5;
                if (ringSumIndex >= 0) {
                    histRingSum[ringSumIndex]->Fill(position,energy);
                }

                if (energy > eventEnergy[id]) {
                    eventEnergy[id] = energy;
                    eventPosition[id] = position;
                }
                ++numFilled[id];
            }
        }

        for (int id=0; id<maxDetectorIndex; ++id) {
            if (eventEnergy[id] > 0 && eventOhmicEnergy[id] > 0 &&
                eventPosition[id] > -900 && histOhmic[id] != nullptr) {
                histOhmic[id]->Fill(eventEnergy[id],eventOhmicEnergy[id]);
                histOhmicPosition[id]->Fill(eventPosition[id],eventOhmicEnergy[id]);

                auto det = mapping.FindDetectorByIndex(id);
                int ringSumIndex = -1;
                if (det != nullptr) {
                    if      (det->ringType == "16")   ringSumIndex = kSum16;
                    else if (det->ringType == "12dE") ringSumIndex = kSum12DE;
                    else if (det->ringType == "12E")  ringSumIndex = kSum12E;
                    else if (det->ringType == "Q")    ringSumIndex = kSumQQQ5;
                }
                if (ringSumIndex >= 0)
                    histOhmicRingSum[ringSumIndex]->Fill(
                        eventPosition[id],eventOhmicEnergy[id]);
            }
        }

        // The current sihit files store each detector as a separate SiHit.
        // Match the 12dE and 12E detectors with the same phi number per event.
        for (int eID=0; eID<maxDetectorIndex; ++eID) {
            int dEID = dEDetectorByE[eID];
            if (histTotalEnergy[eID] == nullptr || dEID < 0)
                continue;
            if (eventEnergy[eID] > 0 && eventEnergy[dEID] > 0) {
                double totalEnergy = eventEnergy[dEID]+eventEnergy[eID];
                histTotalEnergy[eID]->Fill(
                    eventPosition[eID],totalEnergy);
                histTotalVsDE[eID]->Fill(totalEnergy,eventEnergy[dEID]);
                histDEVsE[eID]->Fill(eventEnergy[eID],eventEnergy[dEID]);
                histRingSumTotalEnergy->Fill(eventPosition[eID],totalEnergy);
                histRingSumDEVsE->Fill(eventEnergy[eID],eventEnergy[dEID]);
                histRingSumDEVsTotal->Fill(totalEnergy,eventEnergy[dEID]);
                ++numCoincidence[eID];
            }
            if (eventOhmicEnergy[eID] > 0 && eventOhmicEnergy[dEID] > 0) {
                histOhmicDEVsE[eID]->Fill(
                    eventOhmicEnergy[eID],eventOhmicEnergy[dEID]);
                histRingSumOhmicDEVsE->Fill(
                    eventOhmicEnergy[eID],eventOhmicEnergy[dEID]);
                if (eventPosition[eID] > -900) {
                    double ohmicTotal = eventOhmicEnergy[eID]+eventOhmicEnergy[dEID];
                    histOhmicTotalVsPosition[eID]->Fill(
                        eventPosition[eID],ohmicTotal);
                    histRingSumOhmicTotalVsPosition->Fill(
                        eventPosition[eID],ohmicTotal);
                }
                ++numOhmicCoincidence[eID];
            }
        }

        for (int i=0; i<kNumMult; ++i)
            histMultiplicity[i]->Fill(multiplicity[i]);
        if (event%10000 == 0)
            cout << "event " << event << " / " << numEvents << endl;
    }

    if (outputName.IsNull())
        outputName = Form("data_reco/draw_sihit_%04d.root",runID);
    cout << outputName << endl;
    auto output = new TFile(outputName,"recreate");

    topPosition->Draw();
    position16->GetCanvas()->Write();
    position12DE->GetCanvas()->Write();
    position12E->GetCanvas()->Write();
    positionQQQ5->GetCanvas()->Write();
    position12Sum->GetCanvas()->Write();
    pid12TotalVsDE->GetCanvas()->Write();
    pid12DEVsE->GetCanvas()->Write();

    topOhmic->Draw();
    ohmic16->GetCanvas()->Write();
    ohmic12DE->GetCanvas()->Write();
    ohmic12E->GetCanvas()->Write();
    ohmicQQQ5->GetCanvas()->Write();

    topOhmicPosition->Draw();
    ohmicPosition16->GetCanvas()->Write();
    ohmicPosition12DE->GetCanvas()->Write();
    ohmicPosition12E->GetCanvas()->Write();
    ohmicPositionQQQ5->GetCanvas()->Write();
    pid12OhmicDEVsE->GetCanvas()->Write();
    ohmic12TotalVsPosition->GetCanvas()->Write();

    topRingSum->Draw();
    topRingSum->GetCanvas()->Write();

    topOhmicRingSum->Draw();
    topOhmicRingSum->GetCanvas()->Write();

    topMultiplicity->Draw();
    topMultiplicity->GetCanvas()->Write();

    auto histogramDirectory = output->mkdir("histograms");
    histogramDirectory->cd();
    for (int id=0; id<maxDetectorIndex; ++id) {
        if (histPosition[id] != nullptr) histPosition[id]->Write();
        if (histOhmic[id] != nullptr) histOhmic[id]->Write();
        if (histOhmicPosition[id] != nullptr) histOhmicPosition[id]->Write();
        if (histTotalEnergy[id] != nullptr) histTotalEnergy[id]->Write();
        if (histTotalVsDE[id] != nullptr) histTotalVsDE[id]->Write();
        if (histDEVsE[id] != nullptr) histDEVsE[id]->Write();
        if (histOhmicDEVsE[id] != nullptr) histOhmicDEVsE[id]->Write();
        if (histOhmicTotalVsPosition[id] != nullptr) histOhmicTotalVsPosition[id]->Write();
    }
    for (int i=0; i<kNumMult; ++i)
        histMultiplicity[i]->Write();
    for (int i=0; i<kNumRingSums; ++i) {
        histRingSum[i]->Write();
        histOhmicRingSum[i]->Write();
    }
    histRingSumTotalEnergy->Write();
    histRingSumDEVsE->Write();
    histRingSumDEVsTotal->Write();
    histRingSumOhmicDEVsE->Write();
    histRingSumOhmicTotalVsPosition->Write();
    output->Close();

    cout << endl << "SiHit entries by detector" << endl;
    for (int iDet=0; iDet<mapping.GetNumDetectors(); ++iDet) {
        auto det = mapping.GetDetectorByVectorIndex(iDet);
        if (det != nullptr && det->detIndex >= 0 && det->detIndex < maxDetectorIndex)
            printf("idx=%3d det#=%4d %-5s %-4s entries=%lld\n",
                   det->detIndex,det->detNumber,det->detType.Data(),
                   det->ringType.Data(),numFilled[det->detIndex]);
    }

    cout << endl << "12-ring dE-E coincidences" << endl;
    for (int eID=0; eID<maxDetectorIndex; ++eID) {
        if (histTotalEnergy[eID] == nullptr)
            continue;
        auto eDet = mapping.FindDetectorByIndex(eID);
        auto dEDet = mapping.FindDetectorByIndex(dEDetectorByE[eID]);
        if (eDet != nullptr && dEDet != nullptr)
            printf("phi=%2d dE idx=%3d det#=%4d + E idx=%3d det#=%4d junction=%lld ohmic=%lld\n",
                   eDet->phiNumber,dEDet->detIndex,dEDet->detNumber,
                   eDet->detIndex,eDet->detNumber,
                   numCoincidence[eID],numOhmicCoincidence[eID]);
    }
}
