TString colorHex[] = {"#8000ff","#7019ff","#5e35fe","#4c50fc","#396bf9","#2884f6","#1898f2","#06aeed","#0dc2e8","#1fd3e1","#30e1da","#42edd3","#52f5cb","#64fbc3","#76ffb9","#88ffaf","#9bfba5","#acf59a","#bced8f","#cee184","#e0d377","#f2c26b","#ffae5e","#ff9850","#ff8444","#ff6b37","#ff5029","#ff351b","#ff190d","#ff0000"};

void draw_event_count()
{
    bool useCutGRange = true;
    double y2CutGRange = 0.04;
    double y2ZoomIn = 0.005;
    double ox = -0.5;

    const int numPairs = 28;
    int numEvents[200] = {0};
    int count[200][numPairs] = {0};
    vector<int> runIDExp;
    int runIDMin=INT_MAX, runIDMax=-INT_MAX;
    for (auto runID : runIDExp)
    {
        if (runID<runIDMin) runIDMin = runID;
        if (runID>runIDMax) runIDMax = runID;
    }

    LKParameterContainer parc("config_short_cut.mac");
    vector<int> runIDsE4 = parc.GetParIntRange("e4/RunIDRange");
    vector<int> runIDsE5 = parc.GetParIntRange("e5/RunIDRange");
    vector<int> runIDsE8 = parc.GetParIntRange("e8/RunIDRange");
    if (runIDExp.size()==0) {
        for (auto runID : runIDsE4) runIDExp.push_back(runID);
        for (auto runID : runIDsE5) runIDExp.push_back(runID);
        for (auto runID : runIDsE8) runIDExp.push_back(runID);
    }

    //LKBinning bnnRunID("Event numbers","",120,0,120);
    LKBinning bnnRunID("Event numbers","",120,15,135);
    LKBinning bnnTheta(400,-1,91);
    LKBinning bnnEnergy(400,-1,60);
    LKBinning rngTheta(10,90);
    LKBinning rngEnergy(1.2,50);
    LKBinning bnnRatio(100,0,0.07);
    if (useCutGRange) bnnRatio.SetX2(y2CutGRange);
    bnnRunID.Print();

    TCutG *cutKE4 = (TCutG*) (new TFile("data_x/rough_kinematic_cut_e4.root")) -> Get("CUTG");
    TCutG *cutKE5 = (TCutG*) (new TFile("data_x/rough_kinematic_cut_e5.root")) -> Get("CUTG");
    TCutG *cutKE8 = (TCutG*) (new TFile("data_x/rough_kinematic_cut_e8.root")) -> Get("CUTG");

    auto histETInn = (bnnTheta*bnnEnergy).NewH2("histETInn","Inside the selection");
    auto histETOut = (bnnTheta*bnnEnergy).NewH2("histETOut","Outside the selection");
    auto histRatio = (bnnRunID*bnnRatio).NewH2("histRatio",";runID;#hits/#events"); histRatio -> SetStats(0);
    auto histCount4 = bnnRunID.NewH1("histCount4",";run;count"); histCount4 -> SetStats(0);
    auto histCount5 = bnnRunID.NewH1("histCount5",";run;count"); histCount5 -> SetStats(0);
    auto histCount8 = bnnRunID.NewH1("histCount8",";run;count"); histCount8 -> SetStats(0);
    histRatio -> GetXaxis() -> SetNdivisions(230);
    histCount4 -> GetXaxis() -> SetNdivisions(230);
    histCount5 -> GetXaxis() -> SetNdivisions(230);
    histCount8 -> GetXaxis() -> SetNdivisions(230);
    histCount4 -> SetFillColor(18        ); histCount4 -> SetLineColor(kGray+1);
    histCount5 -> SetFillColor(kCyan-10  ); histCount5 -> SetLineColor(kGray+1);
    histCount8 -> SetFillColor(kYellow-10); histCount8 -> SetLineColor(kGray+1);

    auto top = new LKDrawingGroup();
    auto groupSummary = top -> CreateGroup("summary");
    auto drawETInn = groupSummary -> CreateDrawing("ETInn"); drawETInn -> SetLogz();
    auto drawETout = groupSummary -> CreateDrawing("ETOut"); drawETout -> SetLogz();
    auto drawCount = groupSummary -> CreateDrawing("Count");
    drawETInn -> Add(histETInn);
    drawETout -> Add(histETOut);
    drawCount -> Add(histCount4,"hist","4.41 MeV");
    drawCount -> Add(histCount5,"hist","5.84 MeV");
    drawCount -> Add(histCount8,"hist","8.13 MeV");
    drawCount -> SetCreateLegend();
    //auto drawR1 = groupSummary -> CreateDrawing("drawR1");
    auto drawR1 = top -> CreateGroup("ratio_all") -> CreateDrawing("drawR1");
    drawR1 -> Add(histRatio);
    //auto group1 = top -> CreateGroup("ratio");
    //auto group2 = group1;//top -> CreateGroup("ratio2");

    for (auto idxEnergy : {4,5,8})
    {
        vector<int> runIDArray;
        if (idxEnergy==4) runIDArray = runIDsE4;
        if (idxEnergy==5) runIDArray = runIDsE5;
        if (idxEnergy==8) runIDArray = runIDsE8;
        for (auto runID : runIDArray)
        {
            auto graphShadeRun = new TGraph();
            graphShadeRun -> SetFillColor(18);
            //if (idxEnergy==4) graphShadeRun -> SetFillColor(kYellow-10);
            //if (idxEnergy==5) graphShadeRun -> SetFillColor(kCyan-10);
            //if (idxEnergy==8) graphShadeRun -> SetFillColor(kGreen-10);
            graphShadeRun -> SetFillStyle(3001);
            graphShadeRun -> SetPoint(0,runID  +ox,bnnRatio.x1());
            graphShadeRun -> SetPoint(1,runID  +ox,bnnRatio.x2());
            graphShadeRun -> SetPoint(2,runID+1+ox,bnnRatio.x2());
            graphShadeRun -> SetPoint(3,runID+1+ox,bnnRatio.x1());
            graphShadeRun -> SetPoint(4,runID  +ox,bnnRatio.x1());
            drawR1 -> Add(graphShadeRun,"samef",".");
        }

        auto id1 = runIDArray[0];
        auto graphBoundary1 = new TGraph();
        graphBoundary1 -> SetPoint(0,id1+ox,bnnRatio.x1());
        graphBoundary1 -> SetPoint(1,id1+ox,bnnRatio.x2());
        drawR1 -> Add(graphBoundary1,"samel",".");

        auto id2 = runIDArray.back();
        auto graphBoundary2 = new TGraph();
        graphBoundary2 -> SetPoint(0,id2+1+ox,bnnRatio.x1());
        graphBoundary2 -> SetPoint(1,id2+1+ox,bnnRatio.x2());
        drawR1 -> Add(graphBoundary2,"samel",".");
    }

    //auto lg = new TLegend(0.1,0.1,0.2,0.9);
    auto lg = new TLegend(0.76,0.12,0.89,0.88);
    lg -> SetHeader("Pairs");
    lg -> SetMargin(0.5);
    lg -> SetNColumns(2);
    lg -> SetFillColor(19);
    //lg -> SetFillStyle(0);
    //lg -> SetBorderSize(0);
    TGraph* graphRatio[numPairs] = {0};
    for (auto pairID=0; pairID<numPairs; ++pairID) {
        graphRatio[pairID] = new TGraph();
        if (pairID==22||pairID==23) continue;
        auto colorID = pairID;
        auto mStyle = 4;
        if (colorID>=0&&colorID<4) { colorID = 0; mStyle = pairID+24; }
        else if (colorID>=4&&colorID<12) { colorID = 8; mStyle = pairID-4+24; }
        else mStyle = pairID-12+24;
        graphRatio[pairID] -> SetFillColor(TColor::GetColor(colorHex[colorID]));
        graphRatio[pairID] -> SetLineColor(TColor::GetColor(colorHex[colorID]));
        graphRatio[pairID] -> SetMarkerStyle(mStyle);
        lk_debug << pairID << " " << mStyle << endl;
        //graphRatio[pairID] -> SetMarkerSize(0.6);
        drawR1 -> Add(graphRatio[pairID],"samelp",Form("%d",pairID));
        auto graphRatio2 = (TGraph*) graphRatio[pairID] -> Clone();
        graphRatio2 -> SetMarkerSize(2);
        lg -> AddEntry(graphRatio2,Form("%d",pairID),"fp");
    }
    drawR1 -> Add(lg);
    //drawR1 -> SetCreateLegend(1);

    for (auto idxEnergy : {4,5,8})
    {
        LKDrawingGroup *group = top -> CreateGroup(Form("ratio_%d",idxEnergy));
        group -> SetCanvasDivision(1,2);
        for (auto zz : {1,2})
        {
            //for (auto idxEnergy : {4,5,8})
            vector<int> runIDArray;
            if (idxEnergy==4) runIDArray = runIDsE4;
            if (idxEnergy==5) runIDArray = runIDsE5;
            if (idxEnergy==8) runIDArray = runIDsE8;
            auto drawRE1 = group -> CreateDrawing(Form("drawR%d_e%d",zz,idxEnergy)); drawR1 -> CopyTo(drawRE1);
            double x1 = runIDArray[0];
            double x2 = runIDArray.back();
            double dx = (x2-x1);
            x1 = x1 - 0.05*dx;
            x2 = x2 + 0.28*dx;
            if (zz==1) drawRE1 -> SetRangeUserX(x1,x2);
            if (zz==2) drawRE1 -> SetRangeUserX(x1,x2);
            if (zz==2) drawRE1 -> SetRangeUserY(0,y2ZoomIn);

            auto energy_title = new TLatex();
            energy_title -> SetTextColor(kGray+3);
            energy_title -> SetTextAlign(11);
            energy_title -> SetTextSize(0.05);
            energy_title -> SetTextFont(132);
            TString text = Form("4.41 MeV");
            if (idxEnergy==5) text = Form("5.84 MeV");
            if (idxEnergy==8) text = Form("8.13 MeV");
            energy_title -> SetText(0.5,1.-0.05,text);
            energy_title -> SetNDC();
            energy_title -> SetTextAlign(22);
            if (zz==1) drawRE1 -> Add(energy_title);
            if (zz==1) drawRE1 -> SetCloneReplaceMainHist(Form("+_%d",idxEnergy));
            if (zz==2) drawRE1 -> Add(energy_title);
            if (zz==2) drawRE1 -> SetCloneReplaceMainHist(Form("+_zi_%d",idxEnergy));
        }
    }

    int index = 0;
    int idxEnergyPrev = 0;
    for (auto runID : runIDExp)
    {
        //cout << runID;

        TString inputName = Form("data_reco/stark_%04d.reco7.root",runID);
        auto file = new TFile(inputName);
        auto tree = (TTree*) file -> Get("event");
        TClonesArray* recoHitArray = nullptr;
        tree -> SetBranchAddress("RecoHit",&recoHitArray);
        int numEvents = tree -> GetEntries();

        double beamEnergy = 0;
        int idxEnergy = 4;
        if (LKMisc::ValueIsInArray(runID,runIDsE4)) { idxEnergy = 4; beamEnergy = 4.41;}
        if (LKMisc::ValueIsInArray(runID,runIDsE5)) { idxEnergy = 5; beamEnergy = 5.84;}
        if (LKMisc::ValueIsInArray(runID,runIDsE8)) { idxEnergy = 8; beamEnergy = 8.13;}
        if (idxEnergyPrev==idxEnergy)
            index++;
        else { index = 0; idxEnergyPrev = idxEnergy; }

        //cout << " " << numEvents << endl;

        int countAll = 0;
        int countPair[numPairs] = {0};
        for (auto iEvent=0; iEvent<numEvents; ++iEvent)
        {
            tree -> GetEntry(iEvent);
            auto numHits = recoHitArray -> GetEntries();

            for (auto iHit=0; iHit<numHits; ++iHit)
            {
                auto hit = (EKRecoHit*) recoHitArray -> At(iHit);

                auto pairID     = hit -> GetPairID();
                auto theta_lab  = hit -> GetTheta();
                auto theta_com  = 180. - theta_lab * 2.;
                auto energy     = hit -> GetKeyEnergy();
                auto qvalue     = hit -> GetQValue();

                bool good = false;
                if (useCutGRange) {
                    if (idxEnergy==4 && cutKE4->IsInside(theta_lab,energy)) good = true;
                    if (idxEnergy==5 && cutKE5->IsInside(theta_lab,energy)) good = true;
                    if (idxEnergy==8 && cutKE8->IsInside(theta_lab,energy)) good = true;
                }
                else
                    if (rngEnergy.IsInside(energy) && rngTheta.IsInside(theta_lab))
                        good = true;

                if (good)
                {
                    countAll++;
                    countPair[pairID]++;
                    count[runID][pairID]++;
                    histETInn -> Fill(theta_lab,energy);
                }
                else
                    histETOut -> Fill(theta_lab,energy);
            }
        }
        if (idxEnergy==4) histCount4 -> Fill(runID, countAll);
        if (idxEnergy==5) histCount5 -> Fill(runID, countAll);
        if (idxEnergy==8) histCount8 -> Fill(runID, countAll);
        cout << index << " & " << beamEnergy << " & " << runID << " & " << numEvents << " & " << int(countAll) << " & " << Form("%.3f",double(countAll)/numEvents) << " \\\\\\hline" << endl;
        for (auto pairID=0; pairID<numPairs; ++pairID) {
            graphRatio[pairID] -> SetPoint(graphRatio[pairID]->GetN(),runID+0.5+ox,double(countPair[pairID])/numEvents);
        }
        file -> Close();
    }

    top -> Draw("v");
    //top -> Draw();
}
