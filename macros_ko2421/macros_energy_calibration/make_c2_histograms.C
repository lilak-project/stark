#include "si_calibration.h"

TString make_c2_histograms(bool runViewer=true, int calibration=2, int run=-1, int calRun=-1, bool checkEP=false, TString c3Name="")
{
    MakeRun(run, ((calRun==199&&run==168)?1:-1));
    TString name = Form("run_%04d_d%d",fRun,calibration);
    if (calRun!=run)
        name = Form("check_%d_with_%d_ballistic_parameters",calRun,run);
    auto top = new LKDrawingGroup(name);
    auto groupJ = top -> CreateGroup("HPJ",!checkEP);
    auto groupO = top -> CreateGroup("HPO",!checkEP);
    auto groupFit = top -> CreateGroup("Fit",(calRun<0));
    auto groupA = top -> CreateGroup("HPA",!checkEP);
    auto group1 = top -> CreateGroup("EnergySum",false);
    auto group2 = top -> CreateGroup("LeftRight",false);
    auto group3 = top -> CreateGroup("EnergyPosition",(calibration==3||checkEP));
    auto group4 = top -> CreateGroup("EPO",(!checkEP&&calibration==2));
    auto group0 = top -> CreateGroup("Energy",(!checkEP));


    /////////////////////////////////////////////////////////////////////
    // 1) Get calibration parameter if calibration is 1 (slope correction
    /////////////////////////////////////////////////////////////////////
    TString calName = Form("c%d_",calibration);
    GetC0Parameters(calRun);
    GetC1Parameters(calRun);
    GetC2Parameters(calRun);
    if (c3Name.IsNull()) GetC3Parameters(calRun);
    else                 GetC3Parameters(c3Name);
    //GetCPParameters(calRun);

    /////////////////////////////////////////////////////////////////////
    // 2) Create histograms
    /////////////////////////////////////////////////////////////////////
    for (auto dss : fStripArrayS)
    {
        fHistEnergy[dss.det][dss.side][dss.strip] = MakeHist1(calName+"energy","",dss.det,dss.side,dss.strip,-1,fNBinE,fBinE1,fBinE2);
    }
    for (auto dss : fStripArrayR)
    {
        fHistEnergySum     [dss.det][dss.side][dss.strip] = MakeHist1(calName+"esum",             "",dss.det,dss.side,dss.strip,-1,fNBinE,fBinE1,fBinE2);
        fHistLeftRight     [dss.det][dss.side][dss.strip] = MakeHist2(calName+"left",calName+"right",dss.det,dss.side,dss.strip,-1,fNBinE,fBinE1,fBinE2,fNBinE,fBinE1,fBinE2);
        fHistEnergyPosition[dss.det][dss.side][dss.strip] = MakeHist2(calName+"rpos",calName+"esum", dss.det,dss.side,dss.strip,-1,fNBinX,fBinX1,fBinX2,fNBinE,fBinE1,fBinE2);
        for (auto gate=0; gate<fNumGates; ++gate) {
            for (auto os=0; os<fNumOStrips; ++os) {
                TString calName2 = calName + "os" + os + "_";
                auto hist = MakeHist2(calName2+"rpos",calName2+"esum", dss.det,dss.side,dss.strip,gate,fNBinX,fBinX1,fBinX2,fNBinE,fBinE1,fBinE2);
                fHistEnergyPositionGateOhmic[dss.det][dss.side][dss.strip][gate][os] = hist;
            }
        }
    }
    fHistEnergyPositionAll = MakeHist2(calName+"rpos",calName+"esum", -1,-1,-1,-1,2*fNBinX,fBinX1,fBinX2,2*fNBinE,fBinE1,fBinE2);
    fHistEnergyDetector[0] = MakeHitPatternHist(0,calName,2*fNBinEFix,fBinE1,fBinE2);
    fHistEnergyDetector[1] = MakeHitPatternHist(1,calName,2*fNBinEFix,fBinE1,fBinE2);

    /////////////////////////////////////////////////////////////////////
    // 3) Fill histogram from stark event tree
    /////////////////////////////////////////////////////////////////////
    auto fileIn = new TFile(fRecoFileName);
    auto tree = (TTree*) fileIn -> Get("event");
    TClonesArray *array = nullptr;
    tree -> SetBranchAddress("SiChannel",&array);
    auto numEvents = tree -> GetEntries();
    for (auto iEvent=0; iEvent<numEvents; ++iEvent)
    {
        tree -> GetEntry(iEvent);
        auto numChannels = array -> GetEntries();
        if (iEvent%20000==0) cout << "Filling raw histogram " << iEvent << " / " << numEvents << " (" << 100*iEvent/numEvents << " %)" << endl;
        int os = -1;
        for (auto iChannel=0; iChannel<numChannels; ++iChannel)
        {
            auto channel = (LKSiChannel*) array -> At(iChannel);
            auto det = channel -> GetDetID();
            auto side = channel -> GetSide();
            if (side==1) {
                os = channel -> GetStrip();
                break;
            }
        }
        if (os<0)
            continue;
        for (auto iChannel=0; iChannel<numChannels; ++iChannel)
        {
            auto channel = (LKSiChannel*) array -> At(iChannel);
            auto det = channel -> GetDetID();
            auto side = channel -> GetSide();
            auto strip = channel -> GetStrip();
            auto numStrips = (side==0?fNumJStrips:fNumOStrips);
            if (ContinueRegardingToDataType(det)) continue;
            if (!IsPositionSensitiveStrip(channel))
            {
                auto energy = channel -> GetEnergy();
                CalibrateC0(det,side,strip,energy);
                fHistEnergy[det][side][strip] -> Fill(energy);
                fHistEnergyDetector[side] -> Fill(det*numStrips+strip,energy);
            }
            else
            {
                auto energyR = channel -> GetEnergy();
                auto energyL = channel -> GetEnergy2();
                if (energyL<0)
                    continue;
                CalibrateC1(det,side,strip,0,energyL);
                CalibrateC1(det,side,strip,1,energyR);
                auto sum = energyR + energyL;
                auto pos = (energyR - energyL) / sum;
                CalibrateC2(det,side,strip,sum,pos);
                if (calibration==3)
                    CalibrateC3(det,side,strip,sum);
                fHistEnergyDetector[side] -> Fill(det*numStrips+strip,sum);
                fHistEnergySum     [det][side][strip] -> Fill(sum);
                fHistLeftRight     [det][side][strip] -> Fill(energyR, energyL);
                fHistEnergyPosition[det][side][strip] -> Fill(pos,sum);
                fHistEnergyPositionAll -> Fill(pos,sum);
                for (auto gate=0; gate<fNumGates; ++gate) {
                    double range1 = fEnergyGateRange[gate][0];
                    double range2 = fEnergyGateRange[gate][1];
                    if (sum>range1 && sum<range2)
                        fHistEnergyPositionGateOhmic[det][side][strip][gate][os] -> Fill(pos,sum);
                }
            }
        }
    }

    /////////////////////////////////////////////////////////////////////
    // w
    /////////////////////////////////////////////////////////////////////
    if (calibration==3)
    {
        StartWriteCPParameters();
        for (auto dssGroup : fAllGroupArrayR)
        {
            auto sub = groupFit -> CreateGroup(Form("F%d",dssGroup.GetDet()));
            for (auto dss : dssGroup.array)
            {
                double xCross[2][fNumOStrips+1] = {0};
                for (auto gate=0; gate<fNumGates; ++gate)
                {
                    TH1D* histPosSum;
                    for (auto os=0; os<fNumOStrips; ++os) {
                        auto histEPGO = fHistEnergyPositionGateOhmic[dss.det][dss.side][dss.strip][gate][os];
                        TString calName2 = calName + "os" + os + "_";
                        TString name2 = MakeHistName(calName2+"rpos","",dss.det,dss.side,dss.strip,gate);
                        fHistPositionGateOhmic[dss.det][dss.side][dss.strip][gate][os] = (TH1D*) histEPGO -> ProjectionX(name2);
                        auto histPos = fHistPositionGateOhmic[dss.det][dss.side][dss.strip][gate][os];
                        auto nbins = histPos -> GetXaxis() -> GetNbins();
                        while (nbins>200 && nbins%2==0) {
                            histPos -> Rebin(2);
                            nbins = histPos -> GetXaxis() -> GetNbins();
                        }
                        if (os==0)
                            histPosSum = (TH1D*) histPos -> Clone(MakeHistName(TString("c2_osa_")+"rpos","",dss.det,dss.side,dss.strip,gate));
                        else
                            histPosSum -> Add(histPos);
                    }
                    histPosSum -> SetLineColor(kBlack);
                    double max = 0;
                    for (auto os=0; os<fNumOStrips; ++os) {
                        auto hist = fHistPositionGateOhmic[dss.det][dss.side][dss.strip][gate][os];
                        hist -> SetLineColor(os+1);
                        if (max<hist->GetMaximum()) max = hist -> GetMaximum();
                    }
                    for (auto os=0; os<fNumOStrips-1; ++os) {
                        auto hist0 = fHistPositionGateOhmic[dss.det][dss.side][dss.strip][gate][os];
                        auto hist1 = fHistPositionGateOhmic[dss.det][dss.side][dss.strip][gate][os+1];
                        auto histX = FindHistX(hist0, hist1, 0.05*max);
                        xCross[gate][os+1] = histX -> GetMean();
                        //////////////////////////////////////////////////
                        //if (dss.det==8&&dss.strip==6&&gate==1&&os==2) {
                        //    lk_debug << dss.det << " " << dss.strip << " / " << gate << " " << os << " " << xCross[gate][os+1] << endl;
                        //    auto file = new TFile("dummy.root","recreate");
                        //    hist0 -> Write();
                        //    hist1 -> Write();
                        //    histX -> Write();
                        //    file -> Close();
                        //}
                        //////////////////////////////////////////////////
                    }
                    int nbins = histPosSum -> GetXaxis() -> GetNbins();
                    while (nbins>100) {
                        histPosSum -> Rebin(2);
                        nbins = histPosSum -> GetXaxis() -> GetNbins();
                    }
                    auto drawingStepFit = FitStepHistogram(histPosSum,-0.45,0.45);
                    sub -> Add(drawingStepFit);
                    auto f1 = (TF1*) drawingStepFit -> FindObject("fit");
                    xCross[gate][0]           = f1->GetParameter(0);
                    xCross[gate][fNumOStrips] = f1->GetParameter(1);
                    //////////// XXX ////////////
                    if (abs(xCross[gate][2])>0.2) xCross[gate][2] = 0;
                    if (    xCross[gate][3]==0  ) xCross[gate][3] = xCross[gate][2] + (xCross[gate][2] - xCross[gate][1]);
                    //////////// XXX ////////////
                }
                FillCPParameters(dss.det, dss.side, dss.strip, xCross);
            }
        }
        EndWriteParameters();
    }

    /////////////////////////////////////////////////////////////////////
    // 4) Draw examples
    /////////////////////////////////////////////////////////////////////
    groupA -> AddHist(fHistEnergyPositionAll);
    auto drawingJ = groupJ -> CreateDrawing();
    auto drawingO = groupO -> CreateDrawing();
    drawingJ -> Add(fHistEnergyDetector[0]);
    drawingO -> Add(fHistEnergyDetector[1]);
    for (auto drawing : {drawingJ,drawingO}) {
        auto x2 = (drawing->GetMainHist()) -> GetXaxis() -> GetXmax();
        auto line1 = new TLine(0,fGateEnergy[0],x2,fGateEnergy[0]);
        auto line2 = new TLine(0,fGateEnergy[1],x2,fGateEnergy[1]);
        line1 -> SetLineColor(kRed);
        line2 -> SetLineColor(kRed);
        line1 -> SetLineStyle(2);
        line2 -> SetLineStyle(2);
        drawing -> Add(line1);
        drawing -> Add(line2);
    }
    for (auto dssGroup : fAllGroupArrayS)
    {
        auto sub = group0 -> CreateGroup(Form("E%d",dssGroup.GetDet()));
        for (auto dss : dssGroup.array)
            sub -> AddHist(fHistEnergy[dss.det][dss.side][dss.strip]);
    }
    for (auto dssGroup : fAllGroupArrayR)
    {
        auto sub1 = group1 -> CreateGroup(Form("ES%d",dssGroup.GetDet()));
        auto sub2 = group2 -> CreateGroup(Form("LR%d",dssGroup.GetDet()));
        auto sub3 = group3 -> CreateGroup(Form("EP%d",dssGroup.GetDet()));
        auto sub4 = group4 -> CreateGroup(Form("EPO%d",dssGroup.GetDet()));
        for (auto dss : dssGroup.array)
        {
            sub1 -> AddHist(fHistEnergySum[dss.det][dss.side][dss.strip]);
            sub2 -> AddHist(fHistLeftRight[dss.det][dss.side][dss.strip]);
            auto hist = fHistEnergyPosition[dss.det][dss.side][dss.strip];
            auto drawing3 = sub3 -> CreateDrawing();
            drawing3 -> Add(hist);
            {
                auto drawing4 = sub4 -> CreateDrawing();
                drawing4 -> SetHistCCMode();
                for (auto gate=0; gate<fNumGates; ++gate) {
                    for (auto os=0; os<fNumOStrips; ++os) {
                        auto histEP = fHistEnergyPositionGateOhmic[dss.det][dss.side][dss.strip][gate][os];
                        drawing4 -> Add(histEP);
                    }
                }
                auto line1 = new TLine(-1,fGateEnergy[0],1,fGateEnergy[0]);
                auto line2 = new TLine(-1,fGateEnergy[1],1,fGateEnergy[1]);
                line1 -> SetLineColor(kRed);
                line2 -> SetLineColor(kRed);
                drawing3 -> Add(line1);
                drawing3 -> Add(line2);
                drawing4 -> Add(line1);
                drawing4 -> Add(line2);
                if (!checkEP) {
                    for (auto gate=0; gate<fNumGates; ++gate) {
                        for (auto os=0; os<fNumOStrips+1; ++os) {
                            double energy = fGateEnergy[gate];
                            double y1 = energy - 0.3;
                            double y2 = energy + 0.3;
                            drawing3 -> Add(new TLine(fCPParameters[dss.det][dss.side][dss.strip][gate][os],y1,fCPParameters[dss.det][dss.side][dss.strip][gate][os],y2));
                            drawing4 -> Add(new TLine(fCPParameters[dss.det][dss.side][dss.strip][gate][os],y1,fCPParameters[dss.det][dss.side][dss.strip][gate][os],y2));
                        }
                    }
                    //for (auto os=0; os<fNumOStrips+1; ++os) {
                    for (auto os=1; os<fNumOStrips; ++os) {
                        double x1 = fCPParameters[dss.det][dss.side][dss.strip][0][os];
                        double y1 = fGateEnergy[0];
                        double x2 = fCPParameters[dss.det][dss.side][dss.strip][1][os];
                        double y2 = fGateEnergy[1];
                        if (abs(x2-x1)>0) {
                            auto f1 = new TF1("f1",Form("%f*(x-%f)+%f",(y2-y1)/(x2-x1),x1,y1),-1,1);
                            f1 -> SetLineColor(kBlack);
                            f1 -> SetLineWidth(1);
                            f1 -> SetLineStyle(2);
                            drawing3 -> Add(f1);
                            drawing4 -> Add(f1);
                        }
                        else {
                            auto graph = new TGraph();
                            graph -> SetPoint(0,x1,y1);
                            graph -> SetPoint(1,x2,y2);
                            graph -> SetLineWidth(1);
                            graph -> SetLineStyle(2);
                            drawing3 -> Add(graph);
                            drawing4 -> Add(graph);
                        }
                    }
                }
            }
        }
    }

    /////////////////////////////////////////////////////////////////////
    // 5) Write histograms
    /////////////////////////////////////////////////////////////////////
    TString fileName = fC2HistFileName;
    if (calibration==3)
        fileName = fC3HistFileName;
    if (calRun>=0)
        fileName = MakeCompareFileName(calRun);
    if (calRun!=run)
        fileName = Form("data/compare_%d_with_%d.hist_c%d.root",run,calRun,calibration);
    auto fileHist = new TFile(fileName,"recreate");
    top -> Write();
    cout << fileName << endl;

    if (runViewer)
        top -> Draw("viewer");
    else
        fileHist -> Close();

    return fileName;
}
