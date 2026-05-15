#include "si_calibration.h"

void make_histograms(bool runViewer=true, int calibration=0)
{
    MakeRun();
    int ne = fNBinA;
    double e1 = fBinA1;
    double e2 = fBinA2;

    /////////////////////////////////////////////////////////////////////
    // 1) Get calibration parameter if calibration is 1 (slope correction
    /////////////////////////////////////////////////////////////////////
    TString calName = "";
    if (calibration==1) {
        calName = "c1_";
        GetC1Parameters();
        GetC0Parameters();
        GetCEParameters();
        ne = fNBinE;
        e1 = fBinE1;
        e2 = fBinE2;
    }

    double eGateRange[fNumGates][2] = {0};
    for (auto gate=0; gate<fNumGates; ++gate) {
        double eAlpha = (gate==0?f148GdAlphaEnergy:f241AmAlphaEnergy1);
        eGateRange[gate][1] = eAlpha + 12*eAlpha*fExpectedResolution;
        eGateRange[gate][0] = eAlpha - 12*eAlpha*fExpectedResolution;
        if (gate==0) 
            eGateRange[gate][0] = eAlpha - 24*eAlpha*fExpectedResolution;
        if (gate==fNumGates-1) 
            eGateRange[gate][1] = eAlpha + 24*eAlpha*fExpectedResolution;
    }

    /////////////////////////////////////////////////////////////////////
    // 2) Create histograms
    /////////////////////////////////////////////////////////////////////
    for (auto dss : fStripArrayS)
    {
        fHistEnergy[dss.det][dss.side][dss.strip] = MakeHist1(calName+"energy","",dss.det,dss.side,dss.strip,-1); // 0
    }
    for (auto dss : fStripArrayR)
    {
        fHistEnergySum     [dss.det][dss.side][dss.strip] = MakeHist1(calName+"esum",             "",dss.det,dss.side,dss.strip,-1,ne,e1,e2);
        fHistLeftRight     [dss.det][dss.side][dss.strip] = MakeHist2(calName+"left",calName+"right",dss.det,dss.side,dss.strip,-1,ne,e1,e2,ne,e1,e2); // XXX 0
        fHistEnergyPosition[dss.det][dss.side][dss.strip] = MakeHist2(calName+"rpos",calName+"esum", dss.det,dss.side,dss.strip,-1,fNBinX,fBinX1,fBinX2,ne,e1,1.5*e2); // 1
        for (auto gate=0; gate<fNumGates; ++gate) {
            fHistLeftRightGate     [dss.det][dss.side][dss.strip][gate] = MakeHist2(calName+"left",calName+"right",dss.det,dss.side,dss.strip,gate,ne,e1,e2,ne,e1,e2); // XXX 0
            fHistEnergyPositionGate[dss.det][dss.side][dss.strip][gate] = MakeHist2(calName+"rpos",calName+"esum", dss.det,dss.side,dss.strip,gate,fNBinX,fBinX1,fBinX2,ne,e1,1.5*e2); // 1
        }
    }
    fHistEnergyDetector[0] = MakeHitPatternHist(0,calName,ne,e1,e2);
    fHistEnergyDetector[1] = MakeHitPatternHist(1,calName,ne,e1,e2);

    /////////////////////////////////////////////////////////////////////
    // 3) Fill histogram from stark event tree
    /////////////////////////////////////////////////////////////////////
    auto fileIn = new TFile(fRecoFileName);
    auto tree = (TTree*) fileIn -> Get("event");
    TClonesArray *array = nullptr;
    tree -> SetBranchAddress("SiChannel",&array);
    auto numEvents = tree -> GetEntries();
    //numEvents = 10000;
    for (auto iEvent=0; iEvent<numEvents; ++iEvent)
    {
        tree -> GetEntry(iEvent);
        auto numChannels = array -> GetEntries();
        if (iEvent%20000==0) cout << "Filling raw histogram " << iEvent << " / " << numEvents << " (" << 100*iEvent/numEvents << " %)" << endl;
        for (auto iChannel=0; iChannel<numChannels; ++iChannel)
        {
            auto channel = (LKSiChannel*) array -> At(iChannel);
            auto det = channel -> GetDetID();
            auto side = channel -> GetSide();
            auto strip = channel -> GetStrip();
            if (ContinueRegardingToDataType(det)) continue;
            if (!IsPositionSensitiveStrip(channel))
            {
                auto energy = channel -> GetEnergy();
                if (calibration==1) CalibrateC0(det,side,strip,energy);
                fHistEnergy[det][side][strip] -> Fill(energy);
                fHistEnergyDetector[side] -> Fill(det*(side==0?fNumJStrips:fNumOStrips)+strip,energy);
            }
            else
            {
                auto energyR = channel -> GetEnergy();
                auto energyL = channel -> GetEnergy2();
                if (energyL<0)
                    continue;
                if (calibration==1)
                {
                    CalibrateC1(det,side,strip,0,energyL);
                    CalibrateC1(det,side,strip,1,energyR);
                    if (GetNumOhmicStrips(det)>1)
                    {
                        CalibrateCE(det,side,strip,0,energyL);
                        CalibrateCE(det,side,strip,1,energyR);
                    }
                }
                auto sum = energyR + energyL;
                auto pos = (energyR - energyL) / sum;
                fHistEnergyDetector[side] -> Fill(det*(side==0?fNumJStrips:fNumOStrips)+strip,sum);
                fHistEnergySum[det][side][strip] -> Fill(sum);
                fHistEnergyPosition[det][side][strip] -> Fill(pos,sum);
                if (calibration==0) fHistLeftRight[det][side][strip] -> Fill(energyR, energyL);
                if (0)//calibration==1)
                {
                    for (auto gate=0; gate<fNumGates; ++gate) {
                        double range1 = eGateRange[gate][0];
                        double range2 = eGateRange[gate][1];
                        if (sum>range1 && sum<range2) {
                            fHistLeftRightGate     [det][side][strip][gate] -> Fill(energyR, energyL);
                            fHistEnergyPositionGate[det][side][strip][gate] -> Fill(pos,sum);
                        }
                    }
                }
            }
        }
    }

    /////////////////////////////////////////////////////////////////////
    // 5) Calculate energy gate (0:Gd, 1:Am) for alpha source
    /////////////////////////////////////////////////////////////////////
    double gatingRange[40][8][fNumGates][2] = {0}; // det, strip, gate, range
    //if (calibration==0)
    {
        auto spectrum = new TSpectrum(fNumGates*2);
        TF1 *fitGaus = new TF1("fitGaus","gaus(0)",0,6000);
        for (auto dss : fStripArrayR)
        {
            auto hist = fHistEnergySum[dss.det][dss.side][dss.strip];
            auto numPeaks = spectrum -> Search(hist,5,"goff nodraw");
            double* xPeaks = spectrum -> GetPositionX();
            if (numPeaks<fNumGates) {
                e_warning << hist->GetName() << " #peaks =" << numPeaks << endl;
            }
            if (xPeaks[1]<xPeaks[0]) {
                auto xx = xPeaks[1];
                xPeaks[1] = xPeaks[0];
                xPeaks[0] = xx;
            }
            for (auto iPeak=0; iPeak<fNumGates; ++iPeak)
            {
                if (iPeak+1>numPeaks) {
                    e_warning << hist -> GetName() << " (" << dss.det << ", " << dss.side << ", " << dss.strip << ") : " << " #peaks = " << numPeaks << " #entries = " << hist -> GetEntries() << endl;
                    gatingRange[dss.det][dss.strip][iPeak][0] = 0;
                    gatingRange[dss.det][dss.strip][iPeak][1] = 0;
                    continue;
                }
                double xPeak = xPeaks[iPeak];
                fitGaus -> SetRange(xPeak-5*xPeak*fExpectedResolution,xPeak+5*xPeak*fExpectedResolution);
                fitGaus -> SetParameters(hist->GetBinContent(hist->FindBin(xPeak)),xPeak,xPeak*fExpectedResolution);
                hist -> Fit(fitGaus,"Q0NR");
                gatingRange[dss.det][dss.strip][iPeak][0] = fitGaus -> GetParameter(1) - 5 * fitGaus -> GetParameter(2);
                gatingRange[dss.det][dss.strip][iPeak][1] = fitGaus -> GetParameter(1) + 7 * fitGaus -> GetParameter(2);
            }
        }
    }

    /////////////////////////////////////////////////////////////////////
    // 6) Fill histogram with energy gate
    /////////////////////////////////////////////////////////////////////
    //if (calibration==0)
    {
        for (auto iEvent=0; iEvent<numEvents; ++iEvent)
        {
            tree -> GetEntry(iEvent);
            auto numChannels = array -> GetEntries();
            if (iEvent%20000==0) cout << "Filling gate histogram " << iEvent << " / " << numEvents << " (" << 100*iEvent/numEvents << " %)" << endl;
            if (numChannels!=3) continue;
            for (auto iChannel=0; iChannel<numChannels; ++iChannel)
            {
                auto channel = (LKSiChannel*) array -> At(iChannel);
                auto det = channel -> GetDetID();
                auto side = channel -> GetSide();
                auto strip = channel -> GetStrip();
                if (ContinueRegardingToDataType(det)) continue;
                if (IsPositionSensitiveStrip(channel))
                {
                    auto energyR = channel -> GetEnergy();
                    auto energyL = channel -> GetEnergy2();
                    if (energyL<0)
                        continue;
                    if (calibration==1)
                    {
                        CalibrateC1(det,side,strip,0,energyL);
                        CalibrateC1(det,side,strip,1,energyR);
                        if (GetNumOhmicStrips(det)>1)
                        {
                            CalibrateCE(det,side,strip,0,energyL);
                            CalibrateCE(det,side,strip,1,energyR);
                        }
                    }
                    auto sum = energyR + energyL;
                    auto pos = (energyR - energyL) / sum;
                    for (auto gate=0; gate<fNumGates; ++gate) {
                        //lk_debug << gatingRange[det][strip][gate][0] << " " <<  gatingRange[det][strip][gate][1] << endl;
                        if (sum>gatingRange[det][strip][gate][0] && sum<gatingRange[det][strip][gate][1])
                        {
                            fHistLeftRightGate     [det][side][strip][gate] -> Fill(energyR, energyL);
                            fHistEnergyPositionGate[det][side][strip][gate] -> Fill(pos,sum);
                        }
                    }
                }
            }
        }
    }

    /////////////////////////////////////////////////////////////////////
    // x) Draw examples
    /////////////////////////////////////////////////////////////////////
    auto top = new LKDrawingGroup(Form("run_%04d_d%d",fRun,calibration));
    auto groupJ = top -> CreateGroup("HPJ");
    auto groupO = top -> CreateGroup("HPO");
    if (calibration==1) {
        auto drawingJ = groupJ -> CreateDrawing();
        auto drawingO = groupO -> CreateDrawing();
        drawingJ -> Add(fHistEnergyDetector[0]);
        drawingO -> Add(fHistEnergyDetector[1]);
        for (auto drawing : {drawingJ,drawingO}) {
            auto x2 = (drawing->GetMainHist()) -> GetXaxis() -> GetXmax();
            auto line1 = new TLine(0,f241AmAlphaEnergy1,x2,f241AmAlphaEnergy1);
            auto line2 = new TLine(0,f148GdAlphaEnergy,x2,f148GdAlphaEnergy);
            line1 -> SetLineColor(kRed);
            line2 -> SetLineColor(kRed);
            line1 -> SetLineStyle(2);
            line2 -> SetLineStyle(2);
            drawing -> Add(line1);
            drawing -> Add(line2);
        }
    }
    else {
        groupJ -> AddHist(fHistEnergyDetector[0]);
        groupO -> AddHist(fHistEnergyDetector[1]);
    }
    if (calibration==0)
    {
        auto group0 = top -> CreateGroup("Energy");
        for (auto dssGroup : fAllGroupArrayS)
        {
            auto sub0 = group0 -> CreateGroup(Form("E%d",dssGroup.GetDet()));
            for (auto dss : dssGroup.array) {
                auto drawing0 = sub0 -> CreateDrawing();
                auto hist = fHistEnergy[dss.det][dss.side][dss.strip];
                auto damp = hist -> GetMaximum() * 0.2;
                drawing0 -> Add(hist);
                auto fitArray = FitEnergyResolution(hist, 2);
                auto numFits = fitArray -> GetEntries();
                for (auto iFit=0; iFit<numFits; ++iFit)
                {
                    auto fit = (TF1*) fitArray -> At(iFit);
                    auto amp = fit -> GetParameter(0);
                    auto mean = fit -> GetParameter(1);
                    auto sigma = fit -> GetParameter(2);
                    auto pt = new TPaveText(mean+0.35,damp,mean+1.5,damp*2);
                    pt -> AddText(Form("amp. = %.2f",amp));
                    pt -> AddText(Form("mean = %.2f",mean));
                    pt -> AddText(Form("sigma = %.4f",sigma));
                    pt -> AddText(Form("res. = %.4f",sigma/mean));
                    pt -> SetTextFont(132);
                    pt -> SetTextSize(0.035);
                    pt -> SetTextAlign(13);
                    pt -> SetFillColor(kWhite);
                    drawing0 -> Add(fit,"samel");
                    drawing0 -> Add(pt,"samel");
                }
            }
        }
        auto group2 = top -> CreateGroup("LeftRight");
        auto group4 = top -> CreateGroup("LeftRightGate");
        for (auto dssGroup : fAllGroupArrayR)
        {
            auto sub2 = group2 -> CreateGroup(Form("LR%d",dssGroup.GetDet()));
            auto sub4 = group4 -> CreateGroup(Form("LRG%d",dssGroup.GetDet()));
            for (auto dss : dssGroup.array)
            {
                sub2 -> AddHist(fHistLeftRight[dss.det][dss.side][dss.strip]);
                for (auto gate=0; gate<fNumGates; ++gate)
                    sub4 -> AddHist(fHistLeftRightGate[dss.det][dss.side][dss.strip][gate]);
            }
        }
    }
    else if (calibration==1)
    {
        auto group3 = top -> CreateGroup("EnergyPosition");
        auto group5 = top -> CreateGroup("EnergyPositionGate");
        for (auto dssGroup : fAllGroupArrayR)
        {
            auto sub3 = group3 -> CreateGroup(Form("EP%d",dssGroup.GetDet()));
            auto sub5 = group5 -> CreateGroup(Form("EPG%d",dssGroup.GetDet()));
            for (auto dss : dssGroup.array)
            {
                {
                    auto hist = fHistEnergyPosition[dss.det][dss.side][dss.strip];
                    auto drawing = sub3 -> CreateDrawing();
                    drawing -> Add(hist);
                    auto line1 = new TLine(-1,f241AmAlphaEnergy1,1,f241AmAlphaEnergy1);
                    auto line2 = new TLine(-1,f148GdAlphaEnergy ,1,f148GdAlphaEnergy );
                    line1 -> SetLineColor(kRed);
                    line2 -> SetLineColor(kRed);
                    line1 -> SetLineStyle(2);
                    line2 -> SetLineStyle(2);
                    drawing -> Add(line1);
                    drawing -> Add(line2);
                }
                for (auto gate=0; gate<fNumGates; ++gate)
                {
                    auto drawingEP = sub5 -> CreateDrawing();
                    drawingEP -> Add(fHistEnergyPositionGate[dss.det][dss.side][dss.strip][gate]);
                    if (calibration==1)
                    {
                        for (auto iLU=0; iLU<2; ++iLU) {
                            auto line = new TLine(-1,eGateRange[gate][iLU],1,eGateRange[gate][iLU]);
                            line -> SetLineColor(kRed);
                            line -> SetLineStyle(2);
                            drawingEP -> Add(line);
                        }
                    }
                }
            }
        }
    }

    /////////////////////////////////////////////////////////////////////
    // 7) Write histograms
    /////////////////////////////////////////////////////////////////////
    TString foutName = fHistFileName;
    if (calibration==1)
        foutName = fC1HistFileName;
    auto fileHist = new TFile(foutName,"recreate");
    top -> Write();
    cout << foutName << endl;

    if (runViewer) top -> Draw("viewer");
}
