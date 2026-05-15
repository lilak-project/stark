#include "si_calibration.h"

void test_final_calibration(bool test1=false)
{
    MakeRun();
    auto top = new LKDrawingGroup("test");
    auto seh = new SKEnergyHandler("data/stark_199_168.calibration.dat","data/stark_199_168.position.dat");

    /////////////////////////////////////////////////////////////////////
    // 2) Create histograms
    /////////////////////////////////////////////////////////////////////
    TString calName = "ff_";
    for (auto dss : fStripArrayS)
    {
        fHistEnergy[dss.det][dss.side][dss.strip] = MakeHist1(calName+"energy","",dss.det,dss.side,dss.strip,-1,fNBinE,fBinE1,fBinE2);
    }
    for (auto dss : fStripArrayR)
    {
        fHistEnergySum     [dss.det][dss.side][dss.strip] = MakeHist1(calName+"esum","",dss.det,dss.side,dss.strip,-1,fNBinE,fBinE1,fBinE2);
        fHistEnergyPosition[dss.det][dss.side][dss.strip] = MakeHist2(calName+"rpos",calName+"esum", dss.det,dss.side,dss.strip,-1,fNBinX,fBinX1,fBinX2,fNBinE,fBinE1,fBinE2);
        fHistEnergyPositionCal[dss.det][dss.side][dss.strip] = MakeHist2(calName+"pc_rpos",calName+"pc_esum", dss.det,dss.side,dss.strip,-1,fNBinX,fBinX1,fBinX2,fNBinE,fBinE1,fBinE2);
    }
    fHistEnergyDetector[0] = MakeHitPatternHist(0,calName,2*fNBinEFix,fBinE1,fBinE2);
    fHistEnergyDetector[1] = MakeHitPatternHist(1,calName,2*fNBinEFix,fBinE1,fBinE2);

    /////////////////////////////////////////////////////////////////////
    // 3) Fill histogram from stark event tree
    /////////////////////////////////////////////////////////////////////
    if (test1)
    {
        double dE = 10;
        for (double energyR=0; energyR<4000; energyR+=dE)
        {
            if (int(energyR)%500==0) cout << "Filling raw histogram " << energyR << " / " << 500 << " (" << 100*energyR/500 << " %)" << endl;
            for (double energyL=0; energyL<4000; energyL+=dE)
            {
                auto det = 28;
                auto side = 0;
                auto strip = 0;
                double pos, sum;
                seh -> RestoreEnergyPosition(det, side, strip, energyL, energyR, pos, sum, true);
                fHistEnergySum[det][side][strip] -> Fill(sum);
                fHistEnergyPosition[det][side][strip] -> Fill(pos,sum);
                //
                pos = seh -> RestorePosition(det, side, strip, pos, sum);
                fHistEnergyPositionCal[det][side][strip] -> Fill(pos,sum);
                fHistEnergyDetector[side] -> Fill(det*(side==0?fNumJStrips:fNumOStrips)+strip,sum);
            }
        }
    }
    else {
        auto fileIn = new TFile(fRecoFileName);
        auto tree = (TTree*) fileIn -> Get("event");
        TClonesArray *array = nullptr;
        tree -> SetBranchAddress("SiChannel",&array);
        auto numEvents = tree -> GetEntries();
        //numEvents = 100;
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
                    energy = seh -> RestoreEnergy(det, side, strip, energy);
                    fHistEnergy[det][side][strip] -> Fill(energy);
                    fHistEnergyDetector[side] -> Fill(det*(side==0?fNumJStrips:fNumOStrips)+strip,energy);
                }
                else
                {
                    auto energyR = channel -> GetEnergy();
                    auto energyL = channel -> GetEnergy2();
                    if (energyL<0)
                        continue;
                    auto sum = energyR + energyL;
                    auto pos = (energyR - energyL) / sum;
                    seh -> RestoreEnergyPosition(det, side, strip, energyL, energyR, pos, sum, true);
                    fHistEnergySum[det][side][strip] -> Fill(sum);
                    fHistEnergyPosition[det][side][strip] -> Fill(pos,sum);
                    //
                    pos = seh -> RestorePosition(det, side, strip, pos, sum);
                    fHistEnergyPositionCal[det][side][strip] -> Fill(pos,sum);
                    fHistEnergyDetector[side] -> Fill(det*(side==0?fNumJStrips:fNumOStrips)+strip,sum);
                }
            }
        }
    }

    /////////////////////////////////////////////////////////////////////
    // x) draw
    /////////////////////////////////////////////////////////////////////
    auto group3 = top -> CreateGroup("EP_pc");
    auto group2 = top -> CreateGroup("EP");
    auto group1 = top -> CreateGroup("ECal");
    for (auto dssGroup : fAllGroupArrayS)
    {
        auto sub1 = group1 -> CreateGroup(Form("E%d",dssGroup.GetDet()));
        for (auto dss : dssGroup.array) {
            auto hist = fHistEnergy[dss.det][dss.side][dss.strip];
            auto drawing = sub1 -> CreateDrawing();
            drawing -> Add(hist);
            drawing -> Add(new TLine(fGateEnergy[0],0,fGateEnergy[0],hist->GetMaximum()));
            drawing -> Add(new TLine(fGateEnergy[1],0,fGateEnergy[1],hist->GetMaximum()));
        }
    }
    for (auto dssGroup : fAllGroupArrayR)
    {
        auto sub1 = group1 -> CreateGroup(Form("E%d",dssGroup.GetDet()));
        auto sub2 = group2 -> CreateGroup(Form("E%d",dssGroup.GetDet()));
        auto sub3 = group3 -> CreateGroup(Form("E%d",dssGroup.GetDet()));
        for (auto dss : dssGroup.array)
        {
            auto hist = fHistEnergySum[dss.det][dss.side][dss.strip];
            auto hist2 = fHistEnergyPosition[dss.det][dss.side][dss.strip];
            auto hist3 = fHistEnergyPositionCal[dss.det][dss.side][dss.strip];

            auto drawing2 = sub2 -> CreateDrawing();
            drawing2 -> Add(hist2);
            auto line1 = new TLine(-1,fGateEnergy[0],1,fGateEnergy[0]);
            auto line2 = new TLine(-1,fGateEnergy[1],1,fGateEnergy[1]);
            line1 -> SetLineColor(kRed);
            line2 -> SetLineColor(kRed);
            drawing2 -> Add(line1);
            drawing2 -> Add(line2);
            //auto drawingRP = seh -> RestorePositionAndDraw(dss.det,dss.side,dss.strip,0.2,5);
            //drawing2 -> AddDrawing(drawingRP);
            //auto par = (TParameter<double>*) drawingRP -> FindObject("position");
            //auto lg = new TLegend(0.15,0.7,0.5,0.85);
            //lg -> SetBorderSize(0);
            //lg -> SetFillStyle(0);
            //lg -> AddEntry((TObject*)nullptr,Form("%f",par->GetVal()),"");
            //drawing2 -> Add(lg);

            auto drawing3 = sub3 -> CreateDrawing();
            drawing3 -> Add(hist3);
        }
    }

    top -> Draw("v");
}
