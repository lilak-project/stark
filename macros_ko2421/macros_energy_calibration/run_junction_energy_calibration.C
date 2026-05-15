#include "si_calibration.h"

TString run_junction_energy_calibration(bool runViewer=true, int run=-1, int calRun=-1, TString inputName="")
{
    MakeRun(run, ((calRun==199&&run==168)?1:-1));
    auto top = new LKDrawingGroup(Form("run_%04d_position_calibration",fRun));

    /////////////////////////////////////////////////////////////////////
    // 1)
    /////////////////////////////////////////////////////////////////////
    TString inputFileName = inputName;
    if (inputFileName.IsNull())
        inputFileName = fC2HistFileName;
    auto drawings = new LKDrawingGroup(inputFileName);
    for (auto dss : fStripArrayR)
    {
        for (auto gate=0; gate<fNumGates; ++gate)
        {
            for (auto os=0; os<fNumOStrips; ++os)
            {
                TString calName2 = TString("c2_os") + os + "_";
                TString name = MakeHistName(calName2+"rpos",calName2+"esum", dss.det,dss.side,dss.strip,gate);
                auto hist = (TH2D*) drawings -> FindHist(name);
                if (hist==nullptr) {
                    lk_debug << name << " " << hist << endl;
                    return "";
                }
                fHistEnergyPositionGateOhmic[dss.det][dss.side][dss.strip][gate][os] = hist;
                TString name2 = MakeHistName(calName2+"rpos","",dss.det,dss.side,dss.strip,gate);
                fHistPositionGateOhmic[dss.det][dss.side][dss.strip][gate][os] = (TH1D*) hist -> ProjectionX(name2);
            }
        }
    }

    /////////////////////////////////////////////////////////////////////
    // x) draw
    /////////////////////////////////////////////////////////////////////
    auto groupF = top -> CreateGroup("Fit",false);
    auto group = top -> CreateGroup("EPO");
    auto group0 = top -> CreateGroup("gate0",true);
    auto group1 = top -> CreateGroup("gate1",true);
    auto groupE = top -> CreateGroup("Energy",false);

    TString c3Name;
    if (inputName.IsNull()==false) {
        c3Name = inputName;
        c3Name.ReplaceAll("hist_c2","c3");
        c3Name.ReplaceAll(".root",".dat");
    }
    StartWriteC3Parameters(c3Name);
    for (auto dssGroup : fAllGroupArrayR)
    {
        auto sub = group -> CreateGroup(Form("%d",dssGroup.GetDet()));
        auto subE = groupE -> CreateGroup(Form("%d",dssGroup.GetDet()));
        auto sub0 = group0 -> CreateGroup(Form("%d",dssGroup.GetDet()));
        auto sub1 = group1 -> CreateGroup(Form("%d",dssGroup.GetDet()));
        auto subF = groupF -> CreateGroup(Form("%d",dssGroup.GetDet()));
        for (auto dss : dssGroup.array)
        {

            //////////////////////////////////////////////////////
            double xCross[2][fNumOStrips+1] = {0};
            {
                LKDrawing* drawing[2];
                drawing[0] = sub0 -> CreateDrawing();
                drawing[1] = sub1 -> CreateDrawing();
                for (auto gate=0; gate<fNumGates; ++gate)
                {
                    TH1D* histPosSum;
                    for (auto os=0; os<fNumOStrips; ++os) {
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
                    drawing[gate] -> Add(histPosSum);
                    double max = 0;
                    for (auto os=0; os<fNumOStrips; ++os) {
                        auto hist = fHistPositionGateOhmic[dss.det][dss.side][dss.strip][gate][os];
                        hist -> SetLineColor(os+1);
                        drawing[gate] -> Add(hist);
                        if (max<hist->GetMaximum()) max = hist -> GetMaximum();
                    }
                    for (auto os=0; os<fNumOStrips-1; ++os) {
                        auto hist0 = fHistPositionGateOhmic[dss.det][dss.side][dss.strip][gate][os];
                        auto hist1 = fHistPositionGateOhmic[dss.det][dss.side][dss.strip][gate][os+1];
                        if (hist0->GetEntries()<100 || hist1->GetEntries()<100) {
                            e_warning << hist0->GetName() << " entry is " << hist0->GetEntries() << endl;
                            e_warning << hist1->GetName() << " entry is " << hist1->GetEntries() << endl;
                            break;
                        }
                        auto histX = FindHistX(hist0, hist1, 0.05*max);
                        xCross[gate][os+1] = histX -> GetMean();
                    }
                    if (histPosSum->GetEntries()<100) {
                        e_warning << histPosSum->GetName() << " entry is " << histPosSum->GetEntries() << endl;
                    }
                    else {
                        int nbins = histPosSum -> GetXaxis() -> GetNbins();
                        while (nbins>100) {
                            histPosSum -> Rebin(2);
                            nbins = histPosSum -> GetXaxis() -> GetNbins();
                        }
                        auto drawingFitStep = FitStepHistogram(histPosSum,-0.45,0.45);
                        subF -> Add(drawingFitStep);
                        auto f1 = (TF1*) drawingFitStep -> FindObject("fit");
                        drawing[gate] -> Add(f1);
                        xCross[gate][0]           = f1->GetParameter(0);
                        xCross[gate][fNumOStrips] = f1->GetParameter(1);
                        for (auto os=0; os<fNumOStrips+1; ++os)
                            drawing[gate] -> Add(new TMarker(xCross[gate][os],0.25*max,20));
                    }
                }
            }
            //////////////////////////////////////////////////////

            //////////////////////////////////////////////////////
            {
                auto drawingEP = sub -> CreateDrawing();
                drawingEP -> SetHistCCMode();
                auto drawingE = subE -> CreateDrawing();
                TH1D* histEProj[2];
                histEProj[0] = nullptr;
                histEProj[1] = nullptr;
                double damp = 0;
                for (auto gate=0; gate<fNumGates; ++gate)
                {
                    for (auto os=0; os<fNumOStrips; ++os) {
                        auto histEP = fHistEnergyPositionGateOhmic[dss.det][dss.side][dss.strip][gate][os];
                        drawingEP -> Add(histEP);
                        if (histEProj[gate]==nullptr) {
                            TString name = MakeHistName("c2_osall_esum","",dss.det,dss.side,dss.strip,gate);
                            histEProj[gate] = (TH1D*) histEP -> ProjectionY(name);
                        }
                        else 
                            histEProj[gate] -> Add(histEP->ProjectionY());
                    }
                    drawingE -> Add(histEProj[gate]);
                    double dampt = histEProj[gate] -> GetMaximum() * 0.2;
                    if (damp<dampt) damp = dampt;
                }
                auto graphFit = new TGraph();
                auto fitPol1 = new TF1("fitPol1","pol1",0,10);
                for (auto gate=0; gate<fNumGates; ++gate)
                {
                    auto fitArray = FitEnergyResolution(histEProj[gate], 1, 2.5, 2.5);
                    auto numPeaks = fitArray -> GetEntries();
                    if (numPeaks<1)
                        e_warning << endl;
                    for (auto iPeak=0; iPeak<numPeaks; ++iPeak)
                    {
                        auto fit = (TF1*) fitArray -> At(iPeak);
                        auto amp = fit -> GetParameter(0);
                        auto mean = fit -> GetParameter(1);
                        auto sigma = fit -> GetParameter(2);
                        auto pt = new TPaveText(mean+0.35,damp,mean+1.5,damp*2);
                        pt -> AddText(Form("a = %.2f",amp));
                        pt -> AddText(Form("m = %.2f",mean));
                        pt -> AddText(Form("s = %.4f",sigma));
                        pt -> AddText(Form("r = %.4f",sigma/mean));
                        pt -> SetTextFont(132);
                        pt -> SetTextSize(0.030);
                        pt -> SetTextAlign(13);
                        pt -> SetFillColor(kWhite);
                        drawingE -> Add(fit,"samel");
                        drawingE -> Add(pt,"samel");
                        graphFit -> SetPoint(gate,mean,fGateEnergy[gate]);
                    }
                    //////////// XXX ////////////
                    if (xCross[gate][3]==0) {
                        xCross[gate][3] = xCross[gate][2] + (xCross[gate][2] - xCross[gate][1]);
                    }
                    //////////// XXX ////////////
                    //for (auto os=0; os<fNumOStrips+1; ++os) {
                    for (auto os=1; os<fNumOStrips; ++os) {
                        double energy = fGateEnergy[gate];
                        double y1 = energy - 0.3;
                        double y2 = energy + 0.3;
                        //drawingEP -> Add(new TMarker(xCross[gate][os],0.25*max,20));
                        drawingEP -> Add(new TLine(xCross[gate][os],y1,xCross[gate][os],y2));
                    }
                }
                auto line1 = new TLine(-1,fGateEnergy[0],1,fGateEnergy[0]);
                auto line2 = new TLine(-1,fGateEnergy[1],1,fGateEnergy[1]);
                line1 -> SetLineColor(kRed);
                line2 -> SetLineColor(kRed);
                drawingEP -> Add(line1);
                drawingEP -> Add(line2);
                graphFit -> Fit(fitPol1,"RQN0");
                if (0) {
                    auto drawingFit = subE -> CreateDrawing();
                    graphFit -> SetMarkerStyle(20);
                    drawingFit -> Add(graphFit,"apl");
                    drawingFit -> Add(fitPol1,"samel");
                }
                FillC3Parameters(dss.det, dss.side, dss.strip, 2, fitPol1->GetParameter(0), fitPol1->GetParameter(1));
                for (auto os=1; os<fNumOStrips; ++os) {
                //for (auto os=0; os<fNumOStrips+1; ++os) {
                    double x1 = xCross[0][os];
                    double y1 = fGateEnergy[0];
                    double x2 = xCross[1][os];
                    double y2 = fGateEnergy[1];
                    if (abs(x2-x1)>0) {
                        auto f1 = new TF1("f1",Form("%f*(x-%f)+%f",(y2-y1)/(x2-x1),x1,y1),-1,1);
                        f1 -> SetLineColor(kBlack);
                        f1 -> SetLineWidth(1);
                        f1 -> SetLineStyle(2);
                        drawingEP -> Add(f1);
                    }
                }
            }
            //////////////////////////////////////////////////////

        }
    }
    c3Name = EndWriteParameters();

    if (runViewer)
        top -> Draw("viewer");

    return c3Name;
}
