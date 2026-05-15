#include "si_calibration.h"

void run_slope_correction(bool runViewer=true)
{
    MakeRun();
    auto top = new LKDrawingGroup(Form("run_%04d_slope",fRun));

    /////////////////////////////////////////////////////////////////////
    // 1) Get histograms
    /////////////////////////////////////////////////////////////////////
    auto drawings = new LKDrawingGroup(fHistFileName);
    for (auto dss : fStripArrayR)
    {
        fHistLeftRight[dss.det][dss.side][dss.strip] = (TH2D*) drawings -> FindHist(MakeHistName("left","right",dss.det,dss.side,dss.strip));
        for (auto gate=0; gate<fNumGates; ++gate) {
            fHistLeftRightGate[dss.det][dss.side][dss.strip][gate] = (TH2D*) drawings -> FindHist(MakeHistName("left","right",dss.det,dss.side,dss.strip,gate));
            //fHistLRProjGate[dss.det][dss.side][dss.strip][gate][0] = MakeHist1(calName+"energyL","", dss.det,dss.side,dss.strip,gate,nbinsE,e1,e2);
            //fHistLRProjGate[dss.det][dss.side][dss.strip][gate][1] = MakeHist1(calName+"energyR","", dss.det,dss.side,dss.strip,gate,nbinsE,e1,e2);
        }
    }

    /////////////////////////////////////////////////////////////////////
    // 2) Fit 1st order polynomial for left vs right
    /////////////////////////////////////////////////////////////////////
    auto group = top -> CreateGroup("LeftRight");
    auto groupPjL = top -> CreateGroup("PjL",false);
    auto groupPjR = top -> CreateGroup("PjR",false);
    auto groupFit = top -> CreateGroup("Fit",false);
    //top -> Print(); top -> Draw("viewer"); return;
    StartWriteC1Parameters();
    int det, side, strip;
    double entries, slope, b, x1, y1, x2, y2;
    TF1* fitPol = new TF1("fitPol","pol1",0,6000);
    for (auto dssGroup : fAllGroupArrayR)
    {
        auto sub = group -> CreateGroup(Form("%d",dssGroup.GetDet()));
        auto subPjL = groupPjL -> CreateGroup(Form("%d",dssGroup.GetDet()));
        auto subPjR = groupPjR -> CreateGroup(Form("%d",dssGroup.GetDet()));
        auto subFit = groupFit -> CreateGroup(Form("%d",dssGroup.GetDet()));
        for (auto dss : dssGroup.array)
        {
            det = dss.det;
            side = dss.side;
            strip = dss.strip;
            auto graphL = new TGraph(); graphL -> SetMarkerStyle(20);
            auto graphR = new TGraph(); graphR -> SetMarkerStyle(20);
            auto drawing = sub -> CreateDrawing("");
            auto hist = fHistLeftRight[det][side][strip];
            drawing -> Add(hist);

            double xx[2][2];
            for (auto gate=0; gate<fNumGates; ++gate)
            {
                auto graphSlope = new TGraph();
                auto graphSlopeYX = new TGraph();
                graphSlope -> SetMarkerStyle(20);
                drawing -> Add(graphSlope,"","samelp");
                auto histGated = fHistLeftRightGate[det][side][strip][gate];
                auto histGatedR = histGated -> ProjectionY();
                auto histGatedL = histGated -> ProjectionX();
                for (auto hist : {histGatedR,histGatedL})
                {
                    auto nbins = hist -> GetXaxis() -> GetNbins();
                    while (nbins>100 && nbins%2==0) {
                        hist -> Rebin(2);
                        nbins = hist -> GetXaxis() -> GetNbins();
                    }
                }
                auto drawingR = FitStepHistogram(histGatedR);
                auto drawingL = FitStepHistogram(histGatedL);
                auto fitR = (TF1*) drawingR -> FindObject("fit");
                auto fitL = (TF1*) drawingL -> FindObject("fit");
                xx[gate][1] = fitR -> GetParameter(0);
                xx[gate][0] = fitL -> GetParameter(0);
                subPjR -> AddDrawing(drawingR);
                subPjL -> AddDrawing(drawingL);
                entries = histGated->GetEntries();
                if (entries<fEntriesCut) {
                    e_warning << histGated -> GetName() << " entries = " << entries << endl;
                    slope = 0;
                    b = 0;
                    x1 = 0;
                    y1 = 0;
                    x2 = 0;
                    y2 = 0;
                }
                else {
                    histGated -> Fit(fitPol,"Q0N");
                    slope = fitPol -> GetParameter(1);
                    b = fitPol -> GetParameter(0);
                    x1 = 0;
                    y1 = b;
                    x2 = -y1/slope;
                    y2 = 0;
                    graphSlope -> SetPoint(0,x1,y1);
                    graphSlope -> SetPoint(1,x2,y2);
                    graphSlopeYX -> SetPoint(0,y1,x1);
                    graphSlopeYX -> SetPoint(1,y2,x2);
                    graphL -> SetPoint(gate,x2,fGateEnergy[gate]);
                    graphR -> SetPoint(gate,y1,fGateEnergy[gate]);
                    if (0) {
                        x1 = xx[0][0];
                        y1 = graphSlope -> Eval(x1);
                        y2 = xx[0][1];
                        x2 = graphSlopeYX -> Eval(y2);
                        graphSlope -> SetPoint(0,x1,y1);
                        graphSlope -> SetPoint(1,x2,y2);
                        graphL -> SetPoint(gate,x2,fGateEnergy[gate]);
                        graphR -> SetPoint(gate,y1,fGateEnergy[gate]);
                    }
                }
            }

            TF1* fitL = new TF1("fitL","pol1",0,6000);
            graphL -> Fit(fitL,"RQN0");
            double itcptL = fitL -> GetParameter(0);
            double slopeL = fitL -> GetParameter(1);
            auto drawingL = subFit -> CreateDrawing();
            drawingL -> Add(graphL,"apl");
            drawingL -> Add(fitL,"samel");

            TF1* fitR = new TF1("fitR","pol1",0,6000);
            graphR -> Fit(fitR,"RQN0");
            double itcptR = fitR -> GetParameter(0);
            double slopeR = fitR -> GetParameter(1);
            auto drawingR = subFit -> CreateDrawing();
            drawingR -> Add(graphR,"apl");
            drawingR -> Add(fitR,"samel");

            FillC1Parameters(dss.det, dss.side, dss.strip, 2, itcptL, slopeL, itcptR, slopeR);

            auto lg = new TLegend(0.4,0.55,0.9,0.88);
            lg -> SetBorderSize(0);
            lg -> SetFillStyle(0);
            drawing -> Add(lg);
        }
    }
    EndWriteParameters();

    /////////////////////////////////////////////////////////////////////
    // 3) Draw examples
    /////////////////////////////////////////////////////////////////////

    if (runViewer) {
        top -> Print();
        top -> Draw("viewer:save_all");
        //top -> Draw("?:viewer");
    }
}
