#include "si_calibration.h"

void run_ballistic_correction(bool runViewer=true)
{
    auto top = new LKDrawingGroup(Form("run_%04d_ballistic",fRun));
    MakeRun();

    bool usingSymmetricFit = false;

    /////////////////////////////////////////////////////////////////////
    // 1) Get histograms
    /////////////////////////////////////////////////////////////////////
    //auto drawings = new LKDrawingGroup(fC1HistFileName,"EnergyPosition");
    auto drawings = new LKDrawingGroup(fC1HistFileName);
    for (auto dss : fStripArrayR) {
        fHistEnergyPosition[dss.det][dss.side][dss.strip] = (TH2D*) drawings -> FindHist(MakeHistName("c1_rpos","c1_esum",dss.det,dss.side,dss.strip));
        for (auto gate=0; gate<fNumGates; ++gate) {
            fHistEnergyPositionGate[dss.det][dss.side][dss.strip][gate] = (TH2D*) drawings -> FindHist(MakeHistName("c1_rpos","c1_esum",dss.det,dss.side,dss.strip,gate));
        }
    }
    //fHistEnergyDetector[0] = (TH2D*) drawings->FindHist("hist_j_det_j_c1_esum");
    //fHistEnergyDetector[1] = (TH2D*) drawings->FindHist("hist_o_det_o_c1_esum");
    //top -> CreateGroup("junction") -> AddHist(fHistEnergyDetector[0]);
    //top -> CreateGroup("ohmic")    -> AddHist(fHistEnergyDetector[1]);

    /////////////////////////////////////////////////////////////////////
    // 2) Fit 1st order polynomial for left vs right
    /////////////////////////////////////////////////////////////////////
    auto group = top -> CreateGroup("energy");
    TF1* fitPol = new TF1("fitPol","pol2",-1,1);
    StartWriteC2Parameters();
    for (auto dssGroup : fAllGroupArrayR)
    {
        auto sub = group -> CreateGroup(Form("%d",dssGroup.GetDet()));
        for (auto dss : dssGroup.array)
        {
            auto det = dss.det;
            auto side = dss.side;
            auto strip = dss.strip;
            auto drawing = sub -> CreateDrawing("");
            drawing -> Add(fHistEnergyPosition[dss.det][dss.side][dss.strip]);
            for (auto gate=0; gate<fNumGates; ++gate)
            {
                auto hist = fHistEnergyPositionGate[det][side][strip][gate];
                auto entries = hist -> GetEntries();
                double b0 = -1;
                double b1 = -1;
                double b2 = -1;
                if (entries<fEntriesCut) {
                    e_warning << hist -> GetName() << " entries = " << entries << endl;
                }
                else {
                    fitPol -> SetParameter(0,f241AmAlphaEnergy1);
                    fitPol -> SetParameter(2,0.2);
                    fitPol -> SetParameter(1,0.2);
                    if (usingSymmetricFit)
                        fitPol -> FixParameter(1,0);
                    hist -> Fit(fitPol,"Q0N");
                    b0 = fitPol -> GetParameter(0);
                    b1 = fitPol -> GetParameter(1);
                    b2 = fitPol -> GetParameter(2);
                    drawing -> Add(fitPol->Clone());
                    auto pt = new TPaveText(-0.8,b0-1.2,0.8,b0-0.4);
                    pt -> AddText(Form("%.2fx^{2}+%.2fx+%.2f",b2,b1,b0));
                    pt -> SetTextFont(132);
                    pt -> SetTextSize(0.055);
                    pt -> SetFillColor(kWhite);
                    drawing -> Add(pt);
                }
                if (gate==fChooseGate) {
                    FillC2Parameters(det, side, strip, entries, b0, b1, b2);
                    drawing -> SetRangeUserY(0,b0+3.0);
                }
            }
        }
    }
    EndWriteParameters();

    /////////////////////////////////////////////////////////////////////
    // 3) Draw examples
    /////////////////////////////////////////////////////////////////////
    GetC2Parameters();
    if (runViewer)
        top -> Draw("viewer");
}
