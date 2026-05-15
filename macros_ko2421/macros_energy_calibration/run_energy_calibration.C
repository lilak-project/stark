#include "si_calibration.h"

void run_energy_calibration(bool runViewer=true)
{
    MakeRun();
    auto top = new LKDrawingGroup(Form("run_%04d_energy_calibration",fRun));

    /////////////////////////////////////////////////////////////////////
    // 1) Get histograms
    /////////////////////////////////////////////////////////////////////
    auto drawings = new LKDrawingGroup(fHistFileName,"Energy:HPJ:HPO");
    drawings -> Print();
    for (auto dss : fStripArrayS)
    {
        fHistEnergy[dss.det][dss.side][dss.strip] = (TH1D*) drawings -> FindHist(MakeHistName("energy","",dss.det,dss.side,dss.strip));
    }
    top -> AddGroup(drawings -> FindGroup("HPJ"));
    top -> AddGroup(drawings -> FindGroup("HPO"));

    /////////////////////////////////////////////////////////////////////
    // 2) Calculate energy calibration parameter
    /////////////////////////////////////////////////////////////////////
    //int det, side, strip;
    double entries, a1, m1, s1, a2, m2, s2, g0;
    StartWriteC0Parameters();
    double gausParameters[40][2][8][3][3] = {0}; // det, strip, gate, range
    auto spectrum = new TSpectrum(fNumGates);
    TF1 *fitGaus = new TF1("fitGaus","gaus(0)",0,6000);
    for (auto dss : fStripArrayS)
    {
        //det = dss.det;
        //side = dss.side;
        //strip = dss.strip;
        auto hist = fHistEnergy[dss.det][dss.side][dss.strip];
        fGraphEnergyFit[dss.det][dss.side][dss.strip][0] = new TGraphErrors();
        auto graph = fGraphEnergyFit[dss.det][dss.side][dss.strip][0];
        TString graphName = Form("gr_%s",hist -> GetName());
        TString f1Name = Form("f1_%s",hist -> GetName());
        graph -> SetName(graphName);
        fFitEnergy[dss.det][dss.side][dss.strip][0] = new TF1(f1Name,"pol1",hist->GetXaxis()->GetXmin(),hist->GetXaxis()->GetXmax());
        entries = hist -> GetEntries();
        a1 = 0;
        m1 = 0;
        s1 = 0;
        a2 = 0;
        m2 = 0;
        s2 = 0;
        g0 = 0;
        if (entries<fEntriesCut)
        {
            e_warning << hist->GetName() << " entries = " << entries << endl;
            //FillC0Parameters(dss.det, dss.side, dss.strip, entries, a1, m1, s1, a2, m2, s2, g0);
            FillC0Parameters(dss.det, dss.side, dss.strip, entries, 0, 0);
            continue;
        }
        auto numPeaks = spectrum -> Search(hist,5,"goff nodraw");
        double* xPeaks = spectrum -> GetPositionX();
        if (numPeaks<fNumGates) {
            e_warning << hist->GetName() << " #peaks =" << numPeaks << endl;
            //FillC0Parameters(dss.det, dss.side, dss.strip, entries, a1, m1, s1, a2, m2, s2, g0);
            FillC0Parameters(dss.det, dss.side, dss.strip, entries, 0, 0);
            continue;
        }
        if (xPeaks[1]<xPeaks[0]) {
            auto xx = xPeaks[1];
            xPeaks[1] = xPeaks[0];
            xPeaks[0] = xx;
        }
        for (auto iPeak=0; iPeak<fNumGates; ++iPeak)
        {
            double xPeak = xPeaks[iPeak];
            fitGaus -> SetRange(xPeak-5*xPeak*fExpectedResolution,xPeak+5*xPeak*fExpectedResolution);
            fitGaus -> SetParameters(hist->GetBinContent(hist->FindBin(xPeak)),xPeak,xPeak*fExpectedResolution);
            hist -> Fit(fitGaus,"Q0NR");
            auto amp =   fitGaus -> GetParameter(0);
            auto mean =  fitGaus -> GetParameter(1);
            auto sigma = fitGaus -> GetParameter(2);
            fitGaus -> SetRange(mean-1.0*sigma, mean+2.5*sigma);
            hist -> Fit(fitGaus,"Q0NR");
            amp =   fitGaus -> GetParameter(0);
            mean =  fitGaus -> GetParameter(1);
            sigma = fitGaus -> GetParameter(2);
            gausParameters[dss.det][dss.side][dss.strip][iPeak][0] = amp;
            gausParameters[dss.det][dss.side][dss.strip][iPeak][1] = mean;
            gausParameters[dss.det][dss.side][dss.strip][iPeak][2] = sigma;
            if (iPeak==0) {
                a1 = amp;
                m1 = mean;
                s1 = sigma;
            }
            if (iPeak==1) {
                a2 = amp;
                m2 = mean;
                s2 = sigma;
                g0 = f241AmAlphaEnergy1/mean;
            }
            graph -> SetPoint(graph->GetN(),  mean, fGateEnergy[iPeak]);
            graph -> SetPointError(graph->GetN()-1,sigma,fGateEnergy[iPeak]*0.01);
        }

        auto f1 = fFitEnergy[dss.det][dss.side][dss.strip][0];
        graph -> Fit(f1,"RQN0");
        auto itcpt = f1 -> GetParameter(0);
        auto slope = f1 -> GetParameter(1);
        if (dss.det==12&&dss.strip==0) {
            graph -> Print();
            f1 -> Print();
            graph -> Fit(f1,"N0");
        }
        //FillC0Parameters(dss.det, dss.side, dss.strip, entries, a1, m1, s1, a2, m2, s2, g0);
        FillC0Parameters(dss.det, dss.side, dss.strip, entries, itcpt, slope);
    }
    EndWriteParameters();

    /////////////////////////////////////////////////////////////////////
    // 3) Draw examples
    /////////////////////////////////////////////////////////////////////
    auto group = top -> CreateGroup("energy");
    auto group2 = top -> CreateGroup("Fit");
    //for (auto dssGroup : fExampleGroupArrayS)
    for (auto dssGroup : fAllGroupArrayS)
    {
        auto sub = group -> CreateGroup(Form("%d",dssGroup.GetDet()));
        auto sub2 = group2 -> CreateGroup(Form("%d",dssGroup.GetDet()));
        for (auto dss : dssGroup.array)
        {
            auto drawing = sub -> CreateDrawing("");
            drawing -> Add(fHistEnergy[dss.det][dss.side][dss.strip]);
            //sub -> AddHist(fHistEnergy[dss.det][dss.side][dss.strip]);
            //auto lg = new TLegend(0.4,0.55,0.9,0.88);
            auto lg = new TLegend(0.35,0.60,0.85,0.88);
            lg -> SetBorderSize(0);
            lg -> SetFillStyle(0);
            lg -> SetNColumns(fNumGates);
            for (auto iPeak=0; iPeak<fNumGates; ++iPeak)
            {
                auto amp = gausParameters[dss.det][dss.side][dss.strip][iPeak][0];
                auto mean = gausParameters[dss.det][dss.side][dss.strip][iPeak][1];
                auto sigma = gausParameters[dss.det][dss.side][dss.strip][iPeak][2];
                if (amp==0)
                    continue;
                TF1 *fitGaus = new TF1(Form("fit%d%d%d%d",dss.det,dss.side,dss.strip,iPeak),"gaus(0)",0,6000);
                fitGaus -> SetParameter(0,amp);
                fitGaus -> SetParameter(1,mean);
                fitGaus -> SetParameter(2,sigma);
                fitGaus -> SetRange(mean-1.0*sigma, mean+2.5*sigma);
                drawing -> Add(fitGaus,"","samel");
                //fitGaus -> Draw("samel");
                lg -> AddEntry((TObject*)nullptr,Form("a_{%d} = %.1f",iPeak,amp),"");
                lg -> AddEntry((TObject*)nullptr,Form("m_{%d} = %.1f",iPeak,mean),"");
                lg -> AddEntry((TObject*)nullptr,Form("s_{%d} = %.1f",iPeak,sigma),"");
                lg -> AddEntry((TObject*)nullptr,Form("x_{%d} = %.1f",iPeak,fitGaus->GetChisquare()),"");
            }
            //lg -> Draw();
            drawing -> Add(lg);
            //icvs++;

            auto drawing2 = sub2 -> CreateDrawing();
            auto graph = fGraphEnergyFit[dss.det][dss.side][dss.strip][0];
            graph -> SetMarkerStyle(20);
            auto hist = fHistEnergy[dss.det][dss.side][dss.strip];
            auto frame = new TH2D(graph->GetName(),Form("%s;Raw energy;Alpha energy",hist->GetTitle()),50,hist->GetXaxis()->GetXmin(),hist->GetXaxis()->GetXmax(),50,fBinE1,fBinE2);
            drawing2 -> Add(frame);
            drawing2 -> Add(graph,"pl");
            drawing2 -> Add(fFitEnergy[dss.det][dss.side][dss.strip][0],"samel");
        }
    }

    if (runViewer)
        top -> Draw("viewer");
}
