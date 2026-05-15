const double f241AmAlphaEnergy1 = 5.486;
const double f241AmAlphaEnergy2 = 5.443;
const double f241AmAlphaBranchingRatio1 = .852;
const double f241AmAlphaBranchingRatio2 = .128;
const double f148GdAlphaEnergy = 3.1822;
const double f148GdAlphaBranchingRatio = 1;

const double fGateEnergy[] = {f148GdAlphaEnergy,f241AmAlphaEnergy1};

double EvalX(Double_t *xy, Double_t *par)
{
    double x = xy[0];

    double ADCOffset = par[0];
    double energyResolution = par[1];

    int a = 0;
    if (a==0) {
        double mean1 = par[2]; // ADC mean of 3.1822 MeV peak
        double meanP1 = mean1 - ADCOffset; // pure ADC without ADC-offset
        double sigma1 = energyResolution * meanP1;
        double amplitude1 = par[3];

        double mean2 = f241AmAlphaEnergy1 / f148GdAlphaEnergy * meanP1 + ADCOffset;
        double meanP2 = mean2 - ADCOffset;
        double sigma2 = energyResolution * meanP2;
        double amplitude2 = par[4];

        double mean3 = f241AmAlphaEnergy2 / f148GdAlphaEnergy * meanP1 + ADCOffset;
        double meanP3 = mean3 - ADCOffset;
        double sigma3 = energyResolution * meanP3;
        //double amplitude3 = amplitude2 * sigma1 / (f241AmAlphaBranchingRatio1/f241AmAlphaBranchingRatio2) / sigma2;
        double amplitude3 = amplitude2 / (f241AmAlphaBranchingRatio1/f241AmAlphaBranchingRatio2);

        double value1 = amplitude1*exp(-0.5*((x-mean1)*(x-mean1)/sigma1/sigma1));
        double value2 = amplitude2*exp(-0.5*((x-mean2)*(x-mean2)/sigma2/sigma2));
        double value3 = amplitude3*exp(-0.5*((x-mean3)*(x-mean3)/sigma3/sigma3));
        double value  = value1 + value2 + value3;

        return value;

    }
    else if (a==1) {
        double mean2 = par[2];
        double meanP2 = mean2 - ADCOffset;
        double sigma2 = energyResolution * meanP2;
        double amplitude2 = par[4];

        double mean1 =  f148GdAlphaEnergy / f241AmAlphaEnergy1 * meanP2 + ADCOffset; // ADC mean of 3.1822 MeV peak
        double meanP1 = mean1 - ADCOffset; // pure ADC without ADC-offset
        double sigma1 = energyResolution * meanP1;
        double amplitude1 = par[3];

        double mean3 = f241AmAlphaEnergy2 / f241AmAlphaEnergy1 * meanP2 + ADCOffset;
        double meanP3 = mean3 - ADCOffset;
        double sigma3 = energyResolution * meanP3;
        //double amplitude3 = amplitude2 * sigma1 / (f241AmAlphaBranchingRatio1/f241AmAlphaBranchingRatio2) / sigma2;
        double amplitude3 = amplitude2 / (f241AmAlphaBranchingRatio1/f241AmAlphaBranchingRatio2);

        double value1 = amplitude1*exp(-0.5*((x-mean1)*(x-mean1)/sigma1/sigma1));
        double value2 = amplitude2*exp(-0.5*((x-mean2)*(x-mean2)/sigma2/sigma2));
        double value3 = amplitude3*exp(-0.5*((x-mean3)*(x-mean3)/sigma3/sigma3));
        double value  = value1 + value2 + value3;

        return value;
    }
    else if (a==2) {
        double energyCollectedRatio = 1.-ADCOffset;

        double mean1 = par[2]; // ADC mean of 3.1822 MeV peak
        double meanP1 = mean1 - 0; // pure ADC without ADC-offset
        double sigma1 = energyResolution * meanP1;
        double amplitude1 = par[3] * energyCollectedRatio;

        double mean2 = f241AmAlphaEnergy1 / f148GdAlphaEnergy * meanP1 + 0;
        double meanP2 = mean2 - 0;
        double sigma2 = energyResolution * meanP2;
        double amplitude2 = par[4] * energyCollectedRatio;

        double mean3 = f241AmAlphaEnergy2 / f148GdAlphaEnergy * meanP1 + 0;
        double meanP3 = mean3 - 0;
        double sigma3 = energyResolution * meanP3;
        //double amplitude3 = amplitude2 * sigma1 / (f241AmAlphaBranchingRatio1/f241AmAlphaBranchingRatio2) / sigma2;
        double amplitude3 = amplitude2 / (f241AmAlphaBranchingRatio1/f241AmAlphaBranchingRatio2) * energyCollectedRatio;

        double value1 = amplitude1*exp(-0.5*((x-mean1)*(x-mean1)/sigma1/sigma1));
        double value2 = amplitude2*exp(-0.5*((x-mean2)*(x-mean2)/sigma2/sigma2));
        double value3 = amplitude3*exp(-0.5*((x-mean3)*(x-mean3)/sigma3/sigma3));
        double value  = value1 + value2 + value3;

        return value;

    }

    return 0;
}

/*
void run_energy_calibration()
{
    TH1D* histEnergySum[40][8];

    //////////////////////////////////////////////////////////////////////////////////////////////////////
    int selectedDet[40]; for (auto i=0; i<40; ++i) selectedDet[i] = i;

    bool drawEnergy = true;

    TString histName = "../data/stark_0253.si_cal_hist.root";
    TString calName  = "../data/stark_0253.si_calibration.root";
    //////////////////////////////////////////////////////////////////////////////////////////////////////

    auto spectrum = new TSpectrum(2);

    //auto ana = new LKECal_147Gd_241Am();
    auto fTotalFitFunction = new TF1("EvalX",EvalX,0,6000,5);

    auto fileHist = new TFile(histName,"read");
    for (auto det=0; det<40; ++det)
        for (auto strip=0; strip<8; ++strip)
            histEnergySum[det][strip] = (TH1D*) fileHist -> Get(Form("histEnergy_%d_%d",det,strip));

    //for (auto det=0; det<40; ++det)
    //for (auto det : {0,10,19,27})
    //for (auto det : {19})
    for (auto det : {27})
    {
        //auto cvs = LKPainter::GetPainter() -> CanvasResize(Form("cvsEnergy_%d",det),1000,400,0.9);
        auto cvs = LKPainter::GetPainter() -> CanvasFull(Form("cvsEnergy_%d",det),0.95);
        cvs -> Divide(3,3,0.01,0.01);
        //cvs -> SetFillColor(kGray);
        for (auto strip=0; strip<8; ++strip) {
            cvs -> cd(strip+1);
            //cvs -> cd(strip+1) -> SetFillColor(kGray);
            auto hist = histEnergySum[det][strip];
            hist -> Draw();
            double parIn[6] = {0};
            if (1) {
                auto numPeaks = spectrum -> Search(hist,5,"goff nodraw");
                if (numPeaks!=2) {
                    cout << numPeaks << "!!!"<< endl;
                    continue;
                }
                double* xPeaks = spectrum -> GetPositionX();
                double x1 = xPeaks[0];
                double x2 = xPeaks[1];
                if (x1>x2) {
                    double x3 = x2;
                    x2 = x1;
                    x1 = x3;
                }
                auto max1 = hist -> GetBinContent(hist -> FindBin(x1));
                auto max2 = hist -> GetBinContent(hist -> FindBin(x2));
                //cout << "input = " << x1 << " " << max1 << " " << max2 << endl;
                //fTotalFitFunction -> FixParameter(0,0);
                //fTotalFitFunction -> SetParameter(0,0);
                parIn[0] = -150;
                //parIn[0] = 0.02;
                //parIn[0] = 0.05;
                parIn[1] = 0.01;
                //parIn[2] = x2;
                parIn[2] = x1;
                parIn[3] = max1;
                parIn[4] = max2;
                //lk_debug  << x2/x1 << " " << f241AmAlphaEnergy1/f148GdAlphaEnergy << endl;
                fTotalFitFunction -> SetParameter(0,parIn[0]);
                //fTotalFitFunction -> FixParameter(0,0);
                fTotalFitFunction -> SetParameter(1,parIn[1]);
                fTotalFitFunction -> SetParLimits(1,0.01,0.03);
                fTotalFitFunction -> SetParameter(2,parIn[2]);
                fTotalFitFunction -> SetParameter(3,parIn[3]);
                fTotalFitFunction -> SetParameter(4,parIn[4]);
                hist -> Fit(fTotalFitFunction,"QN0");
            }
            else {
                auto bin = hist -> GetMaximumBin();
                auto mean = hist -> GetBinCenter(bin);
                auto max1 = hist -> GetBinContent(bin);
                auto max2 = max1 * 0.5;
                cout << max1 << endl;
                //auto dev = 0.2 * hist -> GetStdDev();
                fTotalFitFunction -> FixParameter(0,0);
                //fTotalFitFunction -> SetParameter(0,0);
                //fTotalFitFunction -> SetParLimit(0,0);
                fTotalFitFunction -> SetParameter(1,0.01);
                fTotalFitFunction -> SetParameter(2,mean);
                fTotalFitFunction -> SetParameter(3,max1);
                fTotalFitFunction -> SetParameter(4,max2);
                //hist -> Fit(fTotalFitFunction,"QN0");
            }
            //auto f1 = ana -> Analyze(histEnergySum[det][strip]);
            fTotalFitFunction -> SetNpx(1000);
            fTotalFitFunction -> DrawClone("samel");
            double* par = fTotalFitFunction -> GetParameters();
            auto lg1 = new TLegend(0.06,0.55,0.4,0.88);
            lg1 -> SetBorderSize(0);
            lg1 -> SetFillStyle(0);
            lg1 -> AddEntry((TObject*)nullptr,"[Input]","");
            lg1 -> AddEntry((TObject*)nullptr,Form("off = %.1f",parIn[0]),"");
            lg1 -> AddEntry((TObject*)nullptr,Form("res = %.3f",parIn[1]),"");
            lg1 -> AddEntry((TObject*)nullptr,Form("x_{1} = %.1f",parIn[2]),"");
            lg1 -> AddEntry((TObject*)nullptr,Form("a_{1} = %.1f",parIn[3]),"");
            lg1 -> AddEntry((TObject*)nullptr,Form("a_{2} = %.1f",parIn[4]),"");
            lg1 -> Draw();
            auto lg2 = new TLegend(0.6,0.55,0.9,0.88);
            lg2 -> SetBorderSize(0);
            lg2 -> SetFillStyle(0);
            lg2 -> AddEntry((TObject*)nullptr,"[Fit]","");
            lg2 -> AddEntry((TObject*)nullptr,Form("off = %.1f",par[0]),"");
            lg2 -> AddEntry((TObject*)nullptr,Form("res = %.3f",par[1]),"");
            lg2 -> AddEntry((TObject*)nullptr,Form("x1 = %.1f",par[2]),"");
            lg2 -> AddEntry((TObject*)nullptr,Form("a1 = %.1f",par[3]),"");
            //lg2 -> AddEntry((TObject*)nullptr,Form("x2 = %.1f",x2),"");
            lg2 -> AddEntry((TObject*)nullptr,Form("a2 = %.1f",par[4]),"");
            lg2 -> Draw();
        }
    }
}
*/
