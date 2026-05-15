//////////////////////////////////////////////////////////////////////////////////
TObjArray* FitEnergyResolution(TH1D* hist, int n=2, double sigr1=-1, double sigr2=-1);
LKDrawing* FitEnergyAndMakeDrawing(TH1D* hist, int n, double sigr1=-1, double sigr2=-1);
double SiAnaFitStep(double *xyz, double* par);
LKDrawing* FitStepHistogram(TH1D* hist,
        double fillX1 = 0,
        double fillX2 = 0,
        bool   fixSlopeToZero = false,
        double widthErrorRatio = 0,
        double thresholdRatio = 0,
        double valueErrorRatio = 0,
        double zeroErrorRatio = 0,
        double boundaryStiffness = 0);
void SetHistColor(TH2D* hist, int color, int max);
TH1D* FindHistX(TH1D* hist0, TH1D* hist1, double threshold=10, bool repeatWithHalfThreshold=true);

//////////////////////////////////////////////////////////////////////////////////
TObjArray* fListOfFits = nullptr;
TSpectrum* fSpectrum = nullptr;
const double fExpectedResolution = 0.015;

//////////////////////////////////////////////////////////////////////////////////
TObjArray* FitEnergyResolution(TH1D* hist, int n, double sigr1, double sigr2)
{
    if (fSpectrum==nullptr) {
        fSpectrum = new TSpectrum(n);
        fListOfFits = new TObjArray();
    }
    fListOfFits -> Clear();

    if (sigr1<0) sigr1 = 2.0;
    if (sigr2<0) sigr2 = 2.5;

    auto numPeaks = fSpectrum -> Search(hist,2,"goff nodraw");
    double* xPeaks = fSpectrum -> GetPositionX();
    if (numPeaks<n) {
        e_warning << hist->GetName() << " #peaks=" << numPeaks << endl;
        return fListOfFits;
    }
    if (n==2&&xPeaks[1]<xPeaks[0]) {
        auto xx = xPeaks[1];
        xPeaks[1] = xPeaks[0];
        xPeaks[0] = xx;
    }
    for (auto iPeak=0; iPeak<n; ++iPeak)
    {
        double mean = xPeaks[iPeak];
        double amp = hist->GetBinContent(hist->FindBin(mean));
        double sigma = mean*fExpectedResolution;
        //TF1 *fitGaus = new TF1(Form("fit_%d_%s",iPeak,hist->GetName()),"gaus(0)",0,mean+10*mean*fExpectedResolution);
        TF1 *fitGaus = new TF1(Form("fit%d",iPeak),"gaus(0)",0,mean+10*mean*fExpectedResolution);
        fitGaus -> SetRange(mean-2*sigr1*sigma,mean+2*sigr2*sigma);
        fitGaus -> SetParameters(amp,mean,sigma);
        hist -> Fit(fitGaus,"Q0NR");
        mean =  fitGaus -> GetParameter(1);
        sigma = fitGaus -> GetParameter(2);
        fitGaus -> SetRange(mean-sigr1*sigma, mean+sigr2*sigma);
        hist -> Fit(fitGaus,"Q0NR");
        fListOfFits -> Add(fitGaus);
    }

    return fListOfFits;
}

LKDrawing* FitEnergyAndMakeDrawing(TH1D* hist, int n, double sigr1, double sigr2)
{
    auto drawing = new LKDrawing();
    if (hist==nullptr)
        return drawing;
    TObjArray* fitArray = FitEnergyResolution(hist, n, sigr1, sigr2);
    drawing -> Add(hist);
    auto numPeaks = fitArray -> GetEntries();
    auto lg = new TLegend(0.35,0.50,0.85,0.88);
    lg -> SetBorderSize(0);
    lg -> SetFillStyle(0);
    lg -> SetNColumns(numPeaks);
    //lg -> SetTextSize(0.035);
    double amp[10];
    double mean[10];
    double sigma[10];
    double width[10];
    double resol[10];
    for (auto i=0; i<numPeaks; ++i)
    {
        auto fit = (TF1*) fitArray -> At(i);
        drawing -> Add(fit);
        amp[i] = fit -> GetParameter(0);
        mean[i] = fit -> GetParameter(1);
        sigma[i] = fit -> GetParameter(2);
        double xMin, xMax, half;
        //lk_debug << mean[i] << endl;
        LKSAM::GetSAM() -> FWHM(hist, mean[i], amp[i], 0, xMin, xMax, half);
        //lk_debug << xMin << " " << xMax << endl;
        auto line = new TLine(xMin,half,xMax,half);
        width[i] = xMax - xMin;
        resol[i] = width[i]/mean[i];
        drawing -> Add(line);
        drawing -> Add(new TParameter<double>(Form("fwhm%d",i),width[i]));
    }
    for (auto i=0; i<numPeaks; ++i) lg -> AddEntry((TObject*)nullptr,Form("a_{%d} = %.0f",i,amp  [i]),"");
    for (auto i=0; i<numPeaks; ++i) lg -> AddEntry((TObject*)nullptr,Form("m_{%d} = %.2f",i,mean [i]),"");
    for (auto i=0; i<numPeaks; ++i) lg -> AddEntry((TObject*)nullptr,Form("s_{%d} = %.2f",i,sigma[i]),"");
    for (auto i=0; i<numPeaks; ++i) lg -> AddEntry((TObject*)nullptr,Form("w_{%d} = %.3f",i,width[i]),"");
    for (auto i=0; i<numPeaks; ++i) lg -> AddEntry((TObject*)nullptr,Form("r_{%d} = %.2f",i,resol[i]),"");
    drawing -> Add(lg);
    drawing -> AddOption("legend_dx",0.4);
    drawing -> AddOption("legend_line_dy",0.06);
    drawing -> AddOption("opt_stat",10);
    return drawing;
}

double SiAnaFitStep(double *xyz, double *par)
{
    double x = xyz[0];

    double x1 = par[0];
    double x2 = par[1];
    double value1 = par[2];
    double slopeM = par[3];
    double slopeE = par[4];
    double value2 = slopeM * (x2-x1) + value1;
    double dx1 = value1/slopeE;
    double dx2 = value2/slopeE;

    double value = 0;
    if      (x> x1      && x<x2      ) value =  slopeM * (x-x1) + value1;
    else if (x>(x1-dx1) && x<x1      ) value =  slopeE * (x-x1) + value1;
    else if (x> x2      && x<(x2+dx2)) value = -slopeE * (x-x2) + value2;
    return value;
}

LKDrawing* FitStepHistogram(TH1D* hist,
        double fillX1,
        double fillX2,
        bool   fixSlopeToZero,
        double widthErrorRatio,
        double thresholdRatio,
        double valueErrorRatio,
        double zeroErrorRatio,
        double boundaryStiffness)
{
    auto drawing = new LKDrawing();
    auto graph = new TGraphErrors();

    // options ////////////////////////////////////////// 
    if (widthErrorRatio  ==0) widthErrorRatio = 0.1;
    if (thresholdRatio   ==0) thresholdRatio = 0.1;
    if (valueErrorRatio  ==0) valueErrorRatio = 0.1;
    if (zeroErrorRatio   ==0) zeroErrorRatio = 0.005;
    if (boundaryStiffness==0) boundaryStiffness = 0.02;

    // input histogram parameter ////////////////////////////////////////// 
    auto hmax = hist->GetMaximum();
    int nbins = hist -> GetXaxis() -> GetNbins();
    double x1 = hist -> GetXaxis() -> GetXmin();
    double x2 = hist -> GetXaxis() -> GetXmax();
    double wx = (x2-x1)/nbins;
    double ex = 0.5*wx;

    // histogram boundary ////////////////////////////////////////// 
    bool found_bx1 = false;
    bool found_bx2 = false;
    double bx1 = 0;
    double bx2 = 0;
    double valThreshold = thresholdRatio*hmax;
    double errorv = valueErrorRatio*hmax;
    double error0 = zeroErrorRatio*hmax;
    for (auto bin=1; bin<=nbins; ++bin) {
        double value = hist -> GetBinContent(bin);
        double x0 = hist -> GetBinCenter(bin);
        if (value>valThreshold && found_bx1==false) {
            found_bx1 = true;
            bx1 = x0;
            break;
        }
    }
    for (auto bin=nbins; bin>=1; --bin) {
        double value = hist -> GetBinContent(bin);
        double x0 = hist -> GetBinCenter(bin);
        if (value>valThreshold && found_bx2==false) {
            found_bx2 = true;
            bx2 = x0;
            break;
        }
    }

    // fit parameters ////////////////////////////////////////// 
    double width = bx2 - bx1;
    double widthError = widthErrorRatio*width;
    double xLE = bx1 + widthError;
    double xHE = bx2 - widthError;
    int binLE = hist -> FindBin(xLE);
    int binHE = hist -> FindBin(xHE);
    double vLE = hist -> GetBinContent(binLE);
    double vHE = hist -> GetBinContent(binHE);
    double slopeM = (vHE-vLE)/(xHE-xLE);
    double slopeE = (vHE>vLE?vHE:vLE)/(boundaryStiffness*width);
    double slopeMMin = -hmax/width;
    double slopeMMax = +hmax/width;
    double slopeEMin = (vHE>vLE?vHE:vLE)/widthError;
    double slopeEMax = 2*slopeE;
    double val1 = 0;
    int bin1 = hist -> FindBin(bx1);
    for (auto bin=bin1; bin<=binLE; ++bin) {
        if (val1<hist->GetBinContent(bin))
            val1 = hist->GetBinContent(bin);
    }
    //lk_debug << bin1 << " " << binLE << " " << val1 << endl;

    // fit ////////////////////////////////////////// 
    TString expression = "(x>[0]&&x<[1])*([3]*(x-[0])+[2])+ (x>([0]-([2]/[4]))&&x<[0])*([4]*(x-[0])+[2])+ (x>[1]&&x<([1]+([3]*([1]-[0])+[2])/[4]))*(-[4]*(x-[1])+([3]*([1]-[0])+[2]))";
    auto fit = new TF1("fit",expression,x1,x2);
    //auto fit = new TF1("fit",SiAnaFitStep,x1,x2,5);
    fit -> SetParLimits(0,bx1-widthError,bx1+widthError);
    fit -> SetParLimits(1,bx2-widthError,bx2+widthError);
    fit -> SetParLimits(2,0,1.2*hmax);
    if (fixSlopeToZero) { fit -> FixParameter(3,0); slopeM = 0; }
    else fit -> SetParLimits(3,slopeMMin,slopeMMax);
    fit -> SetParLimits(4,slopeEMin,slopeEMax);
    fit -> SetParameters(bx1,bx2,val1,slopeM,slopeE);
    fit -> SetNpx(nbins*10);

    // graph ////////////////////////////////////////// 
    for (auto bin=1; bin<=nbins; ++bin) {
        double value = hist -> GetBinContent(bin);
        double x0 = hist -> GetBinCenter(bin);
        if (value>valThreshold) {
            graph -> SetPoint(bin-1,x0,value);
            graph -> SetPointError(bin-1,ex,errorv);
        }
        else {
            if (fillX1!=fillX2 && (x0>fillX1&&x0<fillX2)) {
                value = fit -> Eval(x0);
                graph -> SetPoint(bin-1,x0,value);
                graph -> SetPointError(bin-1,ex,errorv);
            }
            else {
                graph -> SetPoint(bin-1,x0,value);
                graph -> SetPointError(bin-1,ex,error0);
            }
        }
    }
    graph -> Fit(fit,"RQN0");

    // others ////////////////////////////////////////// 
    auto lineThreshold = new TLine(x1,valThreshold,x2,valThreshold);
    lineThreshold -> SetLineColor(kRed);
    lineThreshold -> SetLineStyle(2);
    TLegend* lg = nullptr;
    if (0) {
        lg = new TLegend(0.65,0.30,0.95,0.75);
        lg -> AddEntry((TObject*)nullptr,Form("x1 = %.2f",fit->GetParameter(0)),"");
        lg -> AddEntry((TObject*)nullptr,Form("x2 = %.2f",fit->GetParameter(1)),"");
        lg -> AddEntry((TObject*)nullptr,Form("v1 = %.2f",fit->GetParameter(2)),"");
        lg -> AddEntry((TObject*)nullptr,Form("sl = %.2f",fit->GetParameter(3)),"");
        lg -> AddEntry((TObject*)nullptr,Form("sE = %.2f",fit->GetParameter(4)),"");
        lg -> AddEntry((TObject*)nullptr,Form("w? = %.4f",widthError),"");
    }
    else {
        lg = new TLegend(0.40,0.30,0.95,0.75);
        lg -> AddEntry((TObject*)nullptr,Form("x1 = %.2f (%.2f)",fit->GetParameter(0),bx1),"");
        lg -> AddEntry((TObject*)nullptr,Form("x2 = %.2f (%.2f)",fit->GetParameter(1),bx2),"");
        lg -> AddEntry((TObject*)nullptr,Form("v1 = %.2f (%.2f)",fit->GetParameter(2),val1),"");
        lg -> AddEntry((TObject*)nullptr,Form("sl = %.2f (%.2f)",fit->GetParameter(3),slopeM),"");
        lg -> AddEntry((TObject*)nullptr,Form("sE = %.2f (%.2f)",fit->GetParameter(4),slopeE),"");
        lg -> AddEntry((TObject*)nullptr,Form("w? = %.4f",widthError),"");
    }
    lg -> SetBorderSize(0);
    lg -> SetFillStyle(0);
    lg -> SetTextColor(kBlue);

    // drawings ////////////////////////////////////////// 
    drawing -> Add(hist);
    drawing -> Add(graph,"samepz");
    //drawing -> Add(graph,"apz");
    drawing -> Add(fit,"samel");
    drawing -> Add(lineThreshold,"samel");
    drawing -> Add(lg,"same");

    return drawing;
}

void SetHistColor(TH2D* hist, int color, int max)
{
    auto nx = hist -> GetXaxis() -> GetNbins();
    auto x1 = hist -> GetXaxis() -> GetXmin();
    auto x2 = hist -> GetXaxis() -> GetXmax();
    auto ny = hist -> GetYaxis() -> GetNbins();
    auto y1 = hist -> GetYaxis() -> GetXmin();
    auto y2 = hist -> GetYaxis() -> GetXmax();
    for (auto ix=1; ix<=nx; ++ix) {
        for (auto iy=1; iy<=ny; ++iy) {
            auto value = hist -> GetBinContent(ix,iy);
            if (value>0) {
                hist -> SetBinContent(ix,iy,color+1);
            }
        }
    }
    hist -> SetMaximum(max+1);
}

TH1D* fHistX = nullptr;
TH1D* FindHistX(TH1D* hist0, TH1D* hist1, double threshold, bool repeatWithHalfThreshold)
{
    auto nx = hist0 -> GetXaxis() -> GetNbins();
    if (fHistX==nullptr) {
        auto x1 = hist0 -> GetXaxis() -> GetXmin();
        auto x2 = hist0 -> GetXaxis() -> GetXmax();
        fHistX = new TH1D("fHistX","",nx,x1,x2);
    }
    else
        fHistX -> Reset();
    for (auto ix=1; ix<=nx; ++ix)
    {
        double value = 0;
        auto value0 = hist0 -> GetBinContent(ix);
        auto value1 = hist1 -> GetBinContent(ix);
        if      (value0==0 || value1==0) continue;
        else if (value0<=value1) value = value0;
        else if (value0> value1) value = value1;
        if (value>threshold) {
            fHistX -> SetBinContent(ix,value);
        }
    }
    if (repeatWithHalfThreshold && fHistX->GetEntries()==0)
    {
        return FindHistX(hist0, hist1, 0.5*threshold, false);
    }
    return fHistX;
};
