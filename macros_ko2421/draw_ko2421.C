#include "LKLogger.h"

void DrawDataCalc();
void DrawQSpectrum();
void DrawEVSTheta();
void GetArgonGraph();

const int numEnergies = 5;
TGraphErrors* fGraphData [numEnergies][3] = {0};
TGraph*       fGraphCalc [numEnergies][2] = {0};
TGraphErrors* fGraphRatio[numEnergies][3] = {0};

bool fUseLogScale = true;
bool fOnlyData = true;
int fEnergyIndices[] = {4,5,8};

double fExpEnergies[] = {4.41,  5.838,    8.13,      7.77,    9.36};
int fColorsData[] = {kRed,      kBlue,    kGreen+3 , kBlack,  kMagenta+2};
int fColorsFill[] = {kYellow-9, kCyan-10, kSpring+1, kGray,   kMagenta-10};
int fColorsCalc[] = {kOrange-3, kAzure,   kGreen+1 , kGray+2, kMagenta};
int    fMStyles[] = { 24,      25,        27,        22,      23};
double fMSizes[]  = {1.4,      1.4,       2.3,       1.4,     1.4};

vector<int> index_array = {0,1,2,3,4};
vector<int> index_others = {3,4};

int fFont = 132;
short fLWidthCalc = 2;
double fCvsScale = 1.2;
double fCvsDx = 700;
double fCvsDy = 600;

void draw_ko2421()
{
    DrawDataCalc();
    //DrawQSpectrum();
    //DrawEVSTheta();
}

void DrawDataCalc()
{
    GetArgonGraph();

    //double yoff_pad = -0.20; // legend above pad
    //double yoff_lg = 0; // legend above pad
    double yoff_pad = 0; // pad above legend
    double yoff_lg = -0.48; // pad above legend

    TString axis_title = ";#theta_{CM} (degrees);#it{d}#sigma/#it{d}#Omega (mb/sr)";
    TString axis_title_zoomin = axis_title;
    //TString axis_title_zoomin = "";

    auto DrawDataPoints = []() {
        if (fOnlyData==false) {
            for (auto i : index_array) if (fGraphData[i][1]!=nullptr) fGraphData[i][1] -> Draw("same 2");
            for (auto i : index_array) if (fGraphCalc[i][0]!=nullptr) fGraphCalc[i][0] -> Draw("same l");
            //for (auto i : index_array) fGraphData[i][2] -> Draw("same 5");
            for (auto i : index_array) if (fGraphData[i][1]!=nullptr) fGraphData[i][1] -> Draw("same p");
        }
        else {
            for (auto i : index_array) if (fGraphData[i][1]!=nullptr) fGraphData[i][1] -> Draw("same 2");
            for (auto i : index_array) if (fGraphData[i][1]!=nullptr) fGraphData[i][1] -> Draw("same pz");
        }
    };

    for (auto i : {0,1,2})
    {
        auto energy = fEnergyIndices[i];
        auto file = new TFile(Form("data_x/x%d.short.root",energy),"read");
        fGraphData[i][0]  = (TGraphErrors*) file -> Get(Form("expdata_e%d",energy));
        fGraphData[i][1]  = (TGraphErrors*) file -> Get(Form("expdata_scaled_e%d",energy));
        fGraphData[i][2]  = (TGraphErrors*) fGraphData[i][1] -> Clone(Form("expdata_scaled_e%d_empty_box",energy));
        fGraphCalc[i][0]  = (TGraph*)       file -> Get(Form("fresco_0_e%d",energy));
        fGraphCalc[i][1]  = (TGraph*)       file -> Get(Form("rutherford_e%d",energy));
        fGraphRatio[i][0] = (TGraphErrors*) file -> Get(Form("expdata_ov_rutherford_e%d",energy));
        //fGraphRatio[i][1] = (TGraphErrors*) file -> Get(Form("expdata_ov_fresco_e%d",energy));
        fGraphRatio[i][1] = (TGraphErrors*) file -> Get(Form("expdata_scaled_ov_rutherford_e%d",energy));
        fGraphRatio[i][2] = (TGraphErrors*) file -> Get(Form("fresco_ov_rutherford_e%d",energy));
    }

    auto cvs = new TCanvas("cvs_x","cross_section_all",fCvsScale*fCvsDx,fCvsScale*fCvsDy);
    cvs -> SetMargin(0.17,0.03,0.15,0.05);
    //auto hist = new TH2D("hist_x",axis_title,100,35,125,100,0,1200);
    TH2D* hist = nullptr;
    if (fUseLogScale)
        hist = new TH2D("hist_x",axis_title,100,35,125,100,3,2000);
    else
        hist = new TH2D("hist_x",axis_title,100,35,125,100,0,1200);
    hist -> GetXaxis() -> SetLabelFont(fFont);
    hist -> GetYaxis() -> SetLabelFont(fFont);
    hist -> GetXaxis() -> SetTitleFont(fFont);
    hist -> GetYaxis() -> SetTitleFont(fFont);
    hist -> GetXaxis() -> SetLabelSize(0.045);
    hist -> GetYaxis() -> SetLabelSize(0.045);
    hist -> GetXaxis() -> SetTitleSize(0.055);
    hist -> GetYaxis() -> SetTitleSize(0.055);
    hist -> GetXaxis() -> SetTitleOffset(1.12);
    hist -> GetYaxis() -> SetTitleOffset(1.4);
    hist -> SetStats(0);
    hist -> Draw();

    for (auto i : index_array)
    {
        auto graph1 = fGraphData[i][1];
        auto graph3 = fGraphData[i][2];
        auto graph2 = fGraphCalc[i][0];
        auto color_calc = fColorsCalc[i];
        auto color_data = fColorsData[i];
        auto color_fill = fColorsFill[i];
        auto mstyle = fMStyles[i];
        auto msize = fMSizes[i];

        graph1 -> SetMarkerStyle(mstyle);
        graph1 -> SetMarkerSize(msize);
        graph1 -> SetMarkerColor(color_data);
        graph1 -> SetLineColor(color_data);
        graph1 -> SetFillColor(color_fill);
        graph1 -> SetFillStyle(1001);

        graph3 -> SetLineColor(color_data);
        graph3 -> SetFillStyle(0);

        graph2 -> SetLineColor(color_calc);
        graph2 -> SetLineWidth(fLWidthCalc);
    }

    DrawDataPoints();

    if (fUseLogScale) {
        hist -> GetYaxis() -> SetRangeUser(1,5000);
        cvs -> SetLogy();
        cvs -> Modified();
        cvs -> Update();
    }
    if (0) {
        //auto pad = new TPad("zoomin_x","",0.45,0.45+yoff_pad,0.93,0.93+yoff_pad);
        auto pad = new TPad("zoomin_x","",0.48,0.45+yoff_pad,0.96,0.93+yoff_pad);
        pad -> SetMargin(0.18,0.05,0.18,0.02);
        pad -> SetFillStyle(4000);
        pad -> SetNumber(1);
        pad -> Draw();

        pad -> cd();
        auto hist_zoomin = new TH2D("hist_zoomin_x",axis_title_zoomin,100,98,122,100,0,62);
        hist_zoomin -> GetXaxis() -> SetLabelFont(fFont);
        hist_zoomin -> GetYaxis() -> SetLabelFont(fFont);
        hist_zoomin -> GetXaxis() -> SetTitleFont(fFont);
        hist_zoomin -> GetYaxis() -> SetTitleFont(fFont);
        hist_zoomin -> GetXaxis() -> SetLabelSize(0.070);
        hist_zoomin -> GetYaxis() -> SetLabelSize(0.070);
        hist_zoomin -> GetXaxis() -> SetTitleSize(0.080);
        hist_zoomin -> GetYaxis() -> SetTitleSize(0.080);
        hist_zoomin -> SetStats(0);
        hist_zoomin -> Draw();

        DrawDataPoints();
    }

    cvs -> cd();
    if (1) {
        //auto lg =  new TLegend(0.48,0.75+yoff_lg,0.92,0.93+yoff_lg);
        if (fUseLogScale) yoff_lg = 0;
        double y2 = 0.91;
        double y1 = y2-0.055*index_array.size();
        double x1_off = 0;
        if (fOnlyData) {
            auto lg =  new TLegend(0.55+x1_off,y1+yoff_lg,0.96,y2+yoff_lg);
            lg -> SetBorderSize(0);
            lg -> SetFillStyle(0);
            //lg -> SetMargin(0.3);
            lg -> SetTextFont(fFont);
            for (auto i : index_array) {
                lg -> AddEntry(fGraphData[i][1],Form("%s (%.2f MeV)",((i<3?"KO2421":"Y. Oda et al.")),fExpEnergies[i]),"plfe");
            }
            lg -> Draw();
        }
        else {
            auto lg =  new TLegend(0.45+x1_off,y1+yoff_lg,0.96,y2+yoff_lg);
            lg -> SetBorderSize(0);
            lg -> SetFillStyle(0);
            lg -> SetNColumns(3);
            lg -> SetMargin(0.3);
            lg -> SetTextFont(fFont);
            for (auto i : index_array) {
                lg -> AddEntry((TObject*)0,Form("%.2f MeV :",fExpEnergies[i]),"");
                lg -> AddEntry(fGraphData[i][1],"Data","plfe");
                if (i<3) lg -> AddEntry(fGraphCalc[i][0],"Fresco(KD)","l");
                else lg -> AddEntry(fGraphCalc[i][0],"Rutherford","l");
            }
            lg -> Draw();
        }
    }
    else {
        auto lg =  new TLegend(0.52,0.75+yoff_lg,0.90,0.93+yoff_lg);
        lg -> SetBorderSize(0);
        lg -> SetFillStyle(0);
        lg -> SetNColumns(2);
        for (auto i : index_array) {
            lg -> AddEntry(fGraphData[i][1],Form("Data (%.2f MeV)",fExpEnergies[i]),"plfe");
            lg -> AddEntry(fGraphCalc[i][0],"Fresco(KD)","l");
        }
        lg -> Draw();
    }

    cvs -> SaveAs("figures_ko2421/figure_data_fresco_all_energies.png");
    cvs -> SaveAs("figures_ko2421/figure_data_fresco_all_energies.eps");
    cvs -> SaveAs("figures_ko2421/figure_data_fresco_all_energies.pdf");
    cvs -> SaveAs("figures_ko2421/figure_data_fresco_all_energies.C");
    cvs -> SaveAs("figures_ko2421/figure_data_fresco_all_energies.root");

    auto cvs_ratio = new TCanvas("cvs_ratio","cvs_ratio",fCvsScale*fCvsDx,fCvsScale*fCvsDy);
    cvs_ratio -> SetMargin(0.17,0.03,0.15,0.05);
    auto hist_ratio = new TH2D("hist_ratio",axis_title,100,35,125,100,0,2);
    hist_ratio -> GetXaxis() -> SetLabelFont(fFont);
    hist_ratio -> GetYaxis() -> SetLabelFont(fFont);
    hist_ratio -> GetXaxis() -> SetTitleFont(fFont);
    hist_ratio -> GetYaxis() -> SetTitleFont(fFont);
    hist_ratio -> GetXaxis() -> SetLabelSize(0.045);
    hist_ratio -> GetYaxis() -> SetLabelSize(0.045);
    hist_ratio -> GetXaxis() -> SetTitleSize(0.055);
    hist_ratio -> GetYaxis() -> SetTitleSize(0.055);
    hist_ratio -> GetXaxis() -> SetTitleOffset(1.12);
    hist_ratio -> GetYaxis() -> SetTitleOffset(1.4);
    hist_ratio -> SetStats(0);
    hist_ratio -> Draw();
    for (auto i : index_array)
    {
        auto graphr1 = fGraphRatio[i][1];
        auto graphr2 = fGraphRatio[i][2];
        auto color_calc = fColorsCalc[i];
        auto color_data = fColorsData[i];
        auto color_fill = fColorsFill[i];
        auto mstyle = fMStyles[i];
        auto msize = fMSizes[i];

        graphr1 -> SetMarkerStyle(mstyle);
        graphr1 -> SetMarkerSize(msize);
        graphr1 -> SetMarkerColor(color_data);
        graphr1 -> SetLineColor(color_data);
        graphr1 -> SetFillColor(color_fill);
        graphr1 -> SetFillStyle(1001);

        graphr2 -> SetLineColor(color_calc);
        graphr2 -> SetLineWidth(fLWidthCalc);
    }
    for (auto i : index_array) fGraphRatio[i][2] -> Draw("samel"); //"fresco_ov_rutherford_e%d"
    for (auto i : index_array) fGraphRatio[i][1] -> Draw("samep"); //"expdata_ov_rutherford_e%d"

    if (1) {
        auto lg2 =  new TLegend(0.32,0.68,0.78,0.9);
        lg2 -> SetBorderSize(0);
        lg2 -> SetFillStyle(0);
        lg2 -> SetNColumns(3);
        lg2 -> SetMargin(0.3);
        lg2 -> SetTextFont(fFont);
        for (auto i : index_array) {
            lg2 -> AddEntry((TObject*)0,Form("%.2f MeV :",fExpEnergies[i]),"");
            lg2 -> AddEntry(fGraphData[i][1],"Data","ple");
            lg2 -> AddEntry(fGraphCalc[i][0],"Fresco(KD)","l");
        }
        lg2 -> Draw();
    }

    cvs_ratio -> SaveAs("figures_ko2421/figure_crosssection_over_rutherford.png");
    cvs_ratio -> SaveAs("figures_ko2421/figure_crosssection_over_rutherford.eps");
    cvs_ratio -> SaveAs("figures_ko2421/figure_crosssection_over_rutherford.pdf");
    cvs_ratio -> SaveAs("figures_ko2421/figure_crosssection_over_rutherford.C");
    cvs_ratio -> SaveAs("figures_ko2421/figure_crosssection_over_rutherford.root");
}

void DrawQSpectrum()
{
    double fit_range1 = -1.8;
    double fit_range2 = 6;

    //auto file = new TFile("data_x/histQT_813_pair14.root","read");
    //auto histQT = (TH2D*) file -> Get("histQT_813_pair14");
    //auto top = new LKDrawingGroup("data_x/histQT_813_pair14_2.root");
    ///auto histQT = (TH2D*) top -> FindHist("histQT_89_pair14");
    auto top = new LKDrawingGroup("data_x/histQT_813_pair14_3.root");
    auto histQT = (TH2D*) top -> FindHist("histQT_86_pair14");
    int bin1 = 75;
    int bin2 = 144;
    auto histQ = histQT -> ProjectionY("histQ",bin1,bin2);
    auto x1 = histQT -> GetXaxis() -> GetBinLowEdge(bin1);
    auto x2 = histQT -> GetXaxis() -> GetBinUpEdge(bin2);
    //histQT -> SetTitle(";Q-value (MeV);Count");
    histQT -> SetStats(0);
    //histQT -> Rebin(2);
    //histQT -> GetXaxis() -> SetRangeUser(-2.7,2.7);
    histQT -> GetXaxis() -> SetLabelFont(fFont);
    histQT -> GetYaxis() -> SetLabelFont(fFont);
    histQT -> GetXaxis() -> SetTitleFont(fFont);
    histQT -> GetYaxis() -> SetTitleFont(fFont);
    histQT -> GetXaxis() -> SetLabelSize(0.045);
    histQT -> GetYaxis() -> SetLabelSize(0.045);
    histQT -> GetXaxis() -> SetTitleSize(0.055);
    histQT -> GetYaxis() -> SetTitleSize(0.055);
    histQT -> GetXaxis() -> SetTitleOffset(1.12);
    histQT -> GetYaxis() -> SetTitleOffset(1.4);
    histQ -> SetTitle(";Q-value (MeV);Counts");
    histQ -> SetStats(0);
    histQ -> Rebin(2);
    histQ -> GetXaxis() -> SetRangeUser(-2.7,2.7);
    histQ -> GetXaxis() -> SetLabelFont(fFont);
    histQ -> GetYaxis() -> SetLabelFont(fFont);
    histQ -> GetXaxis() -> SetTitleFont(fFont);
    histQ -> GetYaxis() -> SetTitleFont(fFont);
    histQ -> GetXaxis() -> SetLabelSize(0.045);
    histQ -> GetYaxis() -> SetLabelSize(0.045);
    histQ -> GetXaxis() -> SetTitleSize(0.055);
    histQ -> GetYaxis() -> SetTitleSize(0.055);
    histQ -> GetXaxis() -> SetTitleOffset(1.12);
    histQ -> GetYaxis() -> SetTitleOffset(1.4);

    auto cvs_qt = new TCanvas("cvs_qt","cvs_qt",fCvsScale*fCvsDx,fCvsScale*fCvsDy);
    cvs_qt -> SetMargin(0.17,0.13,0.15,0.10);
    cvs_qt -> SetGridx();
    cvs_qt -> SetGridy();
    histQT -> SetStats(0);
    histQT -> Draw("colz");
    cvs_qt -> SaveAs("figures_ko2421/figure_qt.png");
    cvs_qt -> SaveAs("figures_ko2421/figure_qt.eps");
    cvs_qt -> SaveAs("figures_ko2421/figure_qt.pdf");
    cvs_qt -> SaveAs("figures_ko2421/figure_qt.C");
    cvs_qt -> SaveAs("figures_ko2421/figure_qt.root");

    auto cvs = new TCanvas("cvs_q","cross_section_all",fCvsScale*fCvsDx,fCvsScale*fCvsDy);
    cvs -> SetMargin(0.17,0.03,0.15,0.05);
    histQ -> Draw();

    //auto lg = new TLegend(0.22,0.6,0.65,0.90);
    //auto lg = new TLegend(0.65,0.6,0.98,0.90);
    //auto lg = new TLegend(0.63,0.6,0.96,0.90);
    auto lg = new TLegend(0.63,0.65,0.96,0.9);
    lg -> SetBorderSize(0);
    lg -> SetFillStyle(0);
    //lg -> SetNColumns(3);
    //lg -> SetMargin(0.3);
    lg -> SetTextSize(0.042);
    lg -> SetTextFont(fFont);
    //lg -> SetHeader(Form("#theta_{CM} = %.0f - %.0f deg.",x1,x2));
    //lg -> AddEntry(histQ,"E_{lab} = 8.3 MeV","f");
    

    if (1) {
        auto fit0 = new TF1("fit0","(gaus(0)+gaus(3)+[6]*exp(-(x-[7]))+[8])",fit_range1,fit_range2);
        fit0 -> SetParameter(0,  277.717); // guas0
        fit0 -> SetParameter(1, 0.283929); // guas0
        fit0 -> SetParameter(2,-0.174612); // guas0
        fit0 -> SetParameter(3,  39.7811); // guas3
        fit0 -> SetParameter(4, -1.19785); // guas3
        fit0 -> SetParameter(5, 0.103295); // guas3
        fit0 -> SetParameter(6, 0.101205); // expb_amp
        fit0 -> SetParameter(7,  5.04704); // expb_off
        fit0 -> SetParameter(8,  11.2889); // constb
        histQ -> Fit(fit0,"RN0");
        auto fit1 = new TF1("fit1","gaus(0)",fit_range1,fit_range2); fit1 -> SetParameters(fit0->GetParameter(0),fit0->GetParameter(1),fit0->GetParameter(2));
        auto fit2 = new TF1("fit2","gaus(0)",fit_range1,fit_range2); fit2 -> SetParameters(fit0->GetParameter(3),fit0->GetParameter(4),fit0->GetParameter(5));
        auto fit3 = new TF1("fit3","[0]*exp(-(x-[1]))+[2]",fit_range1,fit_range2); fit3 -> SetParameters(fit0->GetParameter(6),fit0->GetParameter(7),fit0->GetParameter(8));
        fit1 -> SetRange(fit1->GetParameter(1)+5*fit1->GetParameter(2),fit1->GetParameter(1)-5*fit1->GetParameter(2));
        fit2 -> SetRange(fit2->GetParameter(1)+5*fit2->GetParameter(2),fit2->GetParameter(1)-5*fit2->GetParameter(2));
        for (auto fit : {fit0,fit1,fit2,fit3}) {
            fit -> SetLineWidth(2);
            fit -> SetNpx(1000);
        }
        //fit0 -> SetLineColor(fColorsCalc[0]);
        fit0 -> SetLineColor(kRed);
        fit1 -> SetLineColor(fColorsCalc[1]);
        fit2 -> SetLineColor(fColorsCalc[2]);
        fit3 -> SetLineColor(kGreen+2);
        //fit3 -> Draw("samel");
        fit2 -> Draw("samel");
        fit1 -> Draw("samel");
        fit0 -> Draw("samel");

        cout << "guas0 amp   = "<< fit0 -> GetParameter(0) << endl;
        cout << "guas0 mean  = "<< fit0 -> GetParameter(1) << endl;
        cout << "guas0 sigma = "<< fit0 -> GetParameter(2) << endl;
        cout << "guas3 amp   = "<< fit0 -> GetParameter(3) << endl;
        cout << "guas3 mean  = "<< fit0 -> GetParameter(4) << endl;
        cout << "guas3 sigma = "<< fit0 -> GetParameter(5) << endl;
        cout << "expb_amp    = "<< fit0 -> GetParameter(6) << endl;
        cout << "expb_off    = "<< fit0 -> GetParameter(7) << endl;
        cout << "constb      = "<< fit0 -> GetParameter(8) << endl;

        lg -> AddEntry(fit0,"Total","l");
        fit1 -> SetLineStyle(9);
        fit2 -> SetLineStyle(2);
        lg -> AddEntry(fit1,"Ground state","l");
        lg -> AddEntry(fit2,"1^{st} excited state","l");
    }

    lg -> Draw();

    //auto pp = new TLegend(0.65,0.78,0.98,0.84);
    //pp -> SetBorderSize(0);
    //pp -> SetFillStyle(0);
    //pp -> SetTextSize(0.042);
    //pp -> SetTextFont(fFont);
    //pp -> SetHeader("E_{CM} = 8.13 MeV");
    //pp -> Draw("same");

    cvs -> SaveAs("figures_ko2421/figure_q_value_spectrum.png");
    cvs -> SaveAs("figures_ko2421/figure_q_value_spectrum.eps");
    cvs -> SaveAs("figures_ko2421/figure_q_value_spectrum.pdf");
    cvs -> SaveAs("figures_ko2421/figure_q_value_spectrum.C");
    cvs -> SaveAs("figures_ko2421/figure_q_value_spectrum.root");
}

void DrawEVSTheta()
{
    gStyle -> SetNumberContours(100);
    //gStyle -> SetPalette(kBlueYellow);
    //gStyle -> SetPalette(kCool);
    //gStyle -> SetPalette(kLightTerrain);
    //gStyle -> SetPalette(kBeach);
    //gStyle -> SetPalette(kInvertedDarkBodyRadiator);

    //auto top = new LKDrawingGroup("/home/ejungwoo/lilak/ko2421/macros/data_ana/KO2421_e8_cutX_pidAll_NT360.ana_eve.root");
    //auto top = new LKDrawingGroup("data_ana/KO2421_e8_cutX_pidAll_NT360.ana_eve.root");
    auto top = new LKDrawingGroup("data_x/EVSTheta.root");
    auto hist = top -> FindHist("fHist_ESum_vs_thetaLab_pair3");
    auto hist2 = top -> FindHist("fHist_ESin_vs_thetaLab_det14");
    //{
    //    auto top2 = new LKDrawingGroup("top2");
    //    top2 -> AddHist(hist);
    //    top2 -> AddHist(hist2);
    //    auto file2 = new TFile("data_x/EVSTheta.root","recreate");
    //    top2 -> Print();
    //    top2 -> Write();
    //    return;
    //}
    hist -> Add(hist2);

    //fCvsScale = 10;
    auto cvs = new TCanvas("cvs_e_vs_theta","e_vs_theta",fCvsScale*fCvsDx,fCvsScale*fCvsDy);
    cvs -> SetMargin(0.14,0.12,0.15,0.05);
    cvs -> SetLogz();

    //auto hist = (TH2D*) top -> FindHist("fHist_TTALab_vs_EAll");
    hist -> GetXaxis() -> SetRangeUser(24,78);
    hist -> GetYaxis() -> SetRangeUser(0,25);
    hist -> SetMinimum(1);
    hist -> SetMaximum(110);
    hist -> SetMinimum(0.4);
    //hist -> SetMaximum(200);
    hist -> GetXaxis() -> SetLabelFont(fFont);
    hist -> GetYaxis() -> SetLabelFont(fFont);
    hist -> GetZaxis() -> SetLabelFont(fFont);
    hist -> GetXaxis() -> SetTitleFont(fFont);
    hist -> GetYaxis() -> SetTitleFont(fFont);
    hist -> GetZaxis() -> SetTitleFont(fFont);
    hist -> GetXaxis() -> SetLabelSize(0.045);
    hist -> GetYaxis() -> SetLabelSize(0.045);
    hist -> GetZaxis() -> SetLabelSize(0.045);
    hist -> GetXaxis() -> SetTitleSize(0.055);
    hist -> GetYaxis() -> SetTitleSize(0.055);
    hist -> GetZaxis() -> SetTitleSize(0.055);
    hist -> GetXaxis() -> SetTitleOffset(1.12);
    hist -> GetYaxis() -> SetTitleOffset(1.12);
    hist -> SetTitle("");
    hist -> SetStats(0);
    //hist -> SetTitle(";#theta_{lab} (deg.);E_{p} (MeV)");
    //hist -> SetTitle(";#theta_{lab} (deg.);E_{#it{p}} (MeV)");
    hist -> SetTitle(";#theta_{lab} (deg.);E_{total} (MeV)");
    hist -> Draw("colz");

    double x1 = 29.8936;
    double y1 = 22.0267;
    double x2 = 71.3436;
    double y2 = 1.73658;
    double dy2 = 2.8;
    double dy3 = -3.5;
    double dy4 = dy3+2.8;

    LKGeoLine line1(x1, y1, 0, x2, y2, 0);
    auto pos1i = line1.GetPointAtX(28.4);
    auto pos1f = line1.GetPointAtX(71.4);
    LKGeoLine line2(x1, y1+dy2, 0, x2, y2+dy2, 0);
    auto pos2i = line2.GetPointAtX(31);
    auto pos2f = line2.GetPointAtX(74);

    LKGeoLine line3(x1, y1+dy3, 0, x2, y2+dy3, 0);
    auto pos3i = line3.GetPointAtX(28.4);
    auto pos3f = line3.GetPointAtX(66.4);
    LKGeoLine line4(x1, y1+dy4, 0, x2, y2+dy4, 0);
    auto pos4i = line4.GetPointAtX(31);
    auto pos4f = line4.GetPointAtX(69);

    auto graph1 = new TGraph();
    graph1 -> SetPoint(graph1->GetN(), pos1i.x(),pos1i.y());
    graph1 -> SetPoint(graph1->GetN(), pos1f.x(),pos1f.y());
    graph1 -> SetPoint(graph1->GetN(), pos2f.x(),pos2f.y());
    graph1 -> SetPoint(graph1->GetN(), pos2i.x(),pos2i.y());
    graph1 -> SetPoint(graph1->GetN(), pos1i.x(),pos1i.y());
    //auto graph1S = LKSAM::GetSAM() -> SmoothCorners(graph1,1);

    auto graph2 = new TGraph();
    graph2 -> SetPoint(graph2->GetN(), pos3i.x(),pos3i.y());
    graph2 -> SetPoint(graph2->GetN(), pos3f.x(),pos3f.y());
    graph2 -> SetPoint(graph2->GetN(), pos4f.x(),pos4f.y());
    graph2 -> SetPoint(graph2->GetN(), pos4i.x(),pos4i.y());
    graph2 -> SetPoint(graph2->GetN(), pos3i.x(),pos3i.y());

    graph1 -> SetLineColor(kRed);
    graph2 -> SetLineColor(kBlack);
    graph2 -> SetLineStyle(9);

    graph1 -> SetLineWidth(2);
    graph2 -> SetLineWidth(2);

    //graph1 -> SetLineStyle(2);
    //graph2 -> SetLineStyle(2);

    //graph1 -> SetMarkerStyle(20);
    graph1 -> Draw("samel");
    graph2 -> Draw("samel");

    cvs -> SaveAs("figures_ko2421/figure_e_vs_theta.png");
    cvs -> SaveAs("figures_ko2421/figure_e_vs_theta.eps");
    cvs -> SaveAs("figures_ko2421/figure_e_vs_theta.pdf");
    cvs -> SaveAs("figures_ko2421/figure_e_vs_theta.C");
    cvs -> SaveAs("figures_ko2421/figure_e_vs_theta.root");
}


void GetArgonGraph()
{
    double ep;
    double angle;
    double diff_ov_ruth;
    double diff;
    double diff_ruth;

    for (auto i : index_others) {
        fGraphData[i][0]  = new TGraphErrors();
        fGraphData[i][1]  = new TGraphErrors();
        fGraphData[i][2]  = new TGraphErrors();
        fGraphRatio[i][0] = new TGraphErrors();
        fGraphRatio[i][1] = new TGraphErrors();
        fGraphRatio[i][2] = new TGraphErrors();
        fGraphCalc[i][0]  = new TGraph();
    }

    std::ifstream file("data_x/argon_scattering_data.txt");
    std::string line;
    std::getline(file, line);
    while (file >> ep >> angle >> diff_ov_ruth >> diff >> diff_ruth)
    {
        int i = -1;
        for (auto j : index_others) { if (abs(fExpEnergies[j]-ep)<0.01) { i = j; break; } }
        if (i<0) continue;
        cout << ep << " " << i << endl;

        fGraphData[i][0]  -> SetPoint(fGraphData[i][0] ->GetN(),angle,diff);
        fGraphData[i][1]  -> SetPoint(fGraphData[i][1] ->GetN(),angle,diff);
        fGraphData[i][2]  -> SetPoint(fGraphData[i][2] ->GetN(),angle,diff);
        fGraphRatio[i][0] -> SetPoint(fGraphRatio[i][0]->GetN(),angle,diff_ov_ruth);
        fGraphRatio[i][1] -> SetPoint(fGraphRatio[i][1]->GetN(),angle,diff_ov_ruth);
        fGraphRatio[i][2] -> SetPoint(fGraphRatio[i][2]->GetN(),angle,diff_ov_ruth);
        fGraphCalc[i][0]  -> SetPoint(fGraphCalc[i][0] ->GetN(),angle,diff_ruth);
    }

    std::ifstream file2("data_x/argon_scattering_data..8p5.txt");
    while (file2 >> angle >> diff_ov_ruth) {
    }
}
