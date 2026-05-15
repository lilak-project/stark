#include "LKLogger.h"

void draw_projection_histq()
{
    TString title = "All";
    TString histName = "fHist_TTACom_vs_EAll";
    TString dummyName = "data_viewer/dummy1.root";
    auto bin1 = 300;
    auto bin2 = 399;
    if (1) {
        title = "Det16";
        histName = "fHist_ESin_vs_thetaCoM_det16";
        dummyName = "data_viewer/dummy2.root";
     bin1 = 100;
     bin2 = 170;
    }

    TH2D* hist = nullptr;
    if (1) {
        //auto group = new LKDrawingGroup("/home/ejungwoo/lilak/ko2421/macros/data_reco/ko2421_e8.ana_eve.root");
        auto group = new LKDrawingGroup("/home/ejungwoo/lilak/ko2421/macros/data_reco/ko2421_e5.ana_eve.root");
        //auto group = new LKDrawingGroup("/home/ejungwoo/lilak/ko2421/macros/data_reco/ko2421_e4.ana_eve.root");
        group -> Print();
        //auto hist = (TH2D*) group -> FindHist("fHist_TTACom_vs_EAll");
        hist = (TH2D*) group -> FindHist(histName);
        auto file = new TFile(dummyName,"recreate");
        hist -> Write();
    }
    else {
        auto file = new TFile(dummyName,"read");
        hist = (TH2D*) file -> Get(histName);
    }

    auto top = new LKDrawingGroup();
    auto group1 = top -> CreateGroup("pj");
    auto group0 = top -> CreateGroup("hist");
    auto drawing1 = group0 -> CreateDrawing();
    auto drawing0 = group0 -> CreateDrawing();

    drawing0 -> Add(hist);
    //drawing0 -> SetLogz();
    auto x1 = hist -> GetXaxis() -> GetBinLowEdge(bin1);
    auto x2 = hist -> GetXaxis() -> GetBinUpEdge(bin2);
    auto y1 = hist -> GetYaxis() -> GetXmin();
    auto y2 = hist -> GetYaxis() -> GetXmax();
    auto line1 = new TLine(x1,y1,x1,y2); line1 -> SetLineColor(kRed);
    auto line2 = new TLine(x2,y1,x2,y2); line2 -> SetLineColor(kRed);
    drawing0 -> Add(line1);
    drawing0 -> Add(line2);

    auto histQ = (TH1D*) hist -> ProjectionY("histQ",bin1,bin2);
    histQ -> SetTitle(Form("[%s] #theta_{CoM.} = (%.2f, %.2f);Q-value;count",title.Data(),x1,x2));
    //histQ -> GetXaxis() -> SetRangeUser(-3.5,4);
    drawing1 -> Add(histQ);

    auto fit = new TF1("fit","gaus(0) + gaus(3) + pol1(6)",-1.8,4);
    fit -> SetParameter(0,100);
    fit -> SetParameter(1,-1.1);
    fit -> SetParameter(2,0.1);
    fit -> SetParameter(3,1200);
    fit -> SetParameter(4,0.15);
    fit -> SetParameter(5,0.25);
    fit -> SetParameter(6,200);
    fit -> SetParameter(7,-400./8);
    histQ -> Fit(fit,"RQN0");
    auto fit1 = new TF1("fit1","gaus(0)",-1.8,4);
    auto fit2 = new TF1("fit2","gaus(0)",-1.8,4);
    auto fit3 = new TF1("fit3","pol1(0)",-1.8,4);
    double *par = fit->GetParameters();
    fit1 -> SetParameters(par[0],par[1],par[2]);
    fit2 -> SetParameters(par[3],par[4],par[5]);
    fit3 -> SetParameters(par[6],par[7]);
    fit1 -> SetLineStyle(1);
    fit2 -> SetLineStyle(1);
    fit3 -> SetLineStyle(9);
    fit1 -> SetLineWidth(2);
    fit2 -> SetLineWidth(2);
    fit3 -> SetLineWidth(1);
    fit1 -> SetLineColor(9);
    fit2 -> SetLineColor(8);
    fit3 -> SetLineColor(kBlack);
    fit  -> SetNpx(1000);
    fit1 -> SetNpx(1000);
    fit2 -> SetNpx(1000);
    fit3 -> SetNpx(1000);
    drawing1 -> Add(fit ,"samel");
    drawing1 -> Add(fit1,"samel");
    drawing1 -> Add(fit2,"samel");
    drawing1 -> Add(fit3,"samel");

    double px1 = par[1] - abs(par[2]);
    double py1 = fit -> Eval(par[1]) + par[0]*0.5;
    //double pdx = 2.2;
    double pdx = 1.2;
    //double pdy = histQ -> GetMaximum() * 0.28;
    double pdy = histQ -> GetMaximum() * 0.18;
    auto pv1 = new TPaveText(px1,py1,px1+pdx,py1+pdy);
    pv1 -> SetTextFont(132);
    pv1 -> SetTextAlign(12);
    pv1 -> SetFillColor(0);
    pv1 -> SetFillStyle(0);
    pv1 -> SetBorderSize(0);
    //pv1 -> AddText(Form("A_{1} = %.0f",par[0]));
    pv1 -> AddText(Form("#mu_{1} = %.2f",par[1]));
    pv1 -> AddText(Form("#sigma_{1} = %.2f",abs(par[2])));
    //pv1 -> AddText(Form("R_{0} = %.2f",2.3*par[2]/par[1]*100));
    drawing1 -> Add(pv1,"same");

    double px2 = par[4]+2.3*par[5];
    double py2 = 0.4*par[3];
    auto pv2 = new TPaveText(px2,py2,px2+pdx,py2+pdy);
    pv2 -> SetTextFont(132);
    pv2 -> SetTextAlign(12);
    pv2 -> SetFillColor(0);
    pv2 -> SetFillStyle(0);
    pv2 -> SetBorderSize(0);
    //pv2 -> AddText(Form("A_{0} = %.0f",par[3]));
    pv2 -> AddText(Form("#mu_{0} = %.2f",par[4]));
    pv2 -> AddText(Form("#sigma_{0} = %.2f",abs(par[5])));
    //pv2 -> AddText(Form("R_{0} = %.2f",2.3*par[5]/par[4]*100));
    drawing1 -> Add(pv2,"same");

    auto legend = new TLegend();
    legend -> AddEntry(hist,"Data","f");
    legend -> AddEntry(fit2,"Elastic","l");
    legend -> AddEntry(fit1,"1^{st} excited state","l");
    legend -> AddEntry(fit3,"Background","l");
    drawing1 -> Add(legend);

    //auto cvs0 = new TCanvas("cvs0","",1600,1100);
    //auto cvs1 = new TCanvas("cvs1","",1600,1100);

    //drawing0 -> SetCanvas(cvs0);
    //drawing0 -> Draw();
    //cvs0 -> SaveAs("qvalue_theta.C");

    //drawing1 -> SetCanvas(cvs1);
    //drawing1 -> Draw();

    group0 -> Draw();

    //cvs0 -> SaveAs(Form("qvalue_theta.%s.jpg",title.Data()));
    //cvs1 -> SaveAs(Form("qvalue_spectrum_1st_excited_state.%s.jpg",title.Data()));
    //cvs1 -> SaveAs("qvalue_spectrum_1st_excited_state.%s.C");


    //top -> Draw("v:resize=0.7");
    //top -> Draw();
}
