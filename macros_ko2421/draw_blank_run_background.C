#include "FitParInfo.h"
/*
struct FitParInfo
{
    TString name = "par";
    double value = 0;
    double limit1 = 1;
    double limit2 = 2;
    bool fix = false;

    FitParInfo() = default;
    FitParInfo(TString name_, double value_, double limit1_, double limit2_, bool fix_) : name(name_), value(value_), limit1(limit1_), limit2(limit2_), fix(fix_) {}

    void Print() { cout << name << " " << value << " " << limit1 << " " << limit2 << " " << fix << endl; }
    void SetPar(TF1 *f1, int iPar) {
        f1 -> SetParName(iPar,name);
        f1 -> SetParameter(iPar,value);
        f1 -> SetParLimits(iPar,limit1,limit2);
        if (fix) f1 -> FixParameter(iPar,value);
    }
};
*/

TCutG* cutg_remove;
void RemoveInCut(TH2D* hist);

void draw_blank_run_background()
{
    LKDrawingGroup *save;
    TH2D* fHistET[28] = {0};
    bool removeInCut = true;
    auto nn = 5;

    if (1) {
        nn = 4;
        removeInCut = false;
        save = new LKDrawingGroup("data_ana/KO2421_blank_cutX_pidAll_NT360.ana_eve.root");
        for (auto pairID=0;  pairID<12; ++pairID) fHistET[pairID] = (TH2D*) save -> FindObject(Form("fHist_ESum_vs_thetaLab_pair%d",pairID));
        for (auto pairID=12; pairID<28; ++pairID) fHistET[pairID] = (TH2D*) save -> FindObject(Form("fHist_ESin_vs_thetaLab_det%d",pairID));
    }
    else {
        nn = 2;
        save = new LKDrawingGroup("data_x/example_e8_16.root");
        for (auto pairID=12; pairID<28; ++pairID) fHistET[pairID] = (TH2D*) save -> FindObject(Form("fHist_ESin_vs_thetaLab_det%d",pairID));
    }

    //cutg_remove = (TCutG*) (new TFile("data_lilak/fKDataViewer2/public_EmptyDrawing.CUTG.0.root")) -> Get("CUTG");
    cutg_remove = (TCutG*) (new TFile("data_x/example_cut.root")) -> Get("CUTG");

    TString formula = "(x>[0])*(exp(-(x-[1]))*(x>[0]?TMath::Power(x-[0],[2]):0)+[3])";
    TString formula_short = "(x>[0])*(x^{[2]}*e^{-(x-[1])}+[3])";
    LKBinning bnnSplitX12(nn,30,40);
    LKBinning bnnSplitX16(nn,50,70);
    FitParInfo fit_par[4];
    fit_par[0] = FitParInfo("_0_ x_low"   ,  2,  1.5,  2.2,  false); //fit_par[0].Print();
    fit_par[1] = FitParInfo("_1_ exp_pos" ,  5,  0  ,  10,   false); //fit_par[1].Print();
    fit_par[2] = FitParInfo("_2_ x_pow"   ,  1,  0,  4,    false); //fit_par[2].Print();
    fit_par[3] = FitParInfo("_3_ const_bg",  4,  0  ,  20,   false); //fit_par[3].Print();

    auto top = new LKDrawingGroup("blank_run_bg");

    auto groupPar = top -> CreateGroup();
    TH1D *histPar[4] = {0};
    for (auto iPar=0; iPar<4; ++iPar)
    {
        auto draw = groupPar -> CreateDrawing();
        histPar[iPar] = new TH1D(Form("hist_par%d",iPar),fit_par[iPar].fName,40,0,20);
        draw -> Add(histPar[iPar]);
    }

    auto groupAll = top -> CreateGroup();

    for (auto pairID=0; pairID<28; ++pairID)
    //for (auto pairID : {12})
    {
        auto hist2 = fHistET[pairID];
        if (hist2==nullptr) {
            e_warning << "pair " << pairID << " is null" << endl;
            continue;
        }
        if (removeInCut) RemoveInCut(hist2);
        if (hist2->GetEntries()<1000)
            continue;
        hist2 -> SetTitle(Form("Pair-%d",pairID));
        auto group = groupAll;//top -> CreateGroup();
        //auto group = top -> CreateGroup();
        auto draw = group -> CreateDrawing("",true);
        draw -> Add(hist2,"colz",Form("pair %d",pairID));
        draw -> AddLegendLine(Form("nbinsx=%d",hist2 -> GetXaxis() -> GetNbins()));
        draw -> AddLegendLine(Form("nbinsy=%d",hist2 -> GetYaxis() -> GetNbins()));
        draw -> SetCreateLegend();

        LKBinning bnn(hist2);
        LKBinning bnnSplitX = ((pairID<12)?bnnSplitX12:bnnSplitX16);
        bnn.SetProjectionBinningValues(bnnSplitX);
        auto drawGrid = bnn.CreateProjectionYGrid();
        draw -> Add(drawGrid);

        TGraph* graph[4];
        auto drawPar = group -> CreateDrawing("",false);
        drawPar -> Add((bnnSplitX*LKBinning(100,0,20)).NewH2(Form("histPar_pair%d",pairID),formula_short),".");
        for (auto i : {0,1,2,3}) {
            graph[i] = new TGraph();
            graph[i] -> SetLineColor(i+1);
            drawPar -> Add(graph[i],"l",fit_par[i].fName);
        }
        drawPar -> SetCreateLegend(1);
        while (auto hist = bnn.NextProjectionY(hist2))
        {
            hist -> Rebin(4);
            auto fit = new TF1("fit",formula,bnn.y1(),20);
            fit -> SetNpx(1000);
            for (auto i : {0,1,2,3}) fit_par[i].SetPar(fit,i);
            auto draw = group -> CreateDrawing();
            draw -> SetOptFit(1111);
            draw -> Add(hist);
            draw -> Add(fit);
            for (auto i : {0,1,2,3})
                graph[i] -> SetPoint(graph[i]->GetN(),bnn.GetCurrentProjectionCenter(),fit->GetParameter(i));
            draw -> SetFitObjects(hist,fit);
            hist -> Fit(fit,"RQ0");
            for (auto iPar=0; iPar<4; ++iPar)
                histPar[iPar] -> Fill(fit -> GetParameter(iPar));
        }
    }

    top -> Print();
    top -> WriteFitParameterFile();
    if (0) {
    top -> Draw("v");
    }
    else {
        //groupAll -> SetCanvas(TCanvas* pad) { fCvs = pad; }
        groupPar -> Draw();
        groupPar -> GetCanvas() -> SaveAs("data_x/BlankRunFitPar.png");

        groupAll -> SetCanvasDivision(5,8);
        groupAll -> SetCanvasSize(2000,2400);
        groupAll -> Draw();
        groupAll -> GetCanvas() -> SaveAs("data_x/BlankRunFitting.png");
    }
    top -> WriteFile();
}

void RemoveInCut(TH2D* hist2)
{
    auto nx = hist2 -> GetXaxis() -> GetNbins();
    auto ny = hist2 -> GetYaxis() -> GetNbins();
    for (auto ix=1; ix<=nx; ++ix) {
        for (auto iy=1; iy<=ny; ++iy) {
            auto vx = hist2 -> GetXaxis() -> GetBinCenter(ix+1);
            auto vy = hist2 -> GetYaxis() -> GetBinCenter(iy+1);
            if (hist2 -> GetBinContent(ix,iy)>0)
            {
                if (cutg_remove->IsInside(vx,vy)) {
                    hist2 -> SetBinContent(ix,iy,0);
                }
            }
        }
    }
}
