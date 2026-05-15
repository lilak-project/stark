#include "si_calibration.h"

TH2D* histEnergyDetector[2][2];
TH2D* histEnergyDetectorAll[2];

double c0Array[40][2][8][2];    // ohmic energy calibration
double c1Array[40][2][8][2][2]; // junction slope correction
double c2Array[40][2][8][3];    // junction ballistic correction
double c3Array[40][2][8][2];    // junction energy calibration
double cpArray[40][2][8][2][5]; // junction postion parameters

std::vector<double> fEnergyResolutionArray[40][2];

void detector_summary(int det=0, bool drawFigures=true);
void get_parameters();
void write_parameters();

LKDrawingGroup* fTop = nullptr;
LKDrawingGroup* fGroupEE = nullptr;
LKDrawingGroup* fGroupEP = nullptr;
LKDrawingGroup* fGroupOthers = nullptr;
TH2D* fHistER[2] = {0};
TH1D* fHistNE[2] = {0};

int fRun1, fRun2;

void detector_summary(int det, bool drawFigures)
{
    auto subEE = fGroupEE -> CreateGroup(Form("E%d",det));
    auto subEP = fGroupEP -> CreateGroup(Form("EP%d",det));
    /////////////////////////////////////////////////////////////////////
    // create info
    /////////////////////////////////////////////////////////////////////
    strip_group dssArrayR;
    strip_group dssArrayS;
    auto detector = fStark -> GetSiDetector(det);
    auto numJStrips = detector -> GetNumJunctionStrips();
    auto numOStrips = detector -> GetNumOhmicStrips();
    auto detType = detector -> GetDetType();
    for (auto side=0; side<2; ++side)
    {
        if (IsPositionSensitiveStrip(det,side))
        {
            auto numStrips = (side==0?numJStrips:numOStrips);
            for (auto strip=0; strip<numStrips; ++strip) {
                dssArrayR.array.push_back(strip_info(detType,det,side,strip));
            }
        }
        else {
            auto numStrips = (side==0?numJStrips:numOStrips);
            for (auto strip=0; strip<numStrips; ++strip) {
                dssArrayS.array.push_back(strip_info(detType,det,side,strip));
            }
        }
    }

    /////////////////////////////////////////////////////////////////////
    // t
    /////////////////////////////////////////////////////////////////////
    double sr1 = 2.5;
    double sr2 = 2.5;
    int rebin = 1;
    if (det>=32) {
        sr1 = 2.0;
        sr2 = 2.5;
    }
    if (det==34||det==35) {
        sr1 = 0.5;
        sr2 = 2.5;
        rebin = 4;
    }

    for (auto dss : dssArrayR.array)
    {
        auto drawingEP = subEP -> CreateDrawing();
        auto histEP = fHistEnergyPosition[dss.det][dss.side][dss.strip];
        drawingEP -> Add(histEP);
        auto line1 = new TLine(-1,fGateEnergy[0],1,fGateEnergy[0]);
        auto line2 = new TLine(-1,fGateEnergy[1],1,fGateEnergy[1]);
        line1 -> SetLineColor(kRed);
        line2 -> SetLineColor(kRed);
        drawingEP -> Add(line1);
        drawingEP -> Add(line2);
        for (auto gate=0; gate<fNumGates; ++gate) {
            for (auto os=0; os<fNumOStrips+1; ++os) {
                double energy = fGateEnergy[gate];
                double y1 = energy - 0.3;
                double y2 = energy + 0.3;
                drawingEP -> Add(new TLine(fCPParameters[dss.det][dss.side][dss.strip][gate][os],y1,fCPParameters[dss.det][dss.side][dss.strip][gate][os],y2));
            }
        }
        for (auto os=0; os<fNumOStrips+1; ++os) {
            double x1 = fCPParameters[dss.det][dss.side][dss.strip][0][os];
            double y1 = fGateEnergy[0];
            double x2 = fCPParameters[dss.det][dss.side][dss.strip][1][os];
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

    for (auto dss : dssArrayR.array)
    {
        if (fHistEnergyPosition[dss.det][dss.side][dss.strip] ==nullptr)
            continue;
        auto hist = (TH1D*) fHistEnergyPosition[dss.det][dss.side][dss.strip] -> ProjectionY();
        {
            auto histNE = fHistNE[dss.side];
            auto bin = histNE -> FindBin(dss.det*(dss.side==0?fNumJStrips:fNumOStrips)+dss.strip);
            auto ne = hist -> GetEntries();
            histNE -> SetBinContent(bin,ne);
        }
        auto drawing = FitEnergyAndMakeDrawing(hist, 2, sr1, sr2);
        subEE -> AddDrawing(drawing);
        for (auto i=0; i<2; ++i) {
            if (0) {
                auto fit = (TF1*) drawing -> FindObject(Form("fit%d",i));
                if (fit==nullptr) continue;
                auto resolution = 100*fit->GetParameter(2)/fit->GetParameter(1);
                double fwhm_keV = fit -> GetParameter(2) * 1000. * 2.3;
                fEnergyResolutionArray[det][dss.side].push_back(fwhm_keV);
            }
            else {
                auto par = (TParameter<double>*) drawing -> FindObject(Form("fwhm%d",i));
                if (par==nullptr) continue;
                double fwhm_keV = par -> GetVal() * 1000;
                fEnergyResolutionArray[det][dss.side].push_back(fwhm_keV);
            }
        }
    }

    for (auto dss : dssArrayS.array)
    {
        auto hist = fHistEnergy[dss.det][dss.side][dss.strip];
        if (hist==nullptr)
            continue;
        {
            auto histNE = fHistNE[dss.side];
            auto bin = histNE -> FindBin(dss.det*(dss.side==0?fNumJStrips:fNumOStrips)+dss.strip);
            auto ne = hist -> GetEntries();
            histNE -> SetBinContent(bin,ne);
        }
        if (rebin>1) hist -> Rebin(rebin);
        auto drawing = FitEnergyAndMakeDrawing(hist, 2, sr1, sr2);
        subEE -> AddDrawing(drawing);
        for (auto i=0; i<2; ++i) {
            if (0) {
                auto fit = (TF1*) drawing -> FindObject(Form("fit%d",i));
                if (fit==nullptr) continue;
                auto resolution = 100*fit->GetParameter(2)/fit->GetParameter(1);
                double fwhm_keV = fit -> GetParameter(2) * 1000. * 2.3;
                fEnergyResolutionArray[det][dss.side].push_back(fwhm_keV);
            }
            else {
                auto par = (TParameter<double>*) drawing -> FindObject(Form("fwhm%d",i));
                if (par==nullptr) continue;
                double fwhm_keV = par -> GetVal() * 1000;
                fEnergyResolutionArray[det][dss.side].push_back(fwhm_keV);
            }
        }
    }
}

void get_parameters()
{
    GetC0Parameters();
    GetC1Parameters();
    GetC2Parameters();
    GetC3Parameters();
    GetCPParameters();
}

void write_parameters()
{
    TString name1 = Form("data/stark_%d_%d.calibration.dat",fRun1,fRun2);
    e_info << name1 << endl;
    ofstream file1(name1);
    file1
        << "#det"   << "\t"
        << "side"   << "\t"
        << "strip"  << "\t"
        << "c0_itcpt" << "\t"
        << "c0_slope" << "\t"
        << "c1_itcptL" << "\t"
        << "c1_slopeL" << "\t"
        << "c1_itcptR" << "\t"
        << "c1_slopeR" << "\t"
        << "c2_par0"   << "\t" 
        << "c2_par1"   << "\t" 
        << "c2_par2"   << "\t"
        << "c3_itcpt" << "\t"
        << "c3_slope" << endl;

    TString name2 = Form("data/stark_%d_%d.position.dat",fRun1,fRun2);
    e_info << name2 << endl;
    ofstream file2(name2);
    file2
        << "#det"   << "\t"
        << "side"   << "\t"
        << "strip"  << "\t"
        << "x00"    << "\t"
        << "x01"    << "\t"
        << "x02"    << "\t"
        << "x03"    << "\t"
        << "x04"    << "\t"
        << "x10"    << "\t"
        << "x11"    << "\t"
        << "x12"    << "\t"
        << "x13"    << "\t"
        << "x14"    << endl;

    for (auto groupArray : {fAllGroupArrayR, fAllGroupArrayS})
    {
        for (auto dssGroup : groupArray)
        {
            for (auto dss : dssGroup.array)
            {
                int det = dss.det;
                int side = dss.side;
                int strip = dss.strip;
                file1
                    << det   << "\t" 
                    << side  << "\t" 
                    << strip << "\t" 
                    << fC0Parameters[det][side][strip][0] << "\t"
                    << fC0Parameters[det][side][strip][1] << "\t"
                    << fC1Parameters[det][side][strip][0][0] << "\t"
                    << fC1Parameters[det][side][strip][0][1] << "\t"
                    << fC1Parameters[det][side][strip][1][0] << "\t"
                    << fC1Parameters[det][side][strip][1][1] << "\t"
                    << fC2Parameters[det][side][strip][0] << "\t"
                    << fC2Parameters[det][side][strip][1] << "\t"
                    << fC2Parameters[det][side][strip][2] << "\t"
                    << fC3Parameters[det][side][strip][0] << "\t"
                    << fC3Parameters[det][side][strip][1] << endl;
            }
        }
    }

    for (auto dssGroup : fAllGroupArrayR)
    {
        for (auto dss : dssGroup.array)
        {
            int det = dss.det;
            int side = dss.side;
            int strip = dss.strip;
            file2
                << det   << "\t" 
                << side  << "\t" 
                << strip << "\t" 
                << fCPParameters[det][side][strip][0][0] << "\t"
                << fCPParameters[det][side][strip][0][1] << "\t"
                << fCPParameters[det][side][strip][0][2] << "\t"
                << fCPParameters[det][side][strip][0][3] << "\t"
                << fCPParameters[det][side][strip][0][4] << "\t"
                << fCPParameters[det][side][strip][1][0] << "\t"
                << fCPParameters[det][side][strip][1][1] << "\t"
                << fCPParameters[det][side][strip][1][2] << "\t"
                << fCPParameters[det][side][strip][1][3] << "\t"
                << fCPParameters[det][side][strip][1][4] << endl;
        }
    }

    file1.close();
    file2.close();
}

void summary(bool readSave=false)
{
    fRun1 = 199;
    //fRun2 = 7718;
    fRun2 = 168;

    TString nameSummary = Form("stark_%d_%d",fRun1,fRun2);
    if (readSave) {
        fTop = new LKDrawingGroup(nameSummary);
        fTop -> Draw("v");
        return;
    }

    fTop = new LKDrawingGroup(Form("energy_calibration_summary_%d_%d",fRun1,fRun2));
    LKDrawing *drawing;

    TString calName = "c3_";
    for (auto det=0; det<40; ++det)
        for (auto side=0; side<2; ++side)
            for (auto strip=0; strip<8; ++strip) {
                c0Array[det][side][strip][0] = 0;
                c0Array[det][side][strip][1] = 0;
                c1Array[det][side][strip][0][0] = 0;
                c1Array[det][side][strip][0][1] = 0;
                c1Array[det][side][strip][1][0] = 0;
                c1Array[det][side][strip][1][1] = 0;
                c2Array[det][side][strip][0] = 0;
                c2Array[det][side][strip][1] = 0;
                c2Array[det][side][strip][2] = 0;
                c3Array[det][side][strip][0] = 0;
                c3Array[det][side][strip][1] = 0;
                for (auto i=0; i<5; ++i) cpArray[det][side][strip][0][i] = 0;
                for (auto i=0; i<5; ++i) cpArray[det][side][strip][1][i] = 0;
            }

    ////////////////////////////////////////////////////////////////////
    // 1
    /////////////////////////////////////////////////////////////////////
    if (fRun1>0)
    {
        MakeRun(fRun1);
        get_parameters();
        auto top1 = new LKDrawingGroup(Form("data/stark_0%d.hist_c3.root",fRun));
        for (auto dss : fStripArrayS)
        {
            if (ContinueRegardingToDataType(dss.det)) continue;
            auto name = MakeHistName(calName+"energy","",dss.det,dss.side,dss.strip);
            fHistEnergy[dss.det][dss.side][dss.strip] = (TH1D*) top1 -> FindHist(name);
        }
        for (auto dss : fStripArrayR)
        {
            if (ContinueRegardingToDataType(dss.det)) continue;
            fHistEnergySum     [dss.det][dss.side][dss.strip] = (TH1D*) top1 -> FindHist(MakeHistName(calName+"esum",             "",dss.det,dss.side,dss.strip));
            fHistLeftRight     [dss.det][dss.side][dss.strip] = (TH2D*) top1 -> FindHist(MakeHistName(calName+"left",calName+"right",dss.det,dss.side,dss.strip));
            fHistEnergyPosition[dss.det][dss.side][dss.strip] = (TH2D*) top1 -> FindHist(MakeHistName(calName+"rpos",calName+"esum", dss.det,dss.side,dss.strip));
        }
        histEnergyDetector[0][0] = (TH2D*) top1 -> FindHist(MakeHistName("det_j",calName+"esum",-1,0,-1,-1));
        histEnergyDetector[0][1] = (TH2D*) top1 -> FindHist(MakeHistName("det_o",calName+"esum",-1,1,-1,-1));
    }

    /////////////////////////////////////////////////////////////////////
    // 2
    /////////////////////////////////////////////////////////////////////
    if (fRun2>0)
    {
        MakeRun(fRun2);
        get_parameters();
        auto top2 = new LKDrawingGroup(Form("data/stark_0%d.hist_c3.root",fRun));
        for (auto dss : fStripArrayS)
        {
            if (ContinueRegardingToDataType(dss.det)) continue;
            fHistEnergy[dss.det][dss.side][dss.strip] = (TH1D*) top2 -> FindHist(MakeHistName(calName+"energy","",dss.det,dss.side,dss.strip));
        }
        for (auto dss : fStripArrayR)
        {
            if (ContinueRegardingToDataType(dss.det)) continue;
            fHistEnergySum     [dss.det][dss.side][dss.strip] = (TH1D*) top2 -> FindHist(MakeHistName(calName+"esum",             "",dss.det,dss.side,dss.strip));
            fHistEnergyPosition[dss.det][dss.side][dss.strip] = (TH2D*) top2 -> FindHist(MakeHistName(calName+"rpos",calName+"esum", dss.det,dss.side,dss.strip));
        }
        histEnergyDetector[1][0] = (TH2D*) top2 -> FindHist(MakeHistName("det_j",calName+"esum",-1,0,-1,-1));
        histEnergyDetector[1][1] = (TH2D*) top2 -> FindHist(MakeHistName("det_o",calName+"esum",-1,1,-1,-1));
    }

    GetC3Parameters("data/compare_168_with_199.c3.dat");

    /////////////////////////////////////////////////////////////////////
    // merge hit pattern
    /////////////////////////////////////////////////////////////////////
    histEnergyDetectorAll[0] = (TH2D*) histEnergyDetector[0][0] -> Clone("h0");
    histEnergyDetectorAll[1] = (TH2D*) histEnergyDetector[0][1] -> Clone("h1");
    histEnergyDetectorAll[0] -> SetTitle(Form("[%d,%d] Junction",fRun1,fRun2));
    histEnergyDetectorAll[1] -> SetTitle(Form("[%d,%d] Ohmic",fRun1,fRun2));
    histEnergyDetectorAll[0] -> Add(histEnergyDetector[1][0]);
    histEnergyDetectorAll[1] -> Add(histEnergyDetector[1][1]);
    histEnergyDetectorAll[0] -> SetStats(0);
    histEnergyDetectorAll[1] -> SetStats(0);
    histEnergyDetectorAll[0] -> GetXaxis() -> SetTitleOffset(2.7);
    histEnergyDetectorAll[1] -> GetXaxis() -> SetTitleOffset(2.7);
    for (auto det=0; det<fNumDetectors; ++det)
    {
        auto detector = fStark -> GetSiDetector(det);
        auto name = detector -> GetDetTypeName();
        auto ring = detector -> GetLayer();
        TString sring = "dE"; if (ring==1) sring = "E"; if (ring==2) sring = "16E"; 
        histEnergyDetectorAll[0] -> GetXaxis() -> SetBinLabel(det*fNumJStrips+1,Form("%d (%s,%s)",det,name.Data(),sring.Data()));
        histEnergyDetectorAll[1] -> GetXaxis() -> SetBinLabel(det*fNumOStrips+1,Form("%d (%s,%s)",det,name.Data(),sring.Data()));
    }
    for (auto jo : {0,1})
    {
        auto group = fTop -> CreateGroup(Form("HP%d",jo));
        drawing = group -> CreateDrawing();
        drawing -> Add(histEnergyDetectorAll[jo]);
        auto x2 = (drawing->GetMainHist()) -> GetXaxis() -> GetXmax();
        auto line1 = new TLine(0,fGateEnergy[0],x2,fGateEnergy[0]);
        auto line2 = new TLine(0,fGateEnergy[1],x2,fGateEnergy[1]);
        line1 -> SetLineColor(kRed);
        line2 -> SetLineColor(kRed);
        drawing -> Add(line1);
        drawing -> Add(line2);
        drawing -> SetLogz();
    }
    fGroupOthers = fTop -> CreateGroup("Others");
    fGroupEE = fTop -> CreateGroup("Energy",true);
    fGroupEP = fTop -> CreateGroup("EvsP",true);

    /////////////////////////////////////////////////////////////////////
    //
    /////////////////////////////////////////////////////////////////////

    fHistER[0] = (TH2D*) MakeHitPatternHist(2,"resolution0",50,0,500);
    fHistER[0] -> SetTitle(Form("[%d,%d] Energy Resolution Junction",fRun1,fRun2));
    fHistER[0] -> GetYaxis() -> SetTitle("FWHM [keV]");

    fHistER[1] = (TH2D*) MakeHitPatternHist(2,"resolution1",50,0,500);
    fHistER[1] -> SetTitle(Form("[%d,%d] Energy Resolution Ohmic",fRun1,fRun2));
    fHistER[1] -> GetYaxis() -> SetTitle("FWHM [keV]");

    fHistNE[0] = (TH1D*) MakeDet1DHist(0,"mult0");
    fHistNE[0] -> SetTitle(Form("[%d,%d] Multiplicity Junction",fRun1,fRun2));
    fHistNE[0] -> GetYaxis() -> SetTitle("");

    fHistNE[1] = (TH1D*) MakeDet1DHist(1,"mult1");
    fHistNE[1] -> SetTitle(Form("[%d,%d] Multiplicity Ohmic",fRun1,fRun2));
    fHistNE[1] -> GetYaxis() -> SetTitle("");

    for (TH1* hist : {(TH1*)fHistER[0],(TH1*)fHistER[1],(TH1*)fHistNE[0],(TH1*)fHistNE[1]})
    {
        hist -> SetStats(0);
        drawing = fGroupOthers -> CreateDrawing();
        drawing -> Add(hist);
        drawing -> SetBottomMargin(0.15);
        hist -> GetXaxis() -> SetLabelSize(0.04);
    }

    for (auto det=0; det<40; ++det)
    //for (auto det:{20})
        detector_summary(det,false);

    write_parameters();

    for (auto det=0; det<40; ++det)
    {
        for (auto resolution : fEnergyResolutionArray[det][0]) fHistER[0] -> Fill(det,resolution);
        for (auto resolution : fEnergyResolutionArray[det][1]) fHistER[1] -> Fill(det,resolution);
    }

    TString fileNameSummary = Form("data/%s.summary.root",nameSummary.Data());
    auto file = new TFile(fileNameSummary,"recreate");
    fTop -> Write();
    e_info << fileNameSummary << endl;

    fTop -> Draw("v:save_all");
    //fTop -> Draw("v");
}
