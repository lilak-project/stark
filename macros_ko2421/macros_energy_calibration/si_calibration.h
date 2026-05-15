#ifndef SICALIBRATION_H
#define SICALIBRATION_H

#include "LKLogger.h"
#include "si_energy.h"
#include "si_analysis.h"
#include "si_parameters.h"
#include "stark_energy_handler.h"

//int fRun = 7718;
//int fRun = 7728;
//int fRun = 7777;
int fRun = 199;
//int fRun = 253;
//int fRun = 303;
//
//int fRun = 113;
//int fRun = 167;
//int fRun = 168;
//int fRun = 169;
//int fRun = 170;
//int fRun = 171;

bool fTest = true;

//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
// PREDEFINE METHOD
//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
bool ContinueRegardingToDataType(int det);
bool MakeRun(int run=-1, int type=-1);
int GetNumOhmicStrips(int det);
TString MakeCompareFileName(int calRun);
TString MakeHistName(TString xname, TString yname, int det, int side, int strip, int gate=-1);
TString MakeHistTitle(TString xname, TString yname, int det, int side, int strip, int gate=-1);
TH1D* MakeHist1(TString xname, TString yname, int det, int side, int strip, int gate, int nx=0, double x1=0, double x2=0);
TH2D* MakeHist2(TString xname, TString yname, int det, int side, int strip, int gate, int nx, double x1, double x2, int ny, int y1, int y2);
TCanvas *MakeCanvas(TString cvsName, int npads);
bool IsPositionSensitiveStrip(int det, int side=1);
bool IsPositionSensitiveStrip(LKSiChannel* channel);
void GetFinalHistograms();
TH2D* MakeHitPatternHist(int jo, TString calName="", int nn=0, double e1=0, double e2=0); // jo=0: junction, jo=1: ohmic
TH1D* MakeDet1DHist(int jo, TString calName=""); // jo=0: junction, jo=1: ohmic

//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
// VARIABLES BELOW
//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////

int fDataType = -1;  // 0: 12dE+16E, 1: 12E
int fChooseGate = 1; // 0: Gd, 1: Am
int fEntriesCut = 100; // will be updated as a function of total entries

/***/ int    fNBinA = 200;
const double fBinA1 = 0;
const double fBinA2 = 4000;
/***/ int    fNBinE = 200;
const double fBinE1 = 0;
const double fBinE2 = 10;
const int    fNBinEFix = 400;
/***/ int    fNBinX = 200;
const double fBinX1 = -1;
const double fBinX2 = 1;

const int    fNumDetectors = 40;
const int    fNumStrips = 8;
const int    fNumJStrips = 8;
const int    fNumOStrips = 4;
const int    fNumSides = 2;
const int    fNumGates = 2;
const int    fNumLR = 2;

int fDrawingExampleDetectors[] = {20,4,14,34,39};
double fEnergyGateRange[fNumGates][2];

//////////////////////////////////////////////////////////////////////////////////
TH2D* fHistEnergyDetector[fNumSides];
TH2D* fHistEnergyPositionAll;
TH1D* fHistEnergy            [fNumDetectors][fNumSides][fNumStrips];
TH1D* fHistEnergySum         [fNumDetectors][fNumSides][fNumStrips];
TH2D* fHistLeftRight         [fNumDetectors][fNumSides][fNumStrips];
TH2D* fHistEnergyPosition    [fNumDetectors][fNumSides][fNumStrips];
TH2D* fHistEnergyPositionCal [fNumDetectors][fNumSides][fNumStrips];
TH2D* fHistLeftRightGate     [fNumDetectors][fNumSides][fNumStrips][fNumGates];
TH2D* fHistEnergyPositionGate[fNumDetectors][fNumSides][fNumStrips][fNumGates];

//////////////////////////////////////////////////////////////////////////////////
double  fXParameters           [fNumDetectors][fNumSides][fNumStrips][fNumGates][fNumOStrips][fNumLR][3]; // gaussian parameters
TH2D*   fHistEnergyPositionGateOhmic[fNumDetectors][fNumSides][fNumStrips][fNumGates][fNumOStrips];
TH1D*   fHistPositionGateOhmic [fNumDetectors][fNumSides][fNumStrips][fNumGates][fNumOStrips];
TH2D*   fHistLeftRightGateOhmic[fNumDetectors][fNumSides][fNumStrips][fNumGates][fNumOStrips];
TH1D*   fHistLRProjGateOhmic   [fNumDetectors][fNumSides][fNumStrips][fNumGates][fNumOStrips][fNumLR];
TH1D*   fHistLRProjGate        [fNumDetectors][fNumSides][fNumStrips][fNumGates][fNumOStrips];
TH1D*   fHistPositionGate      [fNumDetectors][fNumSides][fNumStrips][fNumGates];
TGraph* fGraphLeftRightX       [fNumDetectors][fNumSides][fNumStrips][fNumGates];
TGraph* fGraphEnergyPositionX  [fNumDetectors][fNumSides][fNumStrips][fNumGates];
TGraphErrors* fGraphEnergyFit  [fNumDetectors][fNumSides][fNumStrips][fNumLR];
TF1*    fFitEnergy             [fNumDetectors][fNumSides][fNumStrips][fNumLR];

//////////////////////////////////////////////////////////////////////////////////
TString fDataDir;
TString fRecoFileName;
TString fCompareFileName;
TString fHistFileName;
TString fC1HistFileName;
TString fCEHistFileName;
TString fC2HistFileName;
TString fC3HistFileName;
TString fDummyFileName;
TString fAllParFileName;

//////////////////////////////////////////////////////////////////////////////////
TCanvas *MakeCanvas(TString name, int size=8);

//////////////////////////////////////////////////////////////////////////////////
class strip_info;
class strip_group;
vector<strip_info> fStripArrayR;
vector<strip_info> fStripArrayS; /// position sensitive strips
vector<strip_group> fExampleGroupArrayR;
vector<strip_group> fExampleGroupArrayS;
vector<strip_group> fAllGroupArrayR;
vector<strip_group> fAllGroupArrayS;

//////////////////////////////////////////////////////////////////////////////////
//SKSiArrayPlane* fStark = nullptr;
ELARK* fStark = nullptr;

//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
// FUNCTIONS BELOW
//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
class strip_info {
    public:
        strip_info(int type_, int det_, int side_, int strip_) { type = type_; det = det_; side = side_; strip = strip_; }
        int type;
        int det;
        int side;
        int strip;
};
class strip_group {
    public:
        strip_group() {}
        vector<strip_info> array;
        int GetDet() {
            if (array.size()==0)
                return -1;
            return array[0].det;
        }
        TCanvas *MakeGroupCanvas(TString name, int npads=-1) {
            TString cvsName = Form("cvs_%s_%d",name.Data(),array.at(0).det);
            if (npads<0) npads = array.size();
            return MakeCanvas(cvsName,npads);
        }
};

//////////////////////////////////////////////////////////////////////////////////
bool ContinueRegardingToDataType(int det)
{
    if (fDataType==7777)
    {
        if (det==12)
            return false;
        return true;
    }
    else if (fDataType==7728)
    {
        if (det==28)
            return false;
        return true;
    }
    else if (fDataType==7718)
    {
        if (det==18)
            return false;
        return true;
    }
    //else if (fDataType==168)
    //{
    //    if (det<12||det>31)
    //        return true;
    //    return false;
    //}
    else if (fDataType==169)
    {
        if (det==0||det==1||det==10||det==11||det==12)
            return false;
        return true;
    }
    else if (fDataType==170)
    {
        if (det==7||det==8||det==9)
            return false;
        return true;
    }
    else if (fDataType==171)
    {
        if (det==1||det==2||det==3)
            return false;
        return true;
    }
    else if (fDataType==2) {
        return false;
    }
    else if (fDataType==0) {
        if (det<12)
            return true;
        return false;
    }
    else if (fDataType==1) {
        if (det>=12)
            return true;
        return false;
    }

    return false;
}

bool MakeRun(int run, int type)
{
    if (run<0) {
        if (gSystem -> Getenv("RUN"))
            run = atoi(gSystem -> Getenv("RUN"));
        else
            run = fRun;
    }
    fRun = run;

    if (gSystem -> Getenv("TEST"))
        fTest = atoi(gSystem -> Getenv("TEST"));

    fDataType = -1;
    if (type<0) {
        if (gSystem -> Getenv("TYPE"))
            type = atoi(gSystem -> Getenv("TYPE"));
        else
            type = fDataType;
    }
    fDataType = type;

    if (fDataType<0) {
        if      (fRun==199)  fDataType = 1;  // 0: 16dE+16E, 1: 12E
        else if (fRun==303)  fDataType = 0;  // 0: 16dE+16E, 1: 12E
        else if (fRun==253)  fDataType = 0;  // 0: 16dE+16E, 1: 12E
        else if (fRun==113)  fDataType = 0;  // 0: 16dE+16E, 1: 12E
        else if (fRun==167)  fDataType = 0;  // 0: 16dE+16E, 1: 12E
        else if (fRun==168)  fDataType = 0;  // 0: 16dE+16E, 1: 12E
        else                 fDataType = fRun;
        //if (fRun==169)  fDataType = 169;  // 0, 1, 10, 11, 12
        //if (fRun==170)  fDataType = 170;  // 7, 8, 9
        //if (fRun==171)  fDataType = 171;  // 1, 2, 3
        //if (fRun==7777) fDataType = 7777;  // 0: 16dE+16E, 1: 12E
        //if (fRun==169) fDataType = 0;  // 0: 16dE+16E, 1: 12E
    }

    if (fDataType<0)
        fDataType = 2;

    fStripArrayR.clear();
    fStripArrayS.clear();
    fExampleGroupArrayR.clear();
    fExampleGroupArrayS.clear();
    //fAllGroupArrayR.clear();
    //fAllGroupArrayS.clear();

    //fNBinA = 500;
    //fNBinE = 500;
    //fNBinX = 500;
    if (fRun==303||fRun==133||fRun==168||fRun==7728) {
        fNBinA = 500;
        fNBinE = 500;
        fNBinX = 500;
    }

    TString recoDir = "data_reco";
    if (fRun>=7700&&fRun<7799)
        recoDir = "data";
    fDataDir = "data";
    //fRecoFileName   = Form("%s/stark_%04d.reco.root",    recoDir.Data(),run);
    fRecoFileName   = Form("%s/stark_%04d.channel.root",    recoDir.Data(),run);
    fHistFileName   = Form("%s/stark_%04d.hist.root",   fDataDir.Data(),run);
    fCompareFileName= Form("%s/stark_%04d.compare.root",fDataDir.Data(),run);
    fC1HistFileName = Form("%s/stark_%04d.hist_c1.root",fDataDir.Data(),run);
    fCEHistFileName = Form("%s/stark_%04d.hist_ce.root",fDataDir.Data(),run);
    fAllParFileName = Form("%s/stark_%04d.par.root",    fDataDir.Data(),run);
    fC2HistFileName = Form("%s/stark_%04d.hist_c2.root",fDataDir.Data(),run);
    fC3HistFileName = Form("%s/stark_%04d.hist_c3.root",fDataDir.Data(),run);
    InitParameter(run, fDataDir, fTest);

    auto file = new TFile(fRecoFileName,"read");
    auto eventTree = (TTree*) file -> Get("event");
    if (eventTree==nullptr) {
        e_error << fRecoFileName << " " << eventTree << endl;
        return false;
    }
    fEntriesCut = eventTree -> GetEntries() / 10000;
    if (fEntriesCut<100) fEntriesCut = 100;
    e_info << "run " << fRun << " " << eventTree -> GetEntries() << " (" << fEntriesCut << ")" << ", Data-type is " << fDataType << endl;
    file -> Close();

    if (fStark==nullptr)
    {
        fStark = new ELARK();
        fStark -> AddPar("config_stark.mac");
        fStark -> Init();
    }

    auto numDetectors = fStark -> GetNumSiDetectors();
    for (auto det=0; det<numDetectors; ++det) {
        auto detector = fStark -> GetSiDetector(det);
        auto detType = detector -> GetDetType();
        auto numJStrips = detector -> GetNumJunctionStrips();
        auto numOStrips = detector -> GetNumOhmicStrips();
        auto numJDirection = detector -> GetNumJunctionDirection();

        if (ContinueRegardingToDataType(det))
            continue;

        {
            strip_group dss_group1;
            for (auto side=0; side<2; ++side)
            {
                if (side==0&&numJDirection==2)
                    continue;
                auto numStrips = (side==0?numJStrips:numOStrips);
                for (auto strip=0; strip<numStrips; ++strip) {
                    fStripArrayS.push_back(strip_info(detType,det,side,strip));
                    dss_group1.array.push_back(strip_info(detType,det,side,strip));
                }
            }
            bool isExample = false; for (auto ex : fDrawingExampleDetectors) if (det==ex) isExample = true;
            if (isExample) fExampleGroupArrayS.push_back(dss_group1);
            fAllGroupArrayS.push_back(dss_group1);
        }

        if (numJDirection==2) {
            strip_group dss_group2;
            int side = 0;
            for (auto strip=0; strip<numJStrips; ++strip) {
                fStripArrayR.push_back(strip_info(detType,det,side,strip));
                dss_group2.array.push_back(strip_info(detType,det,side,strip));
            }
            bool isExample = false; for (auto ex : fDrawingExampleDetectors) if (det==ex) isExample = true;
            if (isExample) fExampleGroupArrayR.push_back(dss_group2);
            fAllGroupArrayR.push_back(dss_group2);
        }
    }

    for (auto gate=0; gate<fNumGates; ++gate)
    {
        double energy = (gate==0?f148GdAlphaEnergy:f241AmAlphaEnergy1);
        fEnergyGateRange[gate][1] = energy + 12*energy*fExpectedResolution;
        fEnergyGateRange[gate][0] = energy - 12*energy*fExpectedResolution;
        if (gate==0) 
            fEnergyGateRange[gate][0] = energy - 24*energy*fExpectedResolution;
        if (gate==fNumGates-1) 
            fEnergyGateRange[gate][1] = energy + 24*energy*fExpectedResolution;
    }

    return true;
}

int GetNumOhmicStrips(int det)
{
    auto detector = fStark -> GetSiDetector(det);
    if (detector==nullptr)
        return 0;
    return detector -> GetNumOhmicStrips();
}

TString MakeCompareFileName(int calRun)
{
    fCompareFileName = Form("%s/stark_%04d.compare_%d.root",fDataDir.Data(),fRun,calRun);
    return fCompareFileName;
}

TString MakeHistName(TString xname, TString yname, int det, int side, int strip, int gate)
{
    const char* mainName  = yname.IsNull() ? xname.Data() : Form("%s_%s",xname.Data(),yname.Data());
    const char* detName   = ( det  <0 ? "" : Form("_d%d",det   ) );
    const char* sideName  = ( side <0 ? "" : (side==0?"_j":"_o") );
    const char* stripName = ( strip<0 ? "" : Form("_s%d",strip ) );
    const char* gateName  = ( gate <0 ? "" : Form("_g%d",gate  ) );
    TString name = Form("hist%s%s%s%s_%s",detName,sideName,stripName,gateName,mainName);
    //TString name = Form("%s%s%s%s%s",name,detName,sideName,stripName,gateName);
    return name;
}

TString MakeHistTitle(TString xname, TString yname, int det, int side, int strip, int gate)
{
    auto detector = fStark -> GetSiDetector(det);
    const char* runTitle   = (fRun <0 ? "" : Form("[%04d] ",fRun));
    const char* detTitle   = (det  <0 ? "" : Form("detector %d (%s)",det, detector->GetDetTypeName().Data()));
    const char* sideTitle  = (side <0 ? "" : (side==0?", junction":", ohmic"));
    const char* stripTitle = (strip<0 ? "" : Form(", strip %d",strip ));
    const char* gateTitle  = (gate <0 ? "" : Form(", gate %d",gate  ));
    if (TString(detTitle).IsNull()) sideTitle = TString(TString(sideTitle)(2,TString(sideTitle).Sizeof()-3)).Data();
    TString title = Form("%s%s%s%s%s;%s;%s",runTitle,detTitle,sideTitle,stripTitle,gateTitle,xname.Data(),yname.Data());
    return title;
}

TH1D* MakeHist1(TString xname, TString yname, int det, int side, int strip, int gate, int nx, double x1, double x2)
{
    if (nx==0) {
        if (xname.Index("c0")==0 || xname.Index("c1")==0 || xname.Index("c2")==0) {
            nx = 400;
            x1 = 0;
            x2 = 10;
        }
        else {
            if      (det<32  && side==1) { nx = 200; x1 = 200; x2 = 1200; }
            else if (det>=32 && side==1) { nx = 200; x1 = 500; x2 = 2200; }
            else if (det==39 && side==0 && (strip==2 || strip==3) ) { nx = 200; x1 = 200; x2 = 2000; }
            else if (det>=32 && side==0) { nx = 200; x1 = 200; x2 = 2800; }
            else                         { nx = 200; x1 = 200; x2 = 4000; }
        }
    }

    auto hist = new TH1D(MakeHistName(xname,yname,det,side,strip,gate),MakeHistTitle(xname,yname,det,side,strip,gate),nx,x1,x2);
    //hist -> SetFillColor(29);
    return hist;
}

TH2D* MakeHist2(TString xname, TString yname, int det, int side, int strip, int gate, int nx, double x1, double x2, int ny, int y1, int y2)
{
    auto hist = new TH2D(MakeHistName(xname,yname,det,side,strip,gate),MakeHistTitle(xname,yname,det,side,strip,gate),nx,x1,x2,ny,y1,y2);
    return hist;
}

TCanvas *MakeCanvas(TString cvsName, int npads)
{
    TCanvas *cvs;
    if (npads==2) {
        cvs = LKPainter::GetPainter() -> CanvasResize(cvsName,600,500,0.9);
        cvs -> Divide(1,2,0.002,0.002);
    }
    else if (npads==4) {
        cvs = LKPainter::GetPainter() -> CanvasResize(cvsName,550,500,0.7);
        cvs -> Divide(2,2,0.002,0.002);
    }
    else if (npads==8||npads==9) {
        cvs = LKPainter::GetPainter() -> CanvasResize(cvsName,550,500,0.95);
        cvs -> Divide(3,3,0.002,0.002);
    }
    else if (npads==20) {
        cvs = LKPainter::GetPainter() -> CanvasResize(cvsName,700,500,1);
        cvs -> Divide(5,4,0.002,0.002);
    }
    else {
        cvs = LKPainter::GetPainter() -> CanvasDefault(cvsName);//,500,500,0.95);
    }
    for (auto i=1; i<=npads; ++i) {
        cvs -> cd(i) -> SetMargin(0.15, 0.12, 0.12, 0.10);
        cvs -> cd(i) -> SetGridx();
        cvs -> cd(i) -> SetGridy();
    }
    cvs -> cd();
    return cvs;
}

bool IsPositionSensitiveStrip(int det, int side)
{
    auto detector = fStark -> GetSiDetector(det);
    if (detector->GetNumJunctionDirection()==2 && side==0)
        return true;
    return false;
}

bool IsPositionSensitiveStrip(LKSiChannel* channel)
{
    return IsPositionSensitiveStrip(channel->GetDetID(),channel->GetSide());
}

void GetFinalHistograms()
{
    TString calName = "c3_";
    auto fileHist = new TFile(fC3HistFileName,"read");
    for (auto dss : fStripArrayS)
    {
        fHistEnergy[dss.det][dss.side][dss.strip] = (TH1D*) fileHist ->Get(MakeHistName(calName+"energy","",dss.det,dss.side,dss.strip,-1));
    }
    for (auto dss : fStripArrayR)
    {
        fHistEnergySum     [dss.det][dss.side][dss.strip] = (TH1D*) fileHist ->Get(MakeHistName(calName+"esum",             "",dss.det,dss.side,dss.strip,-1));
        fHistLeftRight     [dss.det][dss.side][dss.strip] = (TH2D*) fileHist ->Get(MakeHistName(calName+"left",calName+"right",dss.det,dss.side,dss.strip,-1));
        fHistEnergyPosition[dss.det][dss.side][dss.strip] = (TH2D*) fileHist ->Get(MakeHistName(calName+"rpos",calName+"esum", dss.det,dss.side,dss.strip,-1));
    }
    fHistEnergyDetector[0] = (TH2D*) fileHist ->Get(MakeHistName(calName+"esum","det",-1,0,-1,-1));
    fHistEnergyDetector[1] = (TH2D*) fileHist ->Get(MakeHistName(calName+"esum","det",-1,1,-1,-1));

    cout << "fHistEnergy" << endl;
    cout << "fHistEnergySum" << endl;
    cout << "fHistLeftRight" << endl;
    cout << "fHistEnergyPosition" << endl;
    cout << "fHistEnergyDetector" << endl;
    cout << "fHistEnergyDetector" << endl;
}

TH2D* MakeHitPatternHist(int jo, TString calName, int nn, double e1, double e2) // jo=0: junction, jo=1: ohmic
{
    TH2D* hist;
    if (nn==0) {
        nn = 2*fNBinEFix;
        e1 = fBinE1;
        e2 = 1.5*fBinE2;
    }
    if (jo==0) hist = MakeHist2("det_j",calName+"esum",-1,0,-1,-1,fNumDetectors*fNumJStrips,0,fNumDetectors*fNumJStrips,nn,e1,e2);
    if (jo==1) hist = MakeHist2("det_o",calName+"esum",-1,1,-1,-1,fNumDetectors*fNumOStrips,0,fNumDetectors*fNumOStrips,nn,e1,e2);
    if (jo==2) hist = MakeHist2("det",calName,-1,1,-1,-1,fNumDetectors,0,fNumDetectors,nn,e1,e2);
    hist -> GetXaxis() -> SetTitle("");
    hist -> GetXaxis() -> SetNdivisions(400+fNumDetectors);
    if (jo!=2) 
        hist -> GetXaxis() -> SetNdivisions(fNumDetectors);
    //hist -> GetXaxis() -> SetLabelSize(0.023);
    hist -> GetXaxis() -> SetLabelSize(0.030);
    for (auto det=0; det<fNumDetectors; ++det) {
        auto detector = fStark -> GetSiDetector(det);
        auto name = detector -> GetDetTypeName();
        auto ring = detector -> GetLayer();
        TString sring = "dE"; if (ring==1) sring = "E"; if (ring==2) sring = "16E";
        if (jo==0) hist -> GetXaxis() -> SetBinLabel(det*fNumJStrips+1,Form("%d (%s,%s)",det,name.Data(),sring.Data()));
        if (jo==1) hist -> GetXaxis() -> SetBinLabel(det*fNumOStrips+1,Form("%d (%s,%s)",det,name.Data(),sring.Data()));
        if (jo==2)  {
            hist -> GetXaxis() -> SetBinLabel(det*1+1,Form("%d (%s,%s)",det,name.Data(),sring.Data()));
            hist -> GetXaxis() -> ChangeLabel(det*1+1,270);
        }
    }
    return hist;
}

TH1D* MakeDet1DHist(int jo, TString calName)
{
    TH1D* hist;
    if (jo==0) hist = MakeHist1("det_j",calName,-1,0,-1,-1,fNumDetectors*fNumJStrips,0,fNumDetectors*fNumJStrips);
    if (jo==1) hist = MakeHist1("det_o",calName,-1,1,-1,-1,fNumDetectors*fNumOStrips,0,fNumDetectors*fNumOStrips);
    if (jo==2) hist = MakeHist1("det",calName,-1,1,-1,-1,fNumDetectors,0,fNumDetectors);
    hist -> GetXaxis() -> SetTitle("");
    hist -> GetXaxis() -> SetNdivisions(400+fNumDetectors);
    if (jo!=2) 
        hist -> GetXaxis() -> SetNdivisions(fNumDetectors);
    //hist -> GetXaxis() -> SetLabelSize(0.023);
    hist -> GetXaxis() -> SetLabelSize(0.030);
    for (auto det=0; det<fNumDetectors; ++det) {
        auto detector = fStark -> GetSiDetector(det);
        auto name = detector -> GetDetTypeName();
        auto ring = detector -> GetLayer();
        TString sring = "dE"; if (ring==1) sring = "E"; if (ring==2) sring = "16E";
        if (jo==0) hist -> GetXaxis() -> SetBinLabel(det*fNumJStrips+1,Form("%d (%s,%s)",det,name.Data(),sring.Data()));
        if (jo==1) hist -> GetXaxis() -> SetBinLabel(det*fNumOStrips+1,Form("%d (%s,%s)",det,name.Data(),sring.Data()));
        hist -> GetXaxis() -> ChangeLabel(det*fNumOStrips+1,90);
    }
    return hist;
}

#endif
