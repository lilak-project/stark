#include <fstream>
using namespace std;

#include "SKSiArrayPlane.h"
#include "LKSiDetector.h"
#include "LKSiChannel.h"
#include "TStyle.h"
#include "TH2Poly.h"
#include "LKPainter.h"

ClassImp(SKSiArrayPlane)

SKSiArrayPlane::SKSiArrayPlane()
    :SKSiArrayPlane("SKSiArrayPlane","General mapping class for si array detector")
{
}

SKSiArrayPlane::SKSiArrayPlane(const char *name, const char *title)
    :LKEvePlane(name, title)
{
    fDetectorArray = new TObjArray();
    fChannelArray = new TObjArray();
    fGSelJOSideDisplay = new TGraph();
    fUserDrawingArrayCollection = new TObjArray();

    fPaletteNumber = 1;
    fCtrlBinTextSize = 8.0;
    fCtrlLabelSize = 0.25;
    //fCtrlLabelOffset = 0.10;
}

bool SKSiArrayPlane::Init()
{
    LKEvePlane::Init();

    if (fDetName.IsNull())
        fDetName = "stark";

    fUserDrawingArrayIndex = new int***[4];
    for(int i=0; i<4; ++i) {
        fUserDrawingArrayIndex[i] = new int**[8];
        for(int j=0; j<8; ++j) {
            fUserDrawingArrayIndex[i][j] = new int*[fMaxDetectors];
            for(int k=0; k<fMaxDetectors; ++k) {
                fUserDrawingArrayIndex[i][j][k] = new int[2];
                for(int l=0; l<2; ++l) {
                    fUserDrawingArrayIndex[i][j][k][l] = -1;
                }
            }
        }
    }

    fUserDrawingLeastNDraw = new int***[4];
    for(int i=0; i<4; ++i) {
        fUserDrawingLeastNDraw[i] = new int**[8];
        for(int j=0; j<8; ++j) {
            fUserDrawingLeastNDraw[i][j] = new int*[fMaxDetectors];
            for(int k=0; k<fMaxDetectors; ++k) {
                fUserDrawingLeastNDraw[i][j][k] = new int[2];
                for(int l=0; l<2; ++l) {
                    fUserDrawingLeastNDraw[i][j][k][l] = -1;
                }
            }
        }
    }

    for(int i=0; i<4; ++i)
        for(int j=0; j<7; ++j)
            fUserDrawingName[i][j] = "";

    fMapCAACToChannelIndex = new int***[fNumCobo];
    for(int i=0; i<fNumCobo; ++i) {
        fMapCAACToChannelIndex[i] = new int**[fNumAsad];
        for(int j=0; j<fNumAsad; ++j) {
            fMapCAACToChannelIndex[i][j] = new int*[fNumAget];
            for(int k=0; k<fNumAget; ++k) {
                fMapCAACToChannelIndex[i][j][k] = new int[fNumChan];
                for(int l=0; l<fNumChan; ++l) {
                    fMapCAACToChannelIndex[i][j][k][l] = -1;
                }
            }
        }
    }

    fMapCAACToDetectorIndex = new int***[fNumCobo];
    for(int i=0; i<fNumCobo; ++i) {
        fMapCAACToDetectorIndex[i] = new int**[fNumAsad];
        for(int j=0; j<fNumAsad; ++j) {
            fMapCAACToDetectorIndex[i][j] = new int*[fNumAget];
            for(int k=0; k<fNumAget; ++k) {
                fMapCAACToDetectorIndex[i][j][k] = new int[fNumChan];
                for(int l=0; l<fNumChan; ++l) {
                    fMapCAACToDetectorIndex[i][j][k][l] = -1;
                }
            }
        }
    }

    fPar -> Require(fDetName+"/Mapping","{lilak_common}/stark_channel_mapping_RAON2024_40Arp.txt","cobo asad aget chan chan2 det detID JO strip LR theta phi","",1);
    fPar -> Require(fDetName+"/DetectorPar","{lilak_common}/stark_detector_par_RAON2024_40Arp.txt","detName detType(si,...) numSides numJunctionStrips numOhmicStrips useJunctionLR useOhmicLR","",2);

    fPar -> UpdatePar(fMappingFileName,fDetName+"/Mapping");
    fPar -> UpdatePar(fDetectorParName,fDetName+"/DetectorPar");

    ifstream fileDetector(fDetectorParName);
    if (!fileDetector.is_open()) {
        lk_error << "Cannot open " << fDetectorParName << endl;
        return false;
    }
    lk_info << "detector parameter: " << fDetectorParName << endl;

    fMaxLayerIndex = -INT_MAX;
    fMinPhi = DBL_MAX;
    fMaxPhi = -DBL_MAX;

    TString detName;
    int cobo, asad, aget, chan, chan2, detID, strip, side, lr;
    int numSides, numJunctionStrips, numOhmicStrips, layer, ringIndex, useJunctionLR, useOhmicLR;
    double detDistance, detRadius, phi0, detWidth, detHeight, detThickness, tta, phi;
    while (fileDetector >> detName >> detID
            >> numSides >> numJunctionStrips >> numOhmicStrips >> useJunctionLR >> useOhmicLR
            >> detDistance >> layer >> detRadius >> ringIndex >> phi0
            >> detWidth >> detHeight >> detThickness)
    {
        //lk_debug << detName << " " << detID
        //    << " " << numSides << " " << numJunctionStrips << " " << numOhmicStrips << " " << useJunctionLR << " " << useOhmicLR
        //    << " " << detDistance << " " << layer << " " << detRadius << " " << ringIndex << " " << phi0
        //    << " " << detWidth << " " << detHeight << " " << detThickness << endl;
        auto siDetector = new LKSiDetector();
        if (fDetectorTypeArray.CheckPar(detName)==false)
            fDetectorTypeArray.AddPar(detName, 1);
        int detTypeIndex = fDetectorTypeArray.GetParIndex(detName);
        double dphi = TMath::ATan2(detWidth/2.,detRadius)*TMath::RadToDeg();
        double phi1 = phi0 - dphi;
        double phi2 = phi0 + dphi;
        double tta1 = TMath::ATan2(detRadius,detDistance)*TMath::RadToDeg();
        double tta2 = TMath::ATan2(detRadius,detDistance+detHeight)*TMath::RadToDeg();
        int detIndex = fDetectorArray -> GetEntries();
        siDetector -> SetSiType(detName, detTypeIndex, detIndex, detID, numSides, numJunctionStrips, numOhmicStrips, useJunctionLR, useOhmicLR);
        //siDetector -> SetSiPosition(TVector3(), layer, ringIndex, phi0, tta0);
        siDetector -> SetSiPosition(TVector3(), layer, ringIndex, phi1, phi2, tta1, tta2);
        double layer1 = fLayerYScale * (layer + fLayerYOffset);
        double layer2 = fLayerYScale * (layer + 1 - fLayerYOffset);
        TString name = Form("hist_%s_%d_junction",detName.Data(),detID);
        TString title = Form("%s(%d) Junction (%d,%d)",detName.Data(),detID,numJunctionStrips,(useJunctionLR?2:1));
        siDetector -> CreateHistJunction(name, title, phi1, phi2, layer1, layer2, "!axis");
        name = Form("hist_%s_%d_ohmic",detName.Data(),detID);
        title = Form("%s(%d) Ohmic (%d,%d)",detName.Data(),detID,(useOhmicLR?2:1),numOhmicStrips);
        siDetector -> CreateHistOhmic(name, title, phi1, phi2, layer1, layer2, "!axis");
        fDetectorArray -> Add(siDetector);
        //siDetector -> Print();
        auto histJ = siDetector -> GetHistJunction();
        auto histO = siDetector -> GetHistOhmic();
        histJ -> SetMarkerSize(fCtrlBinTextSize*0.3);
        histO -> SetMarkerSize(fCtrlBinTextSize*0.3);

        if (fMaxLayerIndex < layer) fMaxLayerIndex = layer;
        if (fMinPhi > phi1) fMinPhi = phi1;
        if (fMaxPhi < phi2) fMaxPhi = phi2;
    }

    ifstream fileCAACMap(fMappingFileName);
    if (!fileCAACMap.is_open()) {
        lk_error << "Cannot open " << fMappingFileName << endl;
        return false;
    }
    lk_info << "channel mapping: " << fMappingFileName << endl;

    int globalChannelIndex = 0;
    int currentDetID = -1;
    int localJID = 0;
    int localOID = 0;
    while (fileCAACMap >> cobo >> asad >> aget >> chan >> chan2 >> detName >> detID >> side >> strip >> lr >> tta >> phi)
    {
        side = side - 1;
        if (currentDetID==detID) {
            if (side==0) localJID++;
            if (side==1) localOID++;
        }
        else {
            currentDetID = detID;
            localJID = 0;
            localOID = 0;
        }
        fMapCAACToChannelIndex[cobo][asad][aget][chan] = globalChannelIndex;
        int detTypeIndex = fDetectorTypeArray.GetParIndex(detName);
        auto siChannel = new LKSiChannel();
        fChannelArray -> Add(siChannel);
        siChannel -> SetCobo(cobo);
        siChannel -> SetAsad(asad);
        siChannel -> SetAget(aget);
        siChannel -> SetChan(chan);
        siChannel -> SetDetType(detTypeIndex);
        siChannel -> SetDetID(detID);
        siChannel -> SetChannelID(globalChannelIndex);
        if (side==0) siChannel -> SetLocalID(localJID);
        if (side==1) siChannel -> SetLocalID(localOID);
        siChannel -> SetSide(side);
        siChannel -> SetStrip(strip);
        siChannel -> SetDirection(lr);
        siChannel -> SetTheta1(tta);
        siChannel -> SetPhi1(phi);
        fMapCAACToDetectorIndex[cobo][asad][aget][chan] = detID;

        //auto siDetector = (LKSiDetector*) fDetectorArray -> At(detTypeIndex);
        //siDetector -> SetChannel(channel);

        ++globalChannelIndex;
    }

    lk_info << globalChannelIndex << " channels are mapped!" << endl;

    return true;
}

void SKSiArrayPlane::Print(Option_t *option) const
{
    lk_info << "Silicon array " << fName << endl;

    int numDetectors = fDetectorArray -> GetEntries();
    for (auto iDetector=0; iDetector<numDetectors; ++iDetector)
    {
        auto siDetector = (LKSiDetector*) fDetectorArray -> At(iDetector);
        siDetector -> Print();
    }
}

TCanvas *SKSiArrayPlane::GetCanvas(Option_t *option)
{
    if (fCanvas==nullptr)
    {
        double y1 = 0.00;
        double y2 = 0.10;
        double y3 = 0.20;
        double y4 = 0.30;
        double y5 = 0.65;
        double y6 = 1.00;

        double u1 = 0.10;
        double u2 = (1.-u1)/2.+u1;
        double u3 = 0.80;

        double w1 = 0.10;
        double w2 = (1.-w1)*1/3.+w1;
        double w3 = (1.-w1)*2/3.+w1;

        double x2 = 0.45;

        fCanvas = LKPainter::GetPainter() -> CanvasResize("stark",fDXCanvas,fDYCanvas,0.95);

        fCanvas -> cd();
        fPadControlEvent1    = new TPad("SKSiArrayPlanePadCE1", "", 0.0 , y1, x2, y2);
        fPadControlEvent2    = new TPad("SKSiArrayPlanePadCE2", "", 0.0 , y2, x2, y3);
        fPadCtrlUserDrawing  = new TPad("SKSiArrayPlanePadUM",  "", 0.0 , y3, x2, y4);
        fPadJOSideDisplay[0] = new TPad("SKSiArrayPlanePadJD",  "", 0.0 , y5, 0.5*(x2-0.0),y6);
        fPadJOSideDisplay[1] = new TPad("SKSiArrayPlanePadOD",  "", 0.5*(x2-0.0), y5, x2, y6);
        fPadEventDisplay1    = new TPad("SKSiArrayPlanePadED1", "", 0.0 , y4, x2, y5);

        fPadControlDataDP      = new TPad("SKSiArrayPlanePadDDC0","", x2,  0,   1,  u1);
        fPadDataDisplayFull    = new TPad("SKSiArrayPlanePadDDF0","", x2,  u1,  1,  u3);
        fPadDataDisplayTwo[0]  = new TPad("SKSiArrayPlanePadDDT1","", x2,  u1,  1,  u2);
        fPadDataDisplayTwo[1]  = new TPad("SKSiArrayPlanePadDDT2","", x2,  u2,  1,  1 );

        fPadDataDisplayThree[0]  = new TPad("SKSiArrayPlanePadDD31","", x2,  w1,  1,  w2);
        fPadDataDisplayThree[1]  = new TPad("SKSiArrayPlanePadDD32","", x2,  w2,  1,  w3);
        fPadDataDisplayThree[2]  = new TPad("SKSiArrayPlanePadDD32","", x2,  w3,  1,  1 );

        fPadControlEvent1       -> SetMargin(0.02,0.02,0.22,0.02);
        fPadControlEvent2       -> SetMargin(0.02,0.02,0.22,0.02);
        fPadCtrlUserDrawing     -> SetMargin(0.02,0.02,0.22,0.02);
        fPadEventDisplay1       -> SetMargin(0.02,0.02,0.10,0.02);
        fPadJOSideDisplay[0]    -> SetMargin(0.02,0.13,0.02,0.07);
        fPadJOSideDisplay[1]    -> SetMargin(0.02,0.13,0.02,0.07);
        fPadControlDataDP       -> SetMargin(0.02,0.02,0.22,0.02);
        fPadDataDisplayFull     -> SetMargin(0.10,0.12,0.10,0.10);
        fPadDataDisplayTwo[0]   -> SetMargin(0.10,0.12,0.10,0.10);
        fPadDataDisplayTwo[1]   -> SetMargin(0.10,0.12,0.10,0.10);
        fPadDataDisplayThree[0] -> SetMargin(0.10,0.12,0.10,0.10);
        fPadDataDisplayThree[1] -> SetMargin(0.10,0.12,0.10,0.10);
        fPadDataDisplayThree[2] -> SetMargin(0.10,0.12,0.10,0.10);

        fPadControlDataDP    -> SetGrid(1);
        fPadCtrlUserDrawing  -> SetGrid(1);
        fPadJOSideDisplay[0] -> SetGrid(1,1);
        fPadJOSideDisplay[1] -> SetGrid(1,1);

        fPadJOSideDisplay[0] -> Draw();
        fPadJOSideDisplay[1] -> Draw();
        fPadControlEvent1    -> Draw();
        fPadControlEvent2    -> Draw();
        fPadEventDisplay1    -> Draw();
        fPadCtrlUserDrawing  -> Draw();
        fPadControlDataDP    -> Draw();

        int countPad = 0;
        for (auto iy=0; iy<fNumDDY; ++iy) {
            for (auto ix=0; ix<fNumDDX; ++ix) {
                double a1 = x2 + (1-x2)*(ix+0)*(1./fNumDDX);
                double a2 = x2 + (1-x2)*(ix+1)*(1./fNumDDX);
                double b1 = u1 + (1-u1)*(iy+0)*(1./fNumDDY);
                double b2 = u1 + (1-u1)*(iy+1)*(1./fNumDDY);
                fPadDataDisplaySmall[countPad] = new TPad(Form("SKSiArrayPlanePadDDS%d",countPad),"",a1,b1,a2,b2);
                fPadDataDisplaySmall[countPad] -> Draw();
                countPad++;
            }
        }

        fCanvas -> Modified();
        fCanvas -> Update();

        AddInteractivePad(fPadEventDisplay1);
        AddInteractivePad(fPadControlEvent1);
        AddInteractivePad(fPadControlEvent2);
        AddInteractivePad(fPadJOSideDisplay[0]);
        AddInteractivePad(fPadJOSideDisplay[1]);
        AddInteractivePad(fPadCtrlUserDrawing);
        AddInteractivePad(fPadControlDataDP);
    }

    return fCanvas;
}

TH2* SKSiArrayPlane::GetHist(Option_t *option)
{
    LKEvePlane::GetHist(option);
    GetHistUserDrawing();
    GetHistControlDataDP();
    return (TH2*) fHistEventDisplay1;
}

TH2* SKSiArrayPlane::GetHistControlDataDP(Option_t *option)
{
    if (fHistControlDataDP==nullptr)
    {
        gStyle -> SetHistMinimumZero(); // this will draw text even when content is 0

        fHistControlDataDP = new TH2D("LKEvePlane_DataDisplay1","",3,0,3,1,0,1);
        fHistControlDataDP -> SetStats(0);
        fHistControlDataDP -> GetXaxis() -> SetTickSize(0);
        fHistControlDataDP -> GetYaxis() -> SetTickSize(0);
        fHistControlDataDP -> GetYaxis() -> SetBinLabel(1,"");
        fHistControlDataDP -> GetXaxis() -> SetBinLabel(1,"<");
        fHistControlDataDP -> GetXaxis() -> SetBinLabel(2,"/");
        fHistControlDataDP -> GetXaxis() -> SetBinLabel(3,">");
        fHistControlDataDP -> GetXaxis() -> SetLabelSize(fCtrlLabelSize);
        fHistControlDataDP -> SetMinimum(0);
        fHistControlDataDP -> SetMarkerSize(fCtrlBinTextSize);
        if (fPaletteNumber==0)
            fHistControlDataDP -> SetMarkerColor(kBlack);
        fBinCtrlPrevPage = fHistControlDataDP -> GetBin(1,1);
        fBinCtrlCurrPage = fHistControlDataDP -> GetBin(2,1);
        fBinCtrlNextPage = fHistControlDataDP -> GetBin(3,1);
        fHistControlDataDP -> GetBinContent(fBinCtrlPrevPage,-1);
        fHistControlDataDP -> GetBinContent(fBinCtrlCurrPage,0);
        fHistControlDataDP -> GetBinContent(fBinCtrlNextPage,-1);
    }

    return fHistControlDataDP;
}

TH2* SKSiArrayPlane::GetHistUserDrawing(Option_t *option)
{
    if (fHistCtrlUserDrawing==nullptr)
    {
        gStyle -> SetHistMinimumZero(); // this will draw text even when content is 0

        fHistCtrlUserDrawing = new TH2D("SKSiArrayPlane_UserDrawing","",8,0,8,1,0,1);
        fHistCtrlUserDrawing -> SetStats(0);
        fHistCtrlUserDrawing -> GetXaxis() -> SetTickSize(0);
        fHistCtrlUserDrawing -> GetYaxis() -> SetTickSize(0);
        fHistCtrlUserDrawing -> GetYaxis() -> SetBinLabel(1,"");
        fHistCtrlUserDrawing -> GetXaxis() -> SetBinLabel(1,"Menu");
        fHistCtrlUserDrawing -> GetXaxis() -> SetLabelSize(fCtrlLabelSize);
        //fHistCtrlUserDrawing -> GetXaxis() -> SetLabelOffset(fCtrlLabelOffset);
        fHistCtrlUserDrawing -> SetMarkerSize(fCtrlBinTextSize);
        if (fPaletteNumber==0)
            fHistCtrlUserDrawing -> SetMarkerColor(kBlack);
        fHistCtrlUserDrawing -> SetMinimum(0);
        fHistCtrlUserDrawing -> SetMaximum(1);
        UpdateUserDrawing();
    }

    return fHistCtrlUserDrawing;
}

bool SKSiArrayPlane::AddUserDrawings(TString label, int detID, int joID, TObjArray* userDrawingArray, int leastNDraw)
{
    if (detID<0 && joID<0)
        lk_info << "adding: " << label << " det=" << detID << "(" << (joID==0?"Junction":"Ohmic") << ") containing " << userDrawingArray->GetEntries() << " drawings" << endl;
    else if (joID<0)
        lk_info << "adding: " << label << " det=" << detID << " containing " << userDrawingArray->GetEntries() << " drawings" << endl;
    else
        lk_info << "adding: " << label << " containing " << userDrawingArray->GetEntries() << " drawings" << endl;

    if (label.Index(" ")>=0)
        label.ReplaceAll(" ","_");

    bool isNewLabel = false;
    bool breakAll = false;

    int ucTab = 0;
    int ucBin = 0;

    for (ucTab=0; ucTab<4; ++ucTab)
    {
        for (ucBin=0; ucBin<7; ++ucBin)
        {
            if (fUserDrawingName[ucTab][ucBin]=="") {
                isNewLabel = true;
                breakAll = true;
                break;
            }
            if (fUserDrawingName[ucTab][ucBin]==label) {
                breakAll = true;
                break;
            }
        }
        if (breakAll)
            break;
    }
    if (isNewLabel) {
        fUserDrawingName[ucTab][ucBin] = label;
        lk_info << "new label: " << ucTab << " " << ucBin << " " << fUserDrawingName[ucTab][ucBin] << endl;
        if (label=="waveform")
            fSelUDBin = ucBin+2;
        else if (fSelUDBin<0)
            fSelUDBin = ucBin+2;
        fNumUDLabels++;
    }

    auto haID = fUserDrawingArrayCollection -> GetEntries();
    fUserDrawingArrayCollection -> Add(userDrawingArray);

    if (detID<0 && joID<0)
    {
        for(int k=0; k<fMaxDetectors; ++k) {
            for(int l=0; l<2; ++l) {
                fUserDrawingArrayIndex[ucTab][ucBin][k][l] = haID;
                fUserDrawingLeastNDraw[ucTab][ucBin][k][l] = leastNDraw;
            }
        }
    }
    else if (detID>=0 && joID<0) {
        for(int l=0; l<2; ++l) {
            fUserDrawingArrayIndex[ucTab][ucBin][detID][l] = haID;
            fUserDrawingLeastNDraw[ucTab][ucBin][detID][l] = leastNDraw;
        }
    }
    else {
        fUserDrawingArrayIndex[ucTab][ucBin][detID][joID] = haID;
        fUserDrawingLeastNDraw[ucTab][ucBin][detID][joID] = leastNDraw;
    }

    return true;
}

void SKSiArrayPlane::UpdateUserDrawing()
{
    if (UpdateFlag[kUpdateUserDrawing]) return;
    else UpdateFlag[kUpdateUserDrawing] = true;

    if (fSignalUDTabChange)
    {
        fSignalUDTabChange = false;
        fHistCtrlUserDrawing -> SetBinContent(fHistCtrlUserDrawing->GetBin(1,1),fSelUDTab);
        for (int ucTab=0; ucTab<4; ++ucTab) {
            for (int ucBin=0; ucBin<7; ++ucBin) {
                TString label = fUserDrawingName[ucTab][ucBin];
                if (fSelUDTab==ucTab) {
                    fHistCtrlUserDrawing -> GetXaxis() -> SetBinLabel(ucBin+2,label);
                    fHistCtrlUserDrawing -> SetBinContent(fHistCtrlUserDrawing->GetBin(ucBin+2,1),fSelUDPage);
                }
            }
        }
        fPadCtrlUserDrawing -> cd();
        fHistCtrlUserDrawing -> Draw("col text");
    }
    else if (fSignalUDBinChange && fSelUDBin>=2)
    {
        fSignalUDBinChange = false;

        int selJOID = fSelJOID;
        if (selJOID<0)
            selJOID = 0;
        fSelUDArrayID = fUserDrawingArrayIndex[fSelUDTab][fSelUDBin-2][fSelDetID][selJOID];
        fSelUDLeastNDraw = fUserDrawingLeastNDraw[fSelUDTab][fSelUDBin-2][fSelDetID][selJOID];
        auto selUDArray = (TObjArray*) fUserDrawingArrayCollection -> At(fSelUDArrayID);

        if (selUDArray!=nullptr)// && selUDArray!=fSelUDArray)
        {
            fSelUDArray = selUDArray;
            fNumDrawingsInArray = fSelUDArray -> GetEntries();
            fNumSelUDGroup = fNumDrawingsInArray / fNumDataDisplays;
            fSelUDPage = 0;

            {
                for (auto i=0; i<10; ++i)
                    for (auto j=0; j<16; ++j)
                        fUDGoodIndex[i][j] = -1;
                fNumGoodDrawings = 0;
                fNumUDPage = 1;
                for (auto iDrawing=0; iDrawing<fNumDrawingsInArray; ++iDrawing)
                {
                    bool good = false;
                    auto obj = fSelUDArray -> At(iDrawing);
                    if (obj -> InheritsFrom(TH1::Class())) {
                        auto hist = (TH1*) obj;
                        if (hist -> GetEntries()>0)
                            good = true;
                    }
                    if (good) {
                        fUDGoodIndex[fNumUDPage-1][fNumGoodDrawings] = iDrawing;
                        fNumGoodDrawings++;
                        if (fNumGoodDrawings>=fNumDataDisplays) {
                            fNumGoodDrawings = 0;
                            fNumUDPage++;
                        }
                    }
                }
            }

            fHistCtrlUserDrawing -> SetBinContent(fHistCtrlUserDrawing->GetBin(2,1),0);
            fHistCtrlUserDrawing -> SetBinContent(fHistCtrlUserDrawing->GetBin(3,1),0);
            fHistCtrlUserDrawing -> SetBinContent(fHistCtrlUserDrawing->GetBin(4,1),0);
            fHistCtrlUserDrawing -> SetBinContent(fHistCtrlUserDrawing->GetBin(5,1),0);
            fHistCtrlUserDrawing -> SetBinContent(fHistCtrlUserDrawing->GetBin(6,1),0);
            fHistCtrlUserDrawing -> SetBinContent(fHistCtrlUserDrawing->GetBin(7,1),0);
            fHistCtrlUserDrawing -> SetBinContent(fHistCtrlUserDrawing->GetBin(fSelUDBin,1),1);
            fPadCtrlUserDrawing -> Modified();
            fPadCtrlUserDrawing -> Update();

            UpdateDataDisplays();
        }
    }
}

TH2* SKSiArrayPlane::GetHistEventDisplay1(Option_t *option)
{
    if (fHistEventDisplay1==nullptr)
    {
        fLayerYScale = 100;
        auto histDetectors = new TH2Poly("LKSiArrayHistED1","",fMinPhi-0.02*(fMaxPhi-fMinPhi),fMaxPhi+0.02*(fMaxPhi-fMinPhi),0,fLayerYScale*(fMaxLayerIndex+1));
        fHistEventDisplay1 = histDetectors;
        int numDetectors = fDetectorArray -> GetEntries();
        int maxBins = numDetectors + 10;
        fMapBinToDetector = new int[maxBins]; for (auto i=0; i<maxBins; ++i) fMapBinToDetector[i] = -1;
        for (auto iDetector=0; iDetector<numDetectors; ++iDetector)
        {
            auto siDetector = (LKSiDetector*) fDetectorArray -> At(iDetector);
            auto layer = siDetector -> GetLayer();
            auto phi1 = siDetector -> GetPhi1();
            auto phi2 = siDetector -> GetPhi2();
            double layer1 = fLayerYScale * (layer + fLayerYOffset);
            double layer2 = fLayerYScale * (layer + 1 - fLayerYOffset);
            auto bin = histDetectors -> AddBin(phi1, layer1, phi2, layer2);
            histDetectors -> SetBinContent(bin,iDetector);
            fMapBinToDetector[bin] = iDetector;
            //histDetectors -> SetBinContent(bin,iDetector);
        }

        fHistEventDisplay1 -> SetMarkerSize(1.8);
        fHistEventDisplay1 -> SetStats(0);
        fHistEventDisplay1 -> GetXaxis() -> SetNdivisions(510);
        fHistEventDisplay1 -> GetYaxis() -> SetNdivisions(fMaxLayerIndex);
        fHistEventDisplay1 -> GetXaxis() -> SetTickSize(0.01);
        fHistEventDisplay1 -> GetYaxis() -> SetTickSize(0.01);
        fHistEventDisplay1 -> GetYaxis() -> SetBinLabel(1,"");
        //fHistEventDisplay1 -> GetYaxis() -> SetBinLabel(2,"2");
        //fHistEventDisplay1 -> GetYaxis() -> SetBinLabel(3,"3");
        fHistEventDisplay1 -> GetXaxis() -> SetLabelSize(0.04);
        fHistEventDisplay1 -> GetYaxis() -> SetLabelSize(0.04);
        //fHistEventDisplay1 -> SetMarkerSize(fCtrlBinTextSize);
        //if (fPaletteNumber==0)
        //    fHistEventDisplay1 -> SetMarkerColor(kBlack);
        //fHistEventDisplay1 -> SetMinimum(0);
        //fFrameEventDisplay1 = (TH2Poly*) fHistEventDisplay1 -> Clone();
    }

    return fHistEventDisplay1;
}

TH2* SKSiArrayPlane::GetHistEventDisplay2(Option_t *option)
{
    return fHistEventDisplay2;
}

TH1D* SKSiArrayPlane::GetHistChannelBuffer()
{
    if (fHistChannelBuffer==nullptr)
    {
        fHistChannelBuffer = new TH1D(fName+"_Channel",";tb;y",512,0,512);
        fHistChannelBuffer -> SetStats(0);
        fHistChannelBuffer -> SetLineColor(kBlack);
        fHistChannelBuffer -> GetXaxis() -> SetTitleSize(0.06);
        fHistChannelBuffer -> GetYaxis() -> SetTitleSize(0.06);
        fHistChannelBuffer -> GetYaxis() -> SetTitleOffset(1.0);
        fHistChannelBuffer -> GetXaxis() -> SetLabelSize(0.06);
        fHistChannelBuffer -> GetYaxis() -> SetLabelSize(0.06);
    }
    return fHistChannelBuffer;
}

//LKPad* SKSiArrayPlane::FindPad(int cobo, int asad, int aget, int chan)
//{
//    LKPad *pad = nullptr;
//    auto padID = fMapCAACToChannelIndex[cobo][asad][aget][chan];
//    if (padID>=0) {
//        pad = (LKPad*) fChannelArray -> At(padID);
//        return pad;
//    }
//    return (LKPad*) nullptr;
//}
//
int SKSiArrayPlane::FindPadIDFromHistEventDisplay1Bin(int hbin)
{
    return fMapBinToDetector[hbin];
}

//int SKSiArrayPlane::FindZFromHistEventDisplay2Bin(int hbin)
//{
//    Int_t binx, biny, binz;
//    fHistEventDisplay2 -> GetBinXYZ(hbin, binx, biny, binz);
//    return (binx - 1);
//}


TPad* SKSiArrayPlane::Get3DEventPad()
{
    if (fCurrentView==0)
        return GetPadEventDisplay2();
    else {
        return (TPad *) nullptr;
    }
}

void SKSiArrayPlane::UpdateEventDisplay1()
{
    lk_info << endl;
    if (UpdateFlag[kUpdateEventDisplay1]) return;
    else UpdateFlag[kUpdateEventDisplay1] = true;

    if (fHistEventDisplay1==nullptr) 
        return;

    fPadEventDisplay1 -> Clear();
    fPadEventDisplay1 -> cd();
    //if      (fEnergyMax==0) { fHistEventDisplay1 -> SetMinimum(fEnergyMin); fHistEventDisplay1 -> SetMaximum(-1111); }
    //else if (fEnergyMax==1) { fHistEventDisplay1 -> SetMinimum(fEnergyMin); fHistEventDisplay1 -> SetMaximum(4200); }
    //else
    //{
    //    fHistEventDisplay1 -> SetMinimum(fEnergyMin);
    //    fHistEventDisplay1 -> SetMaximum(fEnergyMax);
    //    //if (fEnergyMax>100)
    //    //    gStyle -> SetNumberContours(100);
    //    //else
    //    //    gStyle -> SetNumberContours(fEnergyMax-fEnergyMin);
    //}
    fEventDisplayDrawOption = "";
    fHistEventDisplay1 -> SetTitle("");
    //fFrameEventDisplay1 -> Draw(fEventDisplayDrawOption);
    //fHistEventDisplay1 -> Draw(fEventDisplayDrawOption+"same");
    fHistEventDisplay1 -> SetMinimum(0);
    fHistEventDisplay1 -> Draw("text");
    int numDetectors = fDetectorArray -> GetEntries();
    double zMax = 0.;
    for (auto iDetector=0; iDetector<numDetectors; ++iDetector)
    {
        auto siDetector = (LKSiDetector*) fDetectorArray -> At(iDetector);
        auto histJ = siDetector -> GetHistJunction();
        auto histO = siDetector -> GetHistOhmic();
        double zMaxJ = histJ -> GetMaximum();
        double zMaxO = histO -> GetMaximum();
        if (zMaxJ>zMax) zMax = zMaxJ;
        if (zMaxO>zMax) zMax = zMaxO;
        if (fSelJOID<0) {
            histJ -> Draw("same col");
            histO -> Draw("same col");
        }
        else if (fSelJOID==0) histJ -> Draw("same col");
        else if (fSelJOID==1) histO -> Draw("same col");
        if (fGSelEventDisplay1!=nullptr)
            fGSelEventDisplay1 -> Draw("samel");
    }
    fHistEventDisplay1 -> SetMinimum(1);
    fHistEventDisplay1 -> SetMaximum(zMax);
    fHistEventDisplay1 -> Draw("same text z");

    UpdateUserDrawing();
}

void SKSiArrayPlane::ExecMouseClickEventOnPad(TVirtualPad *pad, double xOnClick, double yOnClick)
{
    for (auto i=0; i<20; ++i) UpdateFlag[i] = false;
    if (pad==fPadControlEvent1)    { fSignalUDBinChange = true; ClickedControlEvent1(xOnClick, yOnClick); fCountChangeOther++; fCountChangeOther = 0; }
    if (pad==fPadControlEvent2)    { ClickedControlEvent2(xOnClick, yOnClick); fCountChangeOther++; }
    if (pad==fPadEventDisplay1)    { ClickedEventDisplay1(xOnClick, yOnClick); fCountChangeOther++; }
    if (pad==fPadJOSideDisplay[0]) { ClickedJOSideDisplay(0); fCountChangeOther++; }
    if (pad==fPadJOSideDisplay[1]) { ClickedJOSideDisplay(1); fCountChangeOther++; }
    if (pad==fPadCtrlUserDrawing)  { ClickedUserDrawing(xOnClick, yOnClick); fCountChangeOther++; }
    if (pad==fPadControlDataDP)    { ClickedControlDataDP(xOnClick, yOnClick); fCountChangeOther++; }
}

void SKSiArrayPlane::ClickedJOSideDisplay(int side)
{
    auto siDetector = (LKSiDetector*) fDetectorArray -> At(fSelDetID);
    if (siDetector==nullptr)
        return;

    int layer = siDetector -> GetLayer();
    double layer1 = fLayerYScale * (layer + fLayerYOffset);
    double layer2 = fLayerYScale * (layer + 1 - fLayerYOffset);
    auto phi1 = siDetector -> GetPhi1();
    auto phi2 = siDetector -> GetPhi2();
    fGSelJOSideDisplay -> Set(0);
    fGSelJOSideDisplay -> SetPoint(0,phi1,layer1);
    if (fSelJOID==side)
        fSelJOID = -1;
    else
        fSelJOID = side;

    UpdateJunctionOhmic();

    UpdateEventDisplay1();
    if (fSignalUDBinChange||fSelJOID>=0) {
        fSignalUDBinChange = true;
        UpdateUserDrawing();
    }
}

void SKSiArrayPlane::ClickedEventDisplay1(double xOnClick, double yOnClick)
{
    if (fHistEventDisplay1==nullptr)
        return;

    int selectedBin = fHistEventDisplay1 -> FindBin(xOnClick, yOnClick);
    if (selectedBin<=0) {
        lk_error << "SelectedBin is " << selectedBin << " (" << xOnClick << ", " << yOnClick << ") " << endl;
        return;
    }
    auto detID = FindPadIDFromHistEventDisplay1Bin(selectedBin);
    if (detID<0) {
        lk_error << "DetectorID is " << detID << ". gbin = " << selectedBin << endl;
        return;
    }

    fSelDetID = detID;
    fSignalUDBinChange = true;

    auto siDetector = (LKSiDetector*) fDetectorArray -> At(fSelDetID);
    if (siDetector==nullptr) {
        lk_error << "detector at " << fSelDetID << " is nullptr" << endl;
        return;
    }

    siDetector -> Print();

    UpdateJunctionOhmic();
    UpdateUserDrawing();
}

void SKSiArrayPlane::UpdateAll()
{
    UpdateEventDisplay1();
    UpdateControlEvent1();
    UpdateControlEvent2();
    UpdateJunctionOhmic();
    UpdateUserDrawing();
    UpdateDataDisplays();
}

void SKSiArrayPlane::UpdateDataDisplays()
{
    if (UpdateFlag[kUpdateDataDisplays]) return;
    else UpdateFlag[kUpdateDataDisplays] = true;

    fPadControlDataDP -> cd();
    fHistControlDataDP -> SetBinContent(fBinCtrlPrevPage,-1);
    fHistControlDataDP -> SetBinContent(fBinCtrlCurrPage,fSelUDPage);
    fHistControlDataDP -> SetBinContent(fBinCtrlNextPage,fNumUDPage);
    fHistControlDataDP -> Draw("text");

    if (fSelUDArray==nullptr)
        return;

    fPadDataDisplayFull   -> GetListOfPrimitives() -> Clear();
    fPadDataDisplayTwo[0] -> GetListOfPrimitives() -> Clear();
    fPadDataDisplayTwo[1] -> GetListOfPrimitives() -> Clear();
    fPadDataDisplayThree[0] -> GetListOfPrimitives() -> Clear();
    fPadDataDisplayThree[1] -> GetListOfPrimitives() -> Clear();
    fPadDataDisplayThree[2] -> GetListOfPrimitives() -> Clear();
    int countPad = 0;
    for (auto iy=0; iy<fNumDDY; ++iy) {
        for (auto ix=0; ix<fNumDDX; ++ix) {
            fPadDataDisplaySmall[countPad] -> GetListOfPrimitives() -> Clear();
            fPadDataDisplaySmall[countPad] -> Modified();
            fPadDataDisplaySmall[countPad] -> Update();
            countPad++;
        }
    }

    int numDataDisplays = fNumDataDisplays;
    if      (fNumGoodDrawings==1 && fSelUDLeastNDraw<=1) { numDataDisplays = 1; }
    else if (fNumGoodDrawings<=2 && fSelUDLeastNDraw<=2) { numDataDisplays = 2; }
    else if (fNumGoodDrawings<=3 && fSelUDLeastNDraw<=3) { numDataDisplays = 3; }

    for (auto iPad=0; iPad<numDataDisplays; ++iPad)
    {
        auto iDrawing = fUDGoodIndex[fSelUDPage][iPad];
        if (iDrawing<0) break;

        TPad *pad;
        if      (numDataDisplays==1) pad = fPadDataDisplayFull;
        else if (numDataDisplays==2) pad = fPadDataDisplayTwo[iPad];
        else if (numDataDisplays==3) pad = fPadDataDisplayThree[iPad];
        else                         pad = fPadDataDisplaySmall[iPad];

        fCanvas -> cd();
        pad -> Draw();
        pad -> cd();

        auto obj = fSelUDArray -> At(iDrawing);
        if (obj==nullptr) break;
        else if (obj -> InheritsFrom(TH2::Class())) { auto hist = (TH2*) obj; hist -> Draw("colz"); }
        else if (obj -> InheritsFrom(TH1::Class())) { auto hist = (TH1*) obj; hist -> Draw();       }
    }
}

void SKSiArrayPlane::UpdateJunctionOhmic()
{
    lk_info << UpdateFlag[kUpdateJunctionOhmic] << endl;
    if (UpdateFlag[kUpdateJunctionOhmic]) return;
    else UpdateFlag[kUpdateJunctionOhmic] = true;

    auto siDetector = (LKSiDetector*) fDetectorArray -> At(fSelDetID);
    if (siDetector==nullptr)
        return;
    siDetector -> Print();
    auto histJ = siDetector -> GetHistJunction();
    auto histO = siDetector -> GetHistOhmic();

    int layer = siDetector -> GetLayer();
    double layer1 = fLayerYScale * (layer + fLayerYOffset);
    double layer2 = fLayerYScale * (layer + 1 - fLayerYOffset);
    auto phi1 = siDetector -> GetPhi1();
    auto phi2 = siDetector -> GetPhi2();
    fGSelEventDisplay1 -> Set(0);
    fGSelEventDisplay1 -> SetPoint(0,phi1,layer1);
    fGSelEventDisplay1 -> SetPoint(1,phi1,layer2);
    fGSelEventDisplay1 -> SetPoint(2,phi2,layer2);
    fGSelEventDisplay1 -> SetPoint(3,phi2,layer1);
    fGSelEventDisplay1 -> SetPoint(4,phi1,layer1);
    fGSelEventDisplay1 -> SetLineColor(kRed);
    fGSelEventDisplay1 -> SetLineWidth(2);
    fPadEventDisplay1 -> cd();
    fGSelEventDisplay1 -> Draw("samel");

    TString drawOption = "colz";
    if (fAccumulateEvents>0)
        drawOption = "text colz";

    fPadJOSideDisplay[0] -> cd();
    histJ -> Draw(drawOption);

    fPadJOSideDisplay[1] -> cd();
    histO -> Draw(drawOption);

    if (fSelJOID<0) {
        fPadJOSideDisplay[0] -> Modified();
        fPadJOSideDisplay[0] -> Update();
    }
    else {
        if (fSelJOID==0) {
            fPadJOSideDisplay[1] -> Modified();
            fPadJOSideDisplay[1] -> Update();
        }
        else if (fSelJOID==1) {
            fPadJOSideDisplay[0] -> Modified();
            fPadJOSideDisplay[0] -> Update();
        }

        fGSelJOSideDisplay -> SetPoint(0,phi1,layer1);
        fGSelJOSideDisplay -> SetPoint(1,phi1,layer2);
        fGSelJOSideDisplay -> SetPoint(2,phi2,layer2);
        fGSelJOSideDisplay -> SetPoint(3,phi2,layer1);
        fGSelJOSideDisplay -> SetPoint(4,phi1,layer1);
        fGSelJOSideDisplay -> SetLineWidth(3);
        fGSelJOSideDisplay -> SetLineColor(kRed);

        fPadJOSideDisplay[fSelJOID] -> cd();
        fGSelJOSideDisplay -> Draw("samel");
    }
}

void SKSiArrayPlane::UpdateEventDisplay2()
{
}

bool SKSiArrayPlane::SetDataFromBranch()
{
    if (!fBranchIsSet && fRun!=nullptr)
    {
        fRawDataArray = fRun -> GetBranchA("RawData");
        fHitArray = fRun -> GetBranchA("Hit");
        fTrackArray = fRun -> GetBranchA("Track");
        fBranchIsSet = true;
    }

    /*
    if (fRawDataArray==nullptr)
        return false;

    auto numChannels = fRawDataArray -> GetEntries();
    for (auto iRawData=0; iRawData<numChannels; ++iRawData)
    {
        auto channel = (GETChannel*) fRawDataArray -> At(iRawData);
        auto cobo = channel -> GetCobo();
        auto asad = channel -> GetAsad();
        auto aget = channel -> GetAget();
        auto chan = channel -> GetChan();
        //auto pad = FindPad(cobo,asad,aget,chan);
        //if (pad==nullptr)
        //{
        //    if (chan!=11&&chan!=22&&chan!=45&&chan!=56)
        //        lk_error << "FPN! CAAC= " << cobo << " " << asad << " " << aget << " " << chan << endl;
        //    continue;
        //}
        //pad -> SetTime(channel->GetTime());
        //pad -> SetEnergy(channel->GetEnergy());
        //pad -> SetPedestal(channel->GetPedestal());
        //pad -> SetDataIndex(iRawData);

        //if (channel->GetEnergy()>selEnergy) {
        //    auto padID = FindPadID(cobo, asad, aget, chan);
        //    fSelPadID = padID;
        //    fSelRawDataID = iRawData;
        //    selEnergy = channel->GetEnergy();
        //}
    }
    */

    return true;
}

void SKSiArrayPlane::FillDataToHistEventDisplay1(Option_t *option)
{
    lk_info << endl;
    if (fCurrentView==0) // 3d
        return;

    Long64_t currentEventID = 0;
    if (fRun!=nullptr)
    {
        currentEventID = fRun -> GetCurrentEventID();
        if (currentEventID==fPrevEventID) {
            lk_warning << "Skipping because event id is same as previous event id!" << currentEventID << endl;
            return;
        }
        fPrevEventID = currentEventID;
    }

    TString optionString(option);
    if (!fFillOptionSelected.IsNull())
        optionString = fFillOptionSelected;
    if (optionString.IsNull())
        optionString = "preview";
    optionString.ToLower();
    lk_info << "Filling " << optionString << " (" << currentEventID << ")" << endl;

    if (fAccumulateEvents==0) {
        TIter nextDetector(fDetectorArray);
        LKSiDetector *siDetector;
        while ((siDetector = (LKSiDetector*) nextDetector()))
            siDetector -> ClearData();
    }

    TString title;

    if (optionString.Index("hit")>=0&&fHitArray!=nullptr)
    {
        /*
        title = "Hit";
        TIter nextHit(fHitArray);
        LKHit* hit = nullptr;
        while (hit = (LKHit*) nextHit())
        {
            auto pos = hit -> GetPosition(fAxisDrift);
            auto i = pos.I();
            auto j = pos.J();
            auto k = pos.K();
            auto energy = hit -> GetCharge();
            auto z = pos.Z();
            auto y = pos.Y();
        }
        */
    }

    else if (fRawDataArray!=nullptr)
    {
        title = "Raw Data";
        TIter nextChannel(fRawDataArray);
        GETChannel* channel = nullptr;
        while ((channel = (GETChannel*) nextChannel()))
        {
            auto cobo = channel -> GetCobo();
            auto asad = channel -> GetAsad();
            auto aget = channel -> GetAget();
            auto chan = channel -> GetChan();
            auto dummyChannel = GetSiChannel(cobo, asad, aget, chan);
            if (dummyChannel==nullptr) {
                lk_error << "channel is not mapped! CAAC = " << cobo << " " << asad << " " << aget << " " << chan << endl;
                continue;
            }
            auto detID = dummyChannel -> GetDetID();
            if (detID<0)
                continue;
            auto side = dummyChannel -> GetSide();
            auto strip = dummyChannel -> GetStrip();
            auto direction = dummyChannel -> GetDirection();
            fSignalUDBinChange = true;
            fSelDetID = detID;
            auto siDetector = (LKSiDetector*) fDetectorArray -> At(detID);
            if (fAccumulateEvents>0)
                siDetector -> AddChannel(channel,side,strip,direction);
            else
                siDetector -> SetChannel(channel,side,strip,direction);
        }
        int numDetectors = fDetectorArray -> GetEntries();
        for (auto iDetector=0; iDetector<numDetectors; ++iDetector)
        {
            auto siDetector = (LKSiDetector*) fDetectorArray -> At(iDetector);
            if (fAccumulateEvents>0) {
                siDetector -> FillHistCount();
            }
            else {
                auto histJ = siDetector -> GetHistJunction();
                auto histO = siDetector -> GetHistOhmic();
                if (fAccumulateEvents>0)
                {
                    histJ -> SetMinimum(0);
                    histO -> SetMinimum(0);
                    histJ -> SetMaximum(-1111);
                    histO -> SetMaximum(-1111);
                }
                else if (fFillOptionSelected=="preview") {
                    siDetector -> FillHistEnergy();
                    histJ -> SetMinimum(0);
                    histO -> SetMinimum(0);
                    histJ -> SetMaximum(4200);
                    histO -> SetMaximum(4200);
                }
                else {
                    siDetector -> FillHistEnergySum();
                    histJ -> SetMinimum(0);
                    histO -> SetMinimum(0);
                    //histJ -> SetMaximum(-1111);
                    //histO -> SetMaximum(-1111);
                    histJ -> SetMaximum(5.e+6);
                    histO -> SetMaximum(5.e+6);
                }
            }
        }
    }
}

void SKSiArrayPlane::ClickedEventDisplay2(double xOnClick, double yOnClick)
{
    //if (fHistEventDisplay2==nullptr)
    //    return;

    //int selectedBin = fHistEventDisplay2 -> FindBin(xOnClick, yOnClick);
    //fSelIZ = FindZFromHistEventDisplay2Bin(selectedBin);

    //fCountChannelGraph = 0;
    //fAccumulateChannel = 3;
    //fHistControlEvent2 -> SetBinContent(fBinCtrlAcmltCh, 3);
    //fHistControlEvent2 -> SetBinContent(fBinCtrlDrawACh, 3);

    //UpdateChannelBuffer();
}

/*
void SKSiArrayPlane::UpdateChannelBuffer()
{
    if (fHistChannelBuffer==nullptr) 
        return;

    if (fSelIZ>=0)
    {
        if (fUsePixelSpace)
        {
            double z1 = fSelIZ;
            double z2 = fSelIZ+1;
            double y1 = 0;
            double y2 = 512;
            fGSelEventDisplay2 -> SetPoint(0,z1,y1);
            fGSelEventDisplay2 -> SetPoint(1,z1,y2);
            fGSelEventDisplay2 -> SetPoint(2,z2,y2);
            fGSelEventDisplay2 -> SetPoint(3,z2,y1);
            fGSelEventDisplay2 -> SetPoint(4,z1,y1);
            fGSelEventDisplay2 -> SetLineColor(kRed);
            if (fPaletteNumber==0)
                fGSelEventDisplay2 -> SetLineColor(kGreen);
        }
        fPadEventDisplay2 -> cd();
        fGSelEventDisplay2 -> Draw("samel");

        fPadChannelBuffer -> cd();
        fHistChannelBuffer -> Reset();
        fHistChannelBuffer -> SetTitle(Form("iz = %d",fSelIZ));
        fHistChannelBuffer -> Draw();
        double yMin=DBL_MAX, yMax=-DBL_MAX;

        auto iz0 = fSelIZ;
        for (auto iPad=0; iPad<fNX; iPad++)
        {
            auto padID = iz0 + fNZ*iPad;
            //auto padID = iz + ix*fNZ;
            auto pad = (LKPad*) fChannelArray -> At(padID);
            auto iz = pad -> GetI();
            auto ix = pad -> GetJ();
            auto idx = pad -> GetDataIndex();
            if (idx<0)
                continue;
            auto channel = (GETChannel*) fRawDataArray -> At(idx);
            auto graph = (TGraph*) fChannelGraphArray -> ConstructedAt(fCountChannelGraph);
            fCountChannelGraph++;
            channel -> FillGraph(graph);
            graph -> Draw("plc samel");
            double x0, y0;
            auto n = graph -> GetN();
            for (auto i=0; i<n; ++i) {
                graph -> GetPoint(i,x0,y0);
                if (yMin>y0) yMin = y0;
                if (yMax<y0) yMax = y0;
            }
        }

        if (fCountChannelGraph>0)
        {
            lk_info << fCountChannelGraph << " channels in iz = " << fSelIZ << endl;
            double dy = 0.1*(yMax - yMin);
            yMin = yMin - dy;
            yMax = yMax + dy;
            if (yMin<0) yMin = 0;
            if (yMax>4200) yMax = 4200;
            fHistChannelBuffer -> SetMinimum(yMin);
            fHistChannelBuffer -> SetMaximum(yMax);
            fHistChannelBuffer -> SetTitle(Form("iz = %d (%d)",fSelIZ,fCountChannelGraph));
        }

        fSelIZ = -1;
    }
    else
        LKEvePlane::UpdateChannelBuffer();
}
*/

void SKSiArrayPlane::ClickedControlDataDP(double xOnClick, double yOnClick)
{
    if (fHistControlDataDP==nullptr)
        return;
    int bing = fHistControlDataDP -> FindBin(xOnClick, yOnClick);
    //int binx, biny, binz;
    //fHistControlDataDP -> GetBinXYZ(bing,binx,biny,binz);

    if (bing==fBinCtrlCurrPage) 
        return;
    else
    {
        if (bing==fBinCtrlPrevPage) 
        {
            if (fSelUDPage!=0) {
                fSelUDPage--;
                fHistControlDataDP -> SetBinContent(fBinCtrlCurrPage,fSelUDPage);
            }
        }
        else if (bing==fBinCtrlNextPage) 
        {
            if (fSelUDPage<fNumSelUDGroup-1) {
                fSelUDPage++;
                fHistControlDataDP -> SetBinContent(fBinCtrlCurrPage,fSelUDPage);
            }
        }
        UpdateDataDisplays();
    }
}

void SKSiArrayPlane::ClickedUserDrawing(double xOnClick, double yOnClick)
{
    if (fHistCtrlUserDrawing==nullptr)
        return;

    int bing = fHistCtrlUserDrawing -> FindBin(xOnClick, yOnClick);
    int selectedBinX, biny, binz;
    fHistCtrlUserDrawing -> GetBinXYZ(bing,selectedBinX,biny,binz);

    if (selectedBinX==1)
    {
        if (fSelUDTab==0 && fSelUDTab==(fNumUDLabels/7)) ;
        else {
            fSignalUDTabChange = true;
            if (fSelUDTab==(fNumUDLabels/7))
                fSelUDTab = 0;
            else
                fSelUDTab++;
            lk_info << "Toggled next menu " << fSelUDTab << endl;
            fSelUDBin = selectedBinX;
        }
    }
    else if (fSelUDBin!=selectedBinX) {
        fSignalUDBinChange = true;
        if (fUserDrawingName[fSelUDTab][selectedBinX-2]=="") {
            return;
        }
        fSelUDBin = selectedBinX;
        lk_info << "Toggled " << fUserDrawingName[fSelUDTab][fSelUDBin-2] << endl;
        UpdateUserDrawing();
    }
}

int SKSiArrayPlane::FindPadID(int cobo, int asad, int aget, int chan)
{
    return fMapCAACToChannelIndex[cobo][asad][aget][chan];
}

LKSiChannel* SKSiArrayPlane::GetSiChannel(int idx)
{
    TObject *channel = nullptr;
    if (idx >= 0 && idx < fChannelArray -> GetEntriesFast())
        channel = fChannelArray -> At(idx);
        return (LKSiChannel *) channel;
    return (LKSiChannel *) nullptr;
}

LKSiChannel* SKSiArrayPlane::GetSiChannel(int cobo, int asad, int aget, int chan)
{
    LKSiChannel *channel = nullptr;
    auto id = FindSiChannelID(cobo, asad, aget, chan);
    if (id<=0)
        return (LKSiChannel*) nullptr;
    channel = GetSiChannel(id);
    return channel;
}

bool SKSiArrayPlane::SetSiChannelData(LKSiChannel* siChannel, GETChannel* channel)
{
    auto cobo = channel -> GetCobo();
    auto asad = channel -> GetAsad();
    auto aget = channel -> GetAget();
    auto chan = channel -> GetChan();
    int detID, side, strip, direction;
    double theta, phi;
    auto dummyChannel = GetSiChannel(cobo, asad, aget, chan);
    if (dummyChannel==nullptr)
        return false;
    siChannel -> SetCobo(cobo);
    siChannel -> SetAsad(asad);
    siChannel -> SetAget(aget);
    siChannel -> SetChan(chan);
    siChannel -> SetDetType   (dummyChannel->GetDetType   ());
    siChannel -> SetDetID     (dummyChannel->GetDetID     ());
    siChannel -> SetChannelID (dummyChannel->GetChannelID ());
    siChannel -> SetLocalID   (dummyChannel->GetLocalID   ());
    siChannel -> SetSide      (dummyChannel->GetSide      ());
    siChannel -> SetStrip     (dummyChannel->GetStrip     ());
    siChannel -> SetDirection (dummyChannel->GetDirection ());
    siChannel -> SetTheta1    (dummyChannel->GetTheta1    ());
    siChannel -> SetTheta2    (dummyChannel->GetTheta2    ());
    siChannel -> SetPhi1      (dummyChannel->GetPhi1      ());
    siChannel -> SetPhi2      (dummyChannel->GetPhi2      ());
    siChannel -> SetTime      (     channel->GetTime      ());
    siChannel -> SetEnergy    (     channel->GetEnergy    ());
    siChannel -> SetPedestal  (     channel->GetPedestal  ());
    siChannel -> SetNoiseScale(     channel->GetNoiseScale());
    siChannel -> SetBuffer    (     channel->GetBuffer    ());
    return true;
}

int SKSiArrayPlane::FindSiDetectorID(int cobo, int asad, int aget, int chan)
{
    return fMapCAACToDetectorIndex[cobo][asad][aget][chan];
}

LKSiDetector* SKSiArrayPlane::GetSiDetector(int idx) {
    TObject *detector = nullptr;
    if (idx >= 0 && idx < fDetectorArray -> GetEntriesFast())
        detector = fDetectorArray -> At(idx);
    return (LKSiDetector *) detector;
}

LKSiDetector* SKSiArrayPlane::GetSiDetector(int cobo, int asad, int aget, int chan) {
    LKSiDetector *detector = nullptr;
    auto id = FindSiDetectorID(cobo, asad, aget, chan);
    if (id>=0)
        detector = GetSiDetector(id);
    return detector;
}
