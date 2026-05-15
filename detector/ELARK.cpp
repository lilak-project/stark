#include <fstream>
using namespace std;

#include "TApplication.h"
#include "ELARK.h"
#include "LKSiDetector.h"
#include "LKSiChannel.h"
#include "TStyle.h"
#include "LKPainter.h"
#include "TPaveText.h"

ClassImp(ELARK)

ELARK::ELARK()
    :ELARK("ELARK","STARK modified for proton elastic scattering experiment (KO-24-21)")
{
}

ELARK::ELARK(const char *name, const char *title)
    :LKDetectorPlane(name, title)
{
    fDetectorArray = new TObjArray();
    fChannelArray = new TObjArray();
}

bool ELARK::Init()
{
    LKDetectorPlane::Init();

    fDetName = "elark";

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

    fPar -> Require(fDetName+"/Mapping1","","","t");
    fPar -> Require(fDetName+"/Mapping2","","","t");
    fPar -> UpdatePar(fDetectorParName,fDetName+"/Mapping1");
    fPar -> UpdatePar(fMappingFileName,fDetName+"/Mapping2");

    fPar -> Require(fDetName+"/ForceMapping","","","t/");
    if (fPar -> CheckPar(fDetName+"/ForceMapping")) {
        fDetectorParName = fPar -> GetParString(fDetName+"/ForceMapping",0);
        fMappingFileName = fPar -> GetParString(fDetName+"/ForceMapping",1);
        lk_note << "YOU FORCED MAPPING FILES BY USING PARAMETER " << fDetName+"/ForceMapping" << endl;
        lk_note << fDetectorParName << endl;
        lk_note << fMappingFileName << endl;
    }

    ifstream fileDetector(fDetectorParName);
    if (!fileDetector.is_open()) {
        lk_error << "Cannot open (" << fDetName+"/Mapping1" << ") " << fDetectorParName << endl;
        return false;
    }
    lk_info << "detector parameter: " << fDetectorParName << endl;

    fMaxLayerIndex = -INT_MAX;
    fMinPhi = DBL_MAX;
    fMaxPhi = -DBL_MAX;

    TString detName;
    int cobo, asad, aget, chan, chan2, detID, strip, side, lr;
    int numSides, numJunctionStrips, numOhmicStrips, layer, ringIndex, useJunctionLR, useOhmicLR;
    double detDistance, detRadius, detDiameter, phi0, detWidth, detHeight, detThickness, tta, phi;
    TString zapJNo, zapONo, markNo;
    int markID, fb;
    int ringType; // 12 or 16
    int dEEType; // 1:E 2:dE
    int polarID; // 0 - 12(16)

    for (auto i=0; i<12; ++i)
        for (auto j=0; j<2; ++j)
            fdEEPairMapping[i][j] = -1;

    detWidth = 40.3;
    detHeight = 75;
    fPar -> UpdatePar(detWidth ,"stark/SiDetectorWidth");
    fPar -> UpdatePar(detHeight,"stark/SiDetectorHeight");
    while (fileDetector >> detName >> detID >> cobo >> asad
            >> zapJNo >> zapONo >> markNo >> markID >> fb
            >> ringType >> dEEType >> polarID >> detDiameter >> detDistance >> phi0 >> ringIndex)
    {
        detRadius = detDiameter*10/2.;
        detDistance = detDistance*10;
        auto siDetector = new LKSiDetector();
        if (fDetectorTypeArray.CheckPar(detName)==false)
            fDetectorTypeArray.AddPar(detName, 1);
        int detTypeIndex = fDetectorTypeArray.GetParIndex(detName);
        double dphi = TMath::ATan2(detWidth/2.,detRadius)*TMath::RadToDeg();
        double phi1 = phi0 - dphi;
        double phi2 = phi0 + dphi;
        double tta1 = TMath::ATan2(detRadius,detDistance-0.5*detHeight)*TMath::RadToDeg();
        double tta2 = TMath::ATan2(detRadius,detDistance+0.5*detHeight)*TMath::RadToDeg();
        int detIndex = fDetectorArray -> GetEntries();
        if (detName=="x6") // TODO
        {
            numJunctionStrips = 8;
            numOhmicStrips = 4;
            useJunctionLR = 1;
            useOhmicLR = 0;
        }
        else if (detName=="csd")
        {
            numJunctionStrips = 8;
            numOhmicStrips = 1;
            useJunctionLR = 0;
            useOhmicLR = 0;
        }
        siDetector -> SetSiType(detName, detTypeIndex, detIndex, detID, 2, numJunctionStrips, numOhmicStrips, useJunctionLR, useOhmicLR);
        layer = ringIndex;
        TVector3 position(detRadius,0,0);
        position.SetPhi(phi0);
        position.SetZ(detDistance);
        siDetector -> SetSiPosition(position, layer, polarID, phi1, phi2, tta1, tta2);
        siDetector -> SetWidth(detWidth);
        siDetector -> SetHeight(detHeight);
        dEEType = dEEType-1;
        siDetector -> SetIsEDetector();
        if (ringType==12)
        {
            fdEEPairMapping[polarID][dEEType] = detID;
            if (dEEType==0) siDetector -> SetIsEDetector();
            else siDetector -> SetIsdEDetector();
        }
        //else siDetector -> SetRow(-1);
        double layer1 = fLayerYScale * (layer + fLayerYOffset);
        double layer2 = fLayerYScale * (layer + 1 - fLayerYOffset);
        TString name = Form("hist_%s_%d_junction",detName.Data(),detID);
        TString title = Form("%s(%d) Junction (%d,%d)",detName.Data(),detID,numJunctionStrips,(useJunctionLR?2:1));
        auto phi0 = 0.5*(phi1+phi2);
        double phi1_ = phi0 - fUserPhiScale*abs(phi0-phi1);
        double phi2_ = phi0 + fUserPhiScale*abs(phi0-phi2);
        siDetector -> CreateHistJunction(name, title, phi1_, phi2_, layer1, layer2, "!axis");
        //siDetector -> CreateHistJunction(name, title, phi0-0.9*(2*TMath::Pi()/12/2), phi0+0.9*(2*TMath::Pi()/12/2), layer1, layer2, "!axis");
        name = Form("hist_%s_%d_ohmic",detName.Data(),detID);
        title = Form("%s(%d) Ohmic (%d,%d)",detName.Data(),detID,(useOhmicLR?2:1),numOhmicStrips);
        siDetector -> CreateHistOhmic(name, title, phi1_, phi2_, layer1, layer2, "!axis");
        fDetectorArray -> Add(siDetector);
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
    int localJID = -1;
    int localOID = -1;
    while (fileCAACMap >> cobo >> asad >> aget >> chan >> chan2 >> detName >> detID >> side >> strip >> lr >> detRadius >> detDistance)
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
        if (detName=="x6"&&side==0) siChannel -> SetIsPairedChannel();
        else                        siChannel -> SetIsStandaloneChannel();
        fMapCAACToDetectorIndex[cobo][asad][aget][chan] = detID;

        auto siDetector = (LKSiDetector*) fDetectorArray -> At(detID);
        siDetector -> RegisterChannel(siChannel);
        TVector3 position0 = siDetector -> GetPosition();
        TVector3 position1 = siDetector -> GetPosition();
        TVector3 position2 = siDetector -> GetPosition();
        double phi = siDetector -> GetPhi0();
        double theta = siDetector -> GetTheta0();
        double width = siDetector -> GetHeight();
        double height = siDetector -> GetWidth();
        if (side==0) {
            auto numStrips = siDetector -> GetNumJunctionStrips();
            position0.SetPhi(0);
            position0.SetY(0.5*width-(localJID+0.5)*(width/numStrips)); position0.RotateZ(phi);
            position1.SetY(0.5*width-(localJID    )*(width/numStrips)); position1.RotateZ(phi);
            position2.SetY(0.5*width-(localJID+1  )*(width/numStrips)); position2.RotateZ(phi);
        }
        if (side==1) {
            auto numStrips = siDetector -> GetNumOhmicStrips();
            position0.SetZ(position0.Z()+0.5*height-(localOID+0.5)*(height/numStrips));
            position1.SetZ(position1.Z()+0.5*height-(localOID    )*(height/numStrips));
            position2.SetZ(position2.Z()+0.5*height-(localOID+1  )*(height/numStrips));
        }
        siChannel -> SetPosition(position0);
        double phi1 = position1.Phi();
        double phi2 = position2.Phi();
        siChannel -> SetPhi1(phi1);
        siChannel -> SetPhi2(phi2);
        siChannel -> SetTheta1(theta);
        siChannel -> SetTheta2(theta);

        ++globalChannelIndex;
    }

    lk_info << globalChannelIndex << " channels are mapped!" << endl;

    return true;
}

bool ELARK::Init2()
{
    return true;
}

void ELARK::Print(Option_t *option) const
{
    lk_info << "Silicon array " << fName << endl;

    int numDetectors = fDetectorArray -> GetEntries();
    for (auto iDetector=0; iDetector<numDetectors; ++iDetector)
    {
        auto siDetector = (LKSiDetector*) fDetectorArray -> At(iDetector);
        siDetector -> Print();
    }
}

LKDrawing* ELARK::NewHistPlaneDrawing(TString name)
{
    auto drawing = new LKDrawing();

    auto hist = NewHistPlane(name);
    drawing -> Add(hist,"l");

    TH2PolyBin* bin = nullptr;
    auto listOfBins = hist -> GetBins();
    TIter nextBin(listOfBins);
    while ((bin=(TH2PolyBin*)nextBin()))
    {
        int no = bin -> GetBinNumber();
        int detID = no - 1;
        double x2 = bin -> GetXMax();
        double x1 = bin -> GetXMin();
        double y2 = bin -> GetYMax();
        double y1 = bin -> GetYMin();
        auto detector = (LKSiDetector*) fDetectorArray -> At(detID);
        auto pairID = FindEPairDetectorID(detID,0);
        auto type = detector -> GetDetTypeName();
        auto pv = new TPaveText(x1,y1,x2,y2);
        pv -> SetBorderSize(0);
        pv -> SetFillStyle(0);
        pv -> SetTextSize(0.032);
        pv -> SetAllWith("","align",22);
        pv -> AddText(Form("D-%d",detID));
        pv -> AddText(Form("P-%d",pairID));
        if (detID==0 || detID==12 || detID==28) pv -> SetTextColor(kRed);
        type.ToUpper();
        pv -> AddText(Form("%s",type.Data()));
        //pv -> SetTextAlign(22);
        drawing -> Add(pv);
    }

    return drawing;
}

TH2Poly* ELARK::NewHistPlane(TString name)
{
    auto hist = new TH2Poly(name,"ELARK;#phi (deg.);Ring",12,fMinPhi-0.02*(fMaxPhi-fMinPhi),fMaxPhi+0.02*(fMaxPhi-fMinPhi),fMaxLayerIndex,0,fLayerYScale*(fMaxLayerIndex+1));
    int numDetectors = fDetectorArray -> GetEntries();
    int maxBins = numDetectors + 10;
    for (auto det=0; det<numDetectors; ++det)
    {
        auto siDetector = (LKSiDetector*) fDetectorArray -> At(det);
        auto layer = siDetector -> GetLayer();
        auto phi1 = siDetector -> GetPhi1();
        auto phi2 = siDetector -> GetPhi2();
        auto phi0 = 0.5*(phi1+phi2);
        phi1 = phi0 - fUserPhiScale*abs(phi0-phi1);
        phi2 = phi0 + fUserPhiScale*abs(phi0-phi2);
        double layer1 = fLayerYScale * (layer + fLayerYOffset);
        double layer2 = fLayerYScale * (layer + 1 - fLayerYOffset);
        auto bin = hist -> AddBin(phi1, layer1, phi2, layer2);
    }
    hist -> SetStats(0);
    hist -> SetMarkerSize(1.8);
    hist -> GetXaxis() -> SetNdivisions(510);
    hist -> GetYaxis() -> SetNdivisions(fMaxLayerIndex);
    hist -> GetXaxis() -> SetTickSize(0.01);
    hist -> GetYaxis() -> SetTickSize(0.01);
    hist -> GetXaxis() -> SetLabelSize(0.04);
    hist -> GetYaxis() -> SetLabelSize(0.04);
    hist -> GetYaxis() -> SetBinLabel(1,"");
    //hist -> GetYaxis() -> SetBinLabel(2,"12R-E");
    //hist -> GetYaxis() -> SetBinLabel(3,"16R");
    return hist;
}

int ELARK::FindPadID(int cobo, int asad, int aget, int chan)
{
    return fMapCAACToChannelIndex[cobo][asad][aget][chan];
}

LKSiChannel* ELARK::GetSiChannel(int idx)
{
    TObject *channel = nullptr;
    if (idx >= 0 && idx < fChannelArray -> GetEntriesFast())
        channel = fChannelArray -> At(idx);
        return (LKSiChannel *) channel;
    return (LKSiChannel *) nullptr;
}

LKSiChannel* ELARK::GetSiChannel(int cobo, int asad, int aget, int chan)
{
    LKSiChannel *channel = nullptr;
    auto id = fMapCAACToChannelIndex[cobo][asad][aget][chan];;///FindSiChannelID(cobo, asad, aget, chan);
    if (id<0)
        return (LKSiChannel*) nullptr;
    channel = GetSiChannel(id);
    return channel;
}

bool ELARK::SetSiChannelData(LKSiChannel* siChannel, GETChannel* channel)
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

    siChannel -> SetInverted      (dummyChannel->GetInverted ());
    siChannel -> SetPairArrayIndex(dummyChannel->GetPairArrayIndex());
    siChannel -> SetEnergy2       (dummyChannel->GetEnergy2());
    siChannel -> SetPosition      (dummyChannel->GetPosition());

    siChannel -> SetTheta1    (dummyChannel->GetTheta1    ());
    siChannel -> SetTheta2    (dummyChannel->GetTheta2    ());
    siChannel -> SetPhi1      (dummyChannel->GetPhi1      ());
    siChannel -> SetPhi2      (dummyChannel->GetPhi2      ());

    if (dummyChannel -> IsStandaloneChannel()) siChannel -> SetIsStandaloneChannel();
    if (dummyChannel -> IsPairedChannel())     siChannel -> SetIsPairedChannel();

    siChannel -> SetTime      (     channel->GetTime      ());
    siChannel -> SetEnergy    (     channel->GetEnergy    ());
    siChannel -> SetPedestal  (     channel->GetPedestal  ());
    siChannel -> SetNoiseScale(     channel->GetNoiseScale());
    siChannel -> SetBuffer    (     channel->GetBuffer    ());
    if (dummyChannel->IsPairedChannel()) siChannel -> SetIsPairedChannel();
    else                                 siChannel -> SetIsStandaloneChannel();
    return true;
}

int ELARK::FindSiDetectorID(int cobo, int asad, int aget, int chan)
{
    return fMapCAACToDetectorIndex[cobo][asad][aget][chan];
}

LKSiDetector* ELARK::GetSiDetector(int idx) {
    TObject *detector = nullptr;
    if (idx >= 0 && idx < fDetectorArray -> GetEntriesFast())
        detector = fDetectorArray -> At(idx);
    return (LKSiDetector *) detector;
}

LKSiDetector* ELARK::GetSiDetector(int cobo, int asad, int aget, int chan) {
    LKSiDetector *detector = nullptr;
    auto id = FindSiDetectorID(cobo, asad, aget, chan);
    if (id>=0)
        detector = GetSiDetector(id);
    return detector;
}

int ELARK::FindEPairDetectorID(int det, int i)
{
    if (det>=12&&det<28) return det;
    auto detector = GetSiDetector(det);
    auto polarID = detector -> GetRow();
    if (detector->GetLayer()>1)
        return -1;
    if (i>=0)
        return fdEEPairMapping[polarID][i];
    if (detector->IsEDetector())
        return fdEEPairMapping[polarID][1];
    return fdEEPairMapping[polarID][0];
}

void ELARK::FireStrip(int det, int side, int strip, double energy)
{
    auto siDetector = (LKSiDetector*) fDetectorArray -> At(det);
    siDetector -> Fire(side, strip, energy);
}

void ELARK::ClearFiredFlags()
{
    int numDetectors = fDetectorArray -> GetEntries();
    for (auto iDetector=0; iDetector<numDetectors; ++iDetector)
    {
        auto siDetector = (LKSiDetector*) fDetectorArray -> At(iDetector);
        siDetector -> ClearFiredFlags();
    }
}

int ELARK::GetNumFiredDetectors()
{
    int countFiredDetectors = 0;
    int numDetectors = fDetectorArray -> GetEntries();
    for (auto iDetector=0; iDetector<numDetectors; ++iDetector)
    {
        auto siDetector = (LKSiDetector*) fDetectorArray -> At(iDetector);
        if (siDetector->GetNumFiredStrips()>=0)
            countFiredDetectors++;
    }
    return countFiredDetectors;
}

LKSiDetector* ELARK::GetFiredDetector(int iFired)
{
    int countFiredDetectors = 0;
    int numDetectors = fDetectorArray -> GetEntries();
    for (auto iDetector=0; iDetector<numDetectors; ++iDetector)
    {
        auto siDetector = (LKSiDetector*) fDetectorArray -> At(iDetector);
        if (siDetector->GetNumFiredStrips()>=0) {
            if (countFiredDetectors==iFired)
                return siDetector;
            countFiredDetectors++;
        }
    }
    return (LKSiDetector*) nullptr;
}

double ELARK::CalculatePairSolidAngle(double theta1, double theta2, int numDetectors)
{
    double dphi = (2*TMath::Pi()/12) * numDetectors;
    double solid_angle = abs( dphi * (TMath::Cos(theta1)-TMath::Cos(theta2)) );
    return solid_angle;
}

double ELARK::Calculate16RingSolidAngle(double theta1, double theta2, int numDetectors)
{
    double dphi = (2*TMath::Pi()/16) * numDetectors;
    double solid_angle = abs( dphi * (TMath::Cos(theta1)-TMath::Cos(theta2)) );
    return solid_angle;
}
