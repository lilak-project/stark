#include "TStyle.h"

#include "LKRun.h"
#include "LKLogger.h"
#include "GETChannel.h"
#include "LKSiChannel.h"
#include "SKDrawEventStatisticsTask.h"

ClassImp(SKDrawEventStatisticsTask)

SKDrawEventStatisticsTask::SKDrawEventStatisticsTask()
    :LKTask("SKDrawEventStatisticsTask","SKDrawEventStatisticsTask")
{
}

bool SKDrawEventStatisticsTask::Init()
{
    fSiChannelArray = fRun -> GetBranchA("SiChannel");

    fStarkPlane = (SKSiArrayPlane*) fRun -> FindDetectorPlane("SKSiArrayPlane");
    if (fStarkPlane==nullptr) {
        lk_error << "SKSiArrayPlane do not exist!!!" << endl;
        lk_error << fName << " must run with SKSiArrayPlane" << endl;
        return false;
    }

    int nsx = 800;
    double sx1 = 0;
    double sx2 = 8;

    int ne = 400;
    double e1 = 0;
    double e2 = 4100;

    int nt = 40;
    double t1 = 60;
    double t2 = 100;

    fPar -> UpdateBinning(fName+"/binning_energy", ne, e1, e2);
    fPar -> UpdateBinning(fName+"/binning_theta", nt, t1, t2);
    fPar -> UpdateBinning(fName+"/binning_stripx", nsx, sx1, sx2);
    lk_info << "Energy binning is : " << ne << " " << e1 << " " << e2 << endl;
    lk_info << "Theta  binning is : " << nt << " " << t1 << " " << t2 << endl;
    lk_info << "StripX binning is : " << nsx << " " << sx1 << " " << sx2 << endl;

    auto numDetectors = fStarkPlane -> GetNumSiDetectors();
    for (int iDetector=0; iDetector<numDetectors; ++iDetector) {
        auto siDetector = fStarkPlane -> GetSiDetector(iDetector);
        int numJunctionChannels = siDetector -> GetNumJunctionChannels();
        if (fMaxChannels[0]<numJunctionChannels)
            fMaxChannels[0] = numJunctionChannels;
        int numOhmicChannels = siDetector -> GetNumOhmicChannels();
        if (fMaxChannels[1]<numOhmicChannels)
            fMaxChannels[1] = numOhmicChannels;
    }

    fHistHP[0] = new TH2D("fHistHP0","Junction;channel;energy",numDetectors*fMaxChannels[0],0,numDetectors*fMaxChannels[0],ne,e1,e2);
    fHistHP[1] = new TH2D("fHistHP1",   "Ohmic;channel;energy",numDetectors*fMaxChannels[1],0,numDetectors*fMaxChannels[1],ne,e1,e2);
    fHistHP[0] -> GetXaxis() -> SetNdivisions(numDetectors,false);
    fHistHP[1] -> GetXaxis() -> SetNdivisions(numDetectors,false);
    for (int iDetector=0; iDetector<numDetectors; ++iDetector) {
        fHistHP[0] -> GetXaxis() -> SetBinLabel(iDetector*fMaxChannels[0]+1,Form("%d",iDetector));
        fHistHP[1] -> GetXaxis() -> SetBinLabel(iDetector*fMaxChannels[1]+1,Form("%d",iDetector));
    }
    auto array = new TObjArray();
    array -> Add(fHistHP[0]);
    array -> Add(fHistHP[1]);
    fStarkPlane -> AddUserDrawings("hit_pattern", -1, -1, array);

    fHistET[0] = new TH2D("fHistET0","Junction;theta;energy",nt,t1,t2,ne,e1,e2);
    fHistET[1] = new TH2D("fHistET1",   "Ohmic;theta;energy",nt,t1,t2,ne,e1,e2);
    auto array2 = new TObjArray();
    array2 -> Add(fHistET[0]);
    array2 -> Add(fHistET[1]);
    fStarkPlane -> AddUserDrawings("e_vs_t", -1, -1, array2);

    for (auto detID=0; detID<40; ++detID) {
        auto array3 = new TObjArray();
        for (auto strip=0; strip<8; ++strip) {
            fHistEP[detID][strip] = new TH2D(Form("fHistEP_%d_%d",detID,strip),Form("%d_%d;position;energy",detID,strip),nsx,sx1,sx2,ne,e1,2*e2);
            fHistLR[detID][strip] = new TH2D(Form("fHistLR_%d_%d",detID,strip),Form("%d_%d;left;right",detID,strip),ne,e1,e2,ne,e1,e2);
            array3 -> Add(fHistEP[detID][strip]);
            array3 -> Add(fHistLR[detID][strip]);
        }
        fStarkPlane -> AddUserDrawings("e_vs_pos", detID, -1, array3);
    }

    return true;
}

void SKDrawEventStatisticsTask::Exec(Option_t*)
{
    if (!fStarkPlane -> GetAccumulateEvents()) {
        //fHistHP[0] -> Reset();
        //fHistHP[1] -> Reset();
        //fHistET[0] -> Reset();
        //fHistET[1] -> Reset();
    }

    auto numChannels = fSiChannelArray -> GetEntries();
    for (auto iChannel=0; iChannel<numChannels; ++iChannel)
    {
        auto siChannel = (LKSiChannel*) fSiChannelArray -> At(iChannel);
        if (siChannel==nullptr)
            continue;
        int localID = siChannel -> GetLocalID();
        int detID = siChannel -> GetDetID();
        int side = siChannel -> GetSide();
        int strip = siChannel -> GetStrip();
        int direction = siChannel -> GetDirection();
        double energy = siChannel -> GetEnergy();
        double theta = siChannel -> GetTheta1();

        fHistHP[side] -> Fill(detID*fMaxChannels[side]+localID, energy);
        fHistET[side] -> Fill(theta,energy);

        auto pairChannel = siChannel -> GetPairChannel();
        if (pairChannel!=nullptr && siChannel->GetSide()==0 && siChannel->GetDirection()==0) {
            double left = siChannel -> GetEnergy();
            double right = pairChannel -> GetEnergy();
            double energySum = left + right;
            double position = 7.5 * 0.5 * ((left-right)/energySum + 1);
            fHistEP[detID][strip] -> Fill(position, energySum);
            fHistLR[detID][strip] -> Fill(left, right);
        }
    }
}
