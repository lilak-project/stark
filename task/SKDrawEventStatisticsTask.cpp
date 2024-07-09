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

    int ne=400;
    double e1=0, e2=4100;

    int nt=40;
    double t1=60, t2=100;

    fPar -> UpdateBinning(fName+"/binning_energy", ne, e1, e2);
    lk_info << "Energy binning is : " << ne << " " << e1 << " " << e2 << endl;
    lk_info << "Theta  binning is : " << nt << " " << t1 << " " << t2 << endl;

    auto numDetectors = fStarkPlane -> GetNumSiDetectors();
    for (int iDetector=0; iDetector<numDetectors; ++iDetector) {
        auto siDetector = fStarkPlane -> GetSiDetector(iDetector);
        int numJunctionChannels = siDetector -> GetNumJunctionChannels();
        if (fMaxJunctionChannels<numJunctionChannels)
            fMaxJunctionChannels = numJunctionChannels;
        int numOhmicChannels = siDetector -> GetNumOhmicChannels();
        if (fMaxOhmicChannels<numOhmicChannels)
            fMaxOhmicChannels = numOhmicChannels;
    }

    fHistHP[0] = new TH2D("fHistHP0","Junction;channel;energy",numDetectors*fMaxJunctionChannels,0,numDetectors*fMaxJunctionChannels,ne,e1,e2);
    fHistHP[1] = new TH2D("fHistHP1",   "Ohmic;channel;energy",numDetectors*fMaxOhmicChannels   ,0,numDetectors*fMaxOhmicChannels   ,ne,e1,e2);
    fHistHP[0] -> GetXaxis() -> SetNdivisions(numDetectors,false);
    fHistHP[1] -> GetXaxis() -> SetNdivisions(numDetectors,false);
    for (int iDetector=0; iDetector<numDetectors; ++iDetector) {
        fHistHP[0] -> GetXaxis() -> SetBinLabel(iDetector*fMaxJunctionChannels+1,Form("%d",iDetector));
        fHistHP[1] -> GetXaxis() -> SetBinLabel(iDetector*fMaxOhmicChannels+1,Form("%d",iDetector));
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

    return true;
}

void SKDrawEventStatisticsTask::Exec(Option_t*)
{
    if (!fStarkPlane -> GetAccumulateEvents()) {
        fHistHP[0] -> Reset();
        fHistHP[1] -> Reset();
        fHistET[0] -> Reset();
        fHistET[1] -> Reset();
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
        fHistHP[0] -> Fill(detID*fMaxJunctionChannels+localID, energy);
        fHistHP[1] -> Fill(detID*fMaxOhmicChannels+localID, energy);
        fHistET[0] -> Fill(theta,energy);
        fHistET[1] -> Fill(theta,energy);
    }
}
