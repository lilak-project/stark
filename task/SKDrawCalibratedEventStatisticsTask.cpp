#include "TStyle.h"

#include "LKRun.h"
#include "LKLogger.h"
#include "GETChannel.h"
#include "SKSiHit.h"
#include "SKDrawCalibratedEventStatisticsTask.h"

ClassImp(SKDrawCalibratedEventStatisticsTask)

SKDrawCalibratedEventStatisticsTask::SKDrawCalibratedEventStatisticsTask()
    :LKTask("SKDrawCalibratedEventStatisticsTask","SKDrawCalibratedEventStatisticsTask")
{
}

bool SKDrawCalibratedEventStatisticsTask::Init()
{
    fSiHitArray = fRun -> GetBranchA("SiHit");
    if (fSiHitArray==nullptr) {
        lk_error << "Branch SiHit do not exist!!!" << endl;
        return false;
    }

    fStarkPlane = (SKSiArrayPlane*) fRun -> FindDetectorPlane("SKSiArrayPlane");
    if (fStarkPlane==nullptr) {
        lk_error << "SKSiArrayPlane do not exist!!!" << endl;
        lk_error << fName << " must run with SKSiArrayPlane" << endl;
        return false;
    }

    int ne = 400;
    double e1 = 0;
    double e2 = 10;
    fPar -> UpdateBinning(fName+"/binning_cal_energy", ne, e1, e2);

    fHistHP[0] = new TH2D("fHistCHP0","Junction;strip;energy (MeV)",40*fNumJStrips,0,40*fNumJStrips,ne,e1,e2);;
    fHistHP[1] = new TH2D("fHistCHP1",   "Ohmic;strip;energy (MeV)",40*fNumOStrips,0,40*fNumOStrips,ne,e1,e2);;
    fHistHP[0] -> SetStats(0);
    fHistHP[1] -> SetStats(0);
    fHistHP[0] -> GetXaxis() -> SetTitleOffset(2.7);
    fHistHP[1] -> GetXaxis() -> SetTitleOffset(2.7);
    for (auto det=0; det<40; ++det)
    {
        auto detector = fStarkPlane -> GetSiDetector(det);
        auto name = detector -> GetDetTypeName();
        auto ring = detector -> GetLayer();
        TString sring = "dE"; if (ring==1) sring = "E"; if (ring==2) sring = "16E"; 
        fHistHP[0] -> GetXaxis() -> SetBinLabel(det*fNumJStrips+1,Form("%d (%s,%s)",det,name.Data(),sring.Data()));
        fHistHP[1] -> GetXaxis() -> SetBinLabel(det*fNumOStrips+1,Form("%d (%s,%s)",det,name.Data(),sring.Data()));
    }

    auto array = new TObjArray();
    array -> Add(fHistHP[0]);
    array -> Add(fHistHP[1]);
    fStarkPlane -> AddUserDrawings("cal_HP", -1, -1, array);

    return true;
}

void SKDrawCalibratedEventStatisticsTask::Exec(Option_t*)
{
    if (!fStarkPlane -> GetAccumulateEvents()) {
        //fHistHP[0] -> Reset();
        //fHistHP[1] -> Reset();
        //fHistET[0] -> Reset();
        //fHistET[1] -> Reset();
    }

    auto numHits = fSiHitArray -> GetEntries();
    for (auto iHit=0; iHit<numHits; ++iHit)
    {
        auto siHit = (SKSiHit*) fSiHitArray -> At(iHit);
        int detID = siHit -> GetDetID();
        int stripJ = siHit -> GetJunctionStrip();
        int stripO = siHit -> GetOhmicStrip();
        double energySum = siHit -> GetEnergy();
        double energyOhmic = siHit -> GetEnergyOhmic();

        fHistHP[0] -> Fill(detID*fNumJStrips+stripJ, energySum);
        fHistHP[1] -> Fill(detID*fNumOStrips+stripO,energyOhmic);
    }
}
