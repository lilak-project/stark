#include "LKRun.h"
#include "LKLogger.h"
#include "SKPairMatchingTask.h"
#include "SKSiHit.h"

ClassImp(SKPairMatchingTask)

SKPairMatchingTask::SKPairMatchingTask()
    :LKTask("SKPairMatchingTask","SKPairMatchingTask")
{
}

bool SKPairMatchingTask::Init()
{
    fStarkPlane = (SKSiArrayPlane*) fRun -> FindDetectorPlane("SKSiArrayPlane");
    fHitArray = fRun -> KeepBranchA("SiHit");

    return true;
}

void SKPairMatchingTask::Exec(Option_t*)
{
    auto numHits = fHitArray -> GetEntries();
    for (auto iHit=0; iHit<numHits; ++iHit)
    {
        auto siHit1 = (SKSiHit*) fHitArray -> At(iHit);
        auto pairID = fStarkPlane -> FindEPairDetectorID(siHit1->GetDetID());
        if (pairID>=0)
        {
            siHit1 -> SetIsEPairDetector(true);
            for (auto jHit=0; jHit<numHits; ++jHit)
            {
                if (iHit==jHit)
                    continue;
                auto siHit2 = (SKSiHit*) fHitArray -> At(jHit);
                if (pairID==siHit2->GetDetID())
                {
                    if (siHit1->IsEDetector()) {
                        siHit1 -> SetdE(siHit2->GetdE());
                        siHit2 -> SetEnergy(-1);
                    }
                    else if (siHit2->IsEDetector())  {
                        siHit2 -> SetdE(siHit1->GetdE());
                        siHit2 -> SetEnergy(-1);
                    }
                    else
                        lk_debug << siHit1 -> GetDetID() << " " << siHit2 -> GetDetID() << endl;
                }
            }
        }
    }

    for (auto iHit=0; iHit<numHits; ++iHit)
    {
        auto siHit = (SKSiHit*) fHitArray -> At(iHit);
        if (siHit->GetdE()>0 && siHit->GetE()<0)
            fHitArray -> Remove(siHit);
    }

    fHitArray -> Compress();

    lk_info << "Number of si-hits = " << fHitArray->GetEntries() << endl;
}
