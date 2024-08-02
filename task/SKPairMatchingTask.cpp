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
            for (auto jHit=0; jHit<numHits; ++jHit)
            {
                if (iHit==jHit)
                    continue;
                auto siHit2 = (SKSiHit*) fHitArray -> At(jHit);
                if (pairID==siHit2->GetDetID())
                {
                    if      (siHit1->IsEDetector()) siHit1 -> SetdE(siHit2->GetEnergy());
                    else if (siHit2->IsEDetector()) siHit2 -> SetdE(siHit1->GetEnergy());
                    else
                        lk_debug << siHit1 -> GetDetID() << " " << siHit2 -> GetDetID() << endl;
                }
            }
        }
    }

    lk_info << endl;
}
