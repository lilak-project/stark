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
        auto det = siHit1 -> GetDetID();

        if ((det>=0&&det<12) || det>=28)
        {

        auto pairID = fStarkPlane -> FindEPairDetectorID(siHit1->GetDetID());
        if (pairID>=0)
        {
            siHit1 -> SetIsEPairDetector(true);
            for (auto jHit=iHit+1; jHit<numHits; ++jHit)
            {
                auto siHit2 = (SKSiHit*) fHitArray -> At(jHit);
                if (siHit1 -> GetDetID() == siHit2 -> GetDetID())
                    continue; // TODO
                if (pairID==siHit2->GetDetID())
                {
                    //lk_debug << iHit << " " << jHit << " .. " << siHit1 -> GetDetID() << "(" << pairID << ")" << ", " <<  siHit2 -> GetDetID() << " :  " << siHit1 -> GetEnergy() << " " << siHit1 -> GetdE() << ", " << siHit2 -> GetEnergy() << " " << siHit2 -> GetdE() << endl;

                    SKSiHit *siHit_dE = nullptr;
                    SKSiHit *siHit_E = nullptr;
                    if (siHit1->IsEDetector()) {
                        siHit_E = siHit1;
                        siHit_dE = siHit2;
                    }
                    else {
                        siHit_E = siHit2;
                        siHit_dE = siHit1;
                    }

                    if (siHit_dE->fGrab) {
                        double energyPrev = siHit_dE -> GetEnergy();
                        double energy = siHit_E -> GetEnergy();
                        if (energy>energyPrev)
                        {
                            siHit_E -> SetdE(siHit_dE->GetdE());
                            siHit_dE -> SetEnergy(siHit_E->GetEnergy());
                            siHit_dE -> fGrab = true;
                        }
                    }
                    else {
                        siHit_E -> SetdE(siHit_dE->GetdE());
                        siHit_dE -> SetEnergy(siHit_E->GetEnergy());
                        siHit_dE -> fGrab = true;
                    }
                }
            }
        }
        }
    }

    for (auto iHit=0; iHit<numHits; ++iHit)
    {
        auto siHit = (SKSiHit*) fHitArray -> At(iHit);
        //if (siHit->GetdE()>0 && siHit->GetE()<0)
        if (siHit->fGrab)
            fHitArray -> Remove(siHit);
    }

    fHitArray -> Compress();

    lk_info << "Number of si-hits = " << fHitArray->GetEntries() << endl;
}
