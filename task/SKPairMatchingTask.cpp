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

    if (numHits==1)
    {
        auto siHit = (SKSiHit*) fHitArray -> At(0);
        if ( siHit->IsEPairDetector() && siHit->IsEDetector() )
            fHitArray -> Remove(siHit);
    }
    else
    {
        for (auto iHit=0; iHit<numHits; ++iHit)
        {
            auto eHit = (SKSiHit*) fHitArray -> At(iHit);
            auto eID = eHit -> GetDetID();
            if (eHit->IsEPairDetector()==false) continue;
            if (eHit->IsEDetector()==false) continue;

            bool dEIsSet = false;
            for (auto jHit=0; jHit<numHits; ++jHit)
            {
                if (iHit==jHit)
                    continue;
                auto deHit = (SKSiHit*) fHitArray -> At(jHit);
                auto deID = deHit -> GetDetID();
                if (deHit->IsdEDetector()==false) continue;
                if (deID!=fStarkPlane->FindEPairDetectorID(eID)) continue;

                eHit -> SetdEDetID(deHit->GetdEDetID());
                eHit -> SetdE(deHit->GetdE());
                eHit -> SetdEOhmic(deHit->GetEnergyOhmic());
                eHit -> SetRelativeZdE(deHit->GetRelativeZ());
                deHit -> Grab();
                dEIsSet = true;
                break;
            }
            if (dEIsSet==false)
                eHit -> Grab();
        }
        for (auto iHit=0; iHit<numHits; ++iHit)
        {
            auto siHit = (SKSiHit*) fHitArray -> At(iHit);
            if (siHit->IsGrabbed())
                fHitArray -> Remove(siHit);
        }
    }

    fHitArray -> Compress();

    auto numHits2 = fHitArray -> GetEntries();
    //for (auto iHit=0; iHit<numHits2; ++iHit)
    //{
    //    auto siHit = (SKSiHit*) fHitArray -> At(iHit);
    //    double deID = siHit->GetdEDetID();
    //    double eID = siHit->GetDetID();
    //    double de = siHit->GetdEOhmic();
    //    double energy = siHit->GetEnergyOhmic();
    //    int deIs = (deID>=0)?(fStarkPlane->GetSiDetector(deID)->IsEDetector()):-1;
    //    int eIs = (eID>=0)?(fStarkPlane->GetSiDetector(eID)->IsEDetector()):-1;
    //    TString deTitle = (deIs>=0)?((deIs==1)?"E":"dE"):"X";
    //    TString eTitle  = ( eIs>=0)?((eIs==1) ?"E":"dE"):"X";
    //    if (siHit->IsEPairDetector())
    //        lk_debug << "[E-PAIR] de_ohmic=(" << de << "|" << deID << "," << deTitle << "), " << " e_ohmic=(" << energy << "|" << eID << "," << eTitle << ") " << endl;
    //    else
    //        lk_debug << "[SINGLE] de_ohmic=(" << de << "|" << deID << "," << deTitle << "), " << " e_ohmic=(" << energy << "|" << eID << "," << eTitle << ") " << endl;
    //}

    lk_info << "Number of si-hits = " << fHitArray->GetEntries() << " (" << numHits << " -> " << numHits2 << ")" << endl;
}
