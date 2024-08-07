#include "LKRun.h"
#include "LKLogger.h"
#include "LKSiChannel.h"
#include "SKEnergyRestorationTask.h"
#include "SKSiHit.h"

ClassImp(SKEnergyRestorationTask)

SKEnergyRestorationTask::SKEnergyRestorationTask()
    :LKTask("SKEnergyRestorationTask","SKEnergyRestorationTask")
{
}

bool SKEnergyRestorationTask::Init()
{
    fStarkPlane = (SKSiArrayPlane*) fRun -> FindDetectorPlane("SKSiArrayPlane");
    fSiChannelArray = fRun -> GetBranchA("SiChannel");
    fHitArray = fRun -> RegisterBranchA("SiHit","SKSiHit",20);

    TString ecalName = "stark/EnergyCalibrationFile";
    if (fPar->CheckPar(ecalName)==false)
    {
        lk_error << "Energy calibration file parameter " << ecalName << " should be set!" << endl;
        return false;
    }

    int det, side, strip;
    double g0, g1, g2, b0, b1, b2;
    TString parameterFileName = fPar -> GetParString(ecalName);
    auto file = new TFile(parameterFileName,"read");
    auto tree = (TTree*) file -> Get("parameters");
    tree -> SetBranchAddress("det"    ,&det    );
    tree -> SetBranchAddress("side"   ,&side   );
    tree -> SetBranchAddress("strip"  ,&strip  );
    tree -> SetBranchAddress("g0"     ,&g0     );
    tree -> SetBranchAddress("g1"     ,&g1     );
    tree -> SetBranchAddress("g2"     ,&g2     );
    tree -> SetBranchAddress("b0"     ,&b0     );
    tree -> SetBranchAddress("b1"     ,&b1     );
    tree -> SetBranchAddress("b2"     ,&b2     );
    auto n = tree -> GetEntries();
    for (auto i=0; i<n; ++i)
    {
        tree -> GetEntry(i);
        fg0Array[det][side][strip] = g0;
        fg1Array[det][side][strip] = g1;
        fg2Array[det][side][strip] = g2;
        fb0Array[det][side][strip] = b0;
        fb1Array[det][side][strip] = b1;
        fb2Array[det][side][strip] = b2;
    }
    file -> Close();

    return true;
}

void SKEnergyRestorationTask::Exec(Option_t*)
{
    fHitArray -> Clear("C");

    fStarkPlane -> ClearFiredFlags();

    int countHits = 0;
    auto numChannels = fSiChannelArray -> GetEntries();
    for (auto iChannel=0; iChannel<numChannels; ++iChannel)
    {
        auto siChannel = (LKSiChannel*) fSiChannelArray -> At(iChannel);
        auto det = siChannel -> GetDetID();
        auto side = siChannel -> GetSide();
        auto siDetector = fStarkPlane -> GetSiDetector(det);
        auto isEDet = siDetector -> IsEDetector();
        auto dEEPairID = siDetector -> GetRow();
        auto radiusRing = siDetector -> GetRadius();
        auto pairID = fStarkPlane -> FindEPairDetectorID(det);
        if (side==1) continue;
        auto strip = siChannel -> GetStrip();
        auto energy1 = siChannel -> GetEnergy();
        auto energy2 = siChannel -> GetEnergy2();
        auto position = siChannel -> GetPosition();
        auto phi = siChannel -> GetPhi0();
        auto theta = siChannel -> GetTheta0();
        if (siChannel -> IsStandaloneChannel()) 
        {
            double g0 = fg0Array[det][side][strip];
            double energy = energy1 * g0;
            if (energy>0) {
                auto siHit = (SKSiHit*) fHitArray -> ConstructedAt(countHits++);
                siHit -> SetDetID(det);
                siHit -> SetdEDetID(pairID);
                if (isEDet) {
                    siHit -> SetEnergy(energy);
                    siHit -> SetIsEDetector(true);
                }
                else {
                    siHit -> SetdE(energy);
                    siHit -> SetIsEDetector(false);
                }
                if (dEEPairID) siHit -> SetIsEPairDetector(true);
                else siHit -> SetIsEPairDetector(false);
                siHit -> SetJunctionStrip(strip);
                siHit -> SetStripPosition(position);
                siHit -> SetPhi(phi);
                siHit -> SetTheta(theta);
            }
        }
        else if (siChannel -> IsPairedChannel() && energy1>0 && energy2>0)
        {
            double g1 = fg1Array[det][side][strip];
            double g2 = fg2Array[det][side][strip];
            double b0 = fb0Array[det][side][strip];
            double b1 = fb1Array[det][side][strip];
            double b2 = fb2Array[det][side][strip];
            energy1 = energy1 * g1;
            energy2 = energy2 * g2;
            double energy = energy1 + energy2;
            double pos = (energy1 - energy2) / energy;
            energy = energy / (b0 + b1*pos + b2*pos*pos) * f241AmAlphaEnergy1;
            pos = pos / 0.8; // TODO @todo
            //pos = pos / 0.7; // TODO @todo
            if (energy>0) {
                auto siHit = (SKSiHit*) fHitArray -> ConstructedAt(countHits++);
                siHit -> SetDetID(det);
                siHit -> SetdEDetID(pairID);
                if (det>=0&&det<12)
                {
                    siHit -> SetEnergy(energy);
                    siHit -> SetIsEDetector(true);
                }
                else if (det>=28)
                {
                    siHit -> SetdE(energy);
                    siHit -> SetIsEDetector(false);
                }
                else {
                    siHit -> SetEnergy(energy);
                    siHit -> SetIsEDetector(true);
                }
                if (dEEPairID) siHit -> SetIsEPairDetector(true);
                else siHit -> SetIsEPairDetector(false);
                siHit -> SetRelativeZ(pos);
                siHit -> SetEnergyLeft(energy1);
                siHit -> SetEnergyRight(energy2);
                siHit -> SetJunctionStrip(strip);
                siHit -> SetStripPosition(position);
                siHit -> SetPhi(phi);
                theta = TMath::RadToDeg()*TMath::ATan(radiusRing/siHit->GetFinalZ(75));
                siHit -> SetTheta(theta);
            }
        }
        else
            continue;
    }

    for (auto iChannel=0; iChannel<numChannels; ++iChannel)
    {
        auto siChannel = (LKSiChannel*) fSiChannelArray -> At(iChannel);
        auto det = siChannel -> GetDetID();
        auto side = siChannel -> GetSide();
        if (side==0) continue;
        auto strip = siChannel -> GetStrip();
        double g0 = fg0Array[det][side][strip];
        double energy = siChannel -> GetEnergy() * g0;
        for (auto iHit=0; iHit<countHits; ++iHit)
        {
            auto siHit = (SKSiHit*) fHitArray -> At(iHit);
            if (siHit->GetDetID()==det)
            {
                siHit -> SetOhmicStrip(strip);
                siHit -> SetEnergyOhmic(energy);
            }
        }
    }

    lk_info << "Number of si-hits = " << countHits << endl;
}
