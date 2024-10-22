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

    if (fUseOldParameterSet)
    {
        int det, side, strip;
        double g0, g1, g2, b0, b1, b2;
        TString parameterFileName = fPar -> GetParString("stark/EnergyCalibrationFile");
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
    }
    else {
        TString energyCalibrationName;
        TString positionCalibrationName;
        fPar -> UpdatePar(energyCalibrationName,"stark/EnergyCalibrationFile");
        fPar -> UpdatePar(positionCalibrationName,"stark/PositionCalibrationFile");
        fEnergyHandler = new SKEnergyHandler(energyCalibrationName,positionCalibrationName);
    }

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
        bool isdEEPairedDetector = (siDetector->GetLayer()<2);
        auto radiusRing = siDetector -> GetRadius();
        auto pairID = fStarkPlane -> FindEPairDetectorID(det);
        if (side==1) continue;
        auto strip = siChannel -> GetStrip();
        auto energyR = siChannel -> GetEnergy();
        auto energyL = siChannel -> GetEnergy2();
        auto position = siChannel -> GetPosition();
        auto phi = siChannel -> GetPhi0();
        auto theta = siChannel -> GetTheta0();
        double energy;
        if (siChannel -> IsStandaloneChannel()) 
        {
            if (fUseOldParameterSet) {
                double g0 = fg0Array[det][side][strip];
                energy = energyR * g0;
            }
            else {
                energy = fEnergyHandler -> RestoreEnergy(det,side,strip,energyR);
            }
            if (energy>0) {
                auto siHit = (SKSiHit*) fHitArray -> ConstructedAt(countHits++);
                if (isEDet) {
                    siHit -> SetDetID(det);
                    siHit -> SetdEDetID(pairID);
                    siHit -> SetEnergy(energy);
                    siHit -> SetIsEDetector(true);
                }
                else {
                    siHit -> SetdEDetID(det);
                    siHit -> SetDetID(pairID);
                    siHit -> SetdE(energy);
                    siHit -> SetIsEDetector(false);
                }
                if (isdEEPairedDetector) siHit -> SetIsEPairDetector(true);
                else siHit -> SetIsEPairDetector(false);
                siHit -> SetJunctionStrip(strip);
                siHit -> SetStripPosition(position);
                siHit -> SetPhi(phi);
                siHit -> SetTheta(theta);
            }
        }
        else if (siChannel -> IsPairedChannel() && energyR>0 && energyL>0)
        {
            double pos, sum;
            if (fUseOldParameterSet) {
                double g1 = fg1Array[det][side][strip];
                double g2 = fg2Array[det][side][strip];
                double b0 = fb0Array[det][side][strip];
                double b1 = fb1Array[det][side][strip];
                double b2 = fb2Array[det][side][strip];
                energyR = energyR * g1;
                energyL = energyL * g2;
                energy = energyR + energyL;
                pos = (energyR - energyL) / energy;
                energy = energy / (b0 + b1*pos + b2*pos*pos) * f241AmAlphaEnergy1;
                pos = pos / 0.8;
            }
            else {
                fEnergyHandler -> RestoreEnergyPosition(det, side, strip, energyL, energyR, pos, sum);
            }
            if (sum>0) {
                auto siHit = (SKSiHit*) fHitArray -> ConstructedAt(countHits++);
                if (det>=0&&det<12)
                {
                    siHit -> SetDetID(det);
                    siHit -> SetdEDetID(pairID);
                    siHit -> SetEnergy(sum);
                    siHit -> SetIsEDetector(true);
                }
                else if (det>=28)
                {
                    siHit -> SetdEDetID(det);
                    siHit -> SetDetID(pairID);
                    siHit -> SetdE(sum);
                    siHit -> SetIsEDetector(false);
                }
                else {
                    siHit -> SetDetID(det);
                    siHit -> SetdEDetID(pairID);
                    siHit -> SetEnergy(sum);
                    siHit -> SetIsEDetector(true);
                }
                if (isdEEPairedDetector) siHit -> SetIsEPairDetector(true);
                else siHit -> SetIsEPairDetector(false);
                siHit -> SetRelativeZ(pos);
                siHit -> SetEnergyLeft(energyR);
                siHit -> SetEnergyRight(energyL);
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
        double energy = siChannel -> GetEnergy();
        if (fUseOldParameterSet) {
            double g0 = fg0Array[det][side][strip];
            energy = energy * g0;
        }
        else {
            energy = fEnergyHandler -> RestoreEnergy(det,side,strip,energy);
        }
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

    //auto numHits2 = fHitArray -> GetEntries();
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

    lk_info << "Number of si-hits = " << countHits << endl;
}
