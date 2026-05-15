#include "EKRecoHit.h"

ClassImp(EKRecoHit);

EKRecoHit::EKRecoHit()
{
    Clear();
}

void EKRecoHit::Clear(Option_t *option)
{
    fdEIndex = -1;
    fEIndex = -1;
    fDetID = -1;
    fPairID = -1;
    fIsInGate = true;
    fIsEPairDetector = false;
    fKeyEnergy = 0;
    fQValue = 0;
    fEJunction = -1;
    fEOhmic = 0;
    fdEJunction = 0;
    fdEOhmic = 0;
    fPhi = -1;
    fTheta = -1;
    fSolidAngle = -1;
    fRelativeZ = -999;
    fRealZ = -999;

    fdEHit = nullptr;
    fEHit = nullptr;
}

double EKRecoHit::GetQValue(double beamEnergy)
{
	double massP = 1;
	double massAr = 40;
	double beamEnergyAr = beamEnergy * massAr;
	double value = fKeyEnergy * (1 + massP/massAr) - beamEnergyAr * (1 - massAr/massAr) - 2 / massAr * TMath::Sqrt(massAr * massP * beamEnergyAr * fKeyEnergy) * TMath::Cos(fTheta*TMath::DegToRad());
	return value;
}

