#include "EKSiHit.h"

ClassImp(EKSiHit);

EKSiHit::EKSiHit()
{
    Clear();
}

void EKSiHit::Clear(Option_t *option)
{
    fGrab = false;
    fInGate = true;
    fIsEPairDetector = false;
    fIsEDetector = true;
    fKeyEnergy = 0;
    fDetID = -1;
    fJunctionStrip = -1;
    fOhmicStrip = -1;
    fEnergy = 0;
    fEnergyLeft = 0;
    fEnergyRight = 0;
    fEnergyOhmic = 0;
    fStripPosition = TVector3();
    fRelativeZ = -999;
    fX = -999;
    fPhi = -999;
    fTheta = -999;
}

void EKSiHit::PrintAll() const
{
    e_info << "[SiHit] " << std::endl;
    e_cout << "    - fDetID           = " << fDetID           << std::endl;
    e_cout << "    - fJunctionStrip   = " << fJunctionStrip   << std::endl;
    e_cout << "    - fOhmicStrip      = " << fOhmicStrip      << std::endl;
    e_cout << "    - fIsEPairDetector = " << fIsEPairDetector << std::endl;
    e_cout << "    - fInGate          = " << fInGate          << std::endl;
    e_cout << "    - fIsEDetector     = " << fIsEDetector     << std::endl;
    e_cout << "    - fKeyEnergy       = " << fKeyEnergy       << std::endl;
    e_cout << "    - fEnergy          = " << fEnergy          << std::endl;
    e_cout << "    - fEnergyLeft      = " << fEnergyLeft      << std::endl;
    e_cout << "    - fEnergyRight     = " << fEnergyRight     << std::endl;
    e_cout << "    - fEnergyOhmic     = " << fEnergyOhmic     << std::endl;
    e_cout << "    - fStripPosition   = (" << fStripPosition.X() << ", " << fStripPosition.Y() << ", " << fStripPosition.Z() << ")" << std::endl;
    e_cout << "    - fRelativeZ       = " << fRelativeZ       << std::endl;
    e_cout << "    - fX               = " << fX               << std::endl;
    e_cout << "    - fPhi             = " << fPhi             << std::endl;
    e_cout << "    - fTheta           = " << fTheta           << std::endl;
}

double EKSiHit::GetQvalue(double massNoProton, double massNoAr, double beamEnergy, double energy, double theta)
{
    if (energy==0) energy = fEnergy;
    if (theta<-999) theta = fTheta;
    double E_beam = beamEnergy * massNoAr;
    double value1 = energy * (1 + massNoProton/massNoAr);
    double value2 = E_beam * (1 - massNoAr/massNoAr);
    double value3 = TMath::Sqrt(massNoAr * massNoProton * E_beam * energy);
    double value4 = TMath::Cos(theta*TMath::DegToRad());
    double value = value1 - value2 - 2/massNoAr * value3 * value4;
    return value;
}
