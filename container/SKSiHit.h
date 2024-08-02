#ifndef SKSIHIT_HH
#define SKSIHIT_HH

#include "LKContainer.h"
#include "LKLogger.h"
#include "TVector3.h"

class SKSiHit : public LKContainer
{
    public:
        SKSiHit();
        virtual ~SKSiHit() { ; }

        virtual void Clear(Option_t *option="") {
            fDetID = -1;
            fJunctionStrip = -1;
            fOhmicStrip = -1;
            fIsEPairDetector = false;
            fIsEDetector = true;
            fdE = -999;
            fEnergy = -999;
            fEnergyLeft = -999;
            fEnergyRight = -999;
            fEnergyOhmic = -999;
            fStripPosition = TVector3();
            fZ = -999;
            fX = -999;
            fPhi = -999;
            fTheta = -999;
        }

        virtual void Print(Option_t *option="") const {
            e_cout << "det,de,e|z,phi,tta = " << fDetID << ", " << fdE << ", " << fEnergy << " | " << fZ << ", " << fPhi << ", " << fTheta << std::endl;
        }

        void SetDetID(int id) { fDetID = id; }
        void SetJunctionStrip(int id) { fJunctionStrip = id; }
        void SetOhmicStrip(int id) { fOhmicStrip = id; }
        void SetIsEPairDetector(int isPair) { fIsEPairDetector = isPair; }
        void SetIsEDetector(int isE) { fIsEDetector = isE; }
        void SetdE(double de) { fdE = de; }
        void SetEnergy(double energy ) { fEnergy = energy; }
        void SetEnergyLeft(double energy ) { fEnergyLeft = energy; }
        void SetEnergyRight(double energy ) { fEnergyRight = energy; }
        void SetEnergyOhmic(double energy ) { fEnergyOhmic = energy; }
        void SetStripPosition(TVector3 position) { fStripPosition = position; }
        void SetZ(double z) { fZ = z; }
        void SetX(double x) { fX = x; }
        void SetPhi(double phi) { fPhi = phi; }
        void SetTheta(double theta) { fTheta = theta; }

        int GetDetID() const { return fDetID; }
        int GetJunctionStrip() const { return fJunctionStrip; }
        int GetOhmicStrip() const { return fOhmicStrip; }
        bool IsEPairDetector() const { return fIsEPairDetector; }
        bool IsEDetector() const { return fIsEDetector; }
        bool IsdEDetector() const { return (!fIsEDetector); }
        double GetdE() const { return fdE; }
        double GetEnergy() const { return fEnergy; }
        double GetEnergyLeft() const { return fEnergyLeft; }
        double GetEnergyRight() const { return fEnergyRight; }
        double GetEnergyOhmic() const { return fEnergyOhmic; }
        TVector3 GetStripPosition() const { return fStripPosition; }
        double GetZ() const { return fZ; }
        double GetX() const { return fX; }
        double GetPhi() const { return fPhi; }
        double GetTheta() const { return fTheta; }

    protected:
        int fDetID; ///< detector id
        int fJunctionStrip; ///< strip number of fired junction
        int fOhmicStrip; ///< strip number of fired ohmic
        bool fIsEPairDetector; ///< true: detector is paired with dE-E detectors, false: no dE-detector is use.
        bool fIsEDetector; ///< true: E-detector, false: dE-detector
        double fdE; ///< dE
        double fEnergy; ///< E
        double fEnergyLeft; ///< energy of downstream side
        double fEnergyRight; ///< energy of upstream side
        double fEnergyOhmic; ///< energy of ohmic side
        TVector3 fStripPosition; ///< position of strip;
        double fZ; ///< global z-position in mm
        double fX; ///< x-position (axis through width direction) within si-detector
        double fPhi; ///< phi
        double fTheta; ///< theta

    ClassDef(SKSiHit,1);
};

#endif
