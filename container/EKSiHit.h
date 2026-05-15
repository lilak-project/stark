#ifndef EKSIHIT_HH
#define EKSIHIT_HH

#include "LKContainer.h"
#include "LKLogger.h"
#include "TVector3.h"

class EKSiHit : public LKContainer
{
    public:
        EKSiHit();
        virtual ~EKSiHit() { ; }

        virtual void Clear(Option_t *option="");
        virtual void Print(Option_t *option="") const {}
        void PrintAll() const;

        void SetDetID(int id) { fDetID = id; }
        void SetJunctionStrip(int id) { fJunctionStrip = id; }
        void SetOhmicStrip(int id) { fOhmicStrip = id; }
        void SetIsEPairDetector(int isPair) { fIsEPairDetector = isPair; }
        void SetIsEDetector(int isE) { fIsEDetector = isE; }
        void SetInGate(int in) { fInGate = in; }
        void SetKeyEnergy(double e) { fKeyEnergy = e; }
        void SetEnergy(double energy ) { fEnergy = energy; }
        void SetEnergyLeft(double energy ) { fEnergyLeft = energy; }
        void SetEnergyRight(double energy ) { fEnergyRight = energy; }
        void SetEnergyOhmic(double energy ) { fEnergyOhmic = energy; }
        void SetStripPosition(TVector3 position) { fStripPosition = position; }
        void SetRelativeZ(double z) { fRelativeZ = z; }
        void SetX(double x) { fX = x; }
        void SetPhi(double phi) { fPhi = phi; }
        void SetTheta(double theta) { fTheta = theta; }

        int GetDetID() const { return fDetID; }
        int GetJunctionStrip() const { return fJunctionStrip; }
        int GetOhmicStrip() const { return fOhmicStrip; }
        bool IsEPairDetector() const { return fIsEPairDetector; }     ///< Hit is in 12-ring detector where dE-E detectors are paired. This doesn’t mean both dE and E energy are found.
        bool IsNotEPairDetector() const { return !fIsEPairDetector; } ///< Hit is in 16-ring detector where only E detectors are used.
        bool InGate() const { return fInGate; }
        bool IsEDetector() const { return fIsEDetector; } ///< Hit belong to the E-detector
        double GetKeyEnergy() const { return fKeyEnergy; } ///< Energy left in dE-detector
        double GetE() const { return fEnergy; } ///< Energy left in E-detector
        double GetEnergy() const { return fEnergy; } ///< Same as GetE()
        double GetEnergyLeft() const { return fEnergyLeft; }
        double GetEnergyRight() const { return fEnergyRight; }
        double GetEnergyOhmic() const { return fEnergyOhmic; }
        TVector3 GetStripPosition() const { return fStripPosition; } ///< center position of the strip from the target position 
        double GetRelativeZ(double length=2) const { return fRelativeZ*length/2.; } ///< Relative z = (left-right)/(left+right) * (length/2). By default, range is from -1 to 1. If detector length is given, this is z-position from the center of center of detector-strip
        double GetDetectorZ() const { return fStripPosition.Z(); } ///< z-position of the detector from the target position.
        double GetFinalZ(double length) const { return GetDetectorZ() + GetRelativeZ(length); } ///< detector-z + relative-z (scaled to the detector length)
        double GetX() const { return fX; }
        double GetPhi() const { return fPhi; } ///< Phi angle of the detector
        double GetTheta() const { return fTheta; } ///< Theta angle of the detector

        double GetQvalue(double massNoProton, double massNoAr, double beamEnergy, double energy=0, double theta=-9999);
        double GetQvalueE(double massNoProton, double massNoAr, double beamEnergy) { return GetQvalue(massNoProton, massNoAr, beamEnergy, 0); }

        void Grab(bool value=true) { fGrab = value; }
        bool IsGrabbed() const { return fGrab; }

    protected:
        bool fGrab = false; //!

        bool fInGate; ///< true: hit is on gate, false: hit is off gate.
        bool fIsEPairDetector; ///< true: detector is paired with dE-E detectors, false: no dE-detector is use.
        bool fIsEDetector; ///< true: E-detector, false: dE-detector

        int fDetID; ///< detector id

        int fJunctionStrip; ///< strip number of fired junction
        int fOhmicStrip; ///< strip number of fired ohmic

        double fKeyEnergy; ///< Representative energy
        double fEnergy; ///< E
        double fEnergyOhmic; ///< energy of ohmic side
        double fEnergyLeft; ///< energy of downstream side
        double fEnergyRight; ///< energy of upstream side

        TVector3 fStripPosition; ///< center position of the strip from the target position
        double fRelativeZ; ///< relative z = (right-left)/(right+left) ranging from = -1 to 1
        double fX; ///< x-position (axis through width direction) within si-detector
        double fPhi; ///< phi
        double fTheta; ///< theta

    ClassDef(EKSiHit,1);
};

#endif
