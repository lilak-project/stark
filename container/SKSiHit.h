#ifndef SKSIHIT_HH
#define SKSIHIT_HH

#include "LKContainer.h"
#include "LKLogger.h"
#include "TVector3.h"

/**
    ## SKSiHit
    Hit left by single particle. The data may cover one or two detectors (dE and E for two).
   
    ## Hit Cases
    - Case1: Particle go through detector in one of 12-ring detectors (dE-E paired detectors). Energy was left in dE-detector and not in E-detector
    - Case2: Particle go through detector in one of 12-ring detectors (dE-E paired detectors). Energy was left in E-detector and not in dE-detector
    - Case3: Particle go through detector in one of 12-ring detectors (dE-E paired detectors). Energy was left both in dE-detector and not in E-detector
    - Case4: Particle go through detector in one of 12-ring detectors (dE-E paired detectors). Energy was left either in dE or E or in both.
    - Case5: Particle go through detector in one of 16-ring detectors (no dE-E pair).
   
    @code{.cpp}
    {
        bool case1 = siHit -> IsEPairAndOnlydE();
        bool case2 = siHit -> IsEPairAndOnlyE();
        bool case3 = siHit -> IsEPairAndBothdEE();
        bool case4 = siHit -> IsEPairDetector();
        bool case5 = siHit -> IsNotEPairDetector();
   
        // above should be same as
   
        case1 = (siHit->IsEPairDetector() && siHit->GetE()==0);
        case2 = (siHit->IsEPairDetector() && siHit->GetdE()==0);
        case3 = (siHit->IsEPairDetector() && siHit->GetdE()>0 && siHit->GetE()>0);
    }
    @endcode
   
    ## Energy
    For any cases, total energy left in hit is GetEnergyTotal().
    This is equal to E for 16-ring detectors, and is equal to dE+E for 12-ring (pair) detectors.
   
    ## z-position
    - GetDetectorZ(): z-position of Detector from the target.
    - GetRelativeZ(), GetRelativeZ(length): z-position of hit from the center of the detector, where length is the
      detector-length and should be given by user. Default length is 2 (z ranging from -1 to 1).
    - GetFinalZ(), GetFinalZ(length): Final z-position of hit from the target.
*/
class SKSiHit : public LKContainer
{
    public:
        SKSiHit();
        virtual ~SKSiHit() { ; }

        virtual void Clear(Option_t *option="");
        virtual void Print(Option_t *option="") const {}
        void PrintAll() const;

        void SetDetID(int id) { fDetID = id; }
        void SetdEDetID(int id) { fdEDetID = id; }
        void SetdEDetAIndex(int idx) { fdEDetAIndex = idx; }
        void SetJunctionStrip(int id) { fJunctionStrip = id; }
        void SetOhmicStrip(int id) { fOhmicStrip = id; }
        void SetIsEPairDetector(int isPair) { fIsEPairDetector = isPair; }
        void SetInGate(int in) { fInGate = in; }
        void SetIsEDetector(int isE) { fIsEDetector = isE; }
        void SetKeyEnergy(double e) { fKeyEnergy = e; }
        void SetdE(double de) { fdE = de; }
        void SetdEOhmic(double de) { fdEOhmic = de; }
        void SetEnergy(double energy ) { fEnergy = energy; }
        void SetEnergyLeft(double energy ) { fEnergyLeft = energy; }
        void SetEnergyRight(double energy ) { fEnergyRight = energy; }
        void SetEnergyOhmic(double energy ) { fEnergyOhmic = energy; }
        void SetStripPosition(TVector3 position) { fStripPosition = position; }
        void SetRelativeZ(double z) { fRelativeZ = z; }
        void SetRelativeZdE(double z) { fRelativeZdE = z; }
        void SetX(double x) { fX = x; }
        void SetPhi(double phi) { fPhi = phi; }
        void SetTheta(double theta) { fTheta = theta; }


        bool IsEPairAndOnlydE()  const { return (fIsEPairDetector && fEnergy==0); } ///< Is dEE-pair detector but only dE energy is left
        bool IsEPairAndOnlyE()   const { return (fIsEPairDetector && fdE==0); } ///< Is dEE-pair detector but only E energy is left
        bool IsEPairAndBothdEE() const { return (fIsEPairDetector && fdE>0 && fEnergy>0); } ///< Is dEE-pair detector and both dE and E energies are left
        bool IsEPairAndOnlyOhmicdE()  const { return (fIsEPairDetector && fEnergyOhmic==0); }
        bool IsEPairAndOnlyOhmicE()   const { return (fIsEPairDetector && fdEOhmic==0); }
        bool IsEPairAndBothOhmicdEE() const { return (fIsEPairDetector && fdEOhmic>0 && fEnergyOhmic>0); }


        int GetDetID() const { return fDetID; }
        int GetdEDetID() const { return fdEDetID; }
        int GetdEDetAIndex() const { return fdEDetAIndex; }
        int GetJunctionStrip() const { return fJunctionStrip; }
        int GetOhmicStrip() const { return fOhmicStrip; }
        bool IsEPairDetector() const { return fIsEPairDetector; }     ///< Hit is in 12-ring detector where dE-E detectors are paired. This doesn’t mean both dE and E energy are found.
        bool IsNotEPairDetector() const { return !fIsEPairDetector; } ///< Hit is in 16-ring detector where only E detectors are used.
        bool InGate() const { return fInGate; }
        bool IsEDetector() const { return fIsEDetector; } ///< Hit belong to the E-detector
        bool IsdEDetector() const { return (!fIsEDetector); } ///< Hit belong to the dE-detector
        double GetKeyEnergy() const { return fKeyEnergy; } ///< Energy left in dE-detector
        double GetdE() const { return fdE; } ///< Energy left in dE-detector
        double GetdEOhmic() const { return fdEOhmic; } ///< Energy left in dE-detector
        double GetE() const { return fEnergy; } ///< Energy left in E-detector
        double GetEnergy() const { return fEnergy; } ///< Same as GetE()
        double GetEnergyTotal() const { return fEnergy+fdE; } ///< dE+E
        double GetEnergyTotalOhmic() const { return fEnergyOhmic+fdEOhmic; } ///< dE+E
        double GetEnergyLeft() const { return fEnergyLeft; }
        double GetEnergyRight() const { return fEnergyRight; }
        double GetEnergyOhmic() const { return fEnergyOhmic; }
        TVector3 GetStripPosition() const { return fStripPosition; } ///< center position of the strip from the target position 
        double GetRelativeZ(double length=2) const { return fRelativeZ*length/2.; } ///< Relative z = (left-right)/(left+right) * (length/2). By default, range is from -1 to 1. If detector length is given, this is z-position from the center of center of detector-strip
        double GetRelativeZdE(double length=2) const { return fRelativeZdE*length/2.; } ///< Relative z = (left-right)/(left+right) * (length/2). By default, range is from -1 to 1. If detector length is given, this is z-position from the center of center of detector-strip
        double GetDetectorZ() const { return fStripPosition.Z(); } ///< z-position of the detector from the target position.
        double GetFinalZ(double length) const { return GetDetectorZ() + GetRelativeZ(length); } ///< detector-z + relative-z (scaled to the detector length)
        double GetFinalZdE(double length) const { return GetDetectorZ() + GetRelativeZdE(length); } ///< detector-z + relative-z (scaled to the detector length)
        double GetX() const { return fX; }
        double GetPhi() const { return fPhi; } ///< Phi angle of the detector
        double GetTheta() const { return fTheta; } ///< Theta angle of the detector

        double GetQvalue(double massNoProton, double massNoAr, double beamEnergy, double energy=0, double theta=-9999);
        double GetQvalueE(double massNoProton, double massNoAr, double beamEnergy) { return GetQvalue(massNoProton, massNoAr, beamEnergy, 0); }
        double GetQvaluedE(double massNoProton, double massNoAr, double beamEnergy) { return GetQvalue(massNoProton, massNoAr, beamEnergy, -1); }

        void Grab(bool value=true) { fGrab = value; }
        bool IsGrabbed() const { return fGrab; }

    protected:
        bool fGrab = false; //!

        bool fInGate; ///< true: hit is on gate, false: hit is off gate.
        bool fIsEPairDetector; ///< true: detector is paired with dE-E detectors, false: no dE-detector is use.
        bool fIsEDetector; ///< true: E-detector, false: dE-detector
        double fKeyEnergy; ///< Representative energy

        int fDetID; ///< detector id
        int fJunctionStrip; ///< strip number of fired junction
        int fOhmicStrip; ///< strip number of fired ohmic
        double fEnergy; ///< E
        double fEnergyLeft; ///< energy of downstream side
        double fEnergyRight; ///< energy of upstream side
        double fEnergyOhmic; ///< energy of ohmic side
        TVector3 fStripPosition; ///< center position of the strip from the target position
        double fRelativeZ; ///< relative z = (right-left)/(right+left) ranging from = -1 to 1
        double fX; ///< x-position (axis through width direction) within si-detector
        double fPhi; ///< phi
        double fTheta; ///< theta

        int fdEDetID; ///< dE detector id
        int fdEDetAIndex; ///< dE hit idx in dEHit clones array
        double fdE; ///< dE
        double fdEOhmic; ///< energy of ohmic side
        double fRelativeZdE;

    ClassDef(SKSiHit,3);
};

#endif
