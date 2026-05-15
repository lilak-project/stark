#ifndef EKRECOHIT_HH
#define EKRECOHIT_HH

#include "LKContainer.h"
#include "LKLogger.h"
#include "TVector3.h"
#include "EKSiHit.h"

class EKRecoHit : public LKContainer
{
    public:
        EKRecoHit();
        virtual ~EKRecoHit() { ; }
        virtual void Clear(Option_t *option="");

        void SetdEIndex   (int value) { fdEIndex = value; }
        void SetEIndex    (int value) { fEIndex = value; }
        void SetDetID     (int value) { fDetID = value; }
        void SetPairID    (int value) { fPairID = value; }
        void SetIsInGate  (bool value) { fIsInGate = value; }
        void SetIsEPairDetector(bool value) { fIsEPairDetector = value; }
        void SetKeyEnergy (double value) { fKeyEnergy = value; }
        void SetQValue    (double value) { fQValue = value; }
        void SetEJunction (double value) { fEJunction = value; }
        void SetEOhmic    (double value) { fEOhmic = value; }
        void SetdEJunction(double value) { fdEJunction = value; }
        void SetdEOhmic   (double value) { fdEOhmic = value; }
        void SetPhi       (double value) { fPhi = value; }
        void SetTheta     (double value) { fTheta = value; }
        void SetSolidAngle(double value) { fSolidAngle = value; }
        void SetRelativeZ (double value) { fRelativeZ = value; }
        void SetRealZ     (double value) { fRealZ = value; }
        void SetJunctionStrip(int value) { fJunctionStrip = value; }
        void SetOhmicStrip(int value)    { fOhmicStrip = value; }

        int GetdEIndex()       const { return fdEIndex; }
        int GetEIndex()        const { return fEIndex; }
        int GetDetID()         const { return fDetID; }
        int GetPairID()        const { return fPairID; }
        bool InGate()          const { return fIsInGate; }
        bool GetIsInGate()     const { return fIsInGate; }
        bool IsEPairDetector() const { return fIsEPairDetector; } ///< Hit is in 12-ring detector where dE-E detectors are paired. This doesn’t mean both dE and E energy are found.
        bool GetIsEPairDetector() const { return fIsEPairDetector; } ///< true: detector is paired with dE-E detectors, false: no dE-detector is use.
        double GetKeyEnergy()  const { return fKeyEnergy; }
        double GetQValue()     const { return fQValue; }
        double GetEJunction()  const { return fEJunction; }
        double GetEOhmic()     const { return fEOhmic; }
        double GetdEJunction() const { return fdEJunction; }
        double GetdEOhmic()    const { return fdEOhmic; }
        double GetPhi()        const { return fPhi; }
        double GetTheta()      const { return fTheta; }
        double GetSolidAngle() const { return fSolidAngle; }
        double GetRelativeZ()  const { return fRelativeZ; }
        double GetRealZ()      const { return fRealZ; }
        int GetJunctionStrip() const { return fJunctionStrip; }
        int GetOhmicStrip   () const { return fOhmicStrip; }

        bool IsEPairAndOnlydE()   const { return (fIsEPairDetector && fEOhmic==0); } ///< Is dEE-pair detector but only dE energy is left
        bool IsEPairAndOnlyE()    const { return (fIsEPairDetector && fdEOhmic==0); } ///< Is dEE-pair detector but only E energy is left
        bool IsEPairAndBothdEE()  const { return (fIsEPairDetector && fdEOhmic>0 && fEOhmic>0); } ///< Is dEE-pair detector and both dE and E energies are left

        void SetdEHit(EKSiHit* dEHit) { fdEHit = dEHit; }
        void SetEHit(EKSiHit* eHit) { fEHit = eHit; }
        EKSiHit* GetdEHit() { return fdEHit; }
        EKSiHit* GetEHit () { return fEHit; }

        double GetQValue(double beamEnergy);

    protected:
        int fdEIndex; ///< dE-Hit index of SiHit array;
        int fEIndex; ///< E-Hit index of SiHit array;
        int fDetID;
        int fPairID;
        bool fIsInGate;
        bool fIsEPairDetector;
        double fKeyEnergy;
        double fQValue;
        double fEJunction;
        double fEOhmic;
        double fdEJunction;
        double fdEOhmic;
        double fPhi;
        double fTheta;
        double fSolidAngle;
        double fRelativeZ;
        double fRealZ;
        int fJunctionStrip;
        int fOhmicStrip;

        EKSiHit* fdEHit = nullptr; //!
        EKSiHit* fEHit = nullptr; //!

    ClassDef(EKRecoHit,2);
};

#endif
