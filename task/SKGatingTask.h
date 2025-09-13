#ifndef SKGATINGTASK_HH
#define SKGATINGTASK_HH

#include "TH1.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TClonesArray.h"
#include "TCutG.h"

#include "LKTask.h"
#include "SKSiArrayPlane.h"

class SKGatingTask : public LKTask
{
    public:
        SKGatingTask();
        virtual ~SKGatingTask() {}

        bool Init();
        void Exec(Option_t*);

    private:
        TClonesArray* fSiHitArray = nullptr;
        TClonesArray* fRecoHeaderArray = nullptr;

  		double fBeamEnergy = 5.838; // [MeV/u]
  		double fBeamRate = 1.9E6; // [pps]
  		double fRunTime = 6540; // [sec] for run20-run25

        bool fUseEPairOnlydE = false;
        bool fUseEPairOnlyE = false;
        bool fUseEPairBoth = false;
        bool fUseEPair = false;
        bool fUseNotEPair = false;
        bool fUseScintGate = false;

        int fDetIDRange[2] = {0};
        int fdEDetIDRange[2] = {0};

        TCutG* fCutG_Theta_ETot = nullptr;
        TCutG* fCutG_ETot_dE = nullptr;
        TCutG* fCutG_E_dE = nullptr;
        int fCutDetID = -1;
        int fCutdEDetID = -1;

        int fKeyEnergyType = 0;

    ClassDef(SKGatingTask, 1)
};

#endif
