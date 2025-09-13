#ifndef SKAnalysisJW_HH
#define SKAnalysisJW_HH

#include "TH1.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TClonesArray.h"

#include "LKTask.h"
#include "SKSiArrayPlane.h"

class SKAnalysisJW : public LKTask
{
    public:
        SKAnalysisJW();
        virtual ~SKAnalysisJW() {}

        bool Init();
        void Exec(Option_t*);

        SKSiArrayPlane* fStarkPlane = nullptr;
        TClonesArray* fSiHitArray = nullptr;

    private:
        TH2D* fHistHP[2];
        TH1D* fHistCount[2];
        const int fNumJStrips = 8;
        const int fNumOStrips = 4;

        double fEThresholdForCount = 0;

    ClassDef(SKAnalysisJW, 1)
};

#endif
