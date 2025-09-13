#ifndef SKANALYSISDK_HH
#define SKANALYSISDK_HH

#include "TH1.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TClonesArray.h"

#include "LKTask.h"
#include "SKSiArrayPlane.h"

class SKAnalysisDK : public LKTask
{
    public:
        SKAnalysisDK();
        virtual ~SKAnalysisDK() {}

        bool Init();
        void Exec(Option_t*);

        SKSiArrayPlane* fStarkPlane = nullptr;
        TClonesArray* fSiHitArray = nullptr;

    private:
        TH2D* fHistdEEAll[2];
        TH2D* fHistHP[2];
	TH2D* fHistZtotE[2];
	TH2D* fHistZE[40][2];
	TH2D* fHistdEE[12][2];
	TH2D* fHistZdet[2];
        const int fNumJStrips = 8;
        const int fNumOStrips = 4;

    ClassDef(SKAnalysisDK, 1)
};

#endif
