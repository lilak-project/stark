#ifndef SKANALYSISTA_HH
#define SKANALYSISTA_HH

#include "TH1.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TClonesArray.h"
#include "TCutG.h"

#include "LKTask.h"
#include "SKSiArrayPlane.h"

class SKAnalysisTA : public LKTask
{
    public:
        SKAnalysisTA();
        virtual ~SKAnalysisTA() {}

        bool Init();
        void Exec(Option_t*);
        double GetQvalue(double beamenergy, double theta, double energy);
       	double GetDiffSig(double ring_idmax, double beamrate, double ntarget, double runtime, double theta);
	void DrawDiffSig();

        SKSiArrayPlane* fStarkPlane = nullptr;
        TClonesArray* fSiHitArray = nullptr;

    private:
        TCutG *cutg;
	TH2D* fHistAvsE[40][2];
	TH2D* fHistAvstotE[2];
	TH2D* fHistdEEangle[12][2];
	TH1D* fHistAcomgonel;
	TH1D* fHistDiffSiggonel;
	TH1D* fHistDiffSigfresco;
        const int fNumJStrips = 8;
        const int fNumOStrips = 4;

        double fBeamEnergy;
        double fBeamRate;
        double fRunTime;
        double fNTarget;

    ClassDef(SKAnalysisTA, 1)
};

#endif
