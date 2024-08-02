#ifndef SKDRAWCALIBRATEDEVENTSTATISTICSTASK_HH
#define SKDRAWCALIBRATEDEVENTSTATISTICSTASK_HH

#include "TH1.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TClonesArray.h"

#include "LKTask.h"
#include "SKSiArrayPlane.h"

class SKDrawCalibratedEventStatisticsTask : public LKTask
{
    public:
        SKDrawCalibratedEventStatisticsTask();
        virtual ~SKDrawCalibratedEventStatisticsTask() {}

        bool Init();
        void Exec(Option_t*);

        SKSiArrayPlane* fStarkPlane = nullptr;
        TClonesArray* fSiHitArray = nullptr;

    private:
        TH2D* fHistHP[2];
        const int fNumJStrips = 8;
        const int fNumOStrips = 4;

    ClassDef(SKDrawCalibratedEventStatisticsTask, 1)
};

#endif
