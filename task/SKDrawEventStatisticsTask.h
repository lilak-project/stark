#ifndef SKDRAWEVENTSTATISTICSTASK_HH
#define SKDRAWEVENTSTATISTICSTASK_HH

#include "TH1.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TClonesArray.h"

#include "LKTask.h"
#include "SKSiArrayPlane.h"

class SKDrawEventStatisticsTask : public LKTask
{
    public:
        SKDrawEventStatisticsTask();
        virtual ~SKDrawEventStatisticsTask() {}

        bool Init();
        void Exec(Option_t*);

        SKSiArrayPlane* fStarkPlane = nullptr;
        TClonesArray* fSiChannelArray = nullptr;

    private:
        TH2D* fHistET[2];
        TH2D* fHistHP[2];
        TH2D* fHistEP[40][8];
        TH2D* fHistLR[40][8];
        // pid
        // raw spectrum
        int fMaxChannels[2] = {0};
        bool fSinglePlot = true;

    ClassDef(SKDrawEventStatisticsTask, 1)
};

#endif
