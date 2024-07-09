#ifndef SKDRAWWAVEFORMTASK_HH
#define SKDRAWWAVEFORMTASK_HH

#include "TH1.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TClonesArray.h"

#include "LKTask.h"
#include "SKSiArrayPlane.h"

class SKDrawWaveformTask : public LKTask
{
    public:
        SKDrawWaveformTask();
        virtual ~SKDrawWaveformTask() {}

        bool Init();
        void Exec(Option_t*);

        SKSiArrayPlane* fStarkPlane = nullptr;
        TClonesArray* fSiChannelArray = nullptr;

    private:
        int fNumDetectors = 0;
        TH1D* fHistEnergy[48][2][16][2];

    ClassDef(SKDrawWaveformTask, 1)
};

#endif
