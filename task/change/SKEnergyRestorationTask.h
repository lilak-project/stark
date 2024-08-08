#ifndef SKENERGYRESTORATIONTASK_HH
#define SKENERGYRESTORATIONTASK_HH

#include "TClonesArray.h"
#include "LKTask.h"
#include "SKSiArrayPlane.h"

class SKEnergyRestorationTask : public LKTask
{
    public:
        SKEnergyRestorationTask();
        virtual ~SKEnergyRestorationTask() {}

        bool Init();
        void Exec(Option_t*);

    private:
        SKSiArrayPlane* fStarkPlane = nullptr;

        TClonesArray *fSiChannelArray = nullptr;
        TClonesArray *fHitArray = nullptr;

        double f241AmAlphaEnergy1 = 5.486;

        double fg0Array[40][2][8];
        double fg1Array[40][2][8];
        double fg2Array[40][2][8];
        double fb0Array[40][2][8];
        double fb1Array[40][2][8];
        double fb2Array[40][2][8];

    ClassDef(SKEnergyRestorationTask, 1)
};

#endif
