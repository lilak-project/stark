#ifndef SKPAIRMATCHINGTASK_HH
#define SKPAIRMATCHINGTASK_HH

#include "TClonesArray.h"
#include "LKTask.h"
#include "SKSiArrayPlane.h"

class SKPairMatchingTask : public LKTask
{
    public:
        SKPairMatchingTask();
        virtual ~SKPairMatchingTask() {}

        bool Init();
        void Exec(Option_t*);

    private:
        SKSiArrayPlane* fStarkPlane = nullptr;

        TClonesArray *fHitArray = nullptr;

    ClassDef(SKPairMatchingTask, 1)
};

#endif
