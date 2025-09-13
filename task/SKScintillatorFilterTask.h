#ifndef SKSCINTILLATORFILTERTASK_HH
#define SKSCINTILLATORFILTERTASK_HH

#include "TClonesArray.h"
#include "LKTask.h"
#include "SKSiArrayPlane.h"
#include "SKRecoHeader.h"

class SKScintillatorFilterTask : public LKTask
{
    public:
        SKScintillatorFilterTask();
        virtual ~SKScintillatorFilterTask() {}

        bool Init();
        void Exec(Option_t*);

    private:
        TClonesArray *fRawDataArray = nullptr;
        TClonesArray *fRecoHeaderArray = nullptr;
        SKRecoHeader *fRecoHeader = nullptr;

        int fTbRange1 = 300;
        int fTbRange2 = 350;
        int fThreshold = 100;

    ClassDef(SKScintillatorFilterTask, 1)
};

#endif
