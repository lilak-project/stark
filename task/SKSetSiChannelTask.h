#ifndef SKSETSICHANNELTASK_HH
#define SKSETSICHANNELTASK_HH

#include "TH1.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TVirtualPad.h"
#include "TClonesArray.h"

#include "LKTask.h"
#include "SKSiArrayPlane.h"

#define fNumX6    36
#define fNumX6Ost 4
#define fNumX6Jst 8
#define fNumSD    12
#define fNumSDJst 8
#define fNP       12

class SKSetSiChannelTask : public LKTask
{
    public:
        SKSetSiChannelTask();
        virtual ~SKSetSiChannelTask() {}

        bool Init();
        void Exec(Option_t*);

    private:
        SKSiArrayPlane* fStarkPlane = nullptr;
        LKChannelAnalyzer* fChannelAnalyzer = nullptr;

    private:
        TClonesArray *fRawDataArray = nullptr;
        TClonesArray *fSiChannelArray = nullptr;

    ClassDef(SKSetSiChannelTask, 1)
};

#endif
