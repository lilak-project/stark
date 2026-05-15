#ifndef SKRECOHEADER_HH
#define SKRECOHEADER_HH

#include "LKContainer.h"

class SKRecoHeader : public LKContainer
{
    public:
        SKRecoHeader();
        virtual ~SKRecoHeader() { ; }
        virtual void Clear(Option_t *option="");
        virtual void Print(Option_t *option="") const {}

        void SetBeamOnScint(bool val) { fBeamOnScint = val; }
        bool BeamOnScint() const { return fBeamOnScint; }

    protected:
        bool fBeamOnScint = true;

    ClassDef(SKRecoHeader,3);
};

#endif
