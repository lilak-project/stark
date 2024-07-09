#ifndef STARK_HH
#define STARK_HH

#include "LKLogger.h"
#include "LKDetector.h"

class STARK : public LKDetector
{
    public:
        STARK();
        virtual ~STARK() { ; }
        virtual bool BuildDetectorPlane();

    ClassDef(STARK,1);
};

#endif
