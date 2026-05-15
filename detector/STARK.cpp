#include "STARK.h"
#include "LKSiliconArray.h"

ClassImp(STARK);

STARK::STARK()
    :LKDetector("stark","")
{
}

bool STARK::BuildDetectorPlane()
{
    AddPlane(new LKSiliconArray);
    return true;
}
