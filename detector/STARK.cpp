#include "STARK.h"
#include "SKSiArrayPlane.h"

ClassImp(STARK);

STARK::STARK()
    :LKDetector("stark","")
{
}

bool STARK::BuildDetectorPlane()
{
    AddPlane(new SKSiArrayPlane);
    return true;
}
