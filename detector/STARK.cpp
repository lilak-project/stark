#include "STARK.h"
#include "SKSiArrayPlane.h"

ClassImp(STARK);

STARK::STARK()
    :LKDetector("stark","Active target TPC for Multiple nuclear eXperiment")
{
}

bool STARK::BuildDetectorPlane()
{
    AddPlane(new SKSiArrayPlane);
    return true;
}
