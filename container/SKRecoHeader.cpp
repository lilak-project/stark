#include "SKRecoHeader.h"

ClassImp(SKRecoHeader);

SKRecoHeader::SKRecoHeader()
{
    Clear();
}

void SKRecoHeader::Clear(Option_t *option)
{
    fBeamOnScint = false;
}
