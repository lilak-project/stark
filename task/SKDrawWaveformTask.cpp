#include "TStyle.h"

#include "LKRun.h"
#include "LKLogger.h"
#include "GETChannel.h"
#include "LKSiChannel.h"
#include "SKDrawWaveformTask.h"

ClassImp(SKDrawWaveformTask)

SKDrawWaveformTask::SKDrawWaveformTask()
    :LKTask("SKDrawWaveformTask","SKDrawWaveformTask")
{
}

bool SKDrawWaveformTask::Init()
{
    fSiChannelArray = fRun -> GetBranchA("SiChannel");

    fStarkPlane = (SKSiArrayPlane*) fRun -> FindDetectorPlane("SKSiArrayPlane");
    if (fStarkPlane==nullptr) {
        lk_error << "SKSiArrayPlane do not exist!!!" << endl;
        lk_error << fName << " must run with SKSiArrayPlane" << endl;
        return false;
    }

    fNumDetectors = fStarkPlane -> GetNumSiDetectors();
    for (int iDetector=0; iDetector<fNumDetectors; ++iDetector)
    {
        auto siDetector = fStarkPlane -> GetSiDetector(iDetector);
        TString nameDetector = siDetector -> GetNameType();
        TString ttleDetector = siDetector -> GetTitleType();
        int detID = siDetector -> GetDetID();
        int numSides = siDetector -> GetNumSides();
        int numJStrips = siDetector -> GetNumJunctionStrips();
        int numOStrips = siDetector -> GetNumOhmicStrips();
        int numJDirection = siDetector -> GetNumJunctionDirection();
        int numODirection = siDetector -> GetNumOhmicDirection();
        auto array = new TObjArray();
        for (int side=0; side<numSides; ++side)
        {
            int numStrips = (side==0?numJStrips:numOStrips);
            int numDirection = (side==0?numJDirection:numODirection);
            for (int strip=0; strip<numStrips; ++strip)
            {
                TString nameStrip = Form("%s_%s%d",  nameDetector.Data(),(side==0?"J"       :"O"    ),strip);
                TString ttleStrip = Form("%s %s(%d)",ttleDetector.Data(),(side==0?"Junction":"Ohmic"),strip);
                for (int direction=0; direction<numDirection; ++direction)
                {
                    TString nameChannel = Form("%s%s" ,nameStrip.Data(),(side==0?(direction==0?"U" :"D"   ):(direction==0?"L"   :"R"    )));
                    TString ttleChannel = Form("%s %s",ttleStrip.Data(),(side==0?(direction==0?"Up":"Down"):(direction==0?"Left":"Right")));
                    if (numDirection==1) {
                        nameChannel = nameStrip;
                        ttleChannel = ttleStrip;
                    }
                    ttleChannel = ttleChannel + ";time-bucket;ADC";
                    fHistEnergy[iDetector][side][strip][direction] = new TH1D(nameChannel,ttleChannel,512,0,512);
                    auto hist = fHistEnergy[iDetector][side][strip][direction];
                    hist -> SetStats(0);
                    array -> Add(hist);
                }
            }
        }
        fStarkPlane -> AddUserDrawings("waveform", iDetector, -1, array, 6);
    }

    return true;
}

void SKDrawWaveformTask::Exec(Option_t*)
{
    if (fStarkPlane -> GetAccumulateEvents())
        return;

    for (int iDetector=0; iDetector<fNumDetectors; ++iDetector)
    {
        auto siDetector = fStarkPlane -> GetSiDetector(iDetector);
        int numSides = siDetector -> GetNumSides();
        int numJStrips = siDetector -> GetNumJunctionStrips();
        int numOStrips = siDetector -> GetNumOhmicStrips();
        int numJDirection = siDetector -> GetNumJunctionDirection();
        int numODirection = siDetector -> GetNumOhmicDirection();
        for (int side=0; side<numSides; ++side)
        {
            int numStrips = (side==0?numJStrips:numOStrips);
            int numDirection = (side==0?numJDirection:numODirection);
            for (int strip=0; strip<numStrips; ++strip)
            {
                for (int direction=0; direction<numDirection; ++direction)
                {
                    fHistEnergy[iDetector][side][strip][direction] -> Reset();
                }
            }
        }
    }

    auto numChannels = fSiChannelArray -> GetEntries();
    for (auto iChannel=0; iChannel<numChannels; ++iChannel)
    {
        auto siChannel = (LKSiChannel*) fSiChannelArray -> At(iChannel);
        if (siChannel==nullptr)
            continue;
        int detID = siChannel -> GetDetID();
        int side = siChannel -> GetSide();
        int strip = siChannel -> GetStrip();
        int direction = siChannel -> GetDirection();
        auto hist = fHistEnergy[detID][side][strip][direction];
        siChannel -> FillHist(hist);
    }
}
