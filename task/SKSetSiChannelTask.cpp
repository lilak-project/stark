#include "TStyle.h"

#include "LKRun.h"
#include "LKLogger.h"
#include "GETChannel.h"
#include "LKSiChannel.h"

#include "SKSetSiChannelTask.h"

ClassImp(SKSetSiChannelTask)

SKSetSiChannelTask::SKSetSiChannelTask()
    :LKTask("SKSetSiChannelTask","SKSetSiChannelTask")
{
}

bool SKSetSiChannelTask::Init()
{
    fStarkPlane = (SKSiArrayPlane*) fRun -> FindDetectorPlane("SKSiArrayPlane");
    fRawDataArray = fRun -> GetBranchA("RawData");
    fSiChannelArray = fRun -> RegisterBranchA("SiChannel","LKSiChannel",20);

    fChannelAnalyzer = new LKChannelAnalyzer();
    fChannelAnalyzer -> Print();

    return true;
}

void SKSetSiChannelTask::Exec(Option_t*)
{
    fSiChannelArray -> Clear("C");

    int channelCount = 0;
    auto numChannels = fRawDataArray -> GetEntries();
    for (auto iChannel=0; iChannel<numChannels; ++iChannel)
    {
        auto channel = (GETChannel*) fRawDataArray -> At(iChannel);
        auto siChannel1 = (LKSiChannel*) fSiChannelArray -> ConstructedAt(channelCount);
        siChannel1 -> SetChannelID(channelCount);
        bool good = fStarkPlane -> SetSiChannelData(siChannel1, channel);
        if (good) {
            for (auto iSiChannel=0; iSiChannel<numChannels; ++iSiChannel)
            {
                if (iSiChannel==channelCount)
                    continue;
                auto siChannel2 = (LKSiChannel*) fSiChannelArray -> ConstructedAt(iSiChannel);
                if (siChannel1 -> IsPair(siChannel2)) {
                    siChannel1 -> SetPairChannel(siChannel2);
                    siChannel2 -> SetPairChannel(siChannel1);
                    siChannel1 -> SetPairArrayIndex(siChannel2 -> GetChannelID());
                    siChannel2 -> SetPairArrayIndex(siChannel1 -> GetChannelID());
                    // TODO
                    //siChannel1 -> Set
                }
            }
            channelCount++;
        }
        else
            fSiChannelArray -> RemoveAt(channelCount);

        if (siChannel1->GetSide()==0)
            fChannelAnalyzer -> SetDataIsInverted(true);
        else
            fChannelAnalyzer -> SetDataIsInverted(false);

        auto data = channel -> GetWaveformY();
        fChannelAnalyzer -> Analyze(data);
        auto numRecoHits = fChannelAnalyzer -> GetNumHits();
        if (numRecoHits>=1) {
            auto energy = fChannelAnalyzer -> GetAmplitude(0);
            siChannel1 -> SetEnergy(energy);
        }
    }

    lk_info << channelCount << " si-channels found" << endl;
}
