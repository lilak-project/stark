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
        if (!good) {
            fSiChannelArray -> RemoveAt(channelCount);
            continue;
        }

        if (siChannel1->GetSide()==0)
            fChannelAnalyzer -> SetDataIsInverted(true);
        else
            fChannelAnalyzer -> SetDataIsInverted(false);

        auto data = channel -> GetWaveformY();
        fChannelAnalyzer -> Analyze(data);
        auto numRecoHits = fChannelAnalyzer -> GetNumHits();
        if (numRecoHits>=1) {
            auto pedestal = fChannelAnalyzer -> GetPedestal();
            auto energy = fChannelAnalyzer -> GetAmplitude(0);
            auto time = fChannelAnalyzer -> GetTbHit(0);
            siChannel1 -> SetPedestal(pedestal);
            siChannel1 -> SetEnergy(energy);
            siChannel1 -> SetTime(time);
            channel -> SetPedestal(pedestal);
            channel -> SetEnergy(energy);
            channel -> SetTime(time);
        }

        channelCount++;
        for (auto iSiChannel2=0; iSiChannel2<channelCount; ++iSiChannel2)
        {
            if (iSiChannel2==channelCount)
                continue;
            auto siChannel2 = (LKSiChannel*) fSiChannelArray -> At(iSiChannel2);
            if (siChannel1 -> IsPair(siChannel2)) {
                if (siChannel1->GetDirection()==0)
                {
                    siChannel1 -> SetPairChannel(siChannel2);
                    //siChannel2 -> SetPairChannel(siChannel1);
                    siChannel1 -> SetPairArrayIndex(siChannel2 -> GetChannelID());
                    //siChannel2 -> SetPairArrayIndex(siChannel1 -> GetChannelID());
                    siChannel1 -> SetEnergy2(siChannel2->GetEnergy());
                }
                else
                {
                    siChannel2 -> SetPairChannel(siChannel1);
                    //siChannel1 -> SetPairChannel(siChannel2);
                    siChannel2 -> SetPairArrayIndex(siChannel1 -> GetChannelID());
                    //siChannel1 -> SetPairArrayIndex(siChannel2 -> GetChannelID());
                    siChannel2 -> SetEnergy2(siChannel1->GetEnergy());
                }
                if (siChannel1->GetEnergy()==siChannel2->GetEnergy())
                    lk_debug << "1=(" << siChannel1->GetLocalID() << ") " << siChannel1->GetEnergy() << " 2=(" << siChannel2->GetLocalID() << ") " << siChannel2->GetEnergy() << endl;
            }
        }
    }

    lk_info << channelCount << " si-channels found" << endl;
}
