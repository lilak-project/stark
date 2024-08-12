#include "TStyle.h"

#include "LKRun.h"
#include "LKLogger.h"
#include "GETChannel.h"
#include "LKSiChannel.h"

#include "SKSetSiChannelTask.h"

#include "LKPulseFitData.h"

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
    fFitDataArray = fRun -> RegisterBranchA("PFData","LKPulseFitData",10,false);

    TString pulseFileName;
    fPar -> UpdatePar(pulseFileName,"stark/pulseFile");

    fChannelAnalyzer = new LKChannelAnalyzer();
    //fChannelAnalyzer -> SetPulse(pulseFileName);
    //fChannelAnalyzer -> Print();

    fChannelAnalyzer2 = new LKChannelAnalyzer();
    fChannelAnalyzer2 -> SetPulse(pulseFileName);
    //fChannelAnalyzer2 -> Print();

    fSlopeFit = new TF1("SlopeFit","[3]+(x>[0]&&x<[1])*([2]/([1]-[0]))*(x-[0])+(x>[1]&&x<[1]+4)*[2]",0,512);

    fHistBuffer = new TH1D("histBufferSSSC","",512,0,512);

    return true;
}

void SKSetSiChannelTask::Exec(Option_t*)
{
    fSiChannelArray -> Clear("C");

    int fitDataCount = 0;
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

        bool isInverted = true;
        if (siChannel1->GetSide()==0) isInverted = true;
        else isInverted = false;

        fChannelAnalyzer -> SetDataIsInverted(isInverted);
        fChannelAnalyzer2 -> SetDataIsInverted(isInverted);

        auto data = channel -> GetWaveformY();
        fChannelAnalyzer -> Analyze(data);
        auto numRecoHits = fChannelAnalyzer -> GetNumHits();
        if (numRecoHits>=1)
        {
            auto pedestal = fChannelAnalyzer -> GetPedestal();
            auto energy = fChannelAnalyzer -> GetAmplitude(0);
            auto time = fChannelAnalyzer -> GetTbHit(0);
            double energy0 = energy;

            bool isSaturated = false;
            if (isInverted)
            {
                for (auto t=0; t<512; ++t) {
                    if (data[t]==0) {
                        isSaturated = true;
                        break;
                    }
                }
            }
            else {
                for (auto t=0; t<512; ++t) {
                    if (data[t]==4095) {
                        isSaturated = true;
                        break;
                    }
                }
            }

            double t1 = 0;
            double t2 = 512;
            double slope = 0;
            double a1 = 0;
            double a0 = 0;

            if (isSaturated)
            {
                t1 = time-10; if (t1<0) t1 = 0;
                t2 = time+5; if (t2>512) t2 = 512;
                fSlopeFit -> SetRange(t1,t2);
                fSlopeFit -> SetParameters(time-5,time,energy,pedestal);
                if (isInverted) {
                    for (auto tb=0; tb<512; ++tb)
                        fHistBuffer -> SetBinContent(tb+1,4095-data[tb]);
                }
                else {
                    for (auto tb=0; tb<512; ++tb)
                        fHistBuffer -> SetBinContent(tb+1,data[tb]);
                }
                fHistBuffer -> Fit(fSlopeFit,"QN0");
                a0 = fSlopeFit -> GetParameter(0);
                a1 = fSlopeFit -> GetParameter(1);
                double a2 = fSlopeFit -> GetParameter(2);
                slope = a2 / (a1-a0);
                double energy2 = -1.90397 + 9.39922*slope;
                //lk_debug << "bad: " << channel -> GetCAAC() << " : " << energy << " -> " << energy2 << endl;
                energy = energy2;
            }
            else {
                //lk_debug << "goo: " << channel -> GetCAAC() << " : " << energy << endl;
            }

            auto fitData = (LKPulseFitData*) fFitDataArray -> ConstructedAt(fitDataCount++);
            fitData -> fHitIndex = channelCount;
            fitData -> fIsSaturated = isSaturated;
            fitData -> fNumHitsInChannel = 1;
            //fitData -> fNDF;
            //fitData -> fChi2;
            fitData -> fFitRange1 = t1;
            fitData -> fFitRange2 = t2;
            fitData -> fTb = time;
            //fitData -> fWidth;
            //fitData -> fIntegral;
            fitData -> fAmplitude = energy0;
            fitData -> fSlope = slope;
            fitData -> fSlopePar0 = a0;
            fitData -> fSlopePar1 = a1;
            fitData -> fSlopeAmplitude = energy;

            siChannel1 -> SetPedestal(pedestal);
            siChannel1 -> SetEnergy(energy);
            siChannel1 -> SetTime(time);
            channel -> SetPedestal(pedestal);
            channel -> SetEnergy(energy);
            channel -> SetTime(time);
        }

        channelCount++;
        if (siChannel1 -> IsPairedChannel())
        {
            for (auto iSiChannel2=0; iSiChannel2<channelCount; ++iSiChannel2)
            {
                if (iSiChannel2==channelCount)
                    continue;
                auto siChannel2 = (LKSiChannel*) fSiChannelArray -> At(iSiChannel2);
                if (siChannel1 -> IsPair(siChannel2))
                {
                    if (siChannel1->GetDirection()==0)
                    {
                        siChannel1 -> SetPairChannel(siChannel2);
                        siChannel1 -> SetPairArrayIndex(siChannel2 -> GetChannelID());
                        siChannel1 -> SetEnergy2(siChannel2->GetEnergy());
                    }
                    else
                    {
                        siChannel2 -> SetPairChannel(siChannel1);
                        siChannel2 -> SetPairArrayIndex(siChannel1 -> GetChannelID());
                        siChannel2 -> SetEnergy2(siChannel1->GetEnergy());
                    }
                    //if (siChannel1->GetEnergy()==siChannel2->GetEnergy())
                        //lk_debug << "1=(" << siChannel1->GetLocalID() << ") " << siChannel1->GetEnergy() << " 2=(" << siChannel2->GetLocalID() << ") " << siChannel2->GetEnergy() << endl;
                }
            }
        }
    }

    lk_info << channelCount << " si-channels found" << endl;
}
