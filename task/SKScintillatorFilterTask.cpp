#include "SKScintillatorFilterTask.h"

ClassImp(SKScintillatorFilterTask)

SKScintillatorFilterTask::SKScintillatorFilterTask()
    :LKTask("SKScintillatorFilterTask","SKScintillatorFilterTask")
{
}

bool SKScintillatorFilterTask::Init()
{
    fRawDataArray = fRun -> GetBranchA("RawData");
    fRecoHeaderArray = fRun -> RegisterBranchA("RecoHeader","SKRecoHeader",1);

    if (fPar->CheckPar(fName+"/tbRange1"))
    {
        fTbRange1 = fPar -> GetParInt(fName+"/tbRange",0);
        fTbRange2 = fPar -> GetParInt(fName+"/tbRange",1);
    }
    fPar -> UpdatePar(fThreshold, fName+"/threshold");

    return true;
}

void SKScintillatorFilterTask::Exec(Option_t*)
{
    fRecoHeaderArray -> Clear("C");
    fRecoHeader = (SKRecoHeader*) fRecoHeaderArray -> ConstructedAt(0);
    fRecoHeader -> SetBeamOnScint(false);

    auto numChannels = fRawDataArray -> GetEntries();
    for (auto iChannel=0; iChannel<numChannels; ++iChannel)
    {
        auto channel = (GETChannel*) fRawDataArray -> At(iChannel);
        if (!(channel->GetCobo()==0 && channel->GetAsad()==2 && channel->GetAget()==0 && channel->GetChan()==65))
            continue;

        double min =  1000000;
        double max = -1000000;
        auto data = channel -> GetWaveformY();
        for (int tb=fTbRange1; tb<fTbRange2; ++tb)
        {
            double value = data[tb];
            if (value<min) min = value;
            if (value>max) max = value;
        }
        double height = max - min;
        if (height>fThreshold)
            fRecoHeader -> SetBeamOnScint(true);
    }

    if (fRecoHeader -> BeamOnScint()) lk_info << "Beam on scintillator" << endl;
    else lk_info << "No beam correlation" << endl;
}
