void check_channel_pulser()
{
    auto file = new TFile("data_conv/stark_0099.conv.root");
    auto tree = (TTree*) file -> Get("event");
    TClonesArray* channelArray = nullptr;//GETChannel
    tree -> SetBranchAddress("RawData", &channelArray);

    auto top = new LKDrawingGroup("la");

    int countEvents = 0;

    auto numEvents = tree -> GetEntries();
    for (auto iEvent=0; iEvent<numEvents; ++iEvent)
    {
        tree -> GetEntry(iEvent);
        auto numChannels = channelArray -> GetEntries();

        bool good = false;
        if (numChannels==3)
        {
            for (auto iChannel=0; iChannel<numChannels; ++iChannel)
            {
                auto channel = (GETChannel*) channelArray -> At(iChannel);
                if (channel -> GetCobo()==0 && channel -> GetAsad()==0) {
                    good = true;
                    break;
                }
            }
        }

        if (good)
        {
            cout << iEvent << " " << numChannels << endl;
            auto group = top -> CreateGroup(Form("e%d",iEvent));
            for (auto iChannel=0; iChannel<numChannels; ++iChannel)
            {
                auto channel = (GETChannel*) channelArray -> At(iChannel);
                if (channel -> GetCobo()==0 && channel -> GetAsad()==0)
                    good = true;
                //cout << channel -> GetCAAC() << " " << channel -> GetEnergy() << endl;
                auto hist = channel -> GetHist(Form("h_%d_%d",iEvent,channel->GetCAAC()));
                cout << hist -> GetName() << " " << hist -> GetMean() << endl;
                group -> AddHist(hist);
            }
            if (countEvents++>40)
                break;
        }
    }

    top -> Draw("v");
}
