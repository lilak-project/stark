void read_conv()
{
    ofstream txtfile("dump_data.txt");

    auto file = new TFile("data/stark_0044.D0.conv.root");
    auto tree = (TTree *) file -> Get("event");

    TClonesArray *channelArray = nullptr;
    tree -> SetBranchAddress("RawData", &channelArray);

    auto numEvents = tree -> GetEntries();
    for (auto iEvent=0; iEvent<numEvents; ++iEvent)
    {
        tree -> GetEntry(iEvent); // channelArray update for event number = iEvent

        auto numChannels = channelArray -> GetEntries();
        for (auto iChannel=0; iChannel<numChannels; ++iChannel)
        {
            auto channel = (GETChannel*) channelArray -> At(iChannel);
            //channel -> Print(); return;
            //channel -> Draw(); return;
            auto cobo = channel -> GetCobo();
            auto asad = channel -> GetAsad();
            auto aget = channel -> GetAget();
            auto chan = channel -> GetChan();
            auto time = channel -> GetTime();
            auto energy = channel -> GetEnergy();
            auto pedestal = channel -> GetPedestal();
            int* buffer = channel -> GetWaveformY();

            cout << iEvent << " " << iChannel << " " << numChannels << " " << cobo << " " << asad << " " << aget << " " << chan << " " << time << " " << energy << " " << pedestal << endl;
            //txtfile << iEvent << " " << iChannel << " " << numChannels << " " << cobo << " " << asad << " " << aget << " " << chan << " " << time << " " << energy << " " << pedestal << endl;
        }
    }

}
