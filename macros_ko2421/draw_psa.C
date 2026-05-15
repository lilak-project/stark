void draw_psa()
{
    auto run = new LKRun();
    run -> SetInputFile("data_conv/stark_0105.conv.root");
    run -> Init();
    run -> Print();

    auto ana = new LKChannelAnalyzer();
    ana -> SetThreshold(400);

    auto top = new LKDrawingGroup("psa");

    auto channelArray = run -> GetBranchA("RawData");
    while (run -> GetNextEvent())
    {
        auto eventID = run->GetCurrentEventID();
        auto numChannels = channelArray -> GetEntries();

        for (auto iChannel=0; iChannel<numChannels; ++iChannel)
        {
            auto channel = (GETChannel*) channelArray -> At(iChannel);

            ana -> Analyze(channel->GetBufferArray());
            auto numHits = ana -> GetNumHits();
            if (numHits!=1)
                continue;

            auto tbReco = ana -> GetTbHit(0);
            auto ampReco = ana -> GetAmplitude(0);
            auto pdReco = ana -> GetPedestal();

            //auto hist = channel -> GetHist(Form("channel_%d_%d",eventID,iChannel));
            //top -> AddHist(hist);

            auto draw = ana -> GetDrawing();
            auto lg0 = draw -> FindObjectNameClass("",TLegend::Class());
            draw -> Remove(lg0);
            //auto lg = new TLegend();
            //lg -> AddEntry((TObject*)0,Form("tb = %d -> %.1f",tbSim,tbReco),"");
            //lg -> AddEntry((TObject*)0,Form("amp = %d -> %.1f",ampSim,ampReco),"");
            //lg -> AddEntry((TObject*)0,Form("pd = %d -> %.1f",pdSim,pdReco),"");
            //draw -> Add(lg);
            //draw -> SetLegendCorner(1,0.45,0.08);
            top -> AddDrawing(draw);
        }
        break;
    }

    top -> Draw();
}
