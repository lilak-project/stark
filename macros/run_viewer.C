void run_viewer()
{
    auto run = new LKRun();
    run -> AddPar("config_viewer.mac");
    run -> AddInputFile("data/stark_0044.D0.conv.root");
    run -> Add(new LKGETChannelViewer);
    run -> Init();
    run -> ExecuteNextEvent();
}

