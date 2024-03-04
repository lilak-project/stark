void run_viewer()
{
    auto run = new LKRun();
    run -> AddPar("config.mac");
    run -> AddInputFile("stark_0044.D0.conv.root");
    run -> Add(new LKGETChannelViewer);

    run -> Init();
    run -> ExecuteNextEvent();
}
