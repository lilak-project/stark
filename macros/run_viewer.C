void run_viewer()
{
    auto run = new LKRun();
    run -> AddPar("config_viewer.mac");
    //run -> AddInputFile("stark_0044.D0.conv.root");
    run -> Add(new LKGETChannelViewer);
    //run -> Add(new SKRunViewer);
    //run -> InitAndCollectParameters(); return;
    run -> Init();
    //run -> Run();
}
