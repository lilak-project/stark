void run_viewer()
{
    auto run = new LKRun();
    run -> AddPar("config_viewer.mac");
    //run -> AddInputFile("stark_0044.D0.conv.root");
    auto viewer = new LKGETChannelViewer;
    viewer -> SelectCAAC(0,2,3,57);
    run -> Add(viewer);
    //run -> Add(new SKRunViewer);
    //run -> InitAndCollectParameters(); return;
    run -> Init();
    //run -> Run();
}
