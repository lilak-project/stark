void run_merger()
{
    auto run = new LKRun();
    run -> AddPar("config_merger.mac");
    run -> Add(new LKMTEMergerTask());
    //run -> InitAndCollectParameters(); return;
    run -> Init();
    run -> Run();
}
