void run_conv()
{
    auto run = new LKRun();
    run -> AddPar("config_conv.mac");
    run -> SetEventTrigger(new LKMFMConversionTask());
    //run -> InitAndCollectParameters();
    run -> Init();
    run -> Print();
    run -> Run();
}
