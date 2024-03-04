void run_conv()
{
    auto run = new LKRun();
    run -> AddPar("config.mac");
    run -> Add(new LKMFMConversionTask());
    run -> Init();
    run -> Run();
}
