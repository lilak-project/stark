void run_reco()
{
    auto run = new LKRun();
    run -> AddDetector(new STARK());
    run -> AddPar("config_reco.mac");
    run -> Add(new SKSetSiChannelTask());
    run -> Add(new SKEnergyRestorationTask());
    run -> Add(new SKPairMatchingTask());
    run -> Init();
    run -> Run();
}
