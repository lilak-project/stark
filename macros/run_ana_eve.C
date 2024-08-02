void run_ana_eve()
{
    auto run = new LKRun();
    run -> AddDetector(new STARK());
    run -> AddPar("config_ana_eve.mac");
    run -> Add(new SKDrawCalibratedEventStatisticsTask);
    run -> Add(new LKEveTask);
    run -> Init();
    run -> ExecuteNextEvent();
}
