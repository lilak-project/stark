void run_eve()
{
    auto run = new LKRun();
    run -> AddDetector(new STARK());
    run -> AddPar("config_eve.mac");
    run -> Add(new SKSetSiChannelTask());
    run -> Add(new SKDrawEventStatisticsTask);
    run -> Add(new SKDrawWaveformTask);
    run -> Add(new LKEveTask);
    run -> Init();
    run -> ExecuteNextEvent();
}
