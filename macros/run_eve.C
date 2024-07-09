void run_eve()
{
    auto run = new LKRun();
    run -> AddDetector(new STARK());
    run -> AddPar("config_eve.mac");
    run -> AddEveTask(new SKSetSiChannelTask);
    run -> AddEveTask(new SKDrawEventStatisticsTask);
    run -> AddEveTask(new SKDrawWaveformTask);
    run -> Add(new LKEveTask);
    run -> Init();
    run -> ExecuteNextEvent();
}
