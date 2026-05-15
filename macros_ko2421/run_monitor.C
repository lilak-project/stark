void run_monitor()
{
    auto run = new LKRun();
    run -> AddDetector(new STARK());
    //run -> AddInputFile("data/stark_0011.reco.root");
    //run -> AddInputFile("data/stark_0012.reco.root");
    //run -> AddInputFile("data/stark_0013.reco.root");
    //run -> AddInputFile("data/stark_0014.reco.root");
    run -> AddPar("config_monitor.mac");
    run -> Add(new SKAnalysisJW);
    run -> Add(new SKAnalysisDK);
    run -> Add(new SKAnalysisTA);
    //run -> Add(new SKDrawWaveformTask);
    //run -> Add(new SKDrawEventStatisticsTask);
    run -> Add(new LKEveTask);
    run -> Init();
    auto plane = (SKSiArrayPlane*) run -> GetDetectorPlane(0);
    //plane -> SetAllowControlLogger(false);
    //run -> SetEntriesLimit(100000);
    run -> ExecuteNextEvent();
    //plane -> RunAllEvents();
}
