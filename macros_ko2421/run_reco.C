void run_reco()
{
    auto run = new LKRun();
    run -> Add(new ELARK());
    run -> AddPar("config_reco.mac");
    run -> Add(new EKSetSiChannelTask);
    run -> Add(new EKEnergyRestorationTask);
    run -> Add(new EKReconstructHitTask);
    run -> Add(new EKShortAnalysisRecoHit);
    run -> Init();
    run -> Run();
}
