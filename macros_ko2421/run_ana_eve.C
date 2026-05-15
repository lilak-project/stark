void run_ana_eve()
{
    auto run = new LKRun();
    run -> Add(new ELARK());
    run -> AddPar("config_ana_eve.mac");
    run -> Add(new EKAnalysisRecoHit);
    run -> Init();

    run -> SetAutoTermination(false);
    //run -> Print();
    //run -> Run(10000);
    run -> Run();
    //run -> Draw("save_all");
    run -> Draw();
}
