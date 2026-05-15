void run_short()
{
    auto run = new LKRun();
    run -> Add(new ELARK());
    run -> AddPar("config_short_cut.mac");
    run -> Add(new EKShortCutTask);
    run -> Init();
    run -> Run();
}
