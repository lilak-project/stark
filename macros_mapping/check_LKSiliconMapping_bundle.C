void check_LKSiliconMapping_bundle(
        TString mappingDir = "/Users/jungwoo/Research/lilak/stark/macros_mapping/mapping_first_bb10_qqq5_test",
        int queryCobo = 0,
        int queryAsad = 1,
        int queryAget = 3,
        int queryChan = 56)
{
    auto mapping = new LKSiliconMapping();
    if (!mapping->Load(mappingDir))
        return;
    mapping -> Print();
    mapping -> PrintDetectors();
    //mapping -> PrintChannels(1000);
    //mapping -> PrintChannelLookup(queryCobo, queryAsad, queryAget, queryChan);
}
