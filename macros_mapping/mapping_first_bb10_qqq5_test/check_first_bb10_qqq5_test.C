void check_LKSiliconMapping_bundle(
        TString mappingDir = "/Users/jungwoo/Research/lilak/stark/macros_mapping/mapping_first_bb10_qqq5_test",
        int queryCobo = 0,
        int queryAsad = 0,
        int queryAget = 3,
        int queryChan = 0);

void check_first_bb10_qqq5_test(
        int queryCobo = 0,
        int queryAsad = 0,
        int queryAget = 3,
        int queryChan = 0)
{
    gROOT->Macro("/Users/jungwoo/Research/lilak/stark/macros_mapping/check_LKSiliconMapping_bundle.C");
    check_LKSiliconMapping_bundle(
        "/Users/jungwoo/Research/lilak/stark/macros_mapping/mapping_first_bb10_qqq5_test",
        queryCobo,
        queryAsad,
        queryAget,
        queryChan
    );
}
