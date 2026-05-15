void make_et_histograms()
{
    //for (auto energyIndex : {4,5,8})
    for (auto energyIndex : {5})
    {
        auto runList = LKParameterContainer("config_short_cut.mac").GetParIntRange(Form("e%d/RunIDRange",energyIndex));
        auto top = new LKDrawingGroup(Form("e%d",energyIndex));
        for (auto id : runList)
            top -> HAdd(Form("data_reco/stark_%04d.hist8.root",id));
        auto file = new TFile(Form("data_x/stark_e%d.hist_et.root",energyIndex),"recreate");
        top -> Write();
        top -> Draw("v");
    }
}
