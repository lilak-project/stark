void create_short_reco_file()
{
    //auto inputName = "../data_reco/stark_0167.reco.root";
    auto inputName = "../data_reco/stark_0168.reco.root";
    auto outputName = "data/stark_7718.reco.root";

    auto file = new TFile(inputName);
    auto tree = (TTree*) file -> Get("event");

    auto file2 = new TFile(outputName,"recreate");
    auto tree2 = tree -> CopyTree("SiChannel.GetDetID()==18");
    file2 -> cd();
    tree2 -> Write();

    //tree2 -> Scan("SiChannel.GetDetID():SiChannel.fStrip:SiChannel.fSide");
}
