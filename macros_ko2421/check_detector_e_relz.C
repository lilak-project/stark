void check_detector_e_relz(int runID = 80)
{
    auto file = new TFile(Form("data/stark_00%d.reco.root",runID));
    auto tree = (TTree*) file -> Get("event");
    auto cvs = new TCanvas("cvs","",1600,900);
    cvs -> Divide(8,4);
    int icvs = 1;
    //for (auto i=0; i<40; ++i) {
    for (auto i=0; i<32; ++i) {
        cout << i << endl;
        cvs -> cd(icvs++);
        auto hist = new TH2D(Form("hist%d",i),Form("run=%d, det=%d",runID,i),200,-1,1,200,0,20000);
        tree -> Draw(Form("SiChannel.GetEnergySum():SiChannel.GetEnergyPos()>>hist%d",i),Form("SiChannel.GetDetID()==%d",i),"colz");
    }
}
