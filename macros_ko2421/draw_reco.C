void draw_reco()
{
    auto run = new LKRun();
    run -> AddInputFile("/home/ejungwoo/lilak/ko2421/macros/data_reco/stark_0095.reco8.root");
    run -> Init();
    run -> PrintEvent(0);
    return;

    auto file = new TFile("/home/ejungwoo/lilak/ko2421/macros/data_reco/stark_0094.reco8.root","read");
    auto event = (TTree*) file -> Get("event");
    event -> Print();

    auto cvs = new TCanvas("cvs","",2400,1800);
    cvs -> Divide(3,3);

    cvs -> cd(1); event -> Draw("SiHit.GetEnergyTotal():SiHit.fDetID>>hist1(40,0,40,200,0,40)","SiHit.IsEPairDetector()","colz");
    cvs -> cd(2); event -> Draw("SiHit.GetEnergyTotalOhmic():SiHit.fDetID>>hist2(40,0,40,200,0,40)","SiHit.IsEPairDetector()","colz");

    cvs -> cd(3); event -> Draw("SiHit.fdEOhmic:SiHit.fdEDetID>>hist3(40,0,40,200,0,40)","SiHit.IsEPairDetector()","colz");
    cvs -> cd(4); event -> Draw("SiHit.fEnergyOhmic:SiHit.fdEDetID>>hist4(40,0,40,200,0,40)","SiHit.IsEPairDetector()","colz");
    cvs -> cd(5); event -> Draw("SiHit.fdEOhmic:SiHit.fDetID>>hist5(40,0,40,200,0,40)","SiHit.IsEPairDetector()","colz");
    cvs -> cd(6); event -> Draw("SiHit.fEnergyOhmic:SiHit.fDetID>>hist6(40,0,40,200,0,40)","SiHit.IsEPairDetector()","colz");

    cvs -> cd(7); event -> Draw("SiHit.fdEOhmic:SiHit.fEnergyOhmic>>hist7(200,0,40,200,0,40)","SiHit.IsEPairDetector()","colz");
    cvs -> cd(8); event -> Draw("SiHit.fdEOhmic:SiHit.fdEOhmic+SiHit.fEnergyOhmic>>hist8(200,0,40,200,0,40)","SiHit.IsEPairDetector()","colz");
}
