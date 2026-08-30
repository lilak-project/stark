void draw_reco()
{
    auto file = new TFile("data_reco/ko2520_0031.reco.root");
    auto tree = (TTree*) file -> Get("event");

    if (1)
    {
        auto top2 = new LKDrawingGroup();
        for (auto i2=0; i2<8; ++i2)
        {
            auto group = top2 -> CreateGroup();
            for (auto i=0; i<6; ++i)
            {
                auto padID = 6*i2+i;
                for (auto side=0; side<2; ++side)
                {
                    TString cut;
                    cut = cut + Form("SiChannel.fPadID==%d",padID);
                    cut = cut + Form("&&SiChannel.fSide==%d",side);
                    TString histName = Form("hist_%d_%d",padID,side);
                    cout << padID << " " << side << " " << histName <<endl;
                    auto hist = new TH2D(histName,";channel;energy",32,0,32,200,0,2000);
                    tree -> Draw(Form("SiChannel.fEnergy:2*SiChannel.fStrip+SiChannel.fDirection>>%s",histName.Data()),cut,"colz goff");
                    group -> CreateDrawing() -> Add(hist);
                }
            }
        }
        top2 -> Draw("viewer");
    }
}
