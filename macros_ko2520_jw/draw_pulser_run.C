void draw_pulser_run()
{
    gStyle -> SetPalette(kRainbow);
    //auto par = new LKParameterContainer("get_common.mac");
    auto runID = 67;//par -> GetParInt("LKRun/RunID");
    TString fileName = Form("data_reco/ko2520_%04d.pulser.root",runID);
    cout << fileName << endl;
    auto file = new TFile(fileName);
    auto tree = (TTree*) file -> Get("event");
    //tree -> Draw("SiChannel.fEnergy:SiChannel.fChan>>%s",histName.Data()),cut,"colz goff");

    TString nameTop = Form("ko2520_%04d_hit_pattern",runID);
    auto top = new LKDrawingGroup();
    auto group = top -> CreateGroup(nameTop);
    for (auto asad=0; asad<4; ++asad)
    {
        for (auto aget=0; aget<4; ++aget)
        {
            TString cut;
            cut = cut + Form("SiChannel.fAsad==%d",asad);
            cut = cut + Form("&&SiChannel.fAget==%d",aget);
            TString histName = Form("hist_%d_%d",asad,aget);
            cout << asad << " " << aget << " " <<histName << endl;
            //auto hist = new TH2D(histName,"",68,0,68,200,0,2000);
            auto hist = new TH2D(histName,"",68,0,68,100,0,8000);
            //tree -> Draw(Form("RawData.fEnergy:RawData.fChan>>%s",histName.Data()),cut,"colz goff");
            tree -> Draw(Form("SiChannel.fEnergy:SiChannel.fChan>>%s",histName.Data()),cut,"colz goff");
            group -> CreateDrawing() -> Add(hist);
        }
    }
    //top -> Draw("viewer");
    top -> Draw();
    //top -> GetCanvas() -> SaveAs(nameTop+".png");
    //top -> GetCanvas() -> SaveAs(nameTop+".pdf");
}
