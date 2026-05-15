void draw_channels()
{
    auto file = new TFile("data_reco/pre_li7_0111.reco.root");
    auto tree = (TTree*) file -> Get("event");
    for (auto cobo : {0}) {
        for (auto asad : {0}) {
            //TString expr = "SiChannel.fEnergy:SiChannel.fChan+68*SiChannel.fAget+68*4*SiChannel.fAsad";
            //TString expr = "SiChannel.fEnergy:SiChannel.fChan+68*SiChannel.fAget";
            TString expr = "SiChannel.fEnergy:SiChannel.fChan";
            TString cvs_name = Form("cvs_%d_%d",cobo,asad);
            TString cut = Form("SiChannel.fCobo==%d && SiChannel.fAsad==%d",cobo,asad);
            //new TH2D(hist_name, ";channels;energy;", 1024,0,1024,2000,0,5000);
            auto cvs = LKPainter::GetPainter() -> CanvasFull(cvs_name,0.8);
            cvs -> SetTitle(file->GetName());
            cvs -> Divide(2,2);
            for (auto aget : {0,1,2,3})
            {
                cvs -> cd(aget+1);
                TString hist_name = Form("hist_%d_%d_%d",cobo,asad,aget);
                TString cut2 = cut + Form(" && SiChannel.fAget==%d",aget);
                new TH2D(hist_name, cut2+";channels;energy;", 80,0,80,2000,0,2000);
                tree -> Draw(expr + ">>" + hist_name ,cut2,"colz");
            }
        }
    }
}
