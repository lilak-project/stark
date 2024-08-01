void draw_sihit()
{
    auto file = new TFile("../data/stark_0199.reco.root");
    //auto file = new TFile("../data/stark_0303.reco.root");
    auto tree = (TTree*) file -> Get("event");

    if (0)
    {
        auto cvsEs = LKPainter::GetPainter() -> CanvasResize("cvsE",900,600,0.6);
        cvsEs -> Divide(3,2);
        cout << 1 << endl; cvsEs -> cd(1); tree -> Draw("SiHit.fEnergy     >>histE1(500,0,10)","SiHit.fDetID>=0 && SiHit.fDetID<12");
        cout << 2 << endl; cvsEs -> cd(4); tree -> Draw("SiHit.fEnergyOhmic>>histE2(500,0,10)","SiHit.fDetID>=0 && SiHit.fDetID<12");
        cout << 3 << endl; cvsEs -> cd(2); tree -> Draw("SiHit.fEnergy     >>histE3(500,0,10)","SiHit.fDetID>=12&& SiHit.fDetID<32");
        cout << 4 << endl; cvsEs -> cd(5); tree -> Draw("SiHit.fEnergyOhmic>>histE4(500,0,10)","SiHit.fDetID>=12&& SiHit.fDetID<32");
        cout << 5 << endl; cvsEs -> cd(3); tree -> Draw("SiHit.fEnergy     >>histE5(500,0,10)","SiHit.fDetID>=32&& SiHit.fDetID<40");
        cout << 6 << endl; cvsEs -> cd(6); tree -> Draw("SiHit.fEnergyOhmic>>histE6(500,0,10)","SiHit.fDetID>=32&& SiHit.fDetID<40");
    }

    if (1)
    {
        auto cvsEZAll = LKPainter::GetPainter() -> CanvasResize("cvsEZAll",600,500,0.8);
        cvsEZAll -> SetMargin(0.1,0.15,0.1,0.1);
        //tree -> Draw("SiHit.fEnergy:SiHit.fZ>>histEZAll(400,-1,1,400,0,10)","SiHit.fDetID>=12&&SiHit.fDetID<32","colz");
        tree -> Draw("SiHit.fEnergy:SiHit.fZ>>histEZAll(400,-1,1,400,0,10)","SiHit.fDetID<12","colz");

        auto cvsJO = LKPainter::GetPainter() -> CanvasResize("cvsJO",1000,500,0.8);
        cvsJO -> Divide(4,2);
        for (auto i=0; i<4; ++i) {
            cvsJO -> cd(i+1);
            tree -> Draw(Form("SiHit.fEnergy:SiHit.fZ>>histJOEX%d(400,-1,1,400,0,10)",i),Form("SiHit.fDetID<12 && SiHit.fOhmicStrip==%d",i),"colz");
        }
        for (auto i=0; i<4; ++i) {
            cvsJO -> cd(4+i+1);
            tree -> Draw(Form("SiHit.fEnergyRight:SiHit.fEnergyLeft>>histJORL%d(400,0,10,400,0,10)",i),Form("SiHit.fDetID<12 && SiHit.fOhmicStrip==%d",i),"colz");
        }
    }

    if (0)
    {
        //auto cvsEZs = LKPainter::GetPainter() -> CanvasResize("cvsEZs",1000,400,1);
        //cvsEZs -> Divide(10,4);
        //for (auto i=0; i<40; ++i)
        auto cvsEZs = LKPainter::GetPainter() -> CanvasResize("cvsEZs",400,300,1);
        cvsEZs -> Divide(4,3);
        for (auto i=0; i<12; ++i)
        {
            cout << i << endl;
            cvsEZs -> cd(i+1);
            tree -> Draw(Form("SiHit.fEnergy:SiHit.fZ>>histEZ%d(200,-1,1,200,0,10)",i),Form("SiHit.fDetID==%d",i),"colz");
        }
    }
}
