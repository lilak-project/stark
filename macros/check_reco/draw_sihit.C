void draw_sihit()
{
    //int run = 199; int det1 = 0;  int det2 = 11;
    int run = 303; int det1 = 12; int det2 = 31;
    //TCut cutActiveDetIDs = Form("SiHit.fDetID>=%d&&SiHit.fDetID<%d",det1,det2);
    TCut cutActiveDetIDs = Form("SiHit.fDetID==13&&SiHit.fJunctionStrip==0",det1,det2);

    auto file = new TFile(Form("../data/stark_%04d.reco.root",run));
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
        //auto cvsEZAll = LKPainter::GetPainter() -> CanvasResize("cvsEZAll",600,500,0.8);
        //cvsEZAll -> SetMargin(0.1,0.15,0.1,0.1);
        //tree -> Draw("SiHit.fEnergy:SiHit.fZ>>histEZAll(400,-1,1,400,0,10)",cutActiveDetIDs,"colz");

        auto cvsJO = LKPainter::GetPainter() -> CanvasResize("cvsJO",1000,500,0.8);
        cvsJO -> Divide(4,2);
        for (auto i=0; i<4; ++i) {
            cout << i << endl;
            cvsJO -> cd(i+1);
            tree -> Draw(Form("SiHit.fEnergy:SiHit.fZ>>histJOEX%d(400,-1,1,400,0,10)",i),cutActiveDetIDs&&Form("SiHit.fOhmicStrip==%d",i),"colz");
        }
        for (auto i=0; i<4; ++i) {
            cout << i << endl;
            cvsJO -> cd(4+i+1);
            tree -> Draw(Form("SiHit.fEnergyRight:SiHit.fEnergyLeft>>histJORL%d(400,0,10,400,0,10)",i),cutActiveDetIDs&&Form("SiHit.fOhmicStrip==%d",i),"colz");
        }
    }

    if (0)
    {
        auto cvsEZs = LKPainter::GetPainter() -> CanvasResize("cvsEZs",400,300,1);
        cvsEZs -> Divide(6,4);
        for (auto i=0; i<(det2-det1+1); ++i)
        {
            cout << i << endl;
            cvsEZs -> cd(i+1);
            tree -> Draw(Form("SiHit.fEnergy:SiHit.fZ>>histEZ%d(200,-1,1,200,0,10)",det1+i),Form("SiHit.fDetID==%d",det1+i),"colz");
        }
    }
}
