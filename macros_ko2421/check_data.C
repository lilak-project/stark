void GetCuts();
TCutG* cutg;

void check_data()
{
    GetCuts();

    TString cutName;
    //cutName = "SiHit.fdEDetID==37";
    //cutName = "pid2"; //cutName = "SiHit.fDetID!=39";
    //cutName = "cut1";
    //cutName = "pid3&&SiHit.fDetID>=28&&SiHit.fDetID<=31";
    //cutName = "SiHit.fdEDetID>=28&&SiHit.fdEDetID<=31";
    //cutName = "etta_lband&&SiHit.fDetID>=32";
    //cutName = "SiHit.fdEDetID<32";
    //cutName = "etta_lband&&SiHit.fdEDetID>=32";
    //cutName = "etta_lband&&SiHit.fDetID>=4&&SiHit.fDetID<11";
    //cutName = "etta_lband&&SiHit.fDetID<4";
    //cutName = "etta_lband";
    //cutName = "etta_lband";
    //cutName = "SiHit.IsEPairAndOnlydE()==0&&SiHit.fdEDetID==36";//&&SiHit.fJunction";
    //cutName = "SiHit.fEnergy==0&&SiHit.fdE>0&&SiHit.fdEDetID>=28&&SiHit.fdEDetID<=31";
    cutName = "SiHit.fdEDetID>=28&&SiHit.fdEDetID<=31";

    //TString cutName0 = "SiHit.IsEPairAndOnlydE()==0";//&&SiHit.fJunction";
    TString cutName0 = "";//&&SiHit.fJunction";

    auto run = new LKRun();
    run -> AddDetector(new STARK());
    run -> AddPar("config_check.mac");
    run -> Init();
    //run -> Print();
    //cout << run -> GetParameterContainer() -> CheckPar("LKRun/RunIDRange") << endl;
    //return;
    run -> SetAutoTermination(false);
    //run -> Run(10000);
    auto runID = run -> GetRunID();
    gStyle -> SetPalette(kRainBow);

    auto tree = run -> GetInputTree();
    auto cvs = new TCanvas("cvs","",1600,900);
    cvs -> Divide(2,2);

    if (0)
    {
        cvs -> cd(1); tree -> Draw("SiHit.fEnergy:SiHit.fRelativeZ>>hist1(200,-1,1,200,0,30)","SiHit.fDetID==8","colz");
        cvs -> cd(2); tree -> Draw("SiHit.fEnergy:SiHit.fRelativeZ>>hist2(200,-1,1,200,0,30)","SiHit.fDetID==9","colz");
        cvs -> cd(3); tree -> Draw("SiHit.fEnergy:SiHit.fRelativeZ>>hist3(200,-1,1,200,0,30)","SiHit.fDetID==10","colz");
        cvs -> cd(4); tree -> Draw("SiHit.fEnergy:SiHit.fRelativeZ>>hist4(200,-1,1,200,0,30)","SiHit.fDetID==11","colz");
    }
    if (1)
    {
        cvs -> cd(1) -> SetLogz(); tree -> Draw("SiHit.fdE:SiHit.fEnergy>>hist1(200,0,22,200,0,15)",cutName,"colz");
        cvs -> cd(2) -> SetLogz(); tree -> Draw("SiHit.GetEnergyTotal():SiHit.fTheta>>hist2(200,20,80,200,0,30)",cutName,"colz"); 
        cvs -> cd(3) -> SetLogz(); tree -> Draw("SiHit.fdE:SiHit.fEnergy>>hist3(200,0,22,200,0,15)",cutName0,"colz");
        cvs -> cd(4) -> SetLogz(); tree -> Draw("SiHit.GetEnergyTotal():SiHit.fTheta>>hist4(200,20,80,200,0,30)",cutName0,"colz"); 
    }
    return;
	TString qValueString = "SiHit.GetEnergyTotal()*(1+1/40.)-(8.04*40)*(1-40./40.)-2/40.*TMath::Sqrt(40.*1*(8.04*40)*SiHit.GetEnergyTotal())*TMath::Cos(SiHit.GetTheta()*TMath::DegToRad())";
    TString thetaCOMString = "180 - SiHit.fTheta * 2";
    cvs -> cd(2) -> SetLogz(); tree -> Draw(qValueString+":"+thetaCOMString+">>(100,30,130,200,-8,8)",cutName,"colz"); 
    cvs -> cd(3) -> SetLogz(); tree -> Draw("SiHit.GetEnergyTotal():SiHit.fTheta>>hist7(150,20,75,200,0,25)",cutName,"colz");
    cvs -> cd(4) -> SetLogz(); tree -> Draw("SiHit.GetEnergyTotal():SiHit.GetFinalZ(75)>>hist6(100,40,180,200,0,25)",cutName,"colz");
}

void GetCuts()
{
    cutg = new TCutG("etta_lband",5);
    cutg->SetVarX("SiHit.fTheta");
    cutg->SetVarY("SiHit.GetEnergyTotal()");
    cutg->SetTitle("Graph");
    cutg->SetFillStyle(1000);
    cutg->SetPoint(0,27.8532,15.597);
    cutg->SetPoint(1,28.1465,12.9531);
    cutg->SetPoint(2,42.7155,7.75041);
    cutg->SetPoint(3,42.6177,11.162);
    cutg->SetPoint(4,27.8532,15.597);

    cutg = new TCutG("pid3",10);
    cutg->SetVarX("SiHit.GetEnergyTotal()");
    cutg->SetVarY("SiHit.fdE");
    cutg->SetTitle("Graph");
    cutg->SetFillStyle(1000);
    cutg->SetPoint(0,13.6351,12.66);
    cutg->SetPoint(1,9.83481,8.65141);
    cutg->SetPoint(2,11.5199,6.47653);
    cutg->SetPoint(3,17.794,3.40611);
    cutg->SetPoint(4,21.1282,4.89868);
    cutg->SetPoint(5,21.1641,6.90297);
    cutg->SetPoint(6,19.0847,7.8838);
    cutg->SetPoint(7,15.5711,9.63223);
    cutg->SetPoint(8,14.2088,12.0203);
    cutg->SetPoint(9,13.6351,12.66);

    cutg = new TCutG("pid1",11);
    cutg->SetVarX("SiHit.GetEnergyTotal()");
    cutg->SetVarY("SiHit.fdE");
    cutg->SetTitle("Graph");
    cutg->SetFillStyle(1000);
    cutg->SetPoint(0,5.05433,4.68098);
    cutg->SetPoint(1,3.6411,3.01611);
    cutg->SetPoint(2,5.55579,2.07118);
    cutg->SetPoint(3,9.61311,0.946274);
    cutg->SetPoint(4,15.1748,0.676296);
    cutg->SetPoint(5,17.6366,1.17126);
    cutg->SetPoint(6,17.3175,2.38616);
    cutg->SetPoint(7,13.9895,2.52115);
    cutg->SetPoint(8,9.20282,2.92612);
    cutg->SetPoint(9,6.87784,4.05103);
    cutg->SetPoint(10,5.05433,4.68098);

    cutg = new TCutG("pid2",8);
    cutg->SetVarX("SiHit.GetEnergyTotal()");
    cutg->SetVarY("SiHit.fdE");
    cutg->SetTitle("Graph");
    cutg->SetFillStyle(1000);
    cutg->SetPoint(0,12.941,12.7803);
    cutg->SetPoint(1,9.11165,8.82064);
    cutg->SetPoint(2,10.9807,5.94087);
    cutg->SetPoint(3,16.816,3.91604);
    cutg->SetPoint(4,23.563,3.73605);
    cutg->SetPoint(5,23.563,6.25585);
    cutg->SetPoint(6,16.7248,9.22561);
    cutg->SetPoint(7,12.941,12.7803);
}
