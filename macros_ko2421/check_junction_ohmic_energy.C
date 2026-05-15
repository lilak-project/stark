void check_junction_ohmic_energy()
{
    auto run = new LKRun();
    run -> AddDetector(new STARK());
    run -> AddPar("config_check.mac");
    run -> Init();
    run -> SetAutoTermination(false);

    auto runID = run -> GetRunID();
    auto tree = run -> GetInputTree();
    auto cvs = new TCanvas(Form("fast_summary_%d",runID),"",1600,940);
    cvs -> Divide(3,2);
    TString cutName = "SiHit.fDetID==25&&SiHit.fJunctionStrip==0";
    cvs -> cd(1) -> SetLogz(); tree -> Draw("SiHit.fEnergy>>hist1(200,0,20)",cutName,"");
    cvs -> cd(2) -> SetLogz(); tree -> Draw("SiHit.fEnergyOhmic>>hist2(200,0,20)",cutName,"");
    cvs -> cd(3) -> SetLogz(); tree -> Draw("SiHit.fEnergy:SiHit.fEnergyOhmic>>hist3(200,0,20,200,0,20)",cutName,"colz");
}
