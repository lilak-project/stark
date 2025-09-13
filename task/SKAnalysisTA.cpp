#include "TStyle.h"

#include "LKRun.h"
#include "LKLogger.h"
#include "GETChannel.h"
#include "SKSiHit.h"
#include "SKAnalysisTA.h"

ClassImp(SKAnalysisTA)

SKAnalysisTA::SKAnalysisTA()
    :LKTask("SKAnalysisTA","SKAnalysisTA")
{
}

bool SKAnalysisTA::Init()
{
    fSiHitArray = fRun -> GetBranchA("SiHit");

    fStarkPlane = (SKSiArrayPlane*) fRun -> FindDetectorPlane("SKSiArrayPlane");
    if (fStarkPlane==nullptr) {
        lk_error << "SKSiArrayPlane do not exist!!!" << endl;
        lk_error << fName << " must run with SKSiArrayPlane" << endl;
        return false;
    }

    fPar -> UpdatePar(fBeamEnergy,"ko2421/beam_energy");
    fPar -> UpdatePar(fBeamRate,"ko2421/beam_rate");
    fPar -> UpdatePar(fRunTime,"ko2421/run_time");
    fPar -> UpdatePar(fNTarget,"ko2421/n_target");

    int ne = 400;
    double e1 = 0;
    double e2 = 10;
    double e1com = -10;
    double e2com = 10;
    fPar -> UpdateBinning("analysis/binning_cal_energy", ne, e1, e2);

    //fHistAvstotE[0] = new TH2D(Form("fHistAlabvsE"),Form("Lab Polar Angle vs Energy;theta (deg);energy (MeV)"),1800,0,180,ne,e1,2*e2);
    fHistAvstotE[0] = new TH2D(Form("fHistAlabvsE"),Form("Lab Polar Angle vs Energy;theta (deg);energy (MeV)"),1800,0,180,ne,e1,e2);
    //fHistAvstotE[1] = new TH2D(Form("fHistAcomvsE"),Form("CoM Polar Angle vs Energy;theta (deg);energy (MeV)"),1800,0,180,ne,e1com,e2com);
    fHistAvstotE[1] = new TH2D(Form("fHistAcomvsE"),Form("CoM Polar Angle vs Energy;theta (deg);energy (MeV)"),1800,0,180,ne,-8,8);
    auto array5 = new TObjArray();
    array5 -> Add(fHistAvstotE[0]);
    array5 -> Add(fHistAvstotE[1]);
    fStarkPlane -> AddUserDrawingArray("AvsE_all", array5);


    auto array7 = new TObjArray();
    fHistAcomgonel = new TH1D(Form("fHistAcomgonel"),Form("CoM Polar Angle;theta (deg)"),180,0,180);
    fHistDiffSiggonel = new TH1D(Form("fHistDiffSiggonel"),Form("Diffrential Crosssection vs CoM Polar Angle;theta (deg);mbarn/sr"),180,0,180);
    fHistDiffSigfresco = new TH1D(Form("fHistDiffSigfresco"),Form("Diffrential Crosssection vs CoM Polar Angle;theta (deg);mbarn/sr"),180,0,180);
    DrawDiffSig();
    array7 -> Add(fHistAcomgonel);
    array7 -> Add(fHistDiffSiggonel);
    fStarkPlane -> AddUserDrawingArray("AcomDiff_all", array7);

    for(int det=0; det<40; det++){
	    auto array6 = new TObjArray();
	    fHistAvsE[det][0] = new TH2D(Form("fHistAlabvsE_%d",det),Form("Lab Polar Angle vs Energy_%d;theta (deg);energy (MeV)",det),60,20,80,ne,e1,e2);
	    fHistAvsE[det][1] = new TH2D(Form("fHistAcomvsE_%d",det),Form("CoM Polar Angle vs Energy_%d;theta (deg);energy (MeV)",det),100,30,130,ne,e1com,e2com);
	    array6 -> Add(fHistAvsE[det][0]);
	    array6 -> Add(fHistAvsE[det][1]);
	    fStarkPlane -> AddUserDrawingArray("AvsE_each", det, array6);
    }

    for(int pairID=0; pairID<12; pairID++){
	    auto array8 = new TObjArray();
	    fHistdEEangle[pairID][0] = new TH2D(Form("fHistdEE_labangle_e_%d",pairID),Form("dE_Eangle_%d;Theta_Lab;Etotal",pairID),1800,0,180,ne,e1,2*e2);
	    fHistdEEangle[pairID][1] = new TH2D(Form("fHistdEE_comangle_e_%d",pairID),Form("dE_Eangle_%d;Theta_CoM;Etotal",pairID),1800,0,180,ne,e1,2*e2);
	    array8 -> Add(fHistdEEangle[pairID][0]);
	    array8 -> Add(fHistdEEangle[pairID][1]);
	    fStarkPlane -> AddUserDrawingArray("dE_E_angle", pairID, array8);
    }

   cutg = new TCutG("elatics",10);
   cutg->SetVarX("Lab Polar Angle vs Energy");
   cutg->SetVarY("");
   cutg->SetTitle("Graph");
   cutg->SetFillStyle(1000);
   cutg->SetPoint(0,28.9392,21.0064);
   cutg->SetPoint(1,39.7499,15.5253);
   cutg->SetPoint(2,48.2376,12.2772);
   cutg->SetPoint(3,51.9007,11.1607);
   cutg->SetPoint(4,53.777,9.94265);
   cutg->SetPoint(5,71.9139,2.63449);
   cutg->SetPoint(6,69.5909,1.61947);
   cutg->SetPoint(7,50.2032,8.01411);
   cutg->SetPoint(8,47.5228,10.6532);
   cutg->SetPoint(9,41.5368,11.9727);
   cutg->SetPoint(10,29.654,16.1343);
   cutg->SetPoint(11,27.5991,18.7733);
   cutg->SetPoint(12,28.9392,21.0064);

    return true;
}

void SKAnalysisTA::Exec(Option_t*)
{
	if (!fStarkPlane -> GetAccumulateEvents()) {
		//fHistHP[0] -> Reset();
		//fHistHP[1] -> Reset();
		//fHistET[0] -> Reset();
		//fHistET[1] -> Reset();
	}

	auto numHits = fSiHitArray -> GetEntries();
	for (auto iHit=0; iHit<numHits; ++iHit)
	{
		auto siHit = (SKSiHit*) fSiHitArray -> At(iHit);
        if (siHit->InGate()==false)
            continue;
		int detID = siHit -> GetDetID();
		int stripJ = siHit -> GetJunctionStrip();
		int stripO = siHit -> GetOhmicStrip();
		double energySum = siHit -> GetEnergy();
		double energyOhmic = siHit -> GetEnergyOhmic();
		double z_relative = siHit -> GetRelativeZ(); 
		double z_install = fStarkPlane -> GetSiDetector(detID) -> GetZ();
		double dE_energy = siHit -> GetdE();
		double E_energy = siHit -> GetE();
		bool check_dEE = siHit -> IsEPairAndBothdEE();
		//bool check_singE = siHit -> IsEDetector();
		bool check_singE = siHit -> IsEPairAndOnlyE();
		bool check_singdE = siHit -> IsEPairAndOnlydE();
		bool check_16E = siHit -> IsNotEPairDetector();
		double ring_radius = fStarkPlane -> GetSiDetector(detID) -> GetRadius();
		double ring_id = fStarkPlane -> GetSiDetector(detID) -> GetRow();

		double z_real = 75*z_relative/2 + z_install;

        int ring_idmax = 0;
        if(check_16E == true){
            ring_idmax = 16;
        } else
            ring_idmax = 12;

        //else if(check_singE == true || check_singdE == true){
		//       	ring_idmax = 12;
		//}else{
		//       	cout << "no check flags." << endl;
		//       	ring_idmax = 360;
		//}

  		double theta_lab = TMath::RadToDeg()*TMath::ATan(ring_radius/z_real);
  		//double phi_lab = ring_id * (360/ring_idmax);
  		//double phi_lab = siHit -> GetPhi();
		double Qvalue = 0;
  		double theta_com = 180 - theta_lab * 2;
  		double beamenergy = fBeamEnergy; // [MeV/u]
  		double beamrate = fBeamRate; // [pps]
  		double ntarget = fNTarget;
  		double runtime = fRunTime; // [sec] for run20-run25
		double diffsig = 0;
        auto detector = fStarkPlane -> GetSiDetector(detID);
        auto ring = detector -> GetLayer();
        //double tot_energy = dE_energy + E_energy;
        //double tot_energy = siHit -> GetEnergyOhmic();
        double tot_energy = siHit -> GetEnergyOhmic() + dE_energy;

        if(check_16E){
            fHistAvstotE[0] -> Fill(theta_lab,tot_energy);
            Qvalue = GetQvalue(beamenergy,theta_lab,tot_energy);
            fHistAvstotE[1] -> Fill(theta_com,Qvalue);
            fHistAcomgonel->Fill(theta_com);
            diffsig = GetDiffSig(ring_idmax-2,beamrate,ntarget,runtime,theta_com);
            fHistDiffSiggonel->Fill(theta_com,diffsig);
            if (ring ==0){
                fHistAvsE[detID][0] -> Fill(theta_lab, tot_energy);
                fHistAvsE[detID][1] -> Fill(theta_com,Qvalue);
            }else{
                fHistAvsE[detID][0] -> Fill(theta_lab, tot_energy);
                fHistAvsE[detID][1] -> Fill(theta_com,Qvalue);
            }
        }
        //if(check_singE==true){
        //    fHistAvstotE[0] -> Fill(theta_lab,E_energy);
        //    Qvalue = GetQvalue(beamenergy,theta_lab,E_energy);
        //    fHistAvstotE[1] -> Fill(theta_com,Qvalue);
        //    fHistAcomgonel->Fill(theta_com);
        //    diffsig = GetDiffSig(ring_idmax-2,beamrate,ntarget,runtime,theta_com);
        //    fHistDiffSiggonel->Fill(theta_com,diffsig);
        //}
        //if(check_singdE==true&&siHit.GetDetID()>=28&&siHit.GetDetID()<=31){
        //    fHistAvstotE[0] -> Fill(theta_lab,dE_energy);
        //    Qvalue = GetQvalue(beamenergy,theta_lab,dE_energy);
        //    fHistAvstotE[1] -> Fill(theta_com,Qvalue);
        //    fHistAcomgonel->Fill(theta_com);
        //    diffsig = GetDiffSig(ring_idmax,beamrate,ntarget,runtime,theta_com);
        //    fHistDiffSiggonel->Fill(theta_com,diffsig);
        //}
        //if(check_dEE==true && detID>3 && detID<12){
        if(check_dEE) {
            fHistAvstotE[0] -> Fill(theta_lab,tot_energy);
            Qvalue = GetQvalue(beamenergy,theta_lab,tot_energy);
            fHistAvstotE[1] -> Fill(theta_com,Qvalue);
            fHistAcomgonel->Fill(theta_com);
            diffsig = GetDiffSig(ring_idmax,beamrate,ntarget,runtime,theta_com);
            fHistDiffSiggonel->Fill(theta_com,diffsig);

            int pairID = fStarkPlane -> GetSiDetector(detID) -> GetRow();
            //cout << "This is dE_E event!" << "\tE energy: " << E_energy << "\tdE energy: " << dE_energy << "\tdetID: " << detID << "\tpair ID " << pairID << endl;
            fHistdEEangle[pairID][0] -> Fill(theta_lab, tot_energy);
            fHistdEEangle[pairID][1] -> Fill(theta_com, Qvalue);
            if (ring ==0){
                fHistAvsE[detID][0] -> Fill(theta_lab, tot_energy);
                fHistAvsE[detID][1] -> Fill(theta_com,Qvalue);
            }else{
                fHistAvsE[detID][0] -> Fill(theta_lab, tot_energy);
                fHistAvsE[detID][1] -> Fill(theta_com,Qvalue);
            }
        }
    }
}

double SKAnalysisTA::GetQvalue(double beamenergy, double theta, double energy)
{
	double massp = 1;
	double massar = 40;
	double E_beam = beamenergy * massar;
	double value = energy * (1 + massp/massar) - E_beam * (1 - massar/massar) - 2 / massar * TMath::Sqrt(massar * massp * E_beam * energy) * TMath::Cos(theta*TMath::DegToRad());
	return value;
}

double SKAnalysisTA::GetDiffSig(double ring_idmax, double beamrate, double ntarget, double runtime, double theta)
{
	double dphi = 0.34; // radian per detector
	double dcostheta = TMath::Cos(int(theta)*TMath::DegToRad()) - TMath::Cos((int(theta)+1)*TMath::DegToRad());
	double domega = dcostheta * dphi * ring_idmax;
	double beamtotal = beamrate * runtime;
	//double ntarget = 1.96E19; //CH2 no1
	//double ntarget = 2.56E19; //CH2 no3
	double cm2tombarn = 1E27;
	double value = cm2tombarn / ntarget / domega / beamtotal;
	return value;
}
void SKAnalysisTA::DrawDiffSig()
{
       //Ebeam = 4.41 MeV/u
       //fHistDiffSigfresco->Fill(30.00,4674.588750);
       //fHistDiffSigfresco->Fill(31.00,4012.411589);
       //fHistDiffSigfresco->Fill(32.00,3452.677122);
       //fHistDiffSigfresco->Fill(33.00,2977.844732);
       //fHistDiffSigfresco->Fill(34.00,2573.727991);
       //fHistDiffSigfresco->Fill(35.00,2228.784069);
       //fHistDiffSigfresco->Fill(36.00,1933.568298);
       //fHistDiffSigfresco->Fill(37.00,1680.312143);
       //fHistDiffSigfresco->Fill(38.00,1462.594219);
       //fHistDiffSigfresco->Fill(39.00,1275.082029);
       //fHistDiffSigfresco->Fill(40.00,1113.327919);
       //fHistDiffSigfresco->Fill(41.00,973.606857);
       //fHistDiffSigfresco->Fill(42.00,852.786742);
       //fHistDiffSigfresco->Fill(43.00,748.224152);
       //fHistDiffSigfresco->Fill(44.00,657.680123);
       //fHistDiffSigfresco->Fill(45.00,579.251781);
       //fHistDiffSigfresco->Fill(46.00,511.316596);
       //fHistDiffSigfresco->Fill(47.00,452.486736);
       //fHistDiffSigfresco->Fill(48.00,401.571541);
       //fHistDiffSigfresco->Fill(49.00,357.546553);
       //fHistDiffSigfresco->Fill(50.00,319.527862);
       //fHistDiffSigfresco->Fill(51.00,286.750791);
       //fHistDiffSigfresco->Fill(52.00,258.552119);
       //fHistDiffSigfresco->Fill(53.00,234.355210);
       //fHistDiffSigfresco->Fill(54.00,213.657536);
       //fHistDiffSigfresco->Fill(55.00,196.020180);
       //fHistDiffSigfresco->Fill(56.00,181.058976);
       //fHistDiffSigfresco->Fill(57.00,168.437010);
       //fHistDiffSigfresco->Fill(58.00,157.858259);
       //fHistDiffSigfresco->Fill(59.00,149.062180);
       //fHistDiffSigfresco->Fill(60.00,141.819092);
       //fHistDiffSigfresco->Fill(61.00,135.926227);
       //fHistDiffSigfresco->Fill(62.00,131.204349);
       //fHistDiffSigfresco->Fill(63.00,127.494840);
       //fHistDiffSigfresco->Fill(64.00,124.657197);
       //fHistDiffSigfresco->Fill(65.00,122.566862);
       //fHistDiffSigfresco->Fill(66.00,121.113349);
       //fHistDiffSigfresco->Fill(67.00,120.198614);
       //fHistDiffSigfresco->Fill(68.00,119.735636);
       //fHistDiffSigfresco->Fill(69.00,119.647182);
       //fHistDiffSigfresco->Fill(70.00,119.864720);
       //fHistDiffSigfresco->Fill(71.00,120.327474);
       //fHistDiffSigfresco->Fill(72.00,120.981585);
       //fHistDiffSigfresco->Fill(73.00,121.779376);
       //fHistDiffSigfresco->Fill(74.00,122.678699);
       //fHistDiffSigfresco->Fill(75.00,123.642366);
       //fHistDiffSigfresco->Fill(76.00,124.637626);
       //fHistDiffSigfresco->Fill(77.00,125.635720);
       //fHistDiffSigfresco->Fill(78.00,126.611470);
       //fHistDiffSigfresco->Fill(79.00,127.542919);
       //fHistDiffSigfresco->Fill(80.00,128.411005);
       //fHistDiffSigfresco->Fill(81.00,129.199272);
       //fHistDiffSigfresco->Fill(82.00,129.893610);
       //fHistDiffSigfresco->Fill(83.00,130.482021);
       //fHistDiffSigfresco->Fill(84.00,130.954410);
       //fHistDiffSigfresco->Fill(85.00,131.302393);
       //fHistDiffSigfresco->Fill(86.00,131.519131);
       //fHistDiffSigfresco->Fill(87.00,131.599176);
       //fHistDiffSigfresco->Fill(88.00,131.538335);
       //fHistDiffSigfresco->Fill(89.00,131.333543);
       //fHistDiffSigfresco->Fill(90.00,130.982759);
       //fHistDiffSigfresco->Fill(91.00,130.484858);
       //fHistDiffSigfresco->Fill(92.00,129.839550);
       //fHistDiffSigfresco->Fill(93.00,129.047298);
       //fHistDiffSigfresco->Fill(94.00,128.109243);
       //fHistDiffSigfresco->Fill(95.00,127.027145);
       //fHistDiffSigfresco->Fill(96.00,125.803327);
       //fHistDiffSigfresco->Fill(97.00,124.440620);
       //fHistDiffSigfresco->Fill(98.00,122.942320);
       //fHistDiffSigfresco->Fill(99.00,121.312152);
       //fHistDiffSigfresco->Fill(100.00,119.55422);
       //fHistDiffSigfresco->Fill(101.00,117.67302);
       //fHistDiffSigfresco->Fill(102.00,115.67333);
       //fHistDiffSigfresco->Fill(103.00,113.56026);
       //fHistDiffSigfresco->Fill(104.00,111.33921);
       //fHistDiffSigfresco->Fill(105.00,109.01581);
       //fHistDiffSigfresco->Fill(106.00,106.59597);
       //fHistDiffSigfresco->Fill(107.00,104.08579);
       //fHistDiffSigfresco->Fill(108.00,101.49159);
       //fHistDiffSigfresco->Fill(109.00,98.819890);
       //fHistDiffSigfresco->Fill(110.00,96.077354);
       //fHistDiffSigfresco->Fill(111.00,93.270828);
       //fHistDiffSigfresco->Fill(112.00,90.407293);
       //fHistDiffSigfresco->Fill(113.00,87.493851);
       //fHistDiffSigfresco->Fill(114.00,84.537716);
       //fHistDiffSigfresco->Fill(115.00,81.546188);
       //fHistDiffSigfresco->Fill(116.00,78.526640);
       //fHistDiffSigfresco->Fill(117.00,75.486496);
       //fHistDiffSigfresco->Fill(118.00,72.433213);
       //fHistDiffSigfresco->Fill(119.00,69.374258);
       //fHistDiffSigfresco->Fill(120.00,66.317085);
       //fHistDiffSigfresco->Fill(121.00,63.269115);
       //fHistDiffSigfresco->Fill(122.00,60.237710);
       //fHistDiffSigfresco->Fill(123.00,57.230147);
       //fHistDiffSigfresco->Fill(124.00,54.253594);
       //fHistDiffSigfresco->Fill(125.00,51.315084);
       //fHistDiffSigfresco->Fill(126.00,48.421488);
       //fHistDiffSigfresco->Fill(127.00,45.579485);
       //fHistDiffSigfresco->Fill(128.00,42.795540);
       //fHistDiffSigfresco->Fill(129.00,40.075872);
       //fHistDiffSigfresco->Fill(130.00,37.426429);
       //fHistDiffSigfresco->Fill(131.00,34.852860);
       //fHistDiffSigfresco->Fill(132.00,32.360492);
       //fHistDiffSigfresco->Fill(133.00,29.954299);
       //fHistDiffSigfresco->Fill(134.00,27.638884);
       //fHistDiffSigfresco->Fill(135.00,25.418451);
       //fHistDiffSigfresco->Fill(136.00,23.296787);
       //fHistDiffSigfresco->Fill(137.00,21.277238);
       //fHistDiffSigfresco->Fill(138.00,19.362697);
       //fHistDiffSigfresco->Fill(139.00,17.555583);
       //fHistDiffSigfresco->Fill(140.00,15.857830);
       //fHistDiffSigfresco->Fill(141.00,14.270874);
       //fHistDiffSigfresco->Fill(142.00,12.795650);
       //fHistDiffSigfresco->Fill(143.00,11.432579);
       //fHistDiffSigfresco->Fill(144.00,10.181570);
       //fHistDiffSigfresco->Fill(145.00,9.042020);
       //fHistDiffSigfresco->Fill(146.00,8.012818);
       //fHistDiffSigfresco->Fill(147.00,7.092349);
       //fHistDiffSigfresco->Fill(148.00,6.278508);
       //fHistDiffSigfresco->Fill(149.00,5.568707);
       //fHistDiffSigfresco->Fill(150.00,4.959900);
       //fHistDiffSigfresco->Fill(151.00,4.448593);
       //fHistDiffSigfresco->Fill(152.00,4.030875);
       //fHistDiffSigfresco->Fill(153.00,3.702437);
       //fHistDiffSigfresco->Fill(154.00,3.458603);
       //fHistDiffSigfresco->Fill(155.00,3.294362);
       //fHistDiffSigfresco->Fill(156.00,3.204398);
       //fHistDiffSigfresco->Fill(157.00,3.183129);
       //fHistDiffSigfresco->Fill(158.00,3.224744);
       //fHistDiffSigfresco->Fill(159.00,3.323243);
       //fHistDiffSigfresco->Fill(160.00,3.472474);
       //fHistDiffSigfresco->Fill(161.00,3.666182);
       //fHistDiffSigfresco->Fill(162.00,3.898048);
       //fHistDiffSigfresco->Fill(163.00,4.161734);
       //fHistDiffSigfresco->Fill(164.00,4.450928);
       //fHistDiffSigfresco->Fill(165.00,4.759388);
       //fHistDiffSigfresco->Fill(166.00,5.080986);
       //fHistDiffSigfresco->Fill(167.00,5.409752);
       //fHistDiffSigfresco->Fill(168.00,5.739916);
       //fHistDiffSigfresco->Fill(169.00,6.065950);
       //fHistDiffSigfresco->Fill(170.00,6.382606);
       //fHistDiffSigfresco->Fill(171.00,6.684955);
       //fHistDiffSigfresco->Fill(172.00,6.968423);
       //fHistDiffSigfresco->Fill(173.00,7.228824);
       //fHistDiffSigfresco->Fill(174.00,7.462389);
       //fHistDiffSigfresco->Fill(175.00,7.665794);
       //fHistDiffSigfresco->Fill(176.00,7.836184);
       //fHistDiffSigfresco->Fill(177.00,7.971195);
       //fHistDiffSigfresco->Fill(178.00,8.068969);
       //fHistDiffSigfresco->Fill(179.00,8.128171);
       //fHistDiffSigfresco->Fill(179.99,8.147993);
       //Ebeam = 5.838 MeV/u
       // GOMP_KD
       fHistDiffSigfresco->Fill(30,2.297E+03);
       fHistDiffSigfresco->Fill(31,2.039E+03);
       fHistDiffSigfresco->Fill(32,1.821E+03);
       fHistDiffSigfresco->Fill(33,1.636E+03);
       fHistDiffSigfresco->Fill(34,1.478E+03);
       fHistDiffSigfresco->Fill(35,1.341E+03);
       fHistDiffSigfresco->Fill(36,1.223E+03);
       fHistDiffSigfresco->Fill(37,1.119E+03);
       fHistDiffSigfresco->Fill(38,1.027E+03);
       fHistDiffSigfresco->Fill(39,9.460E+02);
       fHistDiffSigfresco->Fill(40,8.735E+02);
       fHistDiffSigfresco->Fill(41,8.083E+02);
       fHistDiffSigfresco->Fill(42,7.494E+02);
       fHistDiffSigfresco->Fill(43,6.959E+02);
       fHistDiffSigfresco->Fill(44,6.471E+02);
       fHistDiffSigfresco->Fill(45,6.025E+02);
       fHistDiffSigfresco->Fill(46,5.614E+02);
       fHistDiffSigfresco->Fill(47,5.235E+02);
       fHistDiffSigfresco->Fill(48,4.885E+02);
       fHistDiffSigfresco->Fill(49,4.560E+02);
       fHistDiffSigfresco->Fill(50,4.257E+02);
       fHistDiffSigfresco->Fill(51,3.976E+02);
       fHistDiffSigfresco->Fill(52,3.713E+02);
       fHistDiffSigfresco->Fill(53,3.467E+02);
       fHistDiffSigfresco->Fill(54,3.237E+02);
       fHistDiffSigfresco->Fill(55,3.022E+02);
       fHistDiffSigfresco->Fill(56,2.820E+02);
       fHistDiffSigfresco->Fill(57,2.631E+02);
       fHistDiffSigfresco->Fill(58,2.454E+02);
       fHistDiffSigfresco->Fill(59,2.287E+02);
       fHistDiffSigfresco->Fill(60,2.130E+02);
       fHistDiffSigfresco->Fill(61,1.984E+02);
       fHistDiffSigfresco->Fill(62,1.846E+02);
       fHistDiffSigfresco->Fill(63,1.716E+02);
       fHistDiffSigfresco->Fill(64,1.595E+02);
       fHistDiffSigfresco->Fill(65,1.481E+02);
       fHistDiffSigfresco->Fill(66,1.374E+02);
       fHistDiffSigfresco->Fill(67,1.273E+02);
       fHistDiffSigfresco->Fill(68,1.180E+02);
       fHistDiffSigfresco->Fill(69,1.092E+02);
       fHistDiffSigfresco->Fill(70,1.009E+02);
       fHistDiffSigfresco->Fill(71,9.325E+01);
       fHistDiffSigfresco->Fill(72,8.607E+01);
       fHistDiffSigfresco->Fill(73,7.937E+01);
       fHistDiffSigfresco->Fill(74,7.312E+01);
       fHistDiffSigfresco->Fill(75,6.730E+01);
       fHistDiffSigfresco->Fill(76,6.189E+01);
       fHistDiffSigfresco->Fill(77,5.687E+01);
       fHistDiffSigfresco->Fill(78,5.221E+01);
       fHistDiffSigfresco->Fill(79,4.789E+01);
       fHistDiffSigfresco->Fill(80,4.390E+01);
       fHistDiffSigfresco->Fill(81,4.021E+01);
       fHistDiffSigfresco->Fill(82,3.682E+01);
       fHistDiffSigfresco->Fill(83,3.370E+01);
       fHistDiffSigfresco->Fill(84,3.084E+01);
       fHistDiffSigfresco->Fill(85,2.822E+01);
       fHistDiffSigfresco->Fill(86,2.583E+01);
       fHistDiffSigfresco->Fill(87,2.365E+01);
       fHistDiffSigfresco->Fill(88,2.168E+01);
       fHistDiffSigfresco->Fill(89,1.990E+01);
       fHistDiffSigfresco->Fill(90,1.830E+01);
       fHistDiffSigfresco->Fill(91,1.687E+01);
       fHistDiffSigfresco->Fill(92,1.559E+01);
       fHistDiffSigfresco->Fill(93,1.446E+01);
       fHistDiffSigfresco->Fill(94,1.347E+01);
       fHistDiffSigfresco->Fill(95,1.260E+01);
       fHistDiffSigfresco->Fill(96,1.186E+01);
       fHistDiffSigfresco->Fill(97,1.122E+01);
       fHistDiffSigfresco->Fill(98,1.069E+01);
       fHistDiffSigfresco->Fill(99,1.026E+01);
       fHistDiffSigfresco->Fill(100,9.908E+00);
       fHistDiffSigfresco->Fill(101,9.643E+00);
       fHistDiffSigfresco->Fill(102,9.453E+00);
       fHistDiffSigfresco->Fill(103,9.333E+00);
       fHistDiffSigfresco->Fill(104,9.276E+00);
       fHistDiffSigfresco->Fill(105,9.277E+00);
       fHistDiffSigfresco->Fill(106,9.331E+00);
       fHistDiffSigfresco->Fill(107,9.434E+00);
       fHistDiffSigfresco->Fill(108,9.580E+00);
       fHistDiffSigfresco->Fill(109,9.766E+00);
       fHistDiffSigfresco->Fill(110,9.987E+00);
       fHistDiffSigfresco->Fill(111,1.024E+01);
       fHistDiffSigfresco->Fill(112,1.052E+01);
       fHistDiffSigfresco->Fill(113,1.083E+01);
       fHistDiffSigfresco->Fill(114,1.115E+01);
       fHistDiffSigfresco->Fill(115,1.150E+01);
       fHistDiffSigfresco->Fill(116,1.186E+01);
       fHistDiffSigfresco->Fill(117,1.224E+01);
       fHistDiffSigfresco->Fill(118,1.263E+01);
       fHistDiffSigfresco->Fill(119,1.303E+01);
       fHistDiffSigfresco->Fill(120,1.343E+01);
       fHistDiffSigfresco->Fill(121,1.385E+01);
       fHistDiffSigfresco->Fill(122,1.426E+01);
       fHistDiffSigfresco->Fill(123,1.468E+01);
       fHistDiffSigfresco->Fill(124,1.510E+01);
       fHistDiffSigfresco->Fill(125,1.552E+01);
       fHistDiffSigfresco->Fill(126,1.594E+01);
       fHistDiffSigfresco->Fill(127,1.636E+01);
       fHistDiffSigfresco->Fill(128,1.678E+01);
       fHistDiffSigfresco->Fill(129,1.719E+01);
       fHistDiffSigfresco->Fill(130,1.760E+01);
       // Double foling
       //fHistDiffSigfresco->Fill( 30.00, 2294.219812);
       //fHistDiffSigfresco->Fill( 31.00, 2033.510197);
       //fHistDiffSigfresco->Fill( 32.00, 1812.837846);
       //fHistDiffSigfresco->Fill( 33.00, 1624.432872);
       //fHistDiffSigfresco->Fill( 34.00, 1462.218964);
       //fHistDiffSigfresco->Fill( 35.00, 1321.418419);
       //fHistDiffSigfresco->Fill( 36.00, 1198.254963);
       //fHistDiffSigfresco->Fill( 37.00, 1089.728672);
       //fHistDiffSigfresco->Fill( 38.00,  993.444427);
       //fHistDiffSigfresco->Fill( 39.00,  907.480430);
       //fHistDiffSigfresco->Fill( 40.00,  830.286881);
       //fHistDiffSigfresco->Fill( 41.00,  760.607519);
       //fHistDiffSigfresco->Fill( 42.00,  697.418582);
       //fHistDiffSigfresco->Fill( 43.00,  639.881121);
       //fHistDiffSigfresco->Fill( 44.00,  587.303598);
       //fHistDiffSigfresco->Fill( 45.00,  539.112444);
       //fHistDiffSigfresco->Fill( 46.00,  494.828806);
       //fHistDiffSigfresco->Fill( 47.00,  454.050120);
       //fHistDiffSigfresco->Fill( 48.00,  416.435471);
       //fHistDiffSigfresco->Fill( 49.00,  381.693932);
       //fHistDiffSigfresco->Fill( 50.00,  349.575246);
       //fHistDiffSigfresco->Fill( 51.00,  319.862371);
       //fHistDiffSigfresco->Fill( 52.00,  292.365498);
       //fHistDiffSigfresco->Fill( 53.00,  266.917249);
       //fHistDiffSigfresco->Fill( 54.00,  243.368804);
       //fHistDiffSigfresco->Fill( 55.00,  221.586782);
       //fHistDiffSigfresco->Fill( 56.00,  201.450722);
       //fHistDiffSigfresco->Fill( 57.00,  182.851046);
       //fHistDiffSigfresco->Fill( 58.00,  165.687406);
       //fHistDiffSigfresco->Fill( 59.00,  149.867346);
       //fHistDiffSigfresco->Fill( 60.00,  135.305215);
       //fHistDiffSigfresco->Fill( 61.00,  121.921278);
       //fHistDiffSigfresco->Fill( 62.00,  109.640996);
       //fHistDiffSigfresco->Fill( 63.00,   98.394435);
       //fHistDiffSigfresco->Fill( 64.00,   88.115782);
       //fHistDiffSigfresco->Fill( 65.00,   78.742945);
       //fHistDiffSigfresco->Fill( 66.00,   70.217223);
       //fHistDiffSigfresco->Fill( 67.00,   62.483037);
       //fHistDiffSigfresco->Fill( 68.00,   55.487695);
       //fHistDiffSigfresco->Fill( 69.00,   49.181204);
       //fHistDiffSigfresco->Fill( 70.00,   43.516107);
       //fHistDiffSigfresco->Fill( 71.00,   38.447338);
       //fHistDiffSigfresco->Fill( 72.00,   33.932105);
       //fHistDiffSigfresco->Fill( 73.00,   29.929780);
       //fHistDiffSigfresco->Fill( 74.00,   26.401806);
       //fHistDiffSigfresco->Fill( 75.00,   23.311609);
       //fHistDiffSigfresco->Fill( 76.00,   20.624522);
       //fHistDiffSigfresco->Fill( 77.00,   18.307712);
       //fHistDiffSigfresco->Fill( 78.00,   16.330111);
       //fHistDiffSigfresco->Fill( 79.00,   14.662353);
       //fHistDiffSigfresco->Fill( 80.00,   13.276713);
       //fHistDiffSigfresco->Fill( 81.00,   12.147046);
       //fHistDiffSigfresco->Fill( 82.00,   11.248736);
       //fHistDiffSigfresco->Fill( 83.00,   10.558632);
       //fHistDiffSigfresco->Fill( 84.00,   10.055003);
       //fHistDiffSigfresco->Fill( 85.00,    9.717479);
       //fHistDiffSigfresco->Fill( 86.00,    9.527004);
       //fHistDiffSigfresco->Fill( 87.00,    9.465779);
       //fHistDiffSigfresco->Fill( 88.00,    9.517222);
       //fHistDiffSigfresco->Fill( 89.00,    9.665908);
       //fHistDiffSigfresco->Fill( 90.00,    9.897530);
       //fHistDiffSigfresco->Fill( 91.00,   10.198848);
       //fHistDiffSigfresco->Fill( 92.00,   10.557645);
       //fHistDiffSigfresco->Fill( 93.00,   10.962678);
       //fHistDiffSigfresco->Fill( 94.00,   11.403640);
       //fHistDiffSigfresco->Fill( 95.00,   11.871110);
       //fHistDiffSigfresco->Fill( 96.00,   12.356519);
       //fHistDiffSigfresco->Fill( 97.00,   12.852100);
       //fHistDiffSigfresco->Fill( 98.00,   13.350857);
       //fHistDiffSigfresco->Fill( 99.00,   13.846520);
       //fHistDiffSigfresco->Fill(100.00,   14.333510);
       //fHistDiffSigfresco->Fill(101.00,   14.806899);
       //fHistDiffSigfresco->Fill(102.00,   15.262381);
       //fHistDiffSigfresco->Fill(103.00,   15.696230);
       //fHistDiffSigfresco->Fill(104.00,   16.105272);
       //fHistDiffSigfresco->Fill(105.00,   16.486848);
       //fHistDiffSigfresco->Fill(106.00,   16.838785);
       //fHistDiffSigfresco->Fill(107.00,   17.159367);
       //fHistDiffSigfresco->Fill(108.00,   17.447301);
       //fHistDiffSigfresco->Fill(109.00,   17.701691);
       //fHistDiffSigfresco->Fill(110.00,   17.922010);
       //fHistDiffSigfresco->Fill(111.00,   18.108073);
       //fHistDiffSigfresco->Fill(112.00,   18.260010);
       //fHistDiffSigfresco->Fill(113.00,   18.378240);
       //fHistDiffSigfresco->Fill(114.00,   18.463450);
       //fHistDiffSigfresco->Fill(115.00,   18.516564);
       //fHistDiffSigfresco->Fill(116.00,   18.538728);
       //fHistDiffSigfresco->Fill(117.00,   18.531282);
       //fHistDiffSigfresco->Fill(118.00,   18.495740);
       //fHistDiffSigfresco->Fill(119.00,   18.433767);
       //fHistDiffSigfresco->Fill(120.00,   18.347162);
       //fHistDiffSigfresco->Fill(121.00,   18.237834);
       //fHistDiffSigfresco->Fill(122.00,   18.107786);
       //fHistDiffSigfresco->Fill(123.00,   17.959095);
       //fHistDiffSigfresco->Fill(124.00,   17.793892);
       //fHistDiffSigfresco->Fill(125.00,   17.614347);
       //fHistDiffSigfresco->Fill(126.00,   17.422653);
       //fHistDiffSigfresco->Fill(127.00,   17.221007);
       //fHistDiffSigfresco->Fill(128.00,   17.011595);
       //fHistDiffSigfresco->Fill(129.00,   16.796581);
       //fHistDiffSigfresco->Fill(130.00,   16.578086);
       //fHistDiffSigfresco->Fill(30,2296.93593499999);
}
