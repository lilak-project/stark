#include "TStyle.h"

#include "LKRun.h"
#include "LKLogger.h"
#include "GETChannel.h"
#include "SKSiHit.h"
#include "SKAnalysisDK.h"

ClassImp(SKAnalysisDK)

SKAnalysisDK::SKAnalysisDK()
    :LKTask("SKAnalysisDK","SKAnalysisDK")
{
}

bool SKAnalysisDK::Init()
{
    fSiHitArray = fRun -> GetBranchA("SiHit");

    fStarkPlane = (SKSiArrayPlane*) fRun -> FindDetectorPlane("SKSiArrayPlane");
    if (fStarkPlane==nullptr) {
        lk_error << "SKSiArrayPlane do not exist!!!" << endl;
        lk_error << fName << " must run with SKSiArrayPlane" << endl;
        return false;
    }

    int ne = 400;
    double e1 = 0;
    double e2 = 10;
    fPar -> UpdateBinning("ko2421/binning_cal_energy", ne, e1, e2);

    fHistdEEAll[0] = new TH2D("fHistdEE_X6_all","dE_E_X6_all;Etotal;dE",ne,e1,2*e2,ne,e1,2*e2);
    fHistdEEAll[1] = new TH2D("fHistdEE_CSD_all","dE_E_CSD_all;Etotal;dE",ne,e1,2*e2,ne,e1,2*e2);
    auto array0 = new TObjArray();
    array0 -> Add(fHistdEEAll[0]);
    array0 -> Add(fHistdEEAll[1]);
    fStarkPlane -> AddUserDrawingArray("dE_E_all", array0);
    for(int pairID=0; pairID<12; pairID++){
	    auto array2 = new TObjArray();
	    fHistdEE[pairID][0] = new TH2D(Form("fHistdEE_pos_e_%d",pairID),Form("dE_E_%d;dE;Etotal",pairID),ne,e1,2*e2,ne,e1,2*e2);
	    fHistdEE[pairID][1] = new TH2D(Form("fHistdEE_z_e_%d",pairID),Form("dE_E_%d;Zposition;Etotal",pairID),500,0,180,ne,e1,2*e2);
	    array2 -> Add(fHistdEE[pairID][0]);
	    array2 -> Add(fHistdEE[pairID][1]);
	    fStarkPlane -> AddUserDrawingArray("dE_E", pairID, array2, 2);
    }

    fHistZtotE[0] = new TH2D("frelZ_energy","Relative_Z_SumE;Zposition;Energy (MeV)",100,-1,1,ne,e1,e2);
    fHistZtotE[1] = new TH2D("fZ_energy","Z from target_SumE;Zposition (mm);Energy (MeV)",500,0,180,ne,e1,e2);
    /*for (auto det=0; det<40; ++det)
      {
      auto detector = fStarkPlane -> GetSiDetector(det);
      auto name = detector -> GetDetTypeName();
      auto ring = detector -> GetLayer();
      TString sring = "dE"; if (ring==1) sring = "E"; if (ring==2) sring = "16E"; 
      fHistHP[0] -> GetXaxis() -> SetBinLabel(det*fNumJStrips+1,Form("%d (%s,%s)",det,name.Data(),sring.Data()));
      fHistHP[1] -> GetXaxis() -> SetBinLabel(det*fNumOStrips+1,Form("%d (%s,%s)",det,name.Data(),sring.Data()));
      }*/

    auto array = new TObjArray();
    array -> Add(fHistZtotE[0]);
    array -> Add(fHistZtotE[1]);
    fStarkPlane -> AddUserDrawingArray("Z_E_all", array);

    for(int det=0; det<40; det++){
	    auto array4 = new TObjArray();
	    fHistZE[det][0] = new TH2D(Form("fHistZE_rz_e_%d",det),Form("Relative_Z_energy_%d;position;energy",det),100,-1,1,ne,e1,2*e2);
	    fHistZE[det][1] = new TH2D(Form("fHistZE_z_e_%d",det),Form("Z from target_energy_%d;position;energy",det),500,0,180,ne,e1,2*e2);
	    array4 -> Add(fHistZE[det][0]);
	    array4 -> Add(fHistZE[det][1]);
	    fStarkPlane -> AddUserDrawingArray("Z_E_each", det, array4);
    }

    fHistZdet[0] = new TH2D("fZ_Edetector","Z from target_Edetector;Zposition (mm);Energy (MeV)",500,0,180,ne,e1,e2);
    fHistZdet[1] = new TH2D("fZ_dEdetector","Z from target_dEdetector;Zposition (mm);Energy (MeV)",500,0,180,ne,e1,e2);

    auto array3 = new TObjArray();
    array3 -> Add(fHistZdet[0]);
    array3 -> Add(fHistZdet[1]);
    fStarkPlane -> AddUserDrawingArray("Z_dEE", array3);
    
    return true;
}

void SKAnalysisDK::Exec(Option_t*)
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
		bool check_singE = siHit -> IsEPairAndOnlyE();
		bool check_singdE = siHit -> IsEPairAndOnlydE();
		bool check_16E = siHit -> IsNotEPairDetector();

		double z_real = 75*z_relative/2 + z_install;
        //double tot_energy = dE_energy + E_energy;
        //double tot_energy = siHit -> GetEnergyOhmic();
        double tot_energy = siHit -> GetEnergyOhmic() + dE_energy;

		//cout << z_relative << "\t" << z_install << "\t" << z_real << endl;
		//cout << z_real << endl;

		if(check_singE==true){
			if(detID>11 && detID<28){
				//cout << "This is 16RING event " << "E energy: " << E_energy << "Detector ID: " << detID <<  endl;
				fHistZdet[0] -> Fill(z_real,E_energy);
			}
		}
		if(check_singdE==true){
			//cout << "This is 12RING event " << "dE energy: " << dE_energy << "Detector ID: " << detID <<  endl;
				fHistZdet[1] -> Fill(z_real,dE_energy);
		}
        if(check_dEE==true && dE_energy>0){
            int pairID = fStarkPlane -> GetSiDetector(detID) -> GetRow();
            //cout << "This is dE_E event!" << "\tE energy: " << E_energy << "\tdE energy: " << dE_energy << "\tdetID: " << detID << "\tpair ID " << pairID << endl;
            if(pairID<4) fHistdEEAll[0] -> Fill(tot_energy, dE_energy);
            else fHistdEEAll[1] -> Fill(tot_energy, dE_energy);
            fHistdEE[pairID][0] -> Fill(dE_energy, E_energy);
            fHistdEE[pairID][1] -> Fill(z_real, tot_energy);
        }

      		auto detector = fStarkPlane -> GetSiDetector(detID);
      		auto ring = detector -> GetLayer();
            if (ring ==0){
                fHistZE[detID][0] -> Fill(z_relative, dE_energy);
                fHistZE[detID][1] -> Fill(z_real, dE_energy);
                fHistZtotE[0] -> Fill(z_relative, dE_energy);
                fHistZtotE[1] -> Fill(z_real, dE_energy);
            }
            else{
                fHistZE[detID][0] -> Fill(z_relative, E_energy);
                fHistZE[detID][1] -> Fill(z_real, E_energy);
                fHistZtotE[0] -> Fill(z_relative, E_energy);
                fHistZtotE[1] -> Fill(z_real, tot_energy);
            }
            //cout << "=============================================================================" << endl;
    }
}
