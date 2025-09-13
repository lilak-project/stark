#include "TStyle.h"

#include "LKRun.h"
#include "LKLogger.h"
#include "GETChannel.h"
#include "SKSiHit.h"
#include "SKGatingTask.h"
#include "SKRecoHeader.h"

ClassImp(SKGatingTask)

SKGatingTask::SKGatingTask()
    :LKTask("SKGatingTask","SKGatingTask")
{
}

bool SKGatingTask::Init()
{
    fSiHitArray = fRun -> GetBranchA("SiHit");
    if (fSiHitArray==nullptr) {
        lk_error << "Branch SiHit do not exist!!!" << endl;
        return false;
    }

    fPar -> UpdatePar(fUseEPairOnlydE,"SKGatingTask/useEPairOnlydE");
    fPar -> UpdatePar(fUseEPairOnlyE ,"SKGatingTask/useEPairOnlyE" );
    fPar -> UpdatePar(fUseEPairBoth  ,"SKGatingTask/useEPairBoth"  );
    fPar -> UpdatePar(fUseEPair      ,"SKGatingTask/useEPair"      );
    fPar -> UpdatePar(fUseNotEPair   ,"SKGatingTask/useNotEPair"   );
    fPar -> UpdatePar(fUseScintGate  ,"SKGatingTask/useScintGate"  );

    if (fUseEPairOnlydE) lk_info << "Use EPairOnlydE = " << fUseEPairOnlydE << endl;
    if (fUseEPairOnlyE ) lk_info << "Use EPairOnlyE  = " << fUseEPairOnlyE  << endl;
    if (fUseEPairBoth  ) lk_info << "Use EPairBoth   = " << fUseEPairBoth   << endl;
    if (fUseEPair      ) lk_info << "Use EPair       = " << fUseEPair       << endl;
    if (fUseNotEPair   ) lk_info << "Use NotEPair    = " << fUseNotEPair    << endl;

    if (fUseScintGate)
    {
        fRecoHeaderArray = fRun -> GetBranchA("RecoHeader");
        if (fRecoHeaderArray==nullptr) {
            lk_error << "Trying to use scintillator gate but branch RecoHeader for scintillator gate do not exist!" << endl;
            return false;
        }
        lk_info << "Use ScintGate   = " << fUseScintGate    << endl;
    }

    TString fileName1, fileName2, fileName3;
    fPar -> UpdatePar(fCutDetID,"SKGatingTask/gate_detID");
    fPar -> UpdatePar(fCutdEDetID,"SKGatingTask/gate_deDetID");
    fPar -> UpdatePar(fileName1,"SKGatingTask/file_cutG_Etot_Theta");
    fPar -> UpdatePar(fileName2,"SKGatingTask/file_cutG_dE_Etot");
    fPar -> UpdatePar(fileName3,"SKGatingTask/file_cutG_dE_E");

    fDetIDRange[0] = -1;
    fDetIDRange[1] = -1;
    fdEDetIDRange[0] = -1;
    fdEDetIDRange[1] = -1;
    if (fPar -> CheckPar("SKGatingTask/DetIDRange"))
    {
        fDetIDRange[0] = fPar -> GetParInt("SKGatingTask/DetIDRange",0);
        fDetIDRange[1] = fPar -> GetParInt("SKGatingTask/DetIDRange",1);
    }
    if (fPar -> CheckPar("SKGatingTask/dEDetIDRange"))
    {
        fdEDetIDRange[0] = fPar -> GetParInt("SKGatingTask/dEDetIDRange",0);
        fdEDetIDRange[1] = fPar -> GetParInt("SKGatingTask/dEDetIDRange",1);
    }

    fPar -> Require("SKGatingTask/KeyEnergyType","jj","evaluation of key energy jj: Etot=dEJunction+EJunction, jo: Etot=dEOhmic+EOhmic, oo: Etot=dEOhmic+EOhmic, xj: Etot=EJunction");

    TString keyEnergyType = "jj";
    fPar -> UpdatePar(keyEnergyType,"SKGatingTask/KeyEnergyType");

    if      (keyEnergyType=="jj") fKeyEnergyType = 00;
    else if (keyEnergyType=="jo") fKeyEnergyType = 10;
    else if (keyEnergyType=="oj") fKeyEnergyType = 01;
    else if (keyEnergyType=="oo") fKeyEnergyType = 11;
    else if (keyEnergyType=="xj") fKeyEnergyType = 20;
    else if (keyEnergyType=="xo") fKeyEnergyType = 21;
    else if (keyEnergyType=="ox") fKeyEnergyType = 12;
    else if (keyEnergyType=="jx") fKeyEnergyType = 02;

    if (!fileName1.IsNull()) {
        auto file = new TFile(fileName1);
        if (file->IsOpen()==false)
            lk_error << "File " << fileName1 << " is bad!" << endl;
        else {
            fCutG_Theta_ETot = (TCutG*) file -> Get("CUTG");
            lk_info << "fCutG_Theta_ETot = " << fileName1 << " " << fCutG_Theta_ETot << endl;
        }
    }
    if (!fileName2.IsNull()) {
        auto file = new TFile(fileName2);
        if (file->IsOpen()==false)
            lk_error << "File " << fileName2 << " is bad!" << endl;
        else {
            fCutG_ETot_dE = (TCutG*) file -> Get("CUTG");
            lk_info << "fCutG_ETot_dE = " << fileName2 << " " << fCutG_ETot_dE << endl;
        }
    }
    if (!fileName3.IsNull()) {
        auto file = new TFile(fileName3);
        if (file->IsOpen()==false)
            lk_error << "File " << fileName3 << " is bad!" << endl;
        else {
            fCutG_E_dE = (TCutG*) file -> Get("CUTG");
            lk_info << "fCutG_E_dE = " << fileName2 << " " << fCutG_E_dE << endl;
        }
    }

    return true;
}

void SKGatingTask::Exec(Option_t*)
{
    bool passedScintGate = true;
    if (fUseScintGate) {
        auto recoHeader = (SKRecoHeader*) fRecoHeaderArray -> At(0);
        passedScintGate = recoHeader -> BeamOnScint();
    }

    auto numHits = fSiHitArray -> GetEntries();
    for (auto iHit=0; iHit<numHits; ++iHit)
    {
        auto siHit = (SKSiHit*) fSiHitArray -> At(iHit);

        int    detID = siHit -> GetDetID();
        int    dEDetID = siHit -> GetdEDetID();
        double E_energy_j = siHit -> GetE();
        double E_energy_o = siHit -> GetEnergyOhmic();
        double dE_energy_j = siHit -> GetdE();
        double dE_energy_o = siHit -> GetdEOhmic();
        double relativeZ = siHit -> GetRelativeZ();
        double detZ = siHit -> GetDetectorZ();
        double finalZ = siHit -> GetFinalZ(75);
        double theta = siHit -> GetTheta();

        double key_energy = siHit -> GetEnergyTotal();
        if      (fKeyEnergyType==00) key_energy = E_energy_j + dE_energy_j; // jj
        else if (fKeyEnergyType==10) key_energy = E_energy_j + dE_energy_o; // jo
        else if (fKeyEnergyType==01) key_energy = E_energy_o + dE_energy_j; // oj
        else if (fKeyEnergyType==11) key_energy = E_energy_o + dE_energy_o; // oo
        else if (fKeyEnergyType==20) key_energy =          0 + dE_energy_j; // xj
        else if (fKeyEnergyType==21) key_energy =          0 + dE_energy_o; // xo
        else if (fKeyEnergyType==12) key_energy = E_energy_o +           0; // ox
        else if (fKeyEnergyType==02) key_energy = E_energy_j +           0; // jx
        siHit -> SetKeyEnergy(key_energy);

        bool checkEPairOnlydE = siHit -> IsEPairAndOnlydE();
        bool checkEPairOnlyE = siHit -> IsEPairAndOnlyE();
        bool checkEPairBoth = siHit -> IsEPairAndBothdEE();
        bool checkEPair = siHit -> IsEPairDetector();
        bool checkNotEPair = siHit -> IsNotEPairDetector();

        siHit -> SetInGate(false);

        //if (dEDetID>=28 && dEDetID<=31) lk_debug << checkEPair << " " << fdEDetIDRange[0] << " " << dEDetID << " " << (dEDetID>fdEDetIDRange[0]) << " " <<  (dEDetID<fdEDetIDRange[1]) << endl;

        if (!passedScintGate) continue;

        if (fDetIDRange[0]>=0 && !(detID>fDetIDRange[0]&&detID<fDetIDRange[1])) continue;
        if (fdEDetIDRange[0]>=0 && !(dEDetID>fdEDetIDRange[0]&&dEDetID<fdEDetIDRange[1])) continue;

        if (fUseEPairOnlydE&& !checkEPairOnlydE) continue;
        if (fUseEPairOnlyE && !checkEPairOnlyE ) continue;
        if (fUseEPairBoth  && !checkEPairBoth  ) continue;
        if (fUseEPair      && !checkEPair      ) continue;
        if (fUseNotEPair   && !checkNotEPair   ) continue;

        if (fCutDetID>=0&&fCutDetID!=detID) continue;
        if (fCutdEDetID>=0&&fCutdEDetID!=dEDetID) continue;
        if (fCutG_Theta_ETot!=nullptr&&fCutG_Theta_ETot->IsInside(theta,key_energy)==false) continue;
        if (fCutG_ETot_dE!=nullptr&&fCutG_ETot_dE->IsInside(key_energy,dE_energy_j)==false) continue;
        if (fCutG_E_dE!=nullptr&&fCutG_E_dE->IsInside(E_energy_j,dE_energy_j)==false) continue;

        siHit -> SetInGate(true);
    }
}
