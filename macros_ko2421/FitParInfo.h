struct FitParInfo
{
    TString fName = "par";
    double fValue = 0;
    double fLimit1 = 1;
    double fLimit2 = 2;
    bool fFix = false;

    FitParInfo() = default;
    FitParInfo       (TString name, double value, double limit1, double limit2, bool fix) : fName(name), fValue(value), fLimit1(limit1), fLimit2(limit2), fFix(fix) {}
    void SetNameValue(TString name, double value, double limit1, double limit2, bool fix=false) { fName = name; fValue = value; fLimit1 = limit1; fLimit2 = limit2; fFix = fix; }
    void SetValue                  (double value, double limit1, double limit2, bool fix=false) { fValue = value; fLimit1 = limit1; fLimit2 = limit2; fFix = fix; }
    void FixValue                  (double value) { fValue = value; fLimit1 = value; fLimit2 = value; fFix = true; }

    void Print() {
        if (fFix)
            cout << "Fit parameter: " << fName << " = " << fValue << " FIXED!" << endl;
        else
            cout << "Fit parameter: " << fName << " = " << fValue << " within (" << fLimit1 << "," << fLimit2 << ")" << endl;
    }
    void SetPar(TF1 *f1, int iPar) {
        f1 -> SetParName(iPar,fName);
        f1 -> SetParameter(iPar,fValue);
        f1 -> SetParLimits(iPar,fLimit1,fLimit2);
        if (fFix) f1 -> FixParameter(iPar,fValue);
    }
    void PullPar(TF1 *f1, int iPar) {
        fName = f1 -> GetParName(iPar);
        fValue = f1 -> GetParameter(iPar);
        f1 -> GetParLimits(iPar,fLimit1,fLimit2);
        fFix = false;
        if (fLimit1==1&&fLimit2==1&&fValue==0) fFix = true;
        else if (fLimit1==fLimit2) fFix = true;
    }
    void CopyFrom(FitParInfo fpi) {
        fName = fpi.fName;
        fValue = fpi.fValue;
        fLimit1 = fpi.fLimit1;
        fLimit2 = fpi.fLimit2;
        fFix = fpi.fFix;
    }
};

const int fNumParameters = 9;

TString fFormula = "(x>[0])*(exp(-(x-[1]))*(x>[0]?TMath::Power(x-[0],[2]):0)+[3]) + [4]*exp(-0.5*((x-[5])/([5]*0.01*[6]))**2) + [7]*exp(-0.5*((x-[8])/([8]*0.01*[6]))**2)";
TString fFormula1 = "[0]*exp(-0.5*((x-[1])/([1]*0.01*[2]))**2)";
TString fFormula2 = "[0]*exp(-0.5*((x-[1])/([1]*0.01*[2]))**2)";
TString fFormula3 = "(x>[0])*(exp(-(x-[1]))*(x>[0]?TMath::Power(x-[0],[2]):0)+[3])";

FitParInfo fit_par_energy[10];

void InitFitParameters(int energyIndex, int pairID)
{
    fit_par_energy[0].SetNameValue("v0_x_low", 2, 2, 2, true);
    fit_par_energy[1].SetNameValue("v1_exp_pos", 6, 4, 9, false);
    fit_par_energy[2].SetNameValue("v2_x_pow", 1, 0, 2, false);
    fit_par_energy[3].SetNameValue("v3_const", 4, 0, 20, false);
    if (pairID>=4)
        fit_par_energy[3].SetNameValue("v3_const", 4, 0, 40, false);
    fit_par_energy[4].SetNameValue("v4_gaus0_a", 1000, 0, 1.E+3, false);
    fit_par_energy[5].SetNameValue("v5_gaus0_m", 10, 0, 20, false);
    fit_par_energy[6].SetNameValue("v6_gaus0_s", 7, 3, 15, false);
    fit_par_energy[7].SetNameValue("v7_gaus1_a", 1000, 0, 1.E+2, false);
    fit_par_energy[8].SetNameValue("v8_gaus1_m", 10, 0, 20, false);
    if (pairID>=4&&pairID<12) {
        fit_par_energy[7].SetNameValue("v7_gaus1_a", 0,0,0,true);
        fit_par_energy[8].SetNameValue("v8_gaus1_m", 0,0,0,true);
    }
}
