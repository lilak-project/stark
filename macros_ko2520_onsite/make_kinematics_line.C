#include "TCanvas.h"
#include "TFile.h"
#include "TGraph.h"
#include "TMath.h"
#include "TParameter.h"
#include "TString.h"

#include <fstream>
#include <iomanip>
#include <iostream>

double ko2520_elastic_recoil_energy(double thetaLabDeg,
                                    double beamEnergyPerU=8.0,
                                    double projectileA=21.0,
                                    double targetA=1.0)
{
    if (thetaLabDeg < 0 || thetaLabDeg > 90)
        return 0;

    const double beamEnergy = projectileA*beamEnergyPerU;
    const double massFactor = 4.*projectileA*targetA
                            / TMath::Power(projectileA+targetA,2);
    const double cosTheta = TMath::Cos(thetaLabDeg*TMath::DegToRad());
    return massFactor*beamEnergy*cosTheta*cosTheta;
}

void make_kinematics_line(double beamEnergyPerU=8.0,
                          double projectileA=21.0,
                          double targetA=1.0,
                          TString outputRoot="kinematics_21Na_p_8MeVu.root",
                          TString outputData="kinematics_21Na_p_8MeVu.dat")
{
    constexpr int numSteps = 900;
    auto graph = new TGraph(numSteps+1);
    graph->SetName("graph_kinematics_21Na_p");
    graph->SetTitle("^{21}Na+p elastic recoil proton;#theta_{lab} [deg];proton energy [MeV]");

    std::ofstream data(outputData.Data());
    if (!data.is_open()) {
        std::cerr << "Cannot create " << outputData << std::endl;
        return;
    }
    data << "# theta_lab_deg proton_energy_MeV\n";
    data << std::fixed << std::setprecision(6);

    for (int i=0; i<=numSteps; ++i) {
        const double theta = 90.*i/numSteps;
        const double energy = ko2520_elastic_recoil_energy(
            theta,beamEnergyPerU,projectileA,targetA);
        graph->SetPoint(i,theta,energy);
        data << theta << " " << energy << "\n";
    }
    data.close();

    TFile output(outputRoot,"recreate");
    graph->Write();
    TParameter<double>("beam_energy_per_u_MeV",beamEnergyPerU).Write();
    TParameter<double>("projectile_mass_number",projectileA).Write();
    TParameter<double>("target_mass_number",targetA).Write();
    output.Close();

    auto canvas = new TCanvas("canvas_kinematics_21Na_p",
                              "21Na+p elastic kinematics",1000,700);
    graph->SetLineColorAlpha(kMagenta+1,0.50);
    graph->SetLineStyle(2);
    graph->SetLineWidth(2);
    graph->Draw("AL");

    std::cout << "Created " << outputRoot << " and " << outputData << std::endl;
}
