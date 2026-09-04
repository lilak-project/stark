#include <algorithm>
#include <cmath>

namespace {
// NIST PSTAR CSDA ranges for protons in silicon
// (https://physics.nist.gov/PhysRefData/Star/Text/PSTAR.html).
// Energy is in MeV and
// range is in g/cm^2.  The table is interpolated log-log, matching the
// smooth behavior of the tabulated CSDA range.
constexpr int kNumPstarSi = 52;
constexpr double kPstarSiEnergy[kNumPstarSi] = {
    0.100,0.125,0.150,0.175,0.200,0.225,0.250,0.275,0.300,
    0.350,0.400,0.450,0.500,0.550,0.600,0.650,0.700,0.750,
    0.800,0.850,0.900,0.950,1.000,1.250,1.500,1.750,2.000,
    2.250,2.500,2.750,3.000,3.500,4.000,4.500,5.000,5.500,
    6.000,6.500,7.000,7.500,8.000,8.500,9.000,9.500,10.00,
    12.50,15.00,17.50,20.00,25.00,27.50,30.00
};
constexpr double kPstarSiRange[kNumPstarSi] = {
    2.308e-4,2.821e-4,3.373e-4,3.963e-4,4.588e-4,5.247e-4,
    5.939e-4,6.661e-4,7.413e-4,9.001e-4,1.070e-3,1.249e-3,
    1.439e-3,1.638e-3,1.846e-3,2.063e-3,2.289e-3,2.523e-3,
    2.767e-3,3.020e-3,3.283e-3,3.555e-3,3.835e-3,5.369e-3,
    7.112e-3,9.057e-3,1.120e-2,1.353e-2,1.605e-2,1.874e-2,
    2.162e-2,2.790e-2,3.485e-2,4.247e-2,5.073e-2,5.962e-2,
    6.914e-2,7.926e-2,8.999e-2,1.013e-1,1.132e-1,1.257e-1,
    1.387e-1,1.524e-1,1.665e-1,2.456e-1,3.380e-1,4.432e-1,
    5.608e-1,8.322e-1,9.854e-1,1.150
};

double interpolateLogLog(double x, const double *xs, const double *ys, int n)
{
    if (x <= xs[0]) return ys[0]*x/xs[0];
    int upper = std::lower_bound(xs,xs+n,x)-xs;
    if (upper >= n) upper = n-1;
    int lower = upper-1;
    double fraction = std::log(x/xs[lower])/std::log(xs[upper]/xs[lower]);
    return ys[lower]*std::exp(fraction*std::log(ys[upper]/ys[lower]));
}

double pstarSiRange(double energyMeV)
{
    if (energyMeV <= 0) return 0;
    return interpolateLogLog(
        energyMeV,kPstarSiEnergy,kPstarSiRange,kNumPstarSi);
}

double pstarSiEnergyFromRange(double range)
{
    if (range <= 0) return 0;
    return interpolateLogLog(
        range,kPstarSiRange,kPstarSiEnergy,kNumPstarSi);
}

double qqq5DepositedEnergy(double incidentEnergyMeV, double thetaLabDeg)
{
    constexpr double siliconDensity = 2.329; // g/cm^3
    constexpr double qqq5ThicknessCm = 0.100; // 1 mm
    double cosIncidence = std::abs(
        std::cos(thetaLabDeg*TMath::DegToRad()));
    if (cosIncidence < 1.e-6) return incidentEnergyMeV;
    double traversedArealDensity =
        siliconDensity*qqq5ThicknessCm/cosIncidence;
    double incidentRange = pstarSiRange(incidentEnergyMeV);
    if (incidentRange <= traversedArealDensity)
        return incidentEnergyMeV;
    double exitEnergy = pstarSiEnergyFromRange(
        incidentRange-traversedArealDensity);
    return std::max(0.,incidentEnergyMeV-exitEnergy);
}
}

void draw_sihit(int mode=2, TString inputName="", TString outputName="")
{
    constexpr int maxDetectorIndex = 256;

    if (mode < 0 || mode > 2) {
        cerr << "Invalid mode " << mode
             << ". Use 0 (relative position), 1 (z position), or 2 (theta)."
             << endl;
        return;
    }

    double energyMax = 30;
    Long64_t maxEvents = -1;

    auto par = new LKParameterContainer("get_common.mac");
    auto runID = par->GetParInt("LKRun/RunID");
    cout << "run " << runID << endl;

    auto tree = new TChain("event");
    if (!inputName.IsNull()) {
        cout << inputName << endl;
        tree->Add(inputName);
    }
    else if (runID == 8000) {
        for (int inputRun:{
                83,84,85,86,87,88,89,90,
                97,98,99,100,101,102,105,
                112,113,114,
                134,135,136,139,141,142,146,147,148})
            tree->Add(Form("data_reco/ko2520_%04d.sihit.root",inputRun));
    }
    else if (runID == 8001) {
        // Production runs acquired on 2026-09-02.
        for (int inputRun:{83,84,85,86,87,88,89,90})
            tree->Add(Form("data_reco/ko2520_%04d.sihit.root",inputRun));
    }
    else if (runID == 8002) {
        // Production runs acquired on 2026-09-03.
        for (int inputRun:{97,98,99,100,101,102,105,112,113,114})
            tree->Add(Form("data_reco/ko2520_%04d.sihit.root",inputRun));
    }
    else if (runID == 8003) {
        // Production runs acquired on 2026-09-04.
        for (int inputRun:{134,135,136,139,141,142,146,147,148})
            tree->Add(Form("data_reco/ko2520_%04d.sihit.root",inputRun));
    }
    else if (runID == 8004)
    {
        for (int inputRun:{134, 135, 136, 139, 141, 142})
            tree->Add(Form("data_reco/ko2520_%04d.sihit.root",inputRun));
    }
    else {
        inputName = Form("data_reco/ko2520_%04d.sihit.root",runID);
        cout << inputName << endl;
        tree->Add(inputName);
    }
    if (tree->GetNtrees() == 0) {
        cerr << "No sihit input file was found for run " << runID << endl;
        return;
    }
    if (tree->GetBranch("SiHit") == nullptr) {
        cerr << "SiHit branch is not present in the input file" << endl;
        return;
    }

    LKSiliconMapping mapping;
    if (!mapping.Load("mapping_ko2520_0904")) {
        cerr << "Failed to load mapping_ko2520_0904" << endl;
        return;
    }
    // QQQ5 detIndex 32/33 are on AGET2 (gain x5 in the 8003 runs).
    // detIndex 34/35 are on AGET3 (gain unchanged) and are intentionally
    // omitted from every draw_sihit histogram.
    auto isExcludedQQQ5Aget3 = [](const LKSiliconMapping::DetectorInfo *det) {
        return det != nullptr && det->detType == "QQQ5" &&
               (det->detIndex == 34 || det->detIndex == 35);
    };

    auto x6AxisTitle = [mode]() -> TString {
        if (mode == 1) return "z position [mm]";
        if (mode == 2) return "#theta_{lab} [deg]";
        return "(E_{high}-E_{low})/(E_{high}+E_{low})";
    };
    auto x6Coordinate = [mode](const LKSiliconMapping::DetectorInfo *det,
                               double relativePosition) -> double {
        if (mode == 0 || det == nullptr)
            return relativePosition;
        // ringZ is the detector center and detHeight is the active z length.
        double z = det->ringZ + 0.5*det->detHeight*relativePosition;
        if (mode == 1)
            return z;
        // ringRadius and z are both read from the detector mapping geometry.
        return TMath::ATan2(det->ringRadius,z)*TMath::RadToDeg();
    };
    auto x6AxisRange = [mode, &x6Coordinate](
        const LKSiliconMapping::DetectorInfo *det, double &xMin, double &xMax) {
        if (mode == 0 || det == nullptr) {
            xMin = -1.2;
            xMax = 1.2;
            return;
        }
        double x1 = x6Coordinate(det,-1.0);
        double x2 = x6Coordinate(det, 1.0);
        xMin = TMath::Min(x1,x2);
        xMax = TMath::Max(x1,x2);
    };

    // QQQ5 junction strips are numbered 1--32 in channel_mapping.txt.  The
    // radii below are the measured active inner/outer edges in micrometres,
    // ordered from the innermost strip to the outermost strip.
    const double qqq5InnerRadiusUm[32] = {
        25250,27800,30300,32750,35150,37500,39800,42050,
        44250,46400,48500,50550,52500,54500,56400,58250,
        60050,61800,63500,65150,66750,68300,69800,71250,
        72650,74000,75300,76550,77750,78900,80000,81050
    };
    const double qqq5OuterRadiusUm[32] = {
        27700,30200,32650,35050,37400,39700,41950,44150,
        46300,48400,50450,52450,54400,56300,58150,59950,
        61700,63400,65050,66650,68200,69700,71150,72550,
        73900,75200,76450,77650,78800,79900,80950,81950
    };
    auto qqq5Coordinate = [mode, &qqq5InnerRadiusUm, &qqq5OuterRadiusUm](
        const LKSiliconMapping::DetectorInfo *det, int junctionStrip) -> double {
        if (junctionStrip < 1 || junctionStrip > 32)
            return -999;
        if (mode != 2 || det == nullptr)
            return junctionStrip;
        const int index = junctionStrip-1;
        const double radiusMm = 0.0005*
            (qqq5InnerRadiusUm[index]+qqq5OuterRadiusUm[index]);
        return TMath::ATan2(radiusMm,det->ringZ)*TMath::RadToDeg();
    };
    auto qqq5AxisRange = [mode](const LKSiliconMapping::DetectorInfo *det,
                                double &xMin, double &xMax) {
        if (mode != 2 || det == nullptr) {
            xMin = 0.5;
            xMax = 32.5;
            return;
        }
        // Use the full active radial edges, 25.25--81.95 mm.
        xMin = TMath::ATan2(25.25,det->ringZ)*TMath::RadToDeg();
        xMax = TMath::ATan2(81.95,det->ringZ)*TMath::RadToDeg();
    };

    cout << "X6 horizontal-axis mode " << mode << ": " << x6AxisTitle() << endl;

    TString nameTop = Form("sihit_energy_position_%04d",runID);
    auto topPosition = new LKDrawingGroup(nameTop);
    auto position16   = topPosition->CreateGroup("X6_16_ring");
    auto position12DE = topPosition->CreateGroup("X6_12dE_ring");
    auto position12E  = topPosition->CreateGroup("X6_12E_ring");
    auto positionQQQ5 = topPosition->CreateGroup("QQQ5");
    auto position12Sum = topPosition->CreateGroup("X6_12_dE_plus_E");
    auto pid12TotalVsDE = topPosition->CreateGroup("X6_12_dE_vs_E_plus_dE");
    auto pid12DEVsE = topPosition->CreateGroup("X6_12_dE_vs_E");

    TString nameTop2 = Form("sihit_junction_ohmic_%04d",runID);
    auto topOhmic = new LKDrawingGroup(nameTop2);
    auto ohmic16   = topOhmic->CreateGroup("X6_16_ring");
    auto ohmic12DE = topOhmic->CreateGroup("X6_12dE_ring");
    auto ohmic12E  = topOhmic->CreateGroup("X6_12E_ring");
    auto ohmicQQQ5 = topOhmic->CreateGroup("QQQ5");

    TString nameTopOhmicPosition = Form("sihit_ohmic_energy_position_%04d",runID);
    auto topOhmicPosition = new LKDrawingGroup(nameTopOhmicPosition);
    auto ohmicPosition16   = topOhmicPosition->CreateGroup("X6_16_ring");
    auto ohmicPosition12DE = topOhmicPosition->CreateGroup("X6_12dE_ring");
    auto ohmicPosition12E  = topOhmicPosition->CreateGroup("X6_12E_ring");
    auto ohmicPositionQQQ5 = topOhmicPosition->CreateGroup("QQQ5");
    auto pid12OhmicDEVsE = topOhmicPosition->CreateGroup("X6_12_ohmic_dE_vs_E");
    auto ohmic12TotalVsPosition = topOhmicPosition->CreateGroup(
        "X6_12_ohmic_dE_plus_E_vs_position");

    TString nameTopRingSum = Form("sihit_ring_sum_energy_position_%04d",runID);
    auto topRingSum = new LKDrawingGroup(nameTopRingSum);

    TString nameTopOhmicRingSum = Form("sihit_ring_sum_ohmic_energy_position_%04d",runID);
    auto topOhmicRingSum = new LKDrawingGroup(nameTopOhmicRingSum);

    TString nameTop3 = Form("sihit_multiplicity_%04d",runID);
    auto topMultiplicity = new LKDrawingGroup(nameTop3);

    enum { kAll=0, kRing16=1, kRing12DE=2, kRing12E=3, kQQQ5=4, kNumMult=5 };
    enum { kSum16=0, kSum12DE=1, kSum12E=2, kSumQQQ5=3, kNumRingSums=4 };
    TString ringSumName[kNumRingSums] = {"X6_16", "X6_12dE", "X6_12E", "QQQ5"};
    TH2D *histRingSum[kNumRingSums] = {nullptr};
    TH2D *histOhmicRingSum[kNumRingSums] = {nullptr};
    for (int i=0; i<kNumRingSums; ++i) {
        bool isQQQ5 = (i == kSumQQQ5);
        const LKSiliconMapping::DetectorInfo *ringDetector = nullptr;
        if (!isQQQ5) {
            TString wantedRing = ringSumName[i];
            wantedRing.ReplaceAll("X6_","");
            for (int iDet=0; iDet<mapping.GetNumDetectors(); ++iDet) {
                auto candidate = mapping.GetDetectorByVectorIndex(iDet);
                if (candidate != nullptr && candidate->ringType == wantedRing) {
                    ringDetector = candidate;
                    break;
                }
            }
        }
        else {
            for (int iDet=0; iDet<mapping.GetNumDetectors(); ++iDet) {
                auto candidate = mapping.GetDetectorByVectorIndex(iDet);
                if (candidate != nullptr && candidate->detType == "QQQ5") {
                    ringDetector = candidate;
                    break;
                }
            }
        }
        double ringXMin = 0.5;
        double ringXMax = 32.5;
        if (isQQQ5)
            qqq5AxisRange(ringDetector,ringXMin,ringXMax);
        else
            x6AxisRange(ringDetector,ringXMin,ringXMax);
        TString ringAxisTitle = (isQQQ5 && mode != 2)
                              ? "junction strip" : x6AxisTitle();
        int ringXBins = (isQQQ5 && mode != 2) ? 32 : 200;
        histRingSum[i] = new TH2D(
            Form("hist_sihit_ring_sum_energy_position_%s",ringSumName[i].Data()),
            Form("%s detector-summed junction energy vs position;%s;junction energy",
                 ringSumName[i].Data(),ringAxisTitle.Data()),
            ringXBins,ringXMin,ringXMax,
            300,0,energyMax);
        auto drawing = topRingSum->CreateDrawing(ringSumName[i]);
        drawing->Add(histRingSum[i],"colz");

        histOhmicRingSum[i] = new TH2D(
            Form("hist_sihit_ring_sum_ohmic_energy_position_%s",ringSumName[i].Data()),
            Form("%s detector-summed ohmic energy vs junction position;%s;ohmic energy",
                 ringSumName[i].Data(),ringAxisTitle.Data()),
            ringXBins,ringXMin,ringXMax,
            300,0,energyMax);
        auto ohmicDrawing = topOhmicRingSum->CreateDrawing(ringSumName[i]);
        ohmicDrawing->Add(histOhmicRingSum[i],"colz");
    }
    const LKSiliconMapping::DetectorInfo *representative16 = nullptr;
    const LKSiliconMapping::DetectorInfo *representative12E = nullptr;
    const LKSiliconMapping::DetectorInfo *representativeQQQ5 = nullptr;
    for (int iDet=0; iDet<mapping.GetNumDetectors(); ++iDet) {
        auto candidate = mapping.GetDetectorByVectorIndex(iDet);
        if (candidate == nullptr)
            continue;
        if (representative16 == nullptr && candidate->ringType == "16")
            representative16 = candidate;
        if (representative12E == nullptr && candidate->ringType == "12E")
            representative12E = candidate;
        if (representativeQQQ5 == nullptr && candidate->detType == "QQQ5" &&
            !isExcludedQQQ5Aget3(candidate))
            representativeQQQ5 = candidate;
    }
    double pairXMin = -1.2;
    double pairXMax = 1.2;
    x6AxisRange(representative12E,pairXMin,pairXMax);
    auto histRingSumTotalEnergy = new TH2D(
        "hist_sihit_ring_sum_total_energy_position_X6_12",
        Form("12-ring detector-summed coincidence;E-detector %s;dE+E",
             x6AxisTitle().Data()),
        200,pairXMin,pairXMax,600,0,2*energyMax);
    auto ringSumTotalDrawing = topRingSum->CreateDrawing("X6_12_dE_plus_E");
    ringSumTotalDrawing->Add(histRingSumTotalEnergy,"colz");
    auto histRingSumDEVsE = new TH2D(
        "hist_sihit_ring_sum_de_vs_e_X6_12",
        "12-ring detector-summed junction coincidence;E;dE",
        300,0,energyMax,300,0,energyMax);
    auto ringSumDEVsEDrawing = topRingSum->CreateDrawing("X6_12_dE_vs_E");
    ringSumDEVsEDrawing->Add(histRingSumDEVsE,"colz");
    auto histRingSumDEVsTotal = new TH2D(
        "hist_sihit_ring_sum_de_vs_total_X6_12",
        "12-ring detector-summed junction coincidence;E+dE;dE",
        600,0,2*energyMax,300,0,energyMax);
    auto ringSumDEVsTotalDrawing = topRingSum->CreateDrawing("X6_12_dE_vs_E_plus_dE");
    ringSumDEVsTotalDrawing->Add(histRingSumDEVsTotal,"colz");
    auto histRingSumOhmicDEVsE = new TH2D(
        "hist_sihit_ring_sum_ohmic_de_vs_e_X6_12",
        "12-ring detector-summed ohmic coincidence;E_{ohmic};dE_{ohmic}",
        300,0,energyMax,300,0,energyMax);
    auto ringSumOhmicDEVsEDrawing = topOhmicRingSum->CreateDrawing("X6_12_ohmic_dE_vs_E");
    ringSumOhmicDEVsEDrawing->Add(histRingSumOhmicDEVsE,"colz");
    auto histRingSumOhmicTotalVsPosition = new TH2D(
        "hist_sihit_ring_sum_ohmic_total_vs_position_X6_12",
        Form("12-ring detector-summed ohmic coincidence;E-detector %s;dE_{ohmic}+E_{ohmic}",
             x6AxisTitle().Data()),
        200,pairXMin,pairXMax,600,0,2*energyMax);
    auto ringSumOhmicTotalVsPositionDrawing = topOhmicRingSum->CreateDrawing(
        "X6_12_ohmic_dE_plus_E_vs_position");
    ringSumOhmicTotalVsPositionDrawing->Add(histRingSumOhmicTotalVsPosition,"colz");

    // Join the first (16-ring ohmic) and last (12-ring dE+E ohmic)
    // spectra from the ring-sum canvas into one continuous coordinate range.
    double ring16XMin = -1.2;
    double ring16XMax = 1.2;
    x6AxisRange(representative16,ring16XMin,ring16XMax);
    double qqq5XMin = 0.5;
    double qqq5XMax = 32.5;
    qqq5AxisRange(representativeQQQ5,qqq5XMin,qqq5XMax);
    double combinedXMin = TMath::Min(pairXMin,ring16XMin);
    double combinedXMax = TMath::Max(pairXMax,ring16XMax);
    if (mode == 2) {
        combinedXMin = TMath::Min(combinedXMin,qqq5XMin);
        combinedXMax = TMath::Max(combinedXMax,qqq5XMax);
    }
    auto histCombinedOhmic = new TH2D(
        "hist_sihit_combined_ohmic_X6_12_16",
        Form("QQQ5 junction + X6 ohmic combined spectrum;%s;energy [MeV]",
             x6AxisTitle().Data()),
        400,combinedXMin,combinedXMax,600,0,2*energyMax);

    TGraph *graphKinematics = nullptr;
    TGraph *graphKinematicsQQQ5 = nullptr;
    if (mode == 2) {
        TFile kinematicsFile("kinematics_21Na_p_8MeVu.root","read");
        auto storedGraph = dynamic_cast<TGraph*>(
            kinematicsFile.Get("graph_kinematics_21Na_p"));
        if (storedGraph != nullptr)
            graphKinematics = dynamic_cast<TGraph*>(storedGraph->Clone(
                "graph_kinematics_21Na_p_overlay"));
        else
            cerr << "Kinematics graph is missing. Run make_kinematics_line.C first."
                 << endl;
        if (graphKinematics != nullptr) {
            graphKinematicsQQQ5 = new TGraph();
            graphKinematicsQQQ5->SetName(
                "graph_kinematics_21Na_p_qqq5_1mm_si");
            graphKinematicsQQQ5->SetTitle(
                "^{21}Na+p elastic: energy deposited in 1 mm Si;#theta_{lab} [deg];#DeltaE_{Si} [MeV]");
            for (int i=0; i<graphKinematics->GetN(); ++i) {
                double theta = 0;
                double incidentEnergy = 0;
                graphKinematics->GetPoint(i,theta,incidentEnergy);
                if (theta < qqq5XMin || theta > qqq5XMax)
                    continue;
                graphKinematicsQQQ5->SetPoint(
                    graphKinematicsQQQ5->GetN(),theta,
                    qqq5DepositedEnergy(incidentEnergy,theta));
            }
        }
    }

    constexpr double kinematicsCutMeV = 1.0;
    constexpr double projectileMass = 21.0;
    constexpr double protonMass = 1.0;
    constexpr double beamEnergyPerU = 8.0;
    constexpr double beamEnergy = projectileMass*beamEnergyPerU;
    auto histQValueVsThetaCM = new TH2D(
        "hist_sihit_qvalue_vs_theta_cm_X6_12_16",
        "^{21}Na+p candidates (|E-E_{kin}| < 1 MeV);#theta_{CM} [deg];Q-value [MeV]",
        360,0,180,400,-2,2);
    auto histQValue = new TH1D(
        "hist_sihit_qvalue_X6_12_16",
        "^{21}Na+p candidates (|E-E_{kin}| < 1 MeV);Q-value [MeV];counts",
        400,-2,2);
    auto histThetaCM = new TH1D(
        "hist_sihit_theta_cm_X6_12_16",
        "^{21}Na+p candidates (|E-E_{kin}| < 1 MeV);#theta_{CM} [deg];counts",
        360,0,180);
    auto histQValueVsThetaCMAll = new TH2D(
        "hist_sihit_qvalue_vs_theta_cm_all_X6_12_16",
        "All QQQ5 + X6 ring events (no kinematics cut);#theta_{CM} [deg];Q-value [MeV]",
        360,0,180,700,-10,60);
    auto histQValueAll = new TH1D(
        "hist_sihit_qvalue_all_X6_12_16",
        "All QQQ5 + X6 ring events (no kinematics cut);Q-value [MeV];counts",
        700,-10,60);
    auto histThetaCMAll = new TH1D(
        "hist_sihit_theta_cm_all_X6_12_16",
        "All QQQ5 + X6 ring events (no kinematics cut);#theta_{CM} [deg];counts",
        360,0,180);
    Long64_t numQValueCandidates = 0;
    Long64_t numAllQValueEvents = 0;
    auto fillQValue = [&](double thetaLab, double energy, bool isQQQ5=false) {
        if (mode != 2 || graphKinematics == nullptr ||
            thetaLab < 0 || thetaLab > 90 || energy <= 0)
            return;

        const double cosTheta = TMath::Cos(thetaLab*TMath::DegToRad());
        const double qValue = energy*(1.+protonMass/projectileMass)
            - 2./projectileMass
            * TMath::Sqrt(projectileMass*protonMass*beamEnergy*energy)
            * cosTheta;
        const double thetaCM = 180.-2.*thetaLab;
        histQValueVsThetaCMAll->Fill(thetaCM,qValue);
        histQValueAll->Fill(qValue);
        histThetaCMAll->Fill(thetaCM);
        ++numAllQValueEvents;

        const double expectedEnergy =
            (isQQQ5 && graphKinematicsQQQ5 != nullptr)
            ? graphKinematicsQQQ5->Eval(thetaLab)
            : graphKinematics->Eval(thetaLab);
        if (TMath::Abs(energy-expectedEnergy) > kinematicsCutMeV)
            return;
        histQValueVsThetaCM->Fill(thetaCM,qValue);
        histQValue->Fill(qValue);
        histThetaCM->Fill(thetaCM);
        ++numQValueCandidates;
    };

    TString multName[kNumMult] = {"all", "X6_16", "X6_12dE", "X6_12E", "QQQ5"};
    TH1D *histMultiplicity[kNumMult] = {nullptr};
    for (int i=0; i<kNumMult; ++i) {
        histMultiplicity[i] = new TH1D(
            Form("hist_sihit_multiplicity_%s",multName[i].Data()),
            Form("%s SiHit multiplicity;number of detector hits;events",multName[i].Data()),
            50,0,50);
        auto drawing = topMultiplicity->CreateDrawing(multName[i]);
        drawing->Add(histMultiplicity[i],"hist");
    }

    TH2D *histPosition[maxDetectorIndex] = {nullptr};
    TH2D *histOhmic[maxDetectorIndex] = {nullptr};
    TH2D *histOhmicPosition[maxDetectorIndex] = {nullptr};
    TH2D *histTotalEnergy[maxDetectorIndex] = {nullptr};
    TH2D *histTotalVsDE[maxDetectorIndex] = {nullptr};
    TH2D *histDEVsE[maxDetectorIndex] = {nullptr};
    TH2D *histOhmicDEVsE[maxDetectorIndex] = {nullptr};
    TH2D *histOhmicTotalVsPosition[maxDetectorIndex] = {nullptr};
    int dEDetectorByE[maxDetectorIndex];
    Long64_t numFilled[maxDetectorIndex] = {0};
    Long64_t numCoincidence[maxDetectorIndex] = {0};
    Long64_t numOhmicCoincidence[maxDetectorIndex] = {0};

    // Dedicated QQQ5-only view.  Keep the horizontal coordinate as the
    // physical junction strip even when mode=2 changes other plots to theta.
    auto histQQQ5JunctionEnergyVsStrip = new TH2D(
        "hist_sihit_qqq5_junction_energy_vs_junction_strip",
        "QQQ5 junction energy;Junction strip;Junction energy [MeV]",
        32,0.5,32.5,300,0,energyMax);
    auto histQQQ5OhmicEnergyVsJunctionStrip = new TH2D(
        "hist_sihit_qqq5_ohmic_energy_vs_junction_strip",
        "QQQ5 ohmic energy;Junction strip;Ohmic energy [MeV]",
        32,0.5,32.5,300,0,energyMax);
    for (int id=0; id<maxDetectorIndex; ++id)
        dEDetectorByE[id] = -1;

    for (int iDet=0; iDet<mapping.GetNumDetectors(); ++iDet) {
        auto det = mapping.GetDetectorByVectorIndex(iDet);
        if (det == nullptr || det->detIndex < 0 || det->detIndex >= maxDetectorIndex)
            continue;
        if (isExcludedQQQ5Aget3(det))
            continue;

        int id = det->detIndex;
        if (det->detType == "QQQ5") {
            double xMin = 0.5;
            double xMax = 32.5;
            qqq5AxisRange(det,xMin,xMax);
            TString axisTitle = (mode == 2) ? x6AxisTitle() : "junction strip";
            histPosition[id] = new TH2D(
                Form("hist_sihit_energy_position_det_%03d",id),
                Form("QQQ5 det#%d (idx=%d);%s;restored energy",
                     det->detNumber,id,axisTitle.Data()),
                mode == 2 ? 200 : 32,xMin,xMax,300,0,energyMax);
        }
        else {
            double xMin = -1.2;
            double xMax = 1.2;
            x6AxisRange(det,xMin,xMax);
            histPosition[id] = new TH2D(
                Form("hist_sihit_energy_position_det_%03d",id),
                Form("X6 det#%d (idx=%d, %s);%s;E_{high}+E_{low}",
                     det->detNumber,id,det->ringType.Data(),x6AxisTitle().Data()),
                200,xMin,xMax,300,0,energyMax);
        }
        histOhmic[id] = new TH2D(
            Form("hist_sihit_junction_ohmic_det_%03d",id),
            Form("%s det#%d (idx=%d, %s);junction energy;ohmic energy",
                 det->detType.Data(),det->detNumber,id,det->ringType.Data()),
            300,0,energyMax,300,0,energyMax);
        bool isQQQ5 = (det->detType == "QQQ5");
        double xMin = 0.5;
        double xMax = 32.5;
        if (isQQQ5)
            qqq5AxisRange(det,xMin,xMax);
        else
            x6AxisRange(det,xMin,xMax);
        TString positionAxisTitle = (isQQQ5 && mode != 2)
                                  ? "junction strip" : x6AxisTitle();
        histOhmicPosition[id] = new TH2D(
            Form("hist_sihit_ohmic_energy_position_det_%03d",id),
            Form("%s det#%d (idx=%d, %s);%s;ohmic energy",
                 det->detType.Data(),det->detNumber,id,det->ringType.Data(),
                 positionAxisTitle.Data()),
            (isQQQ5 && mode != 2) ? 32 : 200,xMin,xMax,
            300,0,energyMax);

        LKDrawingGroup *positionGroup = nullptr;
        LKDrawingGroup *ohmicGroup = nullptr;
        LKDrawingGroup *ohmicPositionGroup = nullptr;
        if (det->ringType == "16") {
            positionGroup = position16;
            ohmicGroup = ohmic16;
            ohmicPositionGroup = ohmicPosition16;
        }
        else if (det->ringType == "12dE") {
            positionGroup = position12DE;
            ohmicGroup = ohmic12DE;
            ohmicPositionGroup = ohmicPosition12DE;
        }
        else if (det->ringType == "12E") {
            positionGroup = position12E;
            ohmicGroup = ohmic12E;
            ohmicPositionGroup = ohmicPosition12E;
        }
        else if (det->ringType == "Q") {
            positionGroup = positionQQQ5;
            ohmicGroup = ohmicQQQ5;
            ohmicPositionGroup = ohmicPositionQQQ5;
        }

        TString drawingName = Form("idx%d_det%d",id,det->detNumber);
        if (positionGroup != nullptr) {
            auto drawing = positionGroup->CreateDrawing(drawingName);
            drawing->Add(histPosition[id],"colz");
        }
        if (ohmicGroup != nullptr) {
            auto drawing = ohmicGroup->CreateDrawing(drawingName);
            drawing->Add(histOhmic[id],"colz");
        }
        if (ohmicPositionGroup != nullptr) {
            auto drawing = ohmicPositionGroup->CreateDrawing(drawingName);
            drawing->Add(histOhmicPosition[id],"colz");
        }

        if (det->ringType == "12E") {
            const LKSiliconMapping::DetectorInfo *dEDet = nullptr;
            for (int jDet=0; jDet<mapping.GetNumDetectors(); ++jDet) {
                auto candidate = mapping.GetDetectorByVectorIndex(jDet);
                if (candidate != nullptr && candidate->ringType == "12dE" &&
                    candidate->phiNumber == det->phiNumber) {
                    dEDet = candidate;
                    break;
                }
            }
            if (dEDet != nullptr) {
                double coincidenceXMin = -1.2;
                double coincidenceXMax = 1.2;
                x6AxisRange(det,coincidenceXMin,coincidenceXMax);
                dEDetectorByE[id] = dEDet->detIndex;
                histTotalEnergy[id] = new TH2D(
                    Form("hist_sihit_total_energy_pair_%03d_%03d",dEDet->detIndex,id),
                    Form("12-ring coincidence dE det#%d (idx=%d) + E det#%d (idx=%d);E %s;dE+E",
                         dEDet->detNumber,dEDet->detIndex,det->detNumber,id,
                         x6AxisTitle().Data()),
                    200,coincidenceXMin,coincidenceXMax,600,0,2*energyMax);
                auto drawing = position12Sum->CreateDrawing(
                    Form("phi%d_dEidx%d_Eidx%d",det->phiNumber,dEDet->detIndex,id));
                drawing->Add(histTotalEnergy[id],"colz");

                histTotalVsDE[id] = new TH2D(
                    Form("hist_sihit_total_vs_de_pair_%03d_%03d",dEDet->detIndex,id),
                    Form("12-ring coincidence dE det#%d (idx=%d) + E det#%d (idx=%d);E+dE;dE",
                         dEDet->detNumber,dEDet->detIndex,det->detNumber,id),
                    600,0,2*energyMax,300,0,energyMax);
                auto pidDrawing = pid12TotalVsDE->CreateDrawing(
                    Form("phi%d_dEidx%d_Eidx%d",det->phiNumber,dEDet->detIndex,id));
                pidDrawing->Add(histTotalVsDE[id],"colz");

                histDEVsE[id] = new TH2D(
                    Form("hist_sihit_de_vs_e_pair_%03d_%03d",dEDet->detIndex,id),
                    Form("12-ring junction coincidence dE det#%d + E det#%d;E;dE",
                         dEDet->detNumber,det->detNumber),
                    300,0,energyMax,300,0,energyMax);
                auto deVsEDrawing = pid12DEVsE->CreateDrawing(
                    Form("phi%d_dEidx%d_Eidx%d",det->phiNumber,dEDet->detIndex,id));
                deVsEDrawing->Add(histDEVsE[id],"colz");

                histOhmicDEVsE[id] = new TH2D(
                    Form("hist_sihit_ohmic_de_vs_e_pair_%03d_%03d",dEDet->detIndex,id),
                    Form("12-ring ohmic coincidence dE det#%d + E det#%d;E_{ohmic};dE_{ohmic}",
                         dEDet->detNumber,det->detNumber),
                    300,0,energyMax,300,0,energyMax);
                auto ohmicDEVsEDrawing = pid12OhmicDEVsE->CreateDrawing(
                    Form("phi%d_dEidx%d_Eidx%d",det->phiNumber,dEDet->detIndex,id));
                ohmicDEVsEDrawing->Add(histOhmicDEVsE[id],"colz");

                histOhmicTotalVsPosition[id] = new TH2D(
                    Form("hist_sihit_ohmic_total_vs_position_pair_%03d_%03d",dEDet->detIndex,id),
                    Form("12-ring ohmic coincidence dE det#%d + E det#%d;E-detector %s;dE_{ohmic}+E_{ohmic}",
                         dEDet->detNumber,det->detNumber,x6AxisTitle().Data()),
                    200,coincidenceXMin,coincidenceXMax,600,0,2*energyMax);
                auto ohmicTotalVsPositionDrawing = ohmic12TotalVsPosition->CreateDrawing(
                    Form("phi%d_dEidx%d_Eidx%d",det->phiNumber,dEDet->detIndex,id));
                ohmicTotalVsPositionDrawing->Add(histOhmicTotalVsPosition[id],"colz");
            }
        }
    }

    TClonesArray *siHitArray = nullptr;
    tree->SetBranchAddress("SiHit",&siHitArray);
    Long64_t numEvents = tree->GetEntries();
    if (maxEvents >= 0 && maxEvents < numEvents)
        numEvents = maxEvents;

    // A paired SiHit contains both detectors: fDetID is E and fdEDetID is dE.
    // Fill both detector histograms so that the 12dE ring is not omitted.
    auto isQQQ5HighGainRun = [](int sourceRun) {
        switch (sourceRun) {
            case 134: case 135: case 136: case 139: case 141:
            case 142: case 146: case 147: case 148:
                return true;
            default:
                return false;
        }
    };
    int previousTreeNumber = -1;
    int sourceRun = -1;
    bool scaleQQQ5JunctionByFive = false;
    for (Long64_t event=0; event<numEvents; ++event) {
        tree->GetEntry(event);
        if (tree->GetTreeNumber() != previousTreeNumber) {
            previousTreeNumber = tree->GetTreeNumber();
            sourceRun = -1;
            auto currentFile = tree->GetFile();
            TString sourceName = currentFile != nullptr ? currentFile->GetName() : "";
            Ssiz_t runPosition = sourceName.Index("ko2520_");
            if (runPosition != kNPOS)
                sscanf(sourceName.Data()+runPosition,"ko2520_%d",&sourceRun);
            scaleQQQ5JunctionByFive = isQQQ5HighGainRun(sourceRun);
            cout << "source run " << sourceRun << ": QQQ5 AGET2 junction energy "
                 << (scaleQQQ5JunctionByFive ? "scaled by 1/5" : "unchanged")
                 << endl;
        }
        int multiplicity[kNumMult] = {0};
        double eventEnergy[maxDetectorIndex] = {0};
        double eventPosition[maxDetectorIndex];
        double eventOhmicEnergy[maxDetectorIndex] = {0};
        int eventJunctionStrip[maxDetectorIndex];
        for (int id=0; id<maxDetectorIndex; ++id) {
            eventPosition[id] = -999;
            eventJunctionStrip[id] = -1;
        }

        for (int iHit=0; iHit<siHitArray->GetEntriesFast(); ++iHit) {
            auto hit = (SKSiHit*) siHitArray->At(iHit);

            int detIDs[2] = {hit->GetDetID(), hit->GetdEDetID()};
            double energies[2] = {hit->GetE(), hit->GetdE()};
            double ohmicEnergies[2] = {hit->GetEnergyOhmic(), hit->GetdEOhmic()};
            double relativeZ[2] = {hit->GetRelativeZ(), hit->GetRelativeZdE()};

            for (int component=0; component<2; ++component) {
                int id = detIDs[component];
                double energy = energies[component];
                if (id < 0 || id >= maxDetectorIndex)
                    continue;

                auto det = mapping.FindDetectorByIndex(id);
                if (det == nullptr || isExcludedQQQ5Aget3(det))
                    continue;

                double ohmicEnergy = ohmicEnergies[component];

                bool isQQQ5OhmicOnly = det->detType == "QQQ5" &&
                    ohmicEnergy > 0 && TMath::Abs(energy-ohmicEnergy) < 1.e-12;
                if (component == 0 && det->detType == "QQQ5" &&
                    !isQQQ5OhmicOnly && scaleQQQ5JunctionByFive)
                    energy /= 5.;

                if (component == 0 && det->detType == "QQQ5" &&
                    !isQQQ5OhmicOnly && energy > 0) {
                    int junctionStrip = hit->GetJunctionStrip();
                    if (junctionStrip >= 1 && junctionStrip <= 32) {
                        histQQQ5JunctionEnergyVsStrip->Fill(junctionStrip,energy);
                    }
                }

                if (histPosition[id] == nullptr)
                    continue;

                if (ohmicEnergy > eventOhmicEnergy[id])
                    eventOhmicEnergy[id] = ohmicEnergy;

                // The E component keeps its calibrated endpoint energies.  A
                // matched dE component keeps the already reconstructed value
                // in fRelativeZdE because SKSiHit has no dE endpoint fields.
                if (det->detType == "X6") {
                    if (component == 0) {
                        if (hit->GetEnergyLeft() <= 0 || hit->GetEnergyRight() <= 0)
                            continue;
                        double energyLow = hit->GetEnergyLeft();
                        double energyHigh = hit->GetEnergyRight();
                        energy = energyLow+energyHigh;
                        relativeZ[component] =
                            (energyHigh-energyLow)/(energyHigh+energyLow);
                    }
                    else if (relativeZ[component] < -1 || relativeZ[component] > 1) {
                        continue;
                    }
                }
                // In old sihit files an ohmic-only QQQ5 channel is a separate
                // SiHit with identical main and ohmic energies.  A matched
                // file may legitimately have ohmic energy on a junction hit.
                else if (isQQQ5OhmicOnly) {
                    continue;
                }

                if (energy <= 0)
                    continue;

                ++multiplicity[kAll];
                if      (det->ringType == "16")   ++multiplicity[kRing16];
                else if (det->ringType == "12dE") ++multiplicity[kRing12DE];
                else if (det->ringType == "12E")  ++multiplicity[kRing12E];
                else if (det->ringType == "Q")    ++multiplicity[kQQQ5];

                double position = (det->detType == "QQQ5")
                                ? qqq5Coordinate(det,hit->GetJunctionStrip())
                                : x6Coordinate(det,relativeZ[component]);
                if (position < -900)
                    continue;
                histPosition[id]->Fill(position,energy);

                int ringSumIndex = -1;
                if      (det->ringType == "16")   ringSumIndex = kSum16;
                else if (det->ringType == "12dE") ringSumIndex = kSum12DE;
                else if (det->ringType == "12E")  ringSumIndex = kSum12E;
                else if (det->ringType == "Q")    ringSumIndex = kSumQQQ5;
                if (ringSumIndex >= 0) {
                    histRingSum[ringSumIndex]->Fill(position,energy);
                }

                if (energy > eventEnergy[id]) {
                    eventEnergy[id] = energy;
                    eventPosition[id] = position;
                    if (component == 0 && det->detType == "QQQ5")
                        eventJunctionStrip[id] = hit->GetJunctionStrip();
                }
                ++numFilled[id];
            }
        }

        for (int id=0; id<maxDetectorIndex; ++id) {
            auto det = mapping.FindDetectorByIndex(id);
            // QQQ5 contributes its junction spectrum directly.  It does not
            // require a coincident ohmic-side hit in the combined plot.
            if (mode == 2 && det != nullptr && det->ringType == "Q" &&
                histPosition[id] != nullptr && eventEnergy[id] > 0 &&
                eventPosition[id] > -900) {
                histCombinedOhmic->Fill(eventPosition[id],eventEnergy[id]);
                fillQValue(eventPosition[id],eventEnergy[id],true);
            }
            if (eventEnergy[id] > 0 && eventOhmicEnergy[id] > 0 &&
                eventPosition[id] > -900 && histOhmic[id] != nullptr) {
                histOhmic[id]->Fill(eventEnergy[id],eventOhmicEnergy[id]);
                histOhmicPosition[id]->Fill(eventPosition[id],eventOhmicEnergy[id]);
                if (det != nullptr && det->detType == "QQQ5" &&
                    eventJunctionStrip[id] >= 1 && eventJunctionStrip[id] <= 32)
                    histQQQ5OhmicEnergyVsJunctionStrip->Fill(
                        eventJunctionStrip[id],eventOhmicEnergy[id]);

                int ringSumIndex = -1;
                if (det != nullptr) {
                    if      (det->ringType == "16")   ringSumIndex = kSum16;
                    else if (det->ringType == "12dE") ringSumIndex = kSum12DE;
                    else if (det->ringType == "12E")  ringSumIndex = kSum12E;
                    else if (det->ringType == "Q")    ringSumIndex = kSumQQQ5;
                }
                if (ringSumIndex >= 0)
                    histOhmicRingSum[ringSumIndex]->Fill(
                        eventPosition[id],eventOhmicEnergy[id]);
                if (ringSumIndex == kSum16) {
                    histCombinedOhmic->Fill(
                        eventPosition[id],eventOhmicEnergy[id]);
                    fillQValue(eventPosition[id],eventOhmicEnergy[id]);
                }
            }
        }

        // This also accepts old files containing separate detector-local hits;
        // a newly matched hit fills the same event-level arrays above.
        for (int eID=0; eID<maxDetectorIndex; ++eID) {
            int dEID = dEDetectorByE[eID];
            if (histTotalEnergy[eID] == nullptr || dEID < 0)
                continue;
            if (eventEnergy[eID] > 0 && eventEnergy[dEID] > 0) {
                double totalEnergy = eventEnergy[dEID]+eventEnergy[eID];
                histTotalEnergy[eID]->Fill(
                    eventPosition[eID],totalEnergy);
                histTotalVsDE[eID]->Fill(totalEnergy,eventEnergy[dEID]);
                histDEVsE[eID]->Fill(eventEnergy[eID],eventEnergy[dEID]);
                histRingSumTotalEnergy->Fill(eventPosition[eID],totalEnergy);
                histRingSumDEVsE->Fill(eventEnergy[eID],eventEnergy[dEID]);
                histRingSumDEVsTotal->Fill(totalEnergy,eventEnergy[dEID]);
                ++numCoincidence[eID];
            }
            if (eventOhmicEnergy[eID] > 0 && eventOhmicEnergy[dEID] > 0) {
                histOhmicDEVsE[eID]->Fill(
                    eventOhmicEnergy[eID],eventOhmicEnergy[dEID]);
                histRingSumOhmicDEVsE->Fill(
                    eventOhmicEnergy[eID],eventOhmicEnergy[dEID]);
                if (eventPosition[eID] > -900) {
                    double ohmicTotal = eventOhmicEnergy[eID]+eventOhmicEnergy[dEID];
                    histOhmicTotalVsPosition[eID]->Fill(
                        eventPosition[eID],ohmicTotal);
                    histRingSumOhmicTotalVsPosition->Fill(
                        eventPosition[eID],ohmicTotal);
                    histCombinedOhmic->Fill(eventPosition[eID],ohmicTotal);
                    fillQValue(eventPosition[eID],ohmicTotal);
                }
                ++numOhmicCoincidence[eID];
            }
        }

        for (int i=0; i<kNumMult; ++i)
            histMultiplicity[i]->Fill(multiplicity[i]);
        if (event%10000 == 0)
            cout << "event " << event << " / " << numEvents << endl;
    }

    if (outputName.IsNull())
        outputName = Form("data_reco/draw_sihit_%04d.root",runID);
    cout << outputName << endl;
    auto output = new TFile(outputName,"recreate");

    // Keep detector-summed ring views, but omit all per-detector canvases.
    topRingSum->Draw();
    topRingSum->GetCanvas()->Write();

    topOhmicRingSum->Draw();
    topOhmicRingSum->GetCanvas()->Write();

    topMultiplicity->Draw();
    topMultiplicity->GetCanvas()->Write();

    auto canvasQQQ5EnergyVsStrip = new TCanvas(
        "canvas_sihit_qqq5_energy_vs_junction_strip",
        "QQQ5 energy vs junction strip",1800,750);
    canvasQQQ5EnergyVsStrip->Divide(2,1);
    canvasQQQ5EnergyVsStrip->cd(1);
    gPad->SetLogz();
    histQQQ5JunctionEnergyVsStrip->Draw("colz");
    canvasQQQ5EnergyVsStrip->cd(2);
    gPad->SetLogz();
    histQQQ5OhmicEnergyVsJunctionStrip->Draw("colz");
    canvasQQQ5EnergyVsStrip->Write();

    auto canvasCombinedOhmic = new TCanvas(
        "canvas_sihit_combined_ohmic_X6_12_16",
        "QQQ5 junction + X6 ohmic combined spectrum",1300,900);
    canvasCombinedOhmic->SetLogz();
    histCombinedOhmic->Draw("colz");
    auto drawRegionBoundary = [energyMax](double x, Color_t color) {
        auto line = new TLine(x,0,x,2*energyMax);
        line->SetLineColor(color);
        line->SetLineStyle(2);
        line->SetLineWidth(3);
        line->Draw();
    };
    drawRegionBoundary(pairXMin,kBlue+1);
    drawRegionBoundary(pairXMax,kBlue+1);
    drawRegionBoundary(ring16XMin,kRed+1);
    drawRegionBoundary(ring16XMax,kRed+1);
    if (mode == 2) {
        drawRegionBoundary(qqq5XMin,kGreen+2);
        drawRegionBoundary(qqq5XMax,kGreen+2);
    }
    TGraph *graphKinematicsCutUpper = nullptr;
    TGraph *graphKinematicsCutLower = nullptr;
    TGraph *graphKinematicsQQQ5CutUpper = nullptr;
    TGraph *graphKinematicsQQQ5CutLower = nullptr;
    if (graphKinematics != nullptr) {
        graphKinematicsCutUpper = dynamic_cast<TGraph*>(
            graphKinematics->Clone("graph_kinematics_cut_upper"));
        graphKinematicsCutLower = dynamic_cast<TGraph*>(
            graphKinematics->Clone("graph_kinematics_cut_lower"));
        for (int i=0; i<graphKinematics->GetN(); ++i) {
            double theta = 0;
            double energy = 0;
            graphKinematics->GetPoint(i,theta,energy);
            graphKinematicsCutUpper->SetPoint(i,theta,energy+kinematicsCutMeV);
            graphKinematicsCutLower->SetPoint(i,theta,energy-kinematicsCutMeV);
        }
        for (auto cutGraph : {graphKinematicsCutUpper,graphKinematicsCutLower}) {
            cutGraph->SetLineColorAlpha(kMagenta+1,0.35);
            cutGraph->SetLineStyle(3);
            cutGraph->SetLineWidth(1);
            cutGraph->Draw("L same");
        }
        graphKinematics->SetLineColorAlpha(kMagenta+1,0.50);
        graphKinematics->SetLineStyle(2);
        graphKinematics->SetLineWidth(2);
        graphKinematics->Draw("L same");
    }
    if (graphKinematicsQQQ5 != nullptr) {
        graphKinematicsQQQ5CutUpper = dynamic_cast<TGraph*>(
            graphKinematicsQQQ5->Clone("graph_kinematics_qqq5_1mm_cut_upper"));
        graphKinematicsQQQ5CutLower = dynamic_cast<TGraph*>(
            graphKinematicsQQQ5->Clone("graph_kinematics_qqq5_1mm_cut_lower"));
        for (int i=0; i<graphKinematicsQQQ5->GetN(); ++i) {
            double theta = 0;
            double energyLoss = 0;
            graphKinematicsQQQ5->GetPoint(i,theta,energyLoss);
            graphKinematicsQQQ5CutUpper->SetPoint(
                i,theta,energyLoss+kinematicsCutMeV);
            graphKinematicsQQQ5CutLower->SetPoint(
                i,theta,energyLoss-kinematicsCutMeV);
        }
        for (auto cutGraph : {
                 graphKinematicsQQQ5CutUpper,graphKinematicsQQQ5CutLower}) {
            cutGraph->SetLineColorAlpha(kGreen+2,0.45);
            cutGraph->SetLineStyle(3);
            cutGraph->SetLineWidth(1);
            cutGraph->Draw("L same");
        }
        graphKinematicsQQQ5->SetLineColor(kGreen+2);
        graphKinematicsQQQ5->SetLineStyle(1);
        graphKinematicsQQQ5->SetLineWidth(3);
        graphKinematicsQQQ5->Draw("L same");
        output->cd();
        graphKinematicsQQQ5->Write();
    }
    output->cd();
    auto label12Ring = new TLatex(0.5*(pairXMin+pairXMax),1.85*energyMax,"12-ring");
    label12Ring->SetTextAlign(22);
    label12Ring->SetTextSize(0.045);
    label12Ring->SetTextColor(kBlue+1);
    label12Ring->Draw();
    auto label16Ring = new TLatex(0.5*(ring16XMin+ring16XMax),1.85*energyMax,"16-ring");
    label16Ring->SetTextAlign(22);
    label16Ring->SetTextSize(0.045);
    label16Ring->SetTextColor(kRed+1);
    label16Ring->Draw();
    if (mode == 2) {
        auto labelQQQ5 = new TLatex(
            0.5*(qqq5XMin+qqq5XMax),1.85*energyMax,"QQQ5");
        labelQQQ5->SetTextAlign(22);
        labelQQQ5->SetTextSize(0.045);
        labelQQQ5->SetTextColor(kGreen+2);
        labelQQQ5->Draw();
    }
    if (graphKinematics != nullptr) {
        auto legendKinematics = new TLegend(0.54,0.67,0.89,0.85);
        legendKinematics->SetBorderSize(0);
        legendKinematics->SetFillStyle(0);
        legendKinematics->AddEntry(
            graphKinematics,"X6: full proton energy","l");
        if (graphKinematicsQQQ5 != nullptr)
            legendKinematics->AddEntry(
                graphKinematicsQQQ5,"QQQ5: #DeltaE in 1 mm Si (PSTAR)","l");
        legendKinematics->AddEntry(
            graphKinematicsCutUpper,"#pm1 MeV gate","l");
        legendKinematics->Draw();
    }
    canvasCombinedOhmic->Write();

    if (mode == 2) {
        auto canvasQValue = new TCanvas(
            "canvas_sihit_qvalue_X6_12_16",
            "QQQ5 + X6 kinematics-gated Q-value",2100,700);
        canvasQValue->Divide(3,1);

        canvasQValue->cd(1);
        gPad->SetLogz();
        histQValueVsThetaCM->Draw("colz");
        auto zeroLine2D = new TLine(0,0,180,0);
        zeroLine2D->SetLineColor(kRed+1);
        zeroLine2D->SetLineStyle(2);
        zeroLine2D->Draw();

        canvasQValue->cd(2);
        histQValue->SetLineColor(kBlue+1);
        histQValue->SetLineWidth(2);
        histQValue->Draw("hist");
        auto zeroLine1D = new TLine(0,0,0,1.05*histQValue->GetMaximum());
        zeroLine1D->SetLineColor(kRed+1);
        zeroLine1D->SetLineStyle(2);
        zeroLine1D->Draw();

        canvasQValue->cd(3);
        histThetaCM->SetLineColor(kGreen+2);
        histThetaCM->SetLineWidth(2);
        histThetaCM->Draw("hist");

        canvasQValue->Write();
        cout << "Q-value candidates within +/-" << kinematicsCutMeV
             << " MeV: " << numQValueCandidates << endl;

        auto canvasQValueAll = new TCanvas(
            "canvas_sihit_qvalue_all_X6_12_16",
            "QQQ5 + X6 Q-value without kinematics cut",2100,700);
        canvasQValueAll->Divide(3,1);

        canvasQValueAll->cd(1);
        gPad->SetLogz();
        histQValueVsThetaCMAll->Draw("colz");
        auto zeroLine2DAll = new TLine(0,0,180,0);
        zeroLine2DAll->SetLineColor(kRed+1);
        zeroLine2DAll->SetLineStyle(2);
        zeroLine2DAll->Draw();

        canvasQValueAll->cd(2);
        histQValueAll->SetLineColor(kBlue+1);
        histQValueAll->SetLineWidth(2);
        histQValueAll->Draw("hist");
        auto zeroLine1DAll = new TLine(
            0,0,0,1.05*histQValueAll->GetMaximum());
        zeroLine1DAll->SetLineColor(kRed+1);
        zeroLine1DAll->SetLineStyle(2);
        zeroLine1DAll->Draw();

        canvasQValueAll->cd(3);
        histThetaCMAll->SetLineColor(kGreen+2);
        histThetaCMAll->SetLineWidth(2);
        histThetaCMAll->Draw("hist");

        canvasQValueAll->Write();
        cout << "All Q-value events without kinematics cut: "
             << numAllQValueEvents << endl;
    }
    output->Close();

    cout << endl << "SiHit entries by detector" << endl;
    for (int iDet=0; iDet<mapping.GetNumDetectors(); ++iDet) {
        auto det = mapping.GetDetectorByVectorIndex(iDet);
        if (det != nullptr && !isExcludedQQQ5Aget3(det) &&
            det->detIndex >= 0 && det->detIndex < maxDetectorIndex)
            printf("idx=%3d det#=%4d %-5s %-4s entries=%lld\n",
                   det->detIndex,det->detNumber,det->detType.Data(),
                   det->ringType.Data(),numFilled[det->detIndex]);
    }

    cout << endl << "12-ring dE-E coincidences" << endl;
    for (int eID=0; eID<maxDetectorIndex; ++eID) {
        if (histTotalEnergy[eID] == nullptr)
            continue;
        auto eDet = mapping.FindDetectorByIndex(eID);
        auto dEDet = mapping.FindDetectorByIndex(dEDetectorByE[eID]);
        if (eDet != nullptr && dEDet != nullptr)
            printf("phi=%2d dE idx=%3d det#=%4d + E idx=%3d det#=%4d junction=%lld ohmic=%lld\n",
                   eDet->phiNumber,dEDet->detIndex,dEDet->detNumber,
                   eDet->detIndex,eDet->detNumber,
                   numCoincidence[eID],numOhmicCoincidence[eID]);
    }
}

// Keep filename-first calls used by older scripts in relative-position mode.
void draw_sihit(TString inputName, TString outputName="")
{
    draw_sihit(0,inputName,outputName);
}
