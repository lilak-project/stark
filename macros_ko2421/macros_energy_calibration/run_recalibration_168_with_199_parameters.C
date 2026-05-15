#include "make_c2_histograms.C"
#include "run_junction_energy_calibration.C"

void run_recalibration_168_with_199_parameters()
{
    TString fileName, c3Name;

    fileName = make_c2_histograms(false, 2, 168, 199, true); e_test << fileName << endl;

    fileName = "data/compare_168_with_199.hist_c2.root";
    c3Name = run_junction_energy_calibration(true, 168, 199, fileName); e_test << c3Name << endl;

    c3Name = "data/compare_168_with_199.cp.dat";
    make_c2_histograms(true, 3, 168, 199, true, c3Name);
}
