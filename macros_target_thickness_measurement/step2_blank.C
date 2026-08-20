/**
 * step2_blank.C
 *
 * Reads blank SPE files, fits a Gaussian, and computes the energy
 * calibration (channel → MeV) using the pulser p0 from step1.
 *
 * Calibration formula:
 *   E [MeV] = (channel - p0) * conversion
 *   conversion = kAlphaEnergy / (mean_blank - p0)
 *
 * Output:
 *   data_analysis/energy_calibration.txt
 *     columns: date  p0  conversion  conv_err  sigma_ch  sigma_mev
 *
 * Run: root -l -b -q step2_blank.C
 *   (requires data_analysis/pulser_calibration.txt from step1)
 */

#include "analysis_utils.h"
#include "GaussianFitter.h"
#include "GraphGaussianFitter.h"

#include <map>
#include "TLegend.h"
#include "LKDrawing.h"
#include "LKDrawingGroup.h"

TLegend *FindLegend(LKDrawing *drawing)
{
  if (drawing == nullptr) return nullptr;
  for (int i = 0; i < drawing->GetEntries(); ++i) {
    auto obj = drawing->At(i);
    if (obj != nullptr && obj->InheritsFrom(TLegend::Class()))
      return (TLegend *) obj;
  }
  return nullptr;
}

// ── Load pulser p0 per date ───────────────────────────────────────────────────
std::map<std::string, double> LoadPulserP0(const char *calFile)
{
  std::map<std::string, double> p0map;
  std::ifstream f(calFile);
  if (!f.is_open()) {
    std::cerr << "Cannot open " << calFile
              << " — run step1_pulser.C first." << std::endl;
    return p0map;
  }
  std::string line;
  while (std::getline(f, line)) {
    line = Trim(line);
    if (line.empty() || line[0] == '#') continue;
    std::istringstream ss(line);
    std::string date;
    double p0, p0Err, slope, slopeErr;
    if (ss >> date >> p0 >> p0Err >> slope >> slopeErr)
      p0map[date] = p0;
  }
  return p0map;
}

void step2_blank(const char *listFile = "list.txt",
                 bool        useGraphFit = true)
{
  const auto p0map = LoadPulserP0("data_analysis/pulser_calibration.txt");
  if (p0map.empty()) return;

  const auto entries = ReadSpeList(listFile);

  auto top = new LKDrawingGroup("step2_top");
  auto grp = top->CreateGroup("blank_channel");

  std::ofstream outFile("data_analysis/energy_calibration.txt");
  outFile << "# fit_mode=" << (useGraphFit ? "graph" : "histogram") << "\n";
  outFile << "# date  p0  conversion_MeV_per_ch  conv_err  sigma_ch  sigma_mev\n";
  outFile << "# E[MeV] = (channel - p0) * conversion\n";

  for (const auto &e : entries) {
    if (e.type != SpeType::kBlank) continue;

    const std::string date = DateFromPath(e.path);

    // Find p0 for this date; fall back to first available
    double p0 = 0;
    if (p0map.count(date)) {
      p0 = p0map.at(date);
    } else if (!p0map.empty()) {
      p0 = p0map.begin()->second;
      std::cerr << "No pulser p0 for date " << date
                << ", using " << p0map.begin()->first
                << " p0=" << p0 << std::endl;
    }

    auto hist = ReadSpeHistogram(e);
    if (!hist) continue;
    hist->SetLineWidth(2);

    if (hist->GetStdDev() <= 0) { delete hist; continue; }
    SetFiveSigmaRange(hist);

    auto fitResult = useGraphFit
        ? ProceedGraphGaussianFittingFull(
            hist,
            true,   // graph
            false,  // initial Gaussian
            false,  // fit range lines
            true,   // final Gaussian
            true,   // legend
            "AP",
            2.0,    // fixed sd_range
            10,
            2.0)
        : ProceedGaussianFittingFull(
            hist,
            true,   // histogram
            false,  // initial Gaussian
            false,  // fit range lines
            true,   // final Gaussian
            true,   // legend
            "hist",
            2.0,    // fixed sd_range
            10,
            2.0);
    auto drawing = fitResult.drawing;
    auto f_result = fitResult.f_result;
    const double mu = fitResult.mu;
    const double sigmaCh = fitResult.sigma;

    const double muErr   = f_result->GetParError(1);
    const double corrCh  = mu - p0;
    const double conv    = kAlphaEnergy / corrCh;
    const double convErr = kAlphaEnergy * muErr / (corrCh * corrCh);
    const double sigmaMev = sigmaCh * conv;

    outFile << date << "\t"
            << p0 << "\t"
            << conv << "\t" << convErr << "\t"
            << sigmaCh << "\t" << sigmaMev << "\n";

    std::cout << "[" << date << "] " << e.path << "\n"
              << "  mu=" << mu << " ch  p0=" << p0 << " ch"
              << "  conv=" << conv << " MeV/ch"
              << "  sigma=" << sigmaCh << " ch (" << sigmaMev << " MeV)"
              << std::endl;

    // Add calibration info to the Gaussian fitter legend.
    auto legend = FindLegend(drawing);
    char line[256];
    if (legend != nullptr) {
      legend->SetY1NDC(0.46);
      legend->SetTextSize(0.030);
      std::snprintf(line, sizeof(line), "p_{0} = %.4g ch", p0);
      legend->AddEntry((TObject*)nullptr, line, "");
      std::snprintf(line, sizeof(line), "conv = %.6f MeV/ch", conv);
      legend->AddEntry((TObject*)nullptr, line, "");
      std::snprintf(line, sizeof(line), "#sigma = %.2f ch = %.4f MeV", sigmaCh, sigmaMev);
      legend->AddEntry((TObject*)nullptr, line, "");
    }

    grp->AddDrawing(drawing);
  }

  outFile.close();
  std::cout << "Saved: data_analysis/energy_calibration.txt" << std::endl;
  top->Draw();
  top->Save();
}
