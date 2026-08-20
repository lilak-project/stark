/**
 * step1_pulser.C
 *
 * Reads pulser SPE files, fits a Gaussian to each peak,
 * then fits pol1 (mean channel vs Vamp) per date.
 *
 * Channel cuts applied in ReadSpeHistogram (analysis_utils.h):
 *   - 1V and 1.6V amplitude: skip channel < 50
 *
 * Output:
 *   data_analysis/pulser_calibration.txt
 *     columns: date  p0  p0_err  slope  slope_err
 *
 * Run: root -l -b -q step1_pulser.C
 */

#include "analysis_utils.h"
#include "GaussianFitter.h"
#include "GraphGaussianFitter.h"

#include <map>
#include "TF1.h"
#include "TGraphErrors.h"
#include "TH2D.h"
#include "TLegend.h"
#include "LKDrawing.h"
#include "LKDrawingGroup.h"

void step1_pulser(const char *listFile  = "list.txt",
                  bool        exclude1V = true,
                  bool        useGraphFit = true)
{
  const auto entries = ReadSpeList(listFile);

  // ── Collect pulser histograms grouped by date ─────────────────────────────
  struct PulserPoint { double vamp, mean, sigma; };
  std::map<std::string, std::vector<PulserPoint>> byDate;

  auto top    = new LKDrawingGroup("step1_top");
  auto grpSpe = top->CreateGroup("pulser_spectra");

  for (const auto &e : entries) {
    if (e.type != SpeType::kPulser) continue;

    double vamp = 0;
    if (!GetPulserVamp(e.path, vamp)) {
      std::cerr << "Cannot parse Vamp from: " << e.path << std::endl;
      continue;
    }

    if (exclude1V && vamp > 0.9 && vamp < 1.1) continue;

    auto hist = ReadSpeHistogram(e);  // low-ch cut for 1V/1.6V already applied
    if (!hist) continue;
    hist->SetLineWidth(1);

    if (hist->GetStdDev() <= 0) {
      std::cerr << "Zero sigma, skipping: " << e.path << std::endl;
      delete hist; continue;
    }

    SetFiveSigmaRange(hist);

    auto fitResult = useGraphFit
        ? ProceedGraphGaussianFittingFull(hist)
        : ProceedGaussianFittingFull(hist);
    auto drawing = fitResult.drawing;
    const double mu = fitResult.mu;
    const double sg = fitResult.sigma;

    const std::string date = DateFromPath(e.path);
    byDate[date].push_back({vamp, mu, sg});

    char title[256];
    std::snprintf(title, sizeof(title),
                  "%s  [Vamp=%.3gV  #mu=%.1f  #sigma=%.1f]",
                  e.path.c_str(), vamp, mu, sg);
    hist->SetTitle(title);
    grpSpe->AddDrawing(drawing);
  }

  // ── Helper: fit pol1 to a subset of points ───────────────────────────────
  auto FitPol1 = [](const std::vector<PulserPoint> &pts,
                    const std::string &name,
                    double vMax,
                    int markerStyle, int color)
      -> std::pair<TGraphErrors*, TF1*>
  {
    auto g = new TGraphErrors(static_cast<int>(pts.size()));
    g->SetName(name.c_str());
    g->SetMarkerStyle(markerStyle);
    g->SetMarkerSize(1.0);
    g->SetLineWidth(2);
    g->SetMarkerColor(color);
    g->SetLineColor(color);
    for (int i = 0; i < (int)pts.size(); ++i) {
      g->SetPoint(i, pts[i].vamp, pts[i].mean);
      g->SetPointError(i, 0, pts[i].sigma);
    }
    auto f = new TF1(("f_" + name).c_str(), "pol1", 0, vMax);
    f->SetNpx(500);
    f->SetLineWidth(2);
    f->SetLineColor(color);
    g->Fit(f, "RQ");
    return {g, f};
  };

  // ── Per-date linear fit: mean_channel = p0 + slope * Vamp ────────────────
  auto grpCal = top->CreateGroup("pulser_calibration");

  std::ofstream outFile("data_analysis/pulser_calibration.txt");
  outFile << "# exclude1V=" << (exclude1V ? "true" : "false") << "\n";
  outFile << "# fit_mode=" << (useGraphFit ? "graph" : "histogram") << "\n";
  outFile << "# date  p0  p0_err  slope  slope_err\n";
  outFile << "# mean_channel = p0 + slope * Vamp\n";

  for (auto &[date, points] : byDate) {
    if (points.size() < 2) {
      std::cerr << "Not enough pulser points for date " << date << std::endl;
      continue;
    }
    std::sort(points.begin(), points.end(),
              [](const PulserPoint &a, const PulserPoint &b){
                return a.vamp < b.vamp; });

    const double vMax = points.back().vamp * 1.15;

    // Fit with whatever points are in byDate (1.6V already excluded if option is on)
    auto [graphAll, fitAll] = FitPol1(points, "g_all_" + date, vMax,
                                      20, kBlack);
    TF1 *fitSaved = fitAll;

    const double p0       = fitSaved->GetParameter(0);
    const double p0Err    = fitSaved->GetParError(0);
    const double slope    = fitSaved->GetParameter(1);
    const double slopeErr = fitSaved->GetParError(1);
    const double chi2ndf  = fitSaved->GetNDF() > 0
                          ? fitSaved->GetChisquare() / fitSaved->GetNDF() : 0;

    outFile << date << "\t"
            << p0 << "\t" << p0Err << "\t"
            << slope << "\t" << slopeErr << "\n";

    std::cout << "[" << date << "]"
              << (exclude1V ? " (no 1.0V)" : " (all)")
              << "  p0=" << p0 << "±" << p0Err
              << "  slope=" << slope << "±" << slopeErr
              << "  chi2/ndf=" << chi2ndf << std::endl;

    // Frame
    auto frame = new TH2D(("hframe_" + date).c_str(),
                          ("Pulser " + date + ";V_{amp} (V);Mean channel").c_str(),
                          100, 0, vMax,
                          100, 0, points.back().mean * 1.2);
    frame->SetDirectory(nullptr);
    frame->SetStats(0);

    // Legend
    auto legend = new TLegend(0.12, 0.58, 0.62, 0.88);
    legend->SetBorderSize(0);
    legend->SetFillStyle(0);
    legend->SetTextSize(0.032);
    char line[256];

    const char *pointsLabel = exclude1V ? "Points (1.0V excluded)" : "All points";
    legend->AddEntry(graphAll, pointsLabel, "p");
    std::snprintf(line, sizeof(line),
                  "p_{0}=%.4g#pm%.4g  slope=%.4g#pm%.4g  #chi^{2}/NDF=%.2g",
                  fitAll->GetParameter(0), fitAll->GetParError(0),
                  fitAll->GetParameter(1), fitAll->GetParError(1),
                  fitAll->GetNDF() > 0
                      ? fitAll->GetChisquare() / fitAll->GetNDF() : 0);
    legend->AddEntry(fitAll, line, "l");

    auto draw = grpCal->CreateDrawing();
    draw->Add(frame, "");
    draw->Add(graphAll, "P same");
    draw->Add(fitAll, "samel");
    draw->Add(legend, "same");
  }

  outFile.close();
  std::cout << "Saved: data_analysis/pulser_calibration.txt" << std::endl;
  top->Draw();
  top -> Save();
}
