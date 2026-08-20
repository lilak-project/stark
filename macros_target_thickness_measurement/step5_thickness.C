/**
 * step5_thickness.C
 *
 * Converts remaining-energy spectra into target-thickness histograms using
 * the analytic CD2 conversion (ThicknessConverter.h), mirroring the
 * energy → thickness block of draw_all_spe.C.
 *
 * Self-contained: rebuilds energy histograms from list.txt +
 * energy_calibration.txt (does not depend on step3_energy.root structure).
 *
 * Fitting strategy:
 *   blank / ch2 / carbon  → Gaussian fit on thickness
 *   cd2                    → no fit, histogram only
 *
 * Input:
 *   list.txt
 *   data_analysis/energy_calibration.txt   (from step2)
 *
 * Output:
 *   data_analysis/step5_thickness.root
 *
 * Run: root -l -b -q step5_thickness.C
 */

#include "analysis_utils.h"
#include "GaussianFitter.h"
#include "GraphGaussianFitter.h"
#include "ThicknessConverter.h"

#include <map>
#include <limits>
#include <vector>
#include "TH2D.h"
#include "TGraphErrors.h"
#include "LKDrawing.h"
#include "LKDrawingGroup.h"

// ── Calibration record (same as step3) ────────────────────────────────────────
struct EnergyCal {
  double p0 = 0, conv = 0, convErr = 0, sigmaCh = 0, sigmaMev = 0;
};

std::map<std::string, EnergyCal> LoadEnergyCal(const char *calFile)
{
  std::map<std::string, EnergyCal> cal;
  std::ifstream f(calFile);
  if (!f.is_open()) {
    std::cerr << "Cannot open " << calFile
              << " — run step2_blank.C first." << std::endl;
    return cal;
  }
  std::string line;
  while (std::getline(f, line)) {
    line = Trim(line);
    if (line.empty() || line[0] == '#') continue;
    std::istringstream ss(line);
    std::string date;
    EnergyCal c;
    if (ss >> date >> c.p0 >> c.conv >> c.convErr >> c.sigmaCh >> c.sigmaMev)
      cal[date] = c;
  }
  return cal;
}

std::string CalDate(const std::string &fileDate,
                    const std::map<std::string, EnergyCal> &cal)
{
  if (cal.count(fileDate)) return fileDate;
  if (fileDate == "20260529" && cal.count("20260528")) return "20260528";
  if (!cal.empty()) return cal.rbegin()->first;
  return "";
}

// channel → energy histogram
TH1D *MakeEnergyHist(TH1D *hCh, const EnergyCal &c,
                     const std::string &path, int nBins)
{
  const double eMin = (0.0                  - c.p0) * c.conv;
  const double eMax = ((double)kNChannels   - c.p0) * c.conv;
  auto hE = new TH1D(("he_" + HistNameFromPath(path)).c_str(),
                     (path + ";Remaining Energy (MeV);Counts").c_str(),
                     nBins, eMin, eMax);
  hE->SetDirectory(nullptr);
  for (int i = 1; i <= kNChannels; ++i) {
    const double counts = hCh->GetBinContent(i);
    if (counts <= 0) continue;
    const double channel = hCh->GetBinCenter(i);
    const double energy = (channel - c.p0) * c.conv;
    const int bin = hE->FindBin(energy);
    if (bin >= 1 && bin <= hE->GetNbinsX()) {
      hE->Fill(energy, counts);
    }
  }
  return hE;
}

void AddDiagnosticGraph(LKDrawingGroup *grp, TGraph *graph,
                        const char *name, const char *yTitle)
{
  if (grp == nullptr || graph == nullptr || graph->GetN() <= 0)
    return;

  double x = 0, y = 0;
  double xMin = 1.e99, xMax = -1.e99;
  double yMin = 1.e99, yMax = -1.e99;
  for (int i = 0; i < graph->GetN(); ++i) {
    graph->GetPoint(i, x, y);
    xMin = std::min(xMin, x);
    xMax = std::max(xMax, x);
    yMin = std::min(yMin, y);
    yMax = std::max(yMax, y);
  }
  if (xMax <= xMin) { xMin -= 0.5; xMax += 0.5; }
  if (yMax <= yMin) { yMin -= 0.5; yMax += 0.5; }
  const double yPad = 0.12 * (yMax - yMin);

  auto frame = new TH2D((std::string("frame_") + name).c_str(),
                        (std::string(name) + ";Iteration;" + yTitle).c_str(),
                        100, xMin - 0.2, xMax + 0.2,
                        100, yMin - yPad, yMax + yPad);
  frame->SetDirectory(nullptr);
  frame->SetStats(0);

  auto draw = grp->CreateDrawing();
  draw->Add(frame, "");
  draw->Add(graph, "pl same");
}

void AddFitDiagnostics(LKDrawingGroup *grp,
                       const GaussianFitResult &fitResult,
                       const std::string &tag)
{
  AddDiagnosticGraph(grp, fitResult.g_amp,
                     (tag + "_gaus_amplitude").c_str(), "Amplitude");
  AddDiagnosticGraph(grp, fitResult.g_mean,
                     (tag + "_gaus_mean").c_str(), "Mean");
  AddDiagnosticGraph(grp, fitResult.g_sigma,
                     (tag + "_gaus_sigma").c_str(), "Sigma");
  AddDiagnosticGraph(grp, fitResult.g_chi2ndf,
                     (tag + "_gaus_chi2ndf").c_str(), "#chi^{2}/NDF");
  AddDiagnosticGraph(grp, fitResult.g_sd_range,
                     (tag + "_gaus_sd_range").c_str(), "sd_range");
}

ThicknessConv::Material TargetMaterial(SpeType type)
{
  switch (type) {
    case SpeType::kCd2:    return ThicknessConv::Material::kCD2;
    case SpeType::kCh2:    return ThicknessConv::Material::kCH2;
    case SpeType::kCarbon: return ThicknessConv::Material::kCarbon;
    default:               return ThicknessConv::Material::kCD2;
  }
}

void MakeThicknessPoints(const TH1D *hE,
                         ThicknessConv::Material material,
                         std::vector<double> &xs,
                         std::vector<double> &ys)
{
  xs.clear();
  ys.clear();
  for (int i = 1; i <= hE->GetNbinsX(); ++i) {
    const double counts = hE->GetBinContent(i);
    if (counts <= 0) continue;
    const double energy = hE->GetBinCenter(i);
    if (energy <= 0) continue;
    const double thickness =
        EnergyToThickness(material, std::min(energy, ThicknessConv::kE0));
    xs.push_back(thickness);
    ys.push_back(counts);
  }
}

TGraphErrors *MakeThicknessCountsGraph(const std::vector<double> &xs,
                                       const std::vector<double> &ys,
                                       const char *name,
                                       const char *title)
{
  auto graph = new TGraphErrors((int)xs.size());
  graph->SetName(name);
  graph->SetTitle(TString::Format("%s;Target Thickness (#mum);Counts", title));
  graph->SetMarkerStyle(20);
  graph->SetMarkerSize(0.45);
  graph->SetMarkerColor(kBlack);
  graph->SetLineColor(kBlack);
  for (int i = 0; i < (int)xs.size(); ++i) {
    graph->SetPoint(i, xs[i], ys[i]);
    graph->SetPointError(i, 0, std::sqrt(std::max(ys[i], 1.0)));
  }
  return graph;
}

void step5_thickness(const char *listFile      = "list.txt",
                     int         thicknessBinFactor = 1,
                     int         thicknessBinsBase  = 120,
                     double      sdRangeMin = 1.5,
                     int         minPointsInRange = 10,
                     bool        showFitIterationInfo = false,
                     bool        useGraphFit = true)
{
  const int kEnergyBins = kNChannels;
  const int thicknessBins = thicknessBinsBase * std::max(thicknessBinFactor, 1);

  const auto cal = LoadEnergyCal("data_analysis/energy_calibration.txt");
  if (cal.empty()) return;

  const auto entries = ReadSpeList(listFile);

  auto top    = new LKDrawingGroup("step5_top");
  auto grpRef = top->CreateGroup("blank_ch2_carbon");
  auto grpCd2 = top->CreateGroup("cd2");
  auto grpDiag = top->CreateGroup("fit_iteration");

  for (const auto &e : entries) {
    if (e.type == SpeType::kPulser) continue;

    const std::string fileDate = DateFromPath(e.path);
    const std::string calDate  = CalDate(fileDate, cal);
    if (calDate.empty()) {
      std::cerr << "No calibration for: " << e.path << std::endl;
      continue;
    }
    const EnergyCal &c = cal.at(calDate);

    auto hCh = ReadSpeHistogram(e);
    if (!hCh) continue;
    auto hE = MakeEnergyHist(hCh, c, e.path, kEnergyBins);
    delete hCh;

    // energy → thickness
    auto hT = MakeThicknessHist(hE, "ht_" + HistNameFromPath(e.path),
                                thicknessBins,
                                ThicknessConv::kThicknessRange,
                                TargetMaterial(e.type));
    std::vector<double> thicknessX;
    std::vector<double> thicknessY;
    MakeThicknessPoints(hE, TargetMaterial(e.type), thicknessX, thicknessY);
    delete hE;

    if (e.type == SpeType::kBlank ||
        e.type == SpeType::kCh2   ||
        e.type == SpeType::kCarbon)
    {
      const double sdRangeMax = (e.type == SpeType::kBlank) ? 1.5 : 3.0;
      const double initialMeanOverride =
          (e.type == SpeType::kBlank)
              ? 0.0 : std::numeric_limits<double>::quiet_NaN();
      const std::string graphName = "gt_" + HistNameFromPath(e.path);
      auto fitResult = useGraphFit
          ? ProceedGraphGaussianFittingFull(
              graphName.c_str(),
              e.path.c_str(),
              "Target Thickness (#mum)",
              "Counts",
              thicknessX,
              thicknessY,
              true,   // graph
              true,   // initial Gaussian
              true,   // fit range lines
              true,   // final Gaussian
              true,   // fit parameter legend
              "AP",
              sdRangeMin,
              minPointsInRange,
              sdRangeMax,
              showFitIterationInfo,
              initialMeanOverride)
          : ProceedGaussianFittingFull(
              hT,
              true,   // histogram
              true,   // initial Gaussian
              true,   // fit range lines
              true,   // final Gaussian
              true,   // fit parameter legend
              "hist",
              sdRangeMin,
              minPointsInRange,
              sdRangeMax,
              showFitIterationInfo,
              initialMeanOverride);
      fitResult.drawing->SetLegendCorner(1);
      grpRef->AddDrawing(fitResult.drawing);

      if (e.type == SpeType::kBlank || e.type == SpeType::kCh2)
        AddFitDiagnostics(grpDiag, fitResult,
                          SpeTypeName(e.type) + "_" + fileDate);
    }
    else {
      auto draw = grpCd2->CreateDrawing();
      if (useGraphFit) {
        const std::string graphName = "gt_" + HistNameFromPath(e.path);
        auto graph = MakeThicknessCountsGraph(thicknessX, thicknessY,
                                              graphName.c_str(), e.path.c_str());
        graph->GetXaxis()->SetLimits(0, ThicknessConv::kThicknessRange);
        draw->Add(graph, "AP");
      }
      else {
        hT->GetXaxis()->SetRangeUser(0, ThicknessConv::kThicknessRange);
        draw->Add(hT, "hist");
      }
    }
  }

  top->WriteFile("data_analysis/step5_thickness.root");
  std::cout << "Saved: data_analysis/step5_thickness.root" << std::endl;
  top->Draw();
}
