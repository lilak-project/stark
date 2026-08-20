/**
 * step3_energy.C
 *
 * Converts channel → remaining energy for all data in list.txt.
 * Calibration loaded from data_analysis/energy_calibration.txt (step2).
 *
 * Fitting strategy:
 *   blank  → Gaussian
 *   ch2    → Gaussian  + legend shows sigma_pure = sqrt(σ² - σ_blank²)
 *   carbon → Gaussian  + legend shows sigma_pure
 *   cd2    → no fit, histogram only
 *
 * Date mapping:
 *   20260527 → 20260527 calibration
 *   20260528 → 20260528 calibration
 *   20260529 → 20260528 calibration  (no blank that day)
 *
 * Run: root -l -b -q step3_energy.C
 *   (requires data_analysis/energy_calibration.txt from step2)
 */

#include "analysis_utils.h"
#include "GaussianFitter.h"
#include "GraphGaussianFitter.h"
#include "ThicknessConverter.h"

#include <map>
#include "TLegend.h"
#include "TGraph.h"
#include "TH2D.h"
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

struct ThicknessPoint {
  double energy = 0;
  double thickness = 0;
  std::string label;
};

double TargetEnergyToThickness(SpeType type, double meanEnergy)
{
  switch (type) {
    case SpeType::kCd2:    return EnergyToThickness(ThicknessConv::Material::kCD2, meanEnergy);
    case SpeType::kCh2:    return EnergyToThickness(ThicknessConv::Material::kCH2, meanEnergy);
    case SpeType::kCarbon: return EnergyToThickness(ThicknessConv::Material::kCarbon, meanEnergy);
    default:               return 0.0;
  }
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

void AddThicknessLegendEntry(TLegend *legend, SpeType type, double meanEnergy)
{
  if (legend == nullptr)
    return;
  const double thickness = TargetEnergyToThickness(type, meanEnergy);
  char line[256];
  std::snprintf(line, sizeof(line),
                "target thickness = %.4f #mum", thickness);
  legend->AddEntry((TObject*)nullptr, line, "");
}

TGraph *MakeThicknessEnergyGraph(const std::vector<ThicknessPoint> &points,
                                 const char *name,
                                 int color,
                                 int marker)
{
  auto graph = new TGraph((int)points.size());
  graph->SetName(name);
  graph->SetTitle(name);
  graph->SetMarkerStyle(marker);
  graph->SetMarkerSize(1.1);
  graph->SetMarkerColor(color);
  graph->SetLineColor(color);
  graph->SetLineWidth(2);
  for (int i = 0; i < (int)points.size(); ++i)
    graph->SetPoint(i, points[i].energy, points[i].thickness);
  return graph;
}

void AddThicknessEnergyDrawing(LKDrawingGroup *grp,
                               const std::vector<ThicknessPoint> &points,
                               const char *frameName,
                               const char *graphName,
                               const char *legendLabel,
                               ThicknessConv::Material material,
                               int color,
                               int marker)
{
  auto frame = new TH2D(frameName,
                        ";Remaining Energy (MeV);Target Thickness (#mum)",
                        100, 1.6, 3.5,
                        100, 0.0, MaterialThicknessRange(material));
  frame->SetDirectory(nullptr);

  auto graph = MakeThicknessEnergyGraph(points, graphName, color, marker);
  auto tableGraph = MakeAStarThicknessGraph(
      material, TString::Format("%s_astar_table", graphName), 300, 1.6, 3.5);
  tableGraph->SetLineColor(color);
  tableGraph->SetLineWidth(2);
  tableGraph->SetMarkerStyle(1);
  auto legend = new TLegend(0.12, 0.78, 0.45, 0.88);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->SetTextSize(0.032);
  legend->AddEntry(tableGraph, TString::Format("%s ASTAR table", legendLabel), "l");
  legend->AddEntry(graph, TString::Format("%s data", legendLabel), "p");

  auto draw = grp->CreateDrawing();
  draw->Add(frame, "");
  draw->Add(tableGraph, "l same");
  draw->Add(graph, "p same");
  draw->Add(legend, "same");
}

// ── Calibration record ────────────────────────────────────────────────────────
struct EnergyCal {
  double p0        = 0;
  double conv      = 0;   // MeV/ch
  double convErr   = 0;
  double sigmaCh   = 0;   // blank sigma [ch]
  double sigmaMev  = 0;   // blank sigma [MeV]
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

// Return the calibration date to use for a given file date
std::string CalDate(const std::string &fileDate,
                    const std::map<std::string, EnergyCal> &cal)
{
  if (cal.count(fileDate)) return fileDate;
  if (fileDate == "20260529" && cal.count("20260528")) return "20260528";
  if (!cal.empty()) return cal.rbegin()->first;
  return "";
}

// ── Convert channel histogram → energy histogram ──────────────────────────────
TH1D *MakeEnergyHist(TH1D *hCh, const EnergyCal &c,
                     const std::string &path, int nBins)
{
  const double eMin = (0.0             - c.p0) * c.conv;
  const double eMax = ((double)kNChannels - c.p0) * c.conv;
  const std::string name = "he_" + HistNameFromPath(path);
  auto hE = new TH1D(name.c_str(),
                     (path + ";Remaining Energy (MeV);Counts").c_str(),
                     nBins, eMin, eMax);
  hE->SetDirectory(nullptr);
  hE->SetLineWidth(2);
  for (int i = 1; i <= kNChannels; ++i) {
    const double counts = hCh->GetBinContent(i);
    if (counts <= 0) continue;
    const double channel = hCh->GetBinCenter(i);
    const double energy = (channel - c.p0) * c.conv;
    const int bin = hE->FindBin(energy);
    if (bin >= 1 && bin <= hE->GetNbinsX())
      hE->Fill(energy, counts);
  }
  return hE;
}

TH1D *MakeEnergyGraphSourceHist(TH1D *hCh, const EnergyCal &c,
                                const std::string &path)
{
  const double eMin = (0.0             - c.p0) * c.conv;
  const double eMax = ((double)kNChannels - c.p0) * c.conv;
  const std::string name = "heg_" + HistNameFromPath(path);
  auto hE = new TH1D(name.c_str(),
                     (path + ";Remaining Energy (MeV);Counts").c_str(),
                     kNChannels, eMin, eMax);
  hE->SetDirectory(nullptr);
  hE->SetLineWidth(2);
  for (int i = 1; i <= kNChannels; ++i)
    hE->SetBinContent(i, hCh->GetBinContent(i));
  return hE;
}

void SetPeakFiveSigmaRange(TH1D *hist)
{
  const auto peak = GetPeakParams(hist);
  if (peak.sigma <= 0)
    return;

  const double axMin = hist->GetXaxis()->GetXmin();
  const double axMax = hist->GetXaxis()->GetXmax();
  hist->GetXaxis()->SetRangeUser(
      std::max(axMin, peak.mean - 5.0 * peak.sigma),
      std::min(axMax, peak.mean + 5.0 * peak.sigma));
}

void SetPeakRange(TH1D *hist, double sigmaScale)
{
  const auto peak = GetPeakParams(hist);
  if (peak.sigma <= 0)
    return;

  const double axMin = hist->GetXaxis()->GetXmin();
  const double axMax = hist->GetXaxis()->GetXmax();
  hist->GetXaxis()->SetRangeUser(
      std::max(axMin, peak.mean - sigmaScale * peak.sigma),
      std::min(axMax, peak.mean + sigmaScale * peak.sigma));
}

void SetEnergyDisplayRange(TH1D *hist)
{
  hist->GetXaxis()->SetRangeUser(1.6, 3.5);
}


// ── Main ──────────────────────────────────────────────────────────────────────
void step3_energy(const char *listFile      = "list.txt",
                  int         energyBinFactor = 1,
                  bool        useGraphFit = true)
{
  const int kEnergyBins = kNChannels / std::max(energyBinFactor, 1);

  const auto cal = LoadEnergyCal("data_analysis/energy_calibration.txt");
  if (cal.empty()) return;

  const auto entries = ReadSpeList(listFile);

  auto top        = new LKDrawingGroup("step3_top");
  auto grpRef     = top->CreateGroup("blank_ch2_carbon");
  auto grpCd2     = top->CreateGroup("cd2");
  auto grpThick   = top->CreateGroup("target_thickness_vs_energy");

  std::vector<ThicknessPoint> ch2Points;
  std::vector<ThicknessPoint> carbonPoints;
  std::vector<ThicknessPoint> cd2Points;

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
    auto hGraphFit = useGraphFit ? MakeEnergyGraphSourceHist(hCh, c, e.path) : nullptr;
    delete hCh;

    LKDrawingGroup *grp = (e.type == SpeType::kCd2) ? grpCd2 : grpRef;

    const std::string label = SpeTypeName(e.type) + "  " + fileDate
                            + (calDate != fileDate ? "  [cal:" + calDate + "]" : "");

    if (e.type == SpeType::kBlank ||
        e.type == SpeType::kCh2   ||
        e.type == SpeType::kCarbon)
    {
      auto fitResult = useGraphFit
          ? ProceedGraphGaussianFittingFull(
              hGraphFit,
              true,   // graph
              false,  // initial Gaussian
              false,  // fit range lines
              true,   // final Gaussian
              true,   // legend
              "AP",
              (e.type == SpeType::kBlank) ? 2.0 : 1.5,
              10,
              (e.type == SpeType::kBlank) ? 2.0 : 3.0)
          : ProceedGaussianFittingFull(
              hE,
              true,   // histogram
              false,  // initial Gaussian
              false,  // fit range lines
              true,   // final Gaussian
              true,   // legend
              "hist",
              (e.type == SpeType::kBlank) ? 2.0 : 1.5,
              10,
              (e.type == SpeType::kBlank) ? 2.0 : 3.0);
      auto drawing = fitResult.drawing;
      const double mean = fitResult.mu;
      const double sigma = fitResult.sigma;

      // Add pure sigma info for ch2 / carbon
      const double sigmaBlankMev = (e.type == SpeType::kBlank) ? 0.0 : c.sigmaMev;
      if (sigmaBlankMev > 0) {
        const double s2 = sigma*sigma - sigmaBlankMev*sigmaBlankMev;
        const double sigmaPure = s2 > 0 ? std::sqrt(s2) : 0.0;
        auto legend = FindLegend(drawing);
        if (legend != nullptr) {
          legend->SetY1NDC(0.42);
          legend->SetTextSize(0.030);
        }
        char line[256];
        std::snprintf(line, sizeof(line),
                      "#sigma_{blank} = %.4f MeV", sigmaBlankMev);
        if (legend != nullptr) legend->AddEntry((TObject*)nullptr, line, "");
        if (sigmaPure > 0)
          std::snprintf(line, sizeof(line),
                        "#sigma_{pure} = #sqrt{#sigma^{2}-#sigma_{b}^{2}} = %.4f MeV",
                        sigmaPure);
        else
          std::snprintf(line, sizeof(line),
                        "#sigma_{pure} = 0  (#sigma < #sigma_{blank})");
        if (legend != nullptr) legend->AddEntry((TObject*)nullptr, line, "");
        AddThicknessLegendEntry(legend, e.type, mean);
      }

      if (e.type == SpeType::kCh2)
        ch2Points.push_back({mean, TargetEnergyToThickness(e.type, mean), label});
      else if (e.type == SpeType::kCarbon)
        carbonPoints.push_back({mean, TargetEnergyToThickness(e.type, mean), label});
      else if (e.type != SpeType::kBlank) {
        auto legend = FindLegend(drawing);
        AddThicknessLegendEntry(legend, e.type, mean);
      }

      grp->AddDrawing(drawing);
    }
    else {
      // cd2: 피팅 없이 히스토그램 + 레이블만
      const double mean = hE->GetMean();
      const double thickness = TargetEnergyToThickness(e.type, mean);
      cd2Points.push_back({mean, thickness, label});
      SetEnergyDisplayRange(hE);
      auto draw = grp->CreateDrawing();
      draw->Add(hE, "hist");
      auto legend = new TLegend(0.12, 0.72, 0.65, 0.88);
      legend->SetBorderSize(0); legend->SetFillStyle(0); legend->SetMargin(0.12);
      legend->SetTextSize(0.032);
      legend->AddEntry((TObject*)nullptr, label.c_str(), "");
      char line[256];
      std::snprintf(line, sizeof(line), "target thickness = %.4f #mum", thickness);
      legend->AddEntry((TObject*)nullptr, line, "");
      draw->Add(legend, "same");
    }
  }

  AddThicknessEnergyDrawing(grpThick, cd2Points,
                            "frame_cd2_thickness_vs_energy",
                            "g_cd2_thickness_vs_energy",
                            "cd2", ThicknessConv::Material::kCD2,
                            kRed + 1, 22);
  AddThicknessEnergyDrawing(grpThick, ch2Points,
                            "frame_ch2_thickness_vs_energy",
                            "g_ch2_thickness_vs_energy",
                            "ch2", ThicknessConv::Material::kCH2,
                            kBlue + 1, 20);
  AddThicknessEnergyDrawing(grpThick, carbonPoints,
                            "frame_carbon_thickness_vs_energy",
                            "g_carbon_thickness_vs_energy",
                            "carbon", ThicknessConv::Material::kCarbon,
                            kGreen + 2, 21);

  top->WriteFile("data_analysis/step3_energy.root");
  std::cout << "Saved: data_analysis/step3_energy.root" << std::endl;
  top->Draw();
  top->Save();
}
