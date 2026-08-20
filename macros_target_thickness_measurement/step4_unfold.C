/**
 * step4_unfold.C
 *
 * Loads energy histograms from step3_energy.root and applies
 * iterative unfolding (IterativeUnfolder) to remove detector smearing.
 *
 * IRF: Gaussian with sigma = sigma_blank (from energy_calibration.txt)
 *
 * Input:
 *   data_analysis/step3_energy.root       — energy histograms from step3
 *   data_analysis/energy_calibration.txt  — sigma_blank per date
 *
 * Run: root -l -b -q step4_unfold.C
 */

#include "analysis_utils.h"
#include "IterativeUnfolder.h"

#include <map>
#include "TFile.h"
#include "TKey.h"
#include "LKDrawing.h"
#include "LKDrawingGroup.h"

// ── Load sigma_blank per date ─────────────────────────────────────────────────
std::map<std::string, double> LoadSigmaBlank(const char *calFile)
{
  std::map<std::string, double> sigmas;
  std::ifstream f(calFile);
  if (!f.is_open()) {
    std::cerr << "Cannot open " << calFile << std::endl;
    return sigmas;
  }
  std::string line;
  while (std::getline(f, line)) {
    line = Trim(line);
    if (line.empty() || line[0] == '#') continue;
    std::istringstream ss(line);
    std::string date;
    double p0, conv, convErr, sigmaCh, sigmaMev;
    if (ss >> date >> p0 >> conv >> convErr >> sigmaCh >> sigmaMev)
      sigmas[date] = sigmaMev;
  }
  return sigmas;
}

// Return sigma_blank [MeV] for a given file date
double GetSigmaBlank(const std::string &fileDate,
                     const std::map<std::string, double> &sigmas)
{
  if (sigmas.count(fileDate))        return sigmas.at(fileDate);
  if (fileDate == "20260529" &&
      sigmas.count("20260528"))      return sigmas.at("20260528");
  if (!sigmas.empty())               return sigmas.begin()->second;
  return 0.0;
}

void step4_unfold(const char *step3File = "data_analysis/step3_energy.root",
                  int         nIter     = 5)
{
  // ── Load sigma_blank ───────────────────────────────────────────────────────
  const auto sigmaMap = LoadSigmaBlank("data_analysis/energy_calibration.txt");
  if (sigmaMap.empty()) {
    std::cerr << "No calibration data — run step2_blank.C first." << std::endl;
    return;
  }

  // ── Load step3 root file ───────────────────────────────────────────────────
  auto f3 = TFile::Open(step3File, "READ");
  if (!f3 || f3->IsZombie()) {
    std::cerr << "Cannot open " << step3File
              << " — run step3_energy.C first." << std::endl;
    return;
  }

  // ── Drawing groups ─────────────────────────────────────────────────────────
  auto top      = new LKDrawingGroup("step4_top");
  auto grpRef   = top->CreateGroup("blank_ch2_carbon");
  auto grpCd2   = top->CreateGroup("cd2_unfolded");

  // ── Iterate over histograms in step3 root file ────────────────────────────
  TIter next(f3->GetListOfKeys());
  TKey *key = nullptr;
  while ((key = static_cast<TKey*>(next()))) {
    if (!TString(key->GetClassName()).Contains("TH1")) continue;

    auto hMeas = static_cast<TH1D*>(key->ReadObj());
    if (!hMeas) continue;
    hMeas->SetDirectory(nullptr);

    const std::string hname = hMeas->GetName();

    // Determine type and date from histogram name (encoded as "he_h_<path>")
    // Quick classification by name content
    const bool isBlank  = hname.find("blank") != std::string::npos ||
                          hname.find("Blank") != std::string::npos;
    const bool isCh2    = hname.find("ch2")   != std::string::npos;
    const bool isCarbon = hname.find("arbon")  != std::string::npos;
    const bool isCd2    = !isBlank && !isCh2 && !isCarbon;

    // Extract date from histogram name (8-digit number)
    std::string fileDate;
    for (size_t i = 0; i + 8 <= hname.size(); ++i) {
      const std::string sub = hname.substr(i, 8);
      if (sub.substr(0,4) == "2026") { fileDate = sub; break; }
    }

    const double sigmaBlank = GetSigmaBlank(fileDate, sigmaMap);
    if (sigmaBlank <= 0) {
      std::cerr << "No sigma_blank for " << hname << std::endl;
      continue;
    }

    // ── Unfold ──────────────────────────────────────────────────────────────
    IterativeUnfolder unfolder(hMeas, sigmaBlank);
    auto hUnfolded = unfolder.Unfold(nIter,
                                     ("unfolded_" + hname).c_str(),
                                     (hMeas->GetTitle() + std::string(" [unfolded]")).c_str());

    // Style
    hMeas->SetLineColor(kGray + 1);
    hMeas->SetLineWidth(1);
    hUnfolded->SetLineColor(kBlue + 1);
    hUnfolded->SetLineWidth(2);

    SetFiveSigmaRange(hMeas);
    hUnfolded->GetXaxis()->SetRangeUser(hMeas->GetXaxis()->GetFirst(),
                                        hMeas->GetXaxis()->GetLast());

    // Legend
    auto legend = new TLegend(0.12, 0.68, 0.65, 0.88);
    legend->SetBorderSize(0);
    legend->SetFillStyle(0);
    legend->SetTextSize(0.032);
    char line[256];
    std::snprintf(line, sizeof(line),
                  "#sigma_{blank} = %.4f MeV  (IRF)", sigmaBlank);
    legend->AddEntry(static_cast<TObject*>(nullptr), line, "");
    legend->AddEntry(hMeas,     "Measured", "l");
    std::snprintf(line, sizeof(line), "Unfolded (%d iter)", nIter);
    legend->AddEntry(hUnfolded, line, "l");

    LKDrawingGroup *grp = isCd2 ? grpCd2 : grpRef;
    auto draw = grp->CreateDrawing();
    draw->Add(hMeas,     "hist");
    draw->Add(hUnfolded, "hist same");
    draw->Add(legend,    "same");
    draw->SetLegendCorner(1);
  }

  f3->Close();

  top->WriteFile("data_analysis/step4_unfold.root");
  std::cout << "Saved: data_analysis/step4_unfold.root" << std::endl;
  top->Draw();
}
