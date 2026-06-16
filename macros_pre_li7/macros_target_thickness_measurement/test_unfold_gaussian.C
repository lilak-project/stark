/**
 * test_unfold_gaussian.C
 *
 * Iterative unfolding example on a Gaussian distribution.
 *
 * Drawings:
 *   1. True Gaussian
 *   2. Smeared (measured)
 *   3~12. Unfolded iter 1..10  (best iteration highlighted in red)
 *   13. chi2 vs iteration curve
 */

#include "analysis_utils.h"
#include "IterativeUnfolder.h"
#include "GaussianFitter.h"

#include "TGraph.h"
#include "TLine.h"
#include "TRandom3.h"
#include "LKDrawing.h"
#include "LKDrawingGroup.h"

void test_unfold_gaussian()
{
  constexpr int    kNEvents  = 100000;
  constexpr double kMu       = 5.0;
  constexpr double kSigTrue  = 1.0;
  constexpr double kSigIRF   = 0.5;
  constexpr double kXMin     = 0.0;
  constexpr double kXMax     = 10.0;
  constexpr int    kNBins    = 100;
  constexpr int    kMaxIter  = 10;

  TRandom3 rng(42);

  auto hTrue = new TH1D("hTrue", "True Gaussian;x;Counts",      kNBins, kXMin, kXMax);
  auto hMeas = new TH1D("hMeas", "Smeared (measured);x;Counts", kNBins, kXMin, kXMax);
  hTrue->SetDirectory(nullptr);
  hMeas->SetDirectory(nullptr);

  for (int i = 0; i < kNEvents; ++i) {
    const double x = rng.Gaus(kMu, kSigTrue);
    hTrue->Fill(x);
    hMeas->Fill(rng.Gaus(x, kSigIRF));
  }

  hTrue->SetLineColor(kBlack);  hTrue->SetLineWidth(2);
  hMeas->SetLineColor(kGray+1); hMeas->SetLineWidth(2);

  // ── Run auto stopping to find best iteration ──────────────────────────────
  std::vector<double> chi2History;
  IterativeUnfolder unfolder(hMeas, kSigIRF);
  {
    TH1D *tmp = unfolder.UnfoldAuto(kMaxIter, 0.0, &chi2History);
    delete tmp;
  }
  const int bestIter = unfolder.GetBestIter();

  // ── Unfold each iteration individually ───────────────────────────────────
  TH1D *hUnfold[kMaxIter];
  for (int k = 0; k < kMaxIter; ++k) {
    char name[64], title[128];
    std::snprintf(name,  sizeof(name),  "hUnfold_%d", k + 1);
    std::snprintf(title, sizeof(title), "Unfolded iter %d;x;Counts", k + 1);
    hUnfold[k] = unfolder.Unfold(k + 1, name, title);

    const bool isBest = (k + 1 == bestIter);
    hUnfold[k]->SetLineColor(isBest ? kRed+1 : kBlue+1);
    hUnfold[k]->SetLineWidth(isBest ? 2 : 1);
  }

  // ── LKDrawingGroup ────────────────────────────────────────────────────────
  auto top = new LKDrawingGroup("unfold_gaussian");

  // 1. True
  top->AddDrawing(ProceedGaussianFitting(hTrue));

  // 2. Smeared
  top->AddDrawing(ProceedGaussianFitting(hMeas));

  // 3~12. Each iteration
  for (int k = 0; k < kMaxIter; ++k) {
    const bool isBest = (k + 1 == bestIter);
    auto d = ProceedGaussianFitting(hUnfold[k]);

    // Best iteration: add extra label to legend
    if (isBest) {
      auto mark = new TLegend(0.12, 0.54, 0.65, 0.62);
      mark->SetBorderSize(0);
      mark->SetFillStyle(0);
      mark->SetTextSize(0.034);
      mark->SetMargin(0.12);
      mark->AddEntry((TObject*)nullptr,
                     "#bf{#leftarrow Best iteration (chi^{2} stopping)}", "");
      d->Add(mark, "same");
    }

    top->AddDrawing(d);
  }

  // ── chi2 vs iteration drawing ─────────────────────────────────────────────
  {
    const int n = static_cast<int>(chi2History.size());
    auto gChi2 = new TGraph(n);
    gChi2->SetName("g_chi2");
    gChi2->SetTitle("#chi^{2} vs Iteration;Iteration;#chi^{2}");
    gChi2->SetMarkerStyle(20);
    gChi2->SetMarkerSize(1.0);
    gChi2->SetLineWidth(2);
    gChi2->SetMarkerColor(kBlue+1);
    gChi2->SetLineColor(kBlue+1);
    for (int i = 0; i < n; ++i)
      gChi2->SetPoint(i, i + 1, chi2History[i]);

    const double chi2Max = *std::max_element(chi2History.begin(), chi2History.end());
    auto lBest = new TLine(bestIter, 0, bestIter, chi2Max * 1.05);
    lBest->SetLineColor(kRed+1);
    lBest->SetLineStyle(2);
    lBest->SetLineWidth(2);

    auto leg = new TLegend(0.45, 0.72, 0.88, 0.88);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.034);
    leg->SetMargin(0.12);
    leg->AddEntry(gChi2, "#chi^{2} = #Sigma(sim-meas)^{2}/meas", "lp");
    char line[128];
    std::snprintf(line, sizeof(line), "Best iter = %d  (#chi^{2}=%.3g)",
                  bestIter, unfolder.GetBestChi2());
    leg->AddEntry(lBest, line, "l");

    auto drawChi2 = new LKDrawing("draw_chi2");
    drawChi2->Add(gChi2, "apl");
    drawChi2->Add(lBest, "same");
    drawChi2->Add(leg,   "same");
    top->AddDrawing(drawChi2);
  }

  top->Draw();
}
