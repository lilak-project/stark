/**
 * test_unfold_linear.C
 *
 * Example 2: Iterative unfolding on a linearly-sloped distribution.
 *
 * Procedure:
 *   1. Build a "true" distribution with a linear slope: p(x) ∝ (x - xMin)
 *   2. Smear it with a Gaussian IRF (sigma_irf=0.5) → "measured"
 *      (smearing distorts the sharp edges and slope)
 *   3. Unfold and compare at several iterations
 *
 * This case is more challenging than the Gaussian example because:
 *   - The distribution has sharp boundaries at xMin and xMax
 *   - A linear slope is NOT well described by a simple sigma subtraction
 *   - Shows that iterative unfolding works for non-Gaussian true distributions
 */

#include "IterativeUnfolder.h"

#include "TCanvas.h"
#include "TF1.h"
#include "TH1D.h"
#include "TLegend.h"
#include "TLine.h"
#include "TRandom3.h"
#include "TStyle.h"

void test_unfold_linear()
{
  gStyle->SetOptStat(0);

  // ── Parameters ──────────────────────────────────────────────────────────────
  constexpr int    kNEvents = 200000;
  constexpr double kXMin    = 0.0;
  constexpr double kXMax    = 10.0;
  constexpr int    kNBins   = 100;
  constexpr double kSigIRF  = 0.5;   // smearing (IRF) sigma — from blank

  // True distribution: p(x) ∝ x  (linear slope, zero at xMin)
  // Sample via inverse CDF: x = xMax * sqrt(U),  U ~ Uniform[0,1]
  TRandom3 rng(123);

  // ── Build true histogram ─────────────────────────────────────────────────────
  auto hTrue = new TH1D("hTrue", "Linear-slope unfolding example;x;Counts",
                         kNBins, kXMin, kXMax);
  hTrue->SetDirectory(nullptr);

  auto SampleLinear = [&]() -> double {
    // Inverse CDF of p(x)∝x on [0, xMax]: x = xMax*sqrt(U)
    return kXMax * std::sqrt(rng.Uniform());
  };

  for (int i = 0; i < kNEvents; ++i)
    hTrue->Fill(SampleLinear());

  // ── Smear to get "measured" histogram ────────────────────────────────────────
  auto hMeas = new TH1D("hMeas", "", kNBins, kXMin, kXMax);
  hMeas->SetDirectory(nullptr);

  for (int i = 0; i < kNEvents; ++i)
    hMeas->Fill(rng.Gaus(SampleLinear(), kSigIRF));

  // ── Unfold at several iterations ─────────────────────────────────────────────
  IterativeUnfolder unfolder(hMeas, kSigIRF);

  auto hUnfold3  = unfolder.Unfold(3,  "hUnfold3",  "Unfolded (3 iter)");
  auto hUnfold5  = unfolder.Unfold(5,  "hUnfold5",  "Unfolded (5 iter)");
  auto hUnfold10 = unfolder.Unfold(10, "hUnfold10", "Unfolded (10 iter)");
  auto hUnfold20 = unfolder.Unfold(20, "hUnfold20", "Unfolded (20 iter)");

  // ── Styling ──────────────────────────────────────────────────────────────────
  hTrue->SetLineColor(kBlack);       hTrue->SetLineWidth(2);
  hMeas->SetLineColor(kRed+1);       hMeas->SetLineWidth(2); hMeas->SetLineStyle(2);
  hUnfold3->SetLineColor(kAzure+7);  hUnfold3->SetLineWidth(2);
  hUnfold5->SetLineColor(kGreen+2);  hUnfold5->SetLineWidth(2);
  hUnfold10->SetLineColor(kOrange+7);hUnfold10->SetLineWidth(2);
  hUnfold20->SetLineColor(kMagenta+1);hUnfold20->SetLineWidth(2);

  // ── Canvas 1: overlay ────────────────────────────────────────────────────────
  auto c1 = new TCanvas("c1", "Linear unfolding", 900, 600);
  c1->SetLeftMargin(0.12);

  double yMax = hTrue->GetMaximum() * 1.4;
  hTrue->SetMaximum(yMax);
  hTrue->Draw("hist");
  hMeas->Draw("hist same");
  hUnfold3->Draw("hist same");
  hUnfold5->Draw("hist same");
  hUnfold10->Draw("hist same");
  hUnfold20->Draw("hist same");

  auto leg = new TLegend(0.15, 0.60, 0.50, 0.88);
  leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextSize(0.032);
  leg->AddEntry(hTrue,     "True  (p(x) #propto x)",       "l");
  leg->AddEntry(hMeas,     Form("Measured (IRF #sigma=%.1f)", kSigIRF), "l");
  leg->AddEntry(hUnfold3,  "Unfolded  3 iter", "l");
  leg->AddEntry(hUnfold5,  "Unfolded  5 iter", "l");
  leg->AddEntry(hUnfold10, "Unfolded 10 iter", "l");
  leg->AddEntry(hUnfold20, "Unfolded 20 iter", "l");
  leg->Draw();

  c1->SaveAs("figures/test_unfold_linear_overlay.png");

  // ── Canvas 2: residuals ───────────────────────────────────────────────────────
  auto c2 = new TCanvas("c2", "Linear residuals", 900, 500);
  c2->SetLeftMargin(0.12);

  auto MakeResidual = [&](TH1D *hUnfold, const char *name) -> TH1D* {
    auto hRes = static_cast<TH1D*>(hUnfold->Clone(name));
    hRes->SetDirectory(nullptr);
    for (int i = 1; i <= kNBins; ++i) {
      const double t = hTrue->GetBinContent(i);
      if (t <= 0) { hRes->SetBinContent(i, 0); continue; }
      hRes->SetBinContent(i, (hUnfold->GetBinContent(i) - t) / t);
    }
    return hRes;
  };

  auto hRes3  = MakeResidual(hUnfold3,  "hRes3");
  auto hRes5  = MakeResidual(hUnfold5,  "hRes5");
  auto hRes10 = MakeResidual(hUnfold10, "hRes10");
  auto hRes20 = MakeResidual(hUnfold20, "hRes20");

  hRes3->SetTitle("Residual (unfolded - true) / true;x;Relative residual");
  hRes3->SetMaximum( 0.5);
  hRes3->SetMinimum(-0.5);
  hRes3->Draw("hist");
  hRes5->Draw("hist same");
  hRes10->Draw("hist same");
  hRes20->Draw("hist same");

  TF1 *fZero = new TF1("fZero","0", kXMin, kXMax);
  fZero->SetLineColor(kBlack); fZero->SetLineStyle(2);
  fZero->Draw("same");

  auto legR = new TLegend(0.15, 0.70, 0.45, 0.88);
  legR->SetBorderSize(0); legR->SetFillStyle(0); legR->SetTextSize(0.033);
  legR->AddEntry(hRes3,  " 3 iter", "l");
  legR->AddEntry(hRes5,  " 5 iter", "l");
  legR->AddEntry(hRes10, "10 iter", "l");
  legR->AddEntry(hRes20, "20 iter", "l");
  legR->Draw();

  c2->SaveAs("figures/test_unfold_linear_residual.png");

  // ── Canvas 3: convergence check — smeared estimate vs measured ───────────────
  // After unfolding, (h_unfolded ⊗ IRF) should match h_meas closely
  auto c3 = new TCanvas("c3", "Convergence check", 900, 500);
  c3->SetLeftMargin(0.12);

  auto hSmeared5  = unfolder.Unfold(5,  "tmp5",  "tmp");  // rerun to get internal state
  // Re-create unfolder at 5 iter to get smeared version
  IterativeUnfolder uf2(hMeas, kSigIRF);
  uf2.Unfold(5);
  auto hSmearedCheck = uf2.GetSmeared("hSmearedCheck");

  hMeas->SetMaximum(hMeas->GetMaximum() * 1.3);
  hMeas->Draw("hist");
  hSmearedCheck->SetLineColor(kGreen+2);
  hSmearedCheck->SetLineWidth(2);
  hSmearedCheck->Draw("hist same");

  auto legC = new TLegend(0.15, 0.72, 0.55, 0.88);
  legC->SetBorderSize(0); legC->SetFillStyle(0); legC->SetTextSize(0.033);
  legC->AddEntry(hMeas,         "Measured",                     "l");
  legC->AddEntry(hSmearedCheck, "Unfolded (5 iter) #otimes IRF", "l");
  legC->Draw();

  c3->SetTitle("Convergence: smeared unfolded should match measured");
  c3->SaveAs("figures/test_unfold_linear_convergence.png");

  std::cout << "\n=== Linear unfolding summary ===" << std::endl;
  std::cout << "  True distribution: p(x) proportional to x on ["
            << kXMin << ", " << kXMax << "]" << std::endl;
  std::cout << "  IRF sigma: " << kSigIRF << std::endl;
  std::cout << "  Output: figures/test_unfold_linear_overlay.png" << std::endl;
  std::cout << "          figures/test_unfold_linear_residual.png" << std::endl;
  std::cout << "          figures/test_unfold_linear_convergence.png" << std::endl;

  delete hSmeared5;
}
