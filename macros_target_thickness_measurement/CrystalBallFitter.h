#ifndef CRYSTAL_BALL_FITTER_H
#define CRYSTAL_BALL_FITTER_H

/**
 * CrystalBallFitter.h
 *
 * ProceedCrystalBallFitting — left-tail Crystal Ball fit helper.
 * Mirrors GaussianFitter.h (same flags / return pattern / styling).
 *
 * Left-tail Crystal Ball:
 *   f(x) = A * exp(-t^2/2)                       for t > -alpha
 *        = A * a * (b - t)^(-n)                  for t <= -alpha
 *   with t = (x - mean)/sigma,
 *        a = (n/alpha)^n * exp(-alpha^2/2),  b = n/alpha - alpha
 *
 * Initial parameters:
 *   amplitude = GetBinContent(GetMaximumBin())
 *   mean      = GetBinCenter(GetMaximumBin())   (peak, robust against tail)
 *   sigma     = FWHM-based estimate
 *   alpha     = 1.0
 *   n         = 3.0
 *
 * Fit range: [mean - 7*sigma, mean + 4*sigma]   (wide left side for tail)
 * x-axis range after fit: [mu - 8*sigma, mu + 5*sigma]
 *
 * All TF1: SetNpx(500), Fit option "RQN0".
 */

#include <cmath>
#include <cstdio>

#include "TF1.h"
#include "TH1D.h"
#include "TLegend.h"
#include "TLine.h"

#include "LKDrawing.h"

// ── Left-tail Crystal Ball function ──────────────────────────────────────────
inline double CBLeftTail(double *x, double *p)
{
  const double amplitude = p[0];
  const double mean       = p[1];
  const double sigma      = std::abs(p[2]);
  const double alpha      = std::abs(p[3]);
  const double n          = std::abs(p[4]);
  if (sigma <= 0 || alpha <= 0 || n <= 0) return 0;

  const double t = (x[0] - mean) / sigma;
  if (t > -alpha)
    return amplitude * std::exp(-0.5 * t * t);

  const double a = std::pow(n / alpha, n) * std::exp(-0.5 * alpha * alpha);
  const double b = n / alpha - alpha;
  return amplitude * a * std::pow(b - t, -n);
}

// ── Unique name helper ────────────────────────────────────────────────────────
inline TString CBUniqueName(const char *prefix, const TH1D *hist)
{
  static int counter = 0;
  return TString::Format("%s_%s_%d", prefix, hist->GetName(), ++counter);
}

// ── FWHM-based sigma estimate ─────────────────────────────────────────────────
inline double CBEstimateSigma(TH1D *hist)
{
  const int    maxBin  = hist->GetMaximumBin();
  const double halfMax = 0.5 * hist->GetBinContent(maxBin);
  if (halfMax <= 0) return hist->GetStdDev();

  int lb = maxBin, rb = maxBin;
  while (lb > 1                 && hist->GetBinContent(lb) > halfMax) --lb;
  while (rb < hist->GetNbinsX() && hist->GetBinContent(rb) > halfMax) ++rb;
  const double fwhm = hist->GetBinCenter(rb) - hist->GetBinCenter(lb);
  return fwhm > 0 ? fwhm / 2.35482 : hist->GetStdDev();
}

// ── Result struct ─────────────────────────────────────────────────────────────
struct CrystalBallFitResult {
  LKDrawing *drawing  = nullptr;
  TF1       *f_init   = nullptr;
  TF1       *f_result = nullptr;
  double     mu       = 0;
  double     sigma    = 0;
  double     alpha    = 0;
  double     n        = 0;
  double     chi2ndf  = 0;
};

// ── Core implementation ───────────────────────────────────────────────────────
inline CrystalBallFitResult _CrystalBallFitImpl(
    TH1D       *hist,
    bool        add_hist,
    bool        add_f_init,
    bool        add_range,
    bool        add_f_result,
    bool        add_legend,
    const char *hist_option)
{
  // ── Initial parameters ────────────────────────────────────────────────────
  const int    maxBin = hist->GetMaximumBin();
  const double amp0   = hist->GetBinContent(maxBin);
  const double mean0  = hist->GetBinCenter(maxBin);
  const double sigma0 = CBEstimateSigma(hist);
  const double alpha0 = 1.0;
  const double n0     = 3.0;

  const double axMin  = hist->GetXaxis()->GetXmin();
  const double axMax  = hist->GetXaxis()->GetXmax();
  const double fitMin = std::max(axMin, mean0 - 7.0 * sigma0);
  const double fitMax = std::min(axMax, mean0 + 4.0 * sigma0);

  // ── Initial Crystal Ball (dashed) ─────────────────────────────────────────
  auto f_init = new TF1(CBUniqueName("cb_init", hist), CBLeftTail, fitMin, fitMax, 5);
  f_init->SetParameters(amp0, mean0, sigma0, alpha0, n0);
  f_init->SetNpx(500);
  f_init->SetLineColor(kGray + 1);
  f_init->SetLineStyle(2);
  f_init->SetLineWidth(2);

  // ── Fit range vertical lines ──────────────────────────────────────────────
  const double yMax = amp0;
  auto lMin = new TLine(fitMin, 0, fitMin, yMax * 0.90);
  auto lMax = new TLine(fitMax, 0, fitMax, yMax * 0.90);
  for (auto l : {lMin, lMax}) {
    l->SetLineColor(kGray + 2);
    l->SetLineStyle(2);
    l->SetLineWidth(1);
  }

  // ── Crystal Ball fit ──────────────────────────────────────────────────────
  auto f_result = new TF1(CBUniqueName("cb_result", hist), CBLeftTail, fitMin, fitMax, 5);
  f_result->SetParameters(amp0, mean0, sigma0, alpha0, n0);
  f_result->SetParNames("Amplitude", "Mean", "Sigma", "Alpha", "N");
  f_result->SetParLimits(0, 0,          10.0 * std::max(amp0,   1.e-9));
  f_result->SetParLimits(1, fitMin,     fitMax);
  f_result->SetParLimits(2, 0.2*sigma0, 5.0  * std::max(sigma0, 1.e-9));
  f_result->SetParLimits(3, 0.05,       8.0);
  f_result->SetParLimits(4, 1.01,       30.0);
  f_result->SetNpx(500);
  f_result->SetLineColor(kRed + 1);
  f_result->SetLineWidth(2);
  hist->Fit(f_result, "RQN0");

  const double mu      = f_result->GetParameter(1);
  const double sigma   = std::abs(f_result->GetParameter(2));
  const double alpha   = std::abs(f_result->GetParameter(3));
  const double nPar    = std::abs(f_result->GetParameter(4));
  const double chi2ndf = f_result->GetNDF() > 0
                       ? f_result->GetChisquare() / f_result->GetNDF() : 0;

  // ── Set x-axis range using fit result ────────────────────────────────────
  if (sigma > 0)
    hist->GetXaxis()->SetRangeUser(
        std::max(axMin, mu - 8.0 * sigma),
        std::min(axMax, mu + 5.0 * sigma));

  // ── Manual TLegend ────────────────────────────────────────────────────────
  TLegend *legend = nullptr;
  if (add_legend) {
    legend = new TLegend(0.12, 0.58, 0.65, 0.88);
    legend->SetBorderSize(0);
    legend->SetFillStyle(0);
    legend->SetTextSize(0.032);
    legend->SetMargin(0.12);
    char line[256];

    // initial params — comment out to hide from legend
    //if (add_f_init) {
    //  std::snprintf(line, sizeof(line),
    //                "Initial:  #mu_{0}=%.3f  #sigma_{0}=%.3f", mean0, sigma0);
    //  legend->AddEntry(f_init, line, "l");
    //}
    // fit range — comment out to hide from legend
    //if (add_range) {
    //  std::snprintf(line, sizeof(line),
    //                "Fit range: [%.3f, %.3f]", fitMin, fitMax);
    //  legend->AddEntry(lMin, line, "l");
    //}

    legend->AddEntry(f_result, "Crystal Ball fit", "l");
    std::snprintf(line, sizeof(line), "#mu = %.4f", mu);
    legend->AddEntry((TObject*)nullptr, line, "");
    std::snprintf(line, sizeof(line), "#sigma = %.4f", sigma);
    legend->AddEntry((TObject*)nullptr, line, "");
    std::snprintf(line, sizeof(line), "#alpha = %.4f", alpha);
    legend->AddEntry((TObject*)nullptr, line, "");
    std::snprintf(line, sizeof(line), "n = %.4f", nPar);
    legend->AddEntry((TObject*)nullptr, line, "");
    std::snprintf(line, sizeof(line), "#chi^{2}/NDF = %.3g/%d = %.3g",
                  f_result->GetChisquare(), f_result->GetNDF(), chi2ndf);
    legend->AddEntry((TObject*)nullptr, line, "");
  }

  // ── Assemble LKDrawing ────────────────────────────────────────────────────
  auto drawing = new LKDrawing(CBUniqueName("cbdraw", hist));

  if (add_hist)     drawing->Add(hist,     hist_option);
  if (add_f_init)   drawing->Add(f_init,   "samel");
  if (add_range)  { drawing->Add(lMin,     "same");
                    drawing->Add(lMax,     "same"); }
  if (add_f_result) drawing->Add(f_result, "samel");
  if (add_legend)   drawing->Add(legend,   "same");

  return {drawing, f_init, f_result, mu, sigma, alpha, nPar, chi2ndf};
}

// ── Public: returns LKDrawing* ────────────────────────────────────────────────
inline LKDrawing *ProceedCrystalBallFitting(
    TH1D       *hist,
    bool        add_hist     = true,
    bool        add_f_init   = true,
    bool        add_range    = true,
    bool        add_f_result = true,
    bool        add_legend   = true,
    const char *hist_option  = "hist")
{
  return _CrystalBallFitImpl(hist, add_hist, add_f_init, add_range,
                             add_f_result, add_legend, hist_option).drawing;
}

// ── Public: returns full result struct ────────────────────────────────────────
inline CrystalBallFitResult ProceedCrystalBallFittingFull(
    TH1D       *hist,
    bool        add_hist     = true,
    bool        add_f_init   = true,
    bool        add_range    = true,
    bool        add_f_result = true,
    bool        add_legend   = true,
    const char *hist_option  = "hist")
{
  return _CrystalBallFitImpl(hist, add_hist, add_f_init, add_range,
                             add_f_result, add_legend, hist_option);
}

#endif // CRYSTAL_BALL_FITTER_H
