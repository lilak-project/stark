#ifndef GAUSSIAN_FITTER_H
#define GAUSSIAN_FITTER_H

/**
 * GaussianFitter.h
 *
 * ProceedGaussianFitting — centralised Gaussian fit helper.
 *
 * Initial parameters:
 *   amplitude = GetBinContent(GetMaximumBin())
 *   mean      = GetMean()
 *   sigma     = GetStdDev()
 *   if sigma is zero, use one bin width as a finite fallback
 *
 * Fit range is optimized iteratively.
 *   Each iteration tests sd_range = 0.8, 1.0, 1.6 times the previous
 *   sd_range, clamped to [sd_range_min, sd_range_max], and keeps the lowest
 *   chi2/ndf fit among ranges containing at least min_points_in_range bins.
 *   The next iteration starts from the previous best fit parameters.
 *   Iteration stops when amplitude, mean, sigma change by <= 5%.
 * x-axis range after fit: [mu_fit - 5*sigma_fit, mu_fit + 5*sigma_fit]
 *
 * Returns LKDrawing* containing:
 *   hist       (drawn, not in legend)
 *   f_init     (dashed gray, drawn but not in legend)   ← comment out to re-enable
 *   lMin/lMax  (vertical dashed lines for fit range)    ← comment out to re-enable
 *   f_result   (solid red)
 *   TLegend    (manually created — no LKDrawing internal legend)
 *
 * All TF1: SetNpx(500), Fit option "RQN0".
 */

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <limits>
#include <vector>

#include "TF1.h"
#include "TH1D.h"
#include "TLegend.h"
#include "TLine.h"
#include "TGraph.h"

#include "LKDrawing.h"

// ── Unique name helper ────────────────────────────────────────────────────────
inline TString GFUniqueName(const char *prefix, const TH1D *hist)
{
  static int counter = 0;
  return TString::Format("%s_%s_%d", prefix, hist->GetName(), ++counter);
}

// ── Result struct ─────────────────────────────────────────────────────────────
struct GaussianFitResult {
  LKDrawing *drawing  = nullptr;
  TF1       *f_init   = nullptr;
  TF1       *f_result = nullptr;
  double     mu       = 0;
  double     sigma    = 0;
  double     chi2ndf  = 0;
  double     sd_range = 0;
  TGraph    *g_amp    = nullptr;
  TGraph    *g_mean   = nullptr;
  TGraph    *g_sigma  = nullptr;
  TGraph    *g_chi2ndf = nullptr;
  TGraph    *g_sd_range = nullptr;
};

inline double GFClamp(double x, double xMin, double xMax)
{
  return std::max(xMin, std::min(xMax, x));
}

inline bool GFWithinRelativeChange(double before, double after, double frac = 0.05)
{
  const double scale = std::max(std::abs(before), 1.e-12);
  return std::abs(after - before) / scale <= frac;
}

inline int GFCountBinsInRange(const TH1D *hist, double xMin, double xMax)
{
  int nBins = 0;
  for (int i = 1; i <= hist->GetNbinsX(); ++i) {
    const double x = hist->GetBinCenter(i);
    if (x >= xMin && x <= xMax)
      ++nBins;
  }
  return nBins;
}

inline std::pair<double, double> GFMeanStdDevInRange(
    const TH1D *hist, double xMin, double xMax,
    double fallbackMean, double fallbackSigma)
{
  double sumW = 0;
  double sumWX = 0;
  double sumWX2 = 0;
  for (int i = 1; i <= hist->GetNbinsX(); ++i) {
    const double x = hist->GetBinCenter(i);
    if (x < xMin || x > xMax)
      continue;
    const double w = hist->GetBinContent(i);
    if (w <= 0)
      continue;
    sumW += w;
    sumWX += w * x;
    sumWX2 += w * x * x;
  }
  if (sumW <= 0)
    return {fallbackMean, fallbackSigma};

  const double mean = sumWX / sumW;
  const double variance = sumWX2 / sumW - mean * mean;
  if (variance <= 0)
    return {mean, fallbackSigma};
  return {mean, std::sqrt(variance)};
}

inline TGraph *GFMakeGraph(const char *prefix, const TH1D *hist,
                           const std::vector<double> &x,
                           const std::vector<double> &y,
                           int color)
{
  auto graph = new TGraph((int)x.size());
  graph->SetName(GFUniqueName(prefix, hist));
  graph->SetMarkerStyle(20);
  graph->SetMarkerSize(0.9);
  graph->SetMarkerColor(color);
  graph->SetLineColor(color);
  graph->SetLineWidth(2);
  for (int i = 0; i < (int)x.size(); ++i)
    graph->SetPoint(i, x[i], y[i]);
  return graph;
}

// ── Core implementation (shared by both overloads) ────────────────────────────
inline GaussianFitResult _GaussianFitImpl(
    TH1D       *hist,
    bool        add_hist,
    bool        add_f_init,
    bool        add_range,
    bool        add_f_result,
    bool        add_legend,
    const char *hist_option,
    double      sd_range_min,
    int         min_points_in_range,
    double      sd_range_max,
    bool        add_iteration_info,
    double      initial_mean_override)
{
  // ── Initial parameters ────────────────────────────────────────────────────
  const int    maxBin = hist->GetMaximumBin();
  const double amp0   = hist->GetBinContent(maxBin);
  const double mean0  = hist->GetMean();
  const double sigmaStats = hist->GetStdDev();
  const double sigma0 = sigmaStats > 0 ? sigmaStats : hist->GetBinWidth(maxBin);

  const double axMin  = hist->GetXaxis()->GetXmin();
  const double axMax  = hist->GetXaxis()->GetXmax();
  const double sdRangeMax = std::max(sd_range_max, 0.0);
  const double sdRangeMin = std::min(GFClamp(sd_range_min, 0.0, sdRangeMax),
                                     sdRangeMax);
  const int minPointsInRange = std::max(min_points_in_range, 1);
  const double initialSdRange = std::max(1.0, sdRangeMin);
  const double safeAmp0 = std::max(amp0, 1.e-9);
  const double safeSigma0Full = std::max(sigma0, hist->GetBinWidth(maxBin));
  const double initialSigmaMin = std::max(axMin, mean0 - initialSdRange * safeSigma0Full);
  const double initialSigmaMax = std::min(axMax, mean0 + initialSdRange * safeSigma0Full);
  const auto [mean0Refined, sigma0Refined] =
      GFMeanStdDevInRange(hist, initialSigmaMin, initialSigmaMax,
                          mean0, safeSigma0Full);
  const double safeMean0 = std::isfinite(initial_mean_override)
                          ? initial_mean_override : mean0Refined;
  const double safeSigma0 = std::max(sigma0Refined, hist->GetBinWidth(maxBin));

  double currentAmp = safeAmp0;
  double currentMean = safeMean0;
  double currentSigma = safeSigma0;
  double currentSdRange = initialSdRange;

  std::vector<double> iterValues;
  std::vector<double> ampValues;
  std::vector<double> meanValues;
  std::vector<double> sigmaValues;
  std::vector<double> chi2Values;
  std::vector<double> sdRangeValues;
  iterValues.push_back(0.0);
  ampValues.push_back(currentAmp);
  meanValues.push_back(currentMean);
  sigmaValues.push_back(currentSigma);
  chi2Values.push_back(0.0);
  sdRangeValues.push_back(currentSdRange);

  TF1 *bestFit = nullptr;
  double bestFitMin = std::max(axMin, currentMean - currentSdRange * currentSigma);
  double bestFitMax = std::min(axMax, currentMean + currentSdRange * currentSigma);
  double bestChi2ndf = std::numeric_limits<double>::infinity();

  constexpr int kMaxIterations = 10;
  for (int iter = 1; iter <= kMaxIterations; ++iter) {
    const double iterAmp = currentAmp;
    const double iterMean = currentMean;
    const double iterSigma = std::max(currentSigma, 1.e-12);

    const double candidates[3] = {
      GFClamp(0.8 * currentSdRange, sdRangeMin, sdRangeMax),
      GFClamp(1.0 * currentSdRange, sdRangeMin, sdRangeMax),
      GFClamp(1.6 * currentSdRange, sdRangeMin, sdRangeMax)
    };

    TF1 *iterBestFit = nullptr;
    double iterBestChi2ndf = std::numeric_limits<double>::infinity();
    double iterBestSdRange = currentSdRange;
    double iterBestFitMin = bestFitMin;
    double iterBestFitMax = bestFitMax;

    for (double sdRange : candidates) {
      const double fitMinTry = std::max(axMin, iterMean - sdRange * iterSigma);
      const double fitMaxTry = std::min(axMax, iterMean + sdRange * iterSigma);
      if (fitMaxTry <= fitMinTry)
        continue;
      if (GFCountBinsInRange(hist, fitMinTry, fitMaxTry) < minPointsInRange)
        continue;

      auto fitTry = new TF1(GFUniqueName("gf_try", hist), "gaus",
                            fitMinTry, fitMaxTry);
      fitTry->SetParameters(iterAmp, iterMean, iterSigma);
      fitTry->SetParNames("Amplitude", "Mean", "Sigma");
      fitTry->SetParLimits(0, 0, 10.0 * std::max(iterAmp, 1.e-9));
      fitTry->SetParLimits(1, fitMinTry, fitMaxTry);
      fitTry->SetParLimits(2, 0.2 * std::max(iterSigma, 1.e-9),
                              5.0 * std::max(iterSigma, 1.e-9));
      fitTry->SetNpx(500);
      fitTry->SetLineColor(kRed + 1);
      fitTry->SetLineWidth(2);
      hist->Fit(fitTry, "RQN0");

      const double chi2ndfTry = fitTry->GetNDF() > 0
          ? fitTry->GetChisquare() / fitTry->GetNDF()
          : std::numeric_limits<double>::infinity();

      if (chi2ndfTry < iterBestChi2ndf) {
        if (iterBestFit != nullptr)
          delete iterBestFit;
        iterBestFit = fitTry;
        iterBestChi2ndf = chi2ndfTry;
        iterBestSdRange = sdRange;
        iterBestFitMin = fitMinTry;
        iterBestFitMax = fitMaxTry;
      }
      else {
        delete fitTry;
      }
    }

    if (iterBestFit == nullptr)
      break;

    if (bestFit != nullptr)
      delete bestFit;
    bestFit = iterBestFit;
    bestChi2ndf = iterBestChi2ndf;
    bestFitMin = iterBestFitMin;
    bestFitMax = iterBestFitMax;
    currentSdRange = iterBestSdRange;
    currentAmp = bestFit->GetParameter(0);
    currentMean = bestFit->GetParameter(1);
    currentSigma = std::abs(bestFit->GetParameter(2));

    iterValues.push_back((double)iter);
    ampValues.push_back(currentAmp);
    meanValues.push_back(currentMean);
    sigmaValues.push_back(currentSigma);
    chi2Values.push_back(bestChi2ndf);
    sdRangeValues.push_back(currentSdRange);

    const bool converged =
        GFWithinRelativeChange(iterAmp, currentAmp) &&
        GFWithinRelativeChange(iterMean, currentMean) &&
        GFWithinRelativeChange(iterSigma, currentSigma);
    if (converged)
      break;
  }

  if (bestFit == nullptr) {
    bestFitMin = std::max(axMin, safeMean0 - currentSdRange * safeSigma0);
    bestFitMax = std::min(axMax, safeMean0 + currentSdRange * safeSigma0);
    bestFit = new TF1(GFUniqueName("gf_result", hist), "gaus", bestFitMin, bestFitMax);
    bestFit->SetParameters(safeAmp0, safeMean0, safeSigma0);
    bestFit->SetParNames("Amplitude", "Mean", "Sigma");
    bestFit->SetParLimits(0, 0, 10.0 * safeAmp0);
    bestFit->SetParLimits(1, bestFitMin, bestFitMax);
    bestFit->SetParLimits(2, 0.2 * safeSigma0, 5.0 * safeSigma0);
    bestFit->SetNpx(500);
    bestFit->SetLineColor(kRed + 1);
    bestFit->SetLineWidth(2);
    hist->Fit(bestFit, "RQN0");
    currentAmp = bestFit->GetParameter(0);
    currentMean = bestFit->GetParameter(1);
    currentSigma = std::abs(bestFit->GetParameter(2));
    bestChi2ndf = bestFit->GetNDF() > 0
        ? bestFit->GetChisquare() / bestFit->GetNDF() : 0.0;

    iterValues.push_back(1.0);
    ampValues.push_back(currentAmp);
    meanValues.push_back(currentMean);
    sigmaValues.push_back(currentSigma);
    chi2Values.push_back(bestChi2ndf);
    sdRangeValues.push_back(currentSdRange);
  }
  bestFit->SetName(GFUniqueName("gf_result", hist));

  // ── Initial Gaussian (dashed) ─────────────────────────────────────────────
  auto f_init = new TF1(GFUniqueName("gf_init", hist), "gaus", bestFitMin, bestFitMax);
  f_init->SetParameters(safeAmp0, safeMean0, safeSigma0);
  f_init->SetNpx(500);
  f_init->SetLineColor(kGray + 1);
  f_init->SetLineStyle(2);
  f_init->SetLineWidth(2);

  // ── Fit range vertical lines ──────────────────────────────────────────────
  const double yMax = safeAmp0;
  auto lMin = new TLine(bestFitMin, 0, bestFitMin, yMax * 0.90);
  auto lMax = new TLine(bestFitMax, 0, bestFitMax, yMax * 0.90);
  for (auto l : {lMin, lMax}) {
    l->SetLineColor(kGray + 2);
    l->SetLineStyle(2);
    l->SetLineWidth(1);
  }

  auto f_result = bestFit;

  const double mu      = f_result->GetParameter(1);
  const double sigma   = std::abs(f_result->GetParameter(2));
  const double chi2ndf = f_result->GetNDF() > 0
                       ? f_result->GetChisquare() / f_result->GetNDF() : 0;

  // ── Set x-axis range using fit result ────────────────────────────────────
  double drawMin = axMin;
  double drawMax = axMax;
  if (sigma > 0) {
    drawMin = std::max(axMin, mu - 5.0 * sigma);
    drawMax = std::min(axMax, mu + 5.0 * sigma);
    hist->GetXaxis()->SetRangeUser(drawMin, drawMax);
  }

  // The fit is performed only in [bestFitMin, bestFitMax], but draw the
  // resulting Gaussian across the visible histogram range so the extrapolated
  // shape is visible without spreading Npx over the full raw histogram axis.
  f_result->SetRange(drawMin, drawMax);

  // ── Manual TLegend ────────────────────────────────────────────────────────
  TLegend *legend = nullptr;
  if (add_legend) {
    legend = new TLegend(0.12, 0.62, 0.65, 0.88);
    legend->SetBorderSize(0);
    legend->SetFillStyle(0);
    legend->SetTextSize(0.034);
    legend->SetMargin(0.12);   // narrow icon column for text-only entries
    char line[256];

    // initial params — comment out block below to hide from legend
    //if (add_f_init) {
    //  std::snprintf(line, sizeof(line),
    //                "Initial:  #mu_{0}=%.3f  #sigma_{0}=%.3f", mean0, sigma0);
    //  legend->AddEntry(f_init, line, "l");
    //}
    // fit range — comment out block below to hide from legend
    //if (add_range) {
    //  std::snprintf(line, sizeof(line),
    //                "Fit range: [%.3f, %.3f]  (#pm3#sigma_{0})", fitMin, fitMax);
    //  legend->AddEntry(lMin, line, "l");
    //}

    legend->AddEntry(f_result, "Gaussian fit", "l");
    std::snprintf(line, sizeof(line), "#mu = %.4f", mu);
    legend->AddEntry((TObject*)nullptr, line, "");
    std::snprintf(line, sizeof(line), "#sigma = %.4f", sigma);
    legend->AddEntry((TObject*)nullptr, line, "");
    if (add_iteration_info) {
      std::snprintf(line, sizeof(line), "iterations = %d", (int)iterValues.size() - 1);
      legend->AddEntry((TObject*)nullptr, line, "");
      std::snprintf(line, sizeof(line), "sd_range = %.4g", currentSdRange);
      legend->AddEntry((TObject*)nullptr, line, "");
      std::snprintf(line, sizeof(line), "sd_range_{min} = %.4g", sdRangeMin);
      legend->AddEntry((TObject*)nullptr, line, "");
      std::snprintf(line, sizeof(line), "sd_range_{max} = %.4g", sdRangeMax);
      legend->AddEntry((TObject*)nullptr, line, "");
      std::snprintf(line, sizeof(line), "min points = %d", minPointsInRange);
      legend->AddEntry((TObject*)nullptr, line, "");
    }
    std::snprintf(line, sizeof(line), "#chi^{2}/NDF = %.3g/%d = %.3g",
                  f_result->GetChisquare(), f_result->GetNDF(), chi2ndf);
    legend->AddEntry((TObject*)nullptr, line, "");
  }

  // ── Assemble LKDrawing ────────────────────────────────────────────────────
  auto drawing = new LKDrawing(GFUniqueName("gfdraw", hist));

  if (add_hist)     drawing->Add(hist,     hist_option);
  if (add_f_init)   drawing->Add(f_init,   "samel");
  if (add_range)  { drawing->Add(lMin,     "same");
                    drawing->Add(lMax,     "same"); }
  if (add_f_result) drawing->Add(f_result, "samel");
  if (add_legend)   drawing->Add(legend,   "same");

  auto gAmp = GFMakeGraph("gf_g_amp", hist, iterValues, ampValues, kRed + 1);
  auto gMean = GFMakeGraph("gf_g_mean", hist, iterValues, meanValues, kBlue + 1);
  auto gSigma = GFMakeGraph("gf_g_sigma", hist, iterValues, sigmaValues, kGreen + 2);
  auto gChi2 = GFMakeGraph("gf_g_chi2ndf", hist, iterValues, chi2Values, kMagenta + 1);
  auto gSdRange = GFMakeGraph("gf_g_sd_range", hist, iterValues, sdRangeValues, kOrange + 7);

  return {drawing, f_init, f_result, mu, sigma, chi2ndf, currentSdRange,
          gAmp, gMean, gSigma, gChi2, gSdRange};
}

// ── Public: returns LKDrawing* ────────────────────────────────────────────────
inline LKDrawing *ProceedGaussianFitting(
    TH1D       *hist,
    bool        add_hist     = true,
    bool        add_f_init   = true,
    bool        add_range    = true,
    bool        add_f_result = true,
    bool        add_legend   = true,
    const char *hist_option  = "hist",
    double      sd_range_min = 1.5,
    int         min_points_in_range = 10,
    double      sd_range_max = 3.0,
    bool        add_iteration_info = false,
    double      initial_mean_override = std::numeric_limits<double>::quiet_NaN())
{
  return _GaussianFitImpl(hist, add_hist, add_f_init, add_range,
                          add_f_result, add_legend, hist_option,
                          sd_range_min, min_points_in_range,
                          sd_range_max, add_iteration_info,
                          initial_mean_override).drawing;
}

// ── Public: returns full result struct ────────────────────────────────────────
inline GaussianFitResult ProceedGaussianFittingFull(
    TH1D       *hist,
    bool        add_hist     = true,
    bool        add_f_init   = true,
    bool        add_range    = true,
    bool        add_f_result = true,
    bool        add_legend   = true,
    const char *hist_option  = "hist",
    double      sd_range_min = 1.5,
    int         min_points_in_range = 10,
    double      sd_range_max = 3.0,
    bool        add_iteration_info = false,
    double      initial_mean_override = std::numeric_limits<double>::quiet_NaN())
{
  return _GaussianFitImpl(hist, add_hist, add_f_init, add_range,
                          add_f_result, add_legend, hist_option,
                          sd_range_min, min_points_in_range,
                          sd_range_max, add_iteration_info,
                          initial_mean_override);
}

#endif // GAUSSIAN_FITTER_H
