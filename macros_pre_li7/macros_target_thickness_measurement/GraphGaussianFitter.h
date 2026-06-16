#ifndef GRAPH_GAUSSIAN_FITTER_H
#define GRAPH_GAUSSIAN_FITTER_H

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <limits>
#include <vector>

#include "TF1.h"
#include "TGraphErrors.h"
#include "TLegend.h"
#include "TLine.h"

#include "GaussianFitter.h"
#include "LKDrawing.h"

inline TString GGFUniqueName(const char *prefix, const TH1D *hist)
{
  static int counter = 0;
  return TString::Format("%s_%s_%d", prefix, hist->GetName(), ++counter);
}

inline TString GGFUniqueName(const char *prefix, const char *name)
{
  static int counter = 0;
  return TString::Format("%s_%s_%d", prefix, name, ++counter);
}

struct GraphStats {
  double amp = 0;
  double mean = 0;
  double sigma = 0;
};

inline GraphStats GetGraphStatsFromHist(TH1D *hist, double xMin, double xMax)
{
  GraphStats stats;
  double sumW = 0;
  double sumWX = 0;
  double sumWX2 = 0;

  for (int i = 1; i <= hist->GetNbinsX(); ++i) {
    const double x = hist->GetXaxis()->GetBinLowEdge(i);
    if (x < xMin || x > xMax)
      continue;
    const double y = hist->GetBinContent(i);
    if (y <= 0)
      continue;
    stats.amp = std::max(stats.amp, y);
    sumW += y;
    sumWX += y * x;
    sumWX2 += y * x * x;
  }

  if (sumW <= 0)
    return stats;

  stats.mean = sumWX / sumW;
  const double variance = sumWX2 / sumW - stats.mean * stats.mean;
  stats.sigma = variance > 0 ? std::sqrt(variance) : hist->GetBinWidth(1);
  return stats;
}

inline GraphStats GetGraphStatsFromPoints(const std::vector<double> &xs,
                                          const std::vector<double> &ys,
                                          double xMin, double xMax)
{
  GraphStats stats;
  double sumW = 0;
  double sumWX = 0;
  double sumWX2 = 0;

  for (int i = 0; i < (int)xs.size(); ++i) {
    const double x = xs[i];
    if (x < xMin || x > xMax)
      continue;
    const double y = ys[i];
    if (y <= 0)
      continue;
    stats.amp = std::max(stats.amp, y);
    sumW += y;
    sumWX += y * x;
    sumWX2 += y * x * x;
  }

  if (sumW <= 0)
    return stats;

  stats.mean = sumWX / sumW;
  const double variance = sumWX2 / sumW - stats.mean * stats.mean;
  stats.sigma = variance > 0 ? std::sqrt(variance) : 1.0;
  return stats;
}

inline TGraphErrors *MakeChannelGraphFromHist(TH1D *hist, double xMin, double xMax)
{
  int n = 0;
  for (int i = 1; i <= hist->GetNbinsX(); ++i) {
    const double x = hist->GetXaxis()->GetBinLowEdge(i);
    const double y = hist->GetBinContent(i);
    if (x >= xMin && x <= xMax && y > 0)
      ++n;
  }

  auto graph = new TGraphErrors(n);
  graph->SetName(GGFUniqueName("g_channel", hist));
  graph->SetTitle(TString::Format("%s;%s;%s", hist->GetTitle(),
                                  hist->GetXaxis()->GetTitle(),
                                  hist->GetYaxis()->GetTitle()));
  graph->SetMarkerStyle(20);
  graph->SetMarkerSize(0.45);
  graph->SetMarkerColor(hist->GetLineColor());
  graph->SetLineColor(hist->GetLineColor());

  int ip = 0;
  for (int i = 1; i <= hist->GetNbinsX(); ++i) {
    const double x = hist->GetXaxis()->GetBinLowEdge(i);
    const double y = hist->GetBinContent(i);
    if (x < xMin || x > xMax || y <= 0)
      continue;
    graph->SetPoint(ip, x, y);
    graph->SetPointError(ip, 0, std::sqrt(std::max(y, 1.0)));
    ++ip;
  }
  return graph;
}

inline TGraphErrors *MakeGraphFromPoints(const std::vector<double> &xs,
                                         const std::vector<double> &ys,
                                         const char *name,
                                         const char *title,
                                         const char *xTitle,
                                         const char *yTitle,
                                         double xMin,
                                         double xMax)
{
  int n = 0;
  for (int i = 0; i < (int)xs.size(); ++i) {
    if (xs[i] >= xMin && xs[i] <= xMax && ys[i] > 0)
      ++n;
  }

  auto graph = new TGraphErrors(n);
  graph->SetName(GGFUniqueName("g_points", name));
  graph->SetTitle(TString::Format("%s;%s;%s", title, xTitle, yTitle));
  graph->SetMarkerStyle(20);
  graph->SetMarkerSize(0.45);
  graph->SetMarkerColor(kBlack);
  graph->SetLineColor(kBlack);

  int ip = 0;
  for (int i = 0; i < (int)xs.size(); ++i) {
    const double x = xs[i];
    const double y = ys[i];
    if (x < xMin || x > xMax || y <= 0)
      continue;
    graph->SetPoint(ip, x, y);
    graph->SetPointError(ip, 0, std::sqrt(std::max(y, 1.0)));
    ++ip;
  }
  return graph;
}

inline GaussianFitResult ProceedGraphGaussianFittingFull(
    const char *name,
    const char *title,
    const char *xTitle,
    const char *yTitle,
    const std::vector<double> &xs,
    const std::vector<double> &ys,
    bool        add_graph = true,
    bool        add_f_init = false,
    bool        add_range = false,
    bool        add_f_result = true,
    bool        add_legend = true,
    const char *graph_option = "AP",
    double      sd_range_min = 1.5,
    int         min_points_in_range = 10,
    double      sd_range_max = 3.0,
    bool        add_iteration_info = false,
    double      initial_mean_override = std::numeric_limits<double>::quiet_NaN())
{
  if (xs.empty() || xs.size() != ys.size())
    return {};

  auto [minIt, maxIt] = std::minmax_element(xs.begin(), xs.end());
  const double axMin = *minIt;
  const double axMax = *maxIt;

  GraphStats stats0 = GetGraphStatsFromPoints(xs, ys, axMin, axMax);
  if (stats0.sigma <= 0)
    stats0.sigma = std::max((axMax - axMin) / 100.0, 1.e-9);

  const double sdRangeMax = std::max(sd_range_max, 0.0);
  const double sdRangeMin = std::min(GFClamp(sd_range_min, 0.0, sdRangeMax),
                                     sdRangeMax);
  const double sdRange = sdRangeMax > 0 ? sdRangeMax : sdRangeMin;
  const double fitMin0 = std::max(axMin, stats0.mean - sdRange * stats0.sigma);
  const double fitMax0 = std::min(axMax, stats0.mean + sdRange * stats0.sigma);

  GraphStats stats = GetGraphStatsFromPoints(xs, ys, fitMin0, fitMax0);
  if (stats.amp <= 0) stats.amp = stats0.amp;
  if (stats.mean <= 0) stats.mean = stats0.mean;
  if (stats.sigma <= 0) stats.sigma = stats0.sigma;
  if (std::isfinite(initial_mean_override))
    stats.mean = initial_mean_override;

  const double fitMin = std::max(axMin, stats.mean - sdRange * stats.sigma);
  const double fitMax = std::min(axMax, stats.mean + sdRange * stats.sigma);

  auto graph = MakeGraphFromPoints(xs, ys, name, title, xTitle, yTitle, axMin, axMax);

  auto f_init = new TF1(GGFUniqueName("ggf_init", name), "gaus", fitMin, fitMax);
  f_init->SetParameters(stats.amp, stats.mean, stats.sigma);
  f_init->SetNpx(500);
  f_init->SetLineColor(kGray + 2);
  f_init->SetLineStyle(2);
  f_init->SetLineWidth(2);

  auto f_result = new TF1(GGFUniqueName("ggf_fit", name), "gaus", fitMin, fitMax);
  f_result->SetParameters(stats.amp, stats.mean, stats.sigma);
  f_result->SetParLimits(0, 0, std::max(stats.amp * 5.0, 1.0));
  f_result->SetParLimits(1, fitMin, fitMax);
  f_result->SetParLimits(2, std::max(stats.sigma * 0.2, 1.e-9),
                         std::max(stats.sigma * 5.0, 1.e-9));
  f_result->SetNpx(500);
  f_result->SetLineColor(kRed + 1);
  f_result->SetLineWidth(2);

  graph->Fit(f_result, "RQN0");
  graph->Fit(f_result, "RQN0");

  const double mu = f_result->GetParameter(1);
  const double sigma = std::abs(f_result->GetParameter(2));
  const double chi2ndf = f_result->GetNDF() > 0
                       ? f_result->GetChisquare() / f_result->GetNDF() : 0.0;

  double drawMin = axMin;
  double drawMax = axMax;
  if (sigma > 0) {
    drawMin = std::max(axMin, mu - 5.0 * sigma);
    drawMax = std::min(axMax, mu + 5.0 * sigma);
  }
  graph->GetXaxis()->SetLimits(drawMin, drawMax);
  f_result->SetRange(drawMin, drawMax);

  auto lMin = new TLine(fitMin, 0, fitMin, stats.amp);
  auto lMax = new TLine(fitMax, 0, fitMax, stats.amp);
  lMin->SetLineColor(kBlue + 1);
  lMax->SetLineColor(kBlue + 1);
  lMin->SetLineStyle(2);
  lMax->SetLineStyle(2);

  TLegend *legend = nullptr;
  if (add_legend) {
    legend = new TLegend(0.12, 0.62, 0.65, 0.88);
    legend->SetBorderSize(0);
    legend->SetFillStyle(0);
    legend->SetTextSize(0.034);
    legend->SetMargin(0.12);
    char line[256];
    legend->AddEntry(graph, "Graph", "p");
    legend->AddEntry(f_result, "Gaussian fit", "l");
    std::snprintf(line, sizeof(line), "#mu = %.4f", mu);
    legend->AddEntry((TObject*)nullptr, line, "");
    std::snprintf(line, sizeof(line), "#sigma = %.4f", sigma);
    legend->AddEntry((TObject*)nullptr, line, "");
    if (add_iteration_info) {
      std::snprintf(line, sizeof(line), "sd_range = %.4g", sdRange);
      legend->AddEntry((TObject*)nullptr, line, "");
      std::snprintf(line, sizeof(line), "min points = %d", min_points_in_range);
      legend->AddEntry((TObject*)nullptr, line, "");
    }
    std::snprintf(line, sizeof(line), "#chi^{2}/NDF = %.3g/%d = %.3g",
                  f_result->GetChisquare(), f_result->GetNDF(), chi2ndf);
    legend->AddEntry((TObject*)nullptr, line, "");
  }

  auto drawing = new LKDrawing(GGFUniqueName("ggfdraw", name));
  if (add_graph) drawing->Add(graph, graph_option);
  if (add_f_init) drawing->Add(f_init, "l same");
  if (add_range) {
    drawing->Add(lMin, "same");
    drawing->Add(lMax, "same");
  }
  if (add_f_result) drawing->Add(f_result, "l same");
  if (add_legend) drawing->Add(legend, "same");

  return {drawing, f_init, f_result, mu, sigma, chi2ndf, sdRange,
          nullptr, nullptr, nullptr, nullptr, nullptr};
}

inline GaussianFitResult ProceedGraphGaussianFittingFull(
    TH1D       *hist,
    bool        add_graph = true,
    bool        add_f_init = false,
    bool        add_range = false,
    bool        add_f_result = true,
    bool        add_legend = true,
    const char *graph_option = "AP",
    double      sd_range_min = 1.5,
    int         min_points_in_range = 10,
    double      sd_range_max = 3.0,
    bool        add_iteration_info = false)
{
  const double axMin = hist->GetXaxis()->GetXmin();
  const double axMax = hist->GetXaxis()->GetXmax();
  int first = hist->GetXaxis()->GetFirst();
  int last  = hist->GetXaxis()->GetLast();
  if (first < 1) first = 1;
  if (last < first || last > hist->GetNbinsX()) last = hist->GetNbinsX();
  const double visibleMin = std::max(axMin, hist->GetXaxis()->GetBinLowEdge(first));
  const double visibleMax = std::min(axMax, hist->GetXaxis()->GetBinUpEdge(last));

  GraphStats stats0 = GetGraphStatsFromHist(hist, visibleMin, visibleMax);
  if (stats0.sigma <= 0)
    stats0.sigma = hist->GetBinWidth(1);

  const double sdRangeMax = std::max(sd_range_max, 0.0);
  const double sdRangeMin = std::min(GFClamp(sd_range_min, 0.0, sdRangeMax),
                                     sdRangeMax);
  const double sdRange = sdRangeMax > 0 ? sdRangeMax : sdRangeMin;
  const double fitMin0 = std::max(visibleMin, stats0.mean - sdRange * stats0.sigma);
  const double fitMax0 = std::min(visibleMax, stats0.mean + sdRange * stats0.sigma);

  GraphStats stats = GetGraphStatsFromHist(hist, fitMin0, fitMax0);
  if (stats.amp <= 0) stats.amp = stats0.amp;
  if (stats.mean <= 0) stats.mean = stats0.mean;
  if (stats.sigma <= 0) stats.sigma = stats0.sigma;

  const double fitMin = std::max(visibleMin, stats.mean - sdRange * stats.sigma);
  const double fitMax = std::min(visibleMax, stats.mean + sdRange * stats.sigma);

  auto graph = MakeChannelGraphFromHist(hist, visibleMin, visibleMax);
  auto f_init = new TF1(GGFUniqueName("ggf_init", hist), "gaus", fitMin, fitMax);
  f_init->SetParameters(stats.amp, stats.mean, stats.sigma);
  f_init->SetNpx(500);
  f_init->SetLineColor(kGray + 2);
  f_init->SetLineStyle(2);
  f_init->SetLineWidth(2);

  auto f_result = new TF1(GGFUniqueName("ggf_fit", hist), "gaus", fitMin, fitMax);
  f_result->SetParameters(stats.amp, stats.mean, stats.sigma);
  f_result->SetParLimits(0, 0, std::max(stats.amp * 5.0, 1.0));
  f_result->SetParLimits(1, fitMin, fitMax);
  f_result->SetParLimits(2, std::max(stats.sigma * 0.2, 1.e-9),
                         std::max(stats.sigma * 5.0, 1.e-9));
  f_result->SetNpx(500);
  f_result->SetLineColor(kRed + 1);
  f_result->SetLineWidth(2);

  graph->Fit(f_result, "RQN0");
  graph->Fit(f_result, "RQN0");

  const double mu = f_result->GetParameter(1);
  const double sigma = std::abs(f_result->GetParameter(2));
  const double chi2ndf = f_result->GetNDF() > 0
                       ? f_result->GetChisquare() / f_result->GetNDF() : 0.0;

  double drawMin = visibleMin;
  double drawMax = visibleMax;
  if (sigma > 0) {
    drawMin = std::max(axMin, mu - 5.0 * sigma);
    drawMax = std::min(axMax, mu + 5.0 * sigma);
  }
  graph->GetXaxis()->SetLimits(drawMin, drawMax);
  f_result->SetRange(drawMin, drawMax);

  auto lMin = new TLine(fitMin, 0, fitMin, stats.amp);
  auto lMax = new TLine(fitMax, 0, fitMax, stats.amp);
  lMin->SetLineColor(kBlue + 1);
  lMax->SetLineColor(kBlue + 1);
  lMin->SetLineStyle(2);
  lMax->SetLineStyle(2);

  TLegend *legend = nullptr;
  if (add_legend) {
    legend = new TLegend(0.12, 0.62, 0.65, 0.88);
    legend->SetBorderSize(0);
    legend->SetFillStyle(0);
    legend->SetTextSize(0.034);
    legend->SetMargin(0.12);
    char line[256];
    legend->AddEntry(graph, "Graph", "p");
    legend->AddEntry(f_result, "Gaussian fit", "l");
    std::snprintf(line, sizeof(line), "#mu = %.4f", mu);
    legend->AddEntry((TObject*)nullptr, line, "");
    std::snprintf(line, sizeof(line), "#sigma = %.4f", sigma);
    legend->AddEntry((TObject*)nullptr, line, "");
    if (add_iteration_info) {
      std::snprintf(line, sizeof(line), "sd_range = %.4g", sdRange);
      legend->AddEntry((TObject*)nullptr, line, "");
      std::snprintf(line, sizeof(line), "min points = %d", min_points_in_range);
      legend->AddEntry((TObject*)nullptr, line, "");
    }
    std::snprintf(line, sizeof(line), "#chi^{2}/NDF = %.3g/%d = %.3g",
                  f_result->GetChisquare(), f_result->GetNDF(), chi2ndf);
    legend->AddEntry((TObject*)nullptr, line, "");
  }

  auto drawing = new LKDrawing(GGFUniqueName("ggfdraw", hist));
  if (add_graph) drawing->Add(graph, graph_option);
  if (add_f_init) drawing->Add(f_init, "l same");
  if (add_range) {
    drawing->Add(lMin, "same");
    drawing->Add(lMax, "same");
  }
  if (add_f_result) drawing->Add(f_result, "l same");
  if (add_legend) drawing->Add(legend, "same");

  return {drawing, f_init, f_result, mu, sigma, chi2ndf, sdRange,
          nullptr, nullptr, nullptr, nullptr, nullptr};
}

inline LKDrawing *ProceedGraphGaussianFitting(
    TH1D       *hist,
    bool        add_graph = true,
    bool        add_f_init = false,
    bool        add_range = false,
    bool        add_f_result = true,
    bool        add_legend = true,
    const char *graph_option = "AP",
    double      sd_range_min = 1.5,
    int         min_points_in_range = 10,
    double      sd_range_max = 3.0,
    bool        add_iteration_info = false)
{
  return ProceedGraphGaussianFittingFull(hist, add_graph, add_f_init, add_range,
                                         add_f_result, add_legend, graph_option,
                                         sd_range_min, min_points_in_range,
                                         sd_range_max, add_iteration_info).drawing;
}

#endif // GRAPH_GAUSSIAN_FITTER_H
