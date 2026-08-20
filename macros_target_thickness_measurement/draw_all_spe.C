#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <map>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

#include "TFile.h"
#include "TF1.h"
#include "TGraphErrors.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TSystem.h"

namespace {

constexpr int kHeaderLines = 12;
constexpr int kNChannels = 4096;

std::string Trim(const std::string &text)
{
  const auto first = text.find_first_not_of(" \t\r\n");
  if (first == std::string::npos)
    return "";
  const auto last = text.find_last_not_of(" \t\r\n");
  return text.substr(first, last - first + 1);
}

std::string FirstColumn(const std::string &line)
{
  const auto tab = line.find('\t');
  return Trim(tab == std::string::npos ? line : line.substr(0, tab));
}

std::string HistNameFromPath(std::string path)
{
  std::replace(path.begin(), path.end(), '/', '_');
  std::replace(path.begin(), path.end(), '.', '_');
  std::replace(path.begin(), path.end(), '-', '_');
  return "h_" + path;
}

std::string FitNameFromPath(const std::string &path)
{
  auto name = HistNameFromPath(path);
  name.replace(0, 2, "f_");
  return name;
}

std::string EnergyHistNameFromPath(const std::string &path)
{
  auto name = HistNameFromPath(path);
  name.replace(0, 2, "he_");
  return name;
}

std::string EnergyFitNameFromPath(const std::string &path)
{
  auto name = HistNameFromPath(path);
  name.replace(0, 2, "fe_");
  return name;
}

std::string DateFromPath(const std::string &path)
{
  const std::string prefix = "data/";
  const auto begin = path.find(prefix);
  if (begin == std::string::npos)
    return "";

  const auto dateBegin = begin + prefix.size();
  const auto dateEnd = path.find('/', dateBegin);
  if (dateEnd == std::string::npos)
    return "";

  return path.substr(dateBegin, dateEnd - dateBegin);
}

bool IsPulserPath(const std::string &path)
{
  return path.find("/Pulser/") != std::string::npos ||
         path.find("/pulser/") != std::string::npos;
}

bool IsBlankPath(const std::string &path)
{
  return path.find("/Blank/") != std::string::npos ||
         path.find("/blank/") != std::string::npos ||
         path.find("Blank") != std::string::npos ||
         path.find("blank") != std::string::npos;
}

bool IsOldPath(const std::string &path)
{
  return path.find("/old/") != std::string::npos ||
         path.rfind("old/", 0) == 0;
}

bool GetPulserVamp(const std::string &path, double &vamp)
{
  const auto slash = path.find_last_of('/');
  const auto name = slash == std::string::npos ? path : path.substr(slash + 1);
  const auto dash = name.find('-');
  if (dash == std::string::npos)
    return false;

  const char *begin = name.c_str() + dash + 1;
  char *end = nullptr;
  const double value = std::strtod(begin, &end);
  if (end == begin)
    return false;

  std::string unit(end, 2);
  if (unit == "mV") {
    vamp = value * 1.e-3;
    return true;
  }
  if (unit[0] == 'V') {
    vamp = value;
    return true;
  }

  return false;
}

bool ShouldSkipLowPulserChannel(const std::string &path, int channel)
{
  double vamp = 0;
  if (!IsPulserPath(path) || !GetPulserVamp(path, vamp))
    return false;

  return channel < 50 && (std::abs(vamp - 1.0) < 1.e-9 ||
                          std::abs(vamp - 1.6) < 1.e-9);
}

bool ShouldSkipLowTargetChannel(const std::string &path, int channel)
{
  return channel < 300 && !IsBlankPath(path) && !IsPulserPath(path);
}

bool UseGaussianOnlyFit(const std::string &path)
{
  return IsBlankPath(path) ||
         path.find("/ch2/") != std::string::npos ||
         path.find("/Carbon") != std::string::npos ||
         path.find("/carbon") != std::string::npos ||
         path.find("/5.8.2/") != std::string::npos ||
         path.find("/5_8_2/") != std::string::npos;
}

struct PeakInitialParameters {
  double amplitude = 0;
  double mean = 0;
  double sigma = 0;
};

double EstimatePeakSigma(TH1D *hist)
{
  const int maxBin = hist->GetMaximumBin();
  const double halfMax = 0.5 * hist->GetBinContent(maxBin);
  if (halfMax <= 0)
    return hist->GetStdDev();

  int leftBin = maxBin;
  while (leftBin > 1 && hist->GetBinContent(leftBin) > halfMax)
    --leftBin;

  int rightBin = maxBin;
  while (rightBin < hist->GetNbinsX() && hist->GetBinContent(rightBin) > halfMax)
    ++rightBin;

  const double fwhm = hist->GetBinCenter(rightBin) - hist->GetBinCenter(leftBin);
  if (fwhm > 0)
    return fwhm / 2.35482;

  return hist->GetStdDev();
}

PeakInitialParameters GetPeakInitialParameters(TH1D *hist)
{
  const int maxBin = hist->GetMaximumBin();
  PeakInitialParameters parameters;
  parameters.amplitude = hist->GetBinContent(maxBin);
  parameters.mean = hist->GetBinCenter(maxBin);
  parameters.sigma = EstimatePeakSigma(hist);
  return parameters;
}

double LeftTailCrystalBall(double *x, double *p)
{
  const double amplitude = p[0];
  const double mean = p[1];
  const double sigma = std::abs(p[2]);
  const double alpha = std::abs(p[3]);
  const double n = std::abs(p[4]);
  if (sigma <= 0 || alpha <= 0 || n <= 0)
    return 0;

  const double t = (x[0] - mean) / sigma;
  if (t > -alpha)
    return amplitude * std::exp(-0.5 * t * t);

  const double a = std::pow(n / alpha, n) * std::exp(-0.5 * alpha * alpha);
  const double b = n / alpha - alpha;
  return amplitude * a * std::pow(b - t, -n);
}

TF1 *CreatePeakFit(const std::string &name, const std::string &path,
                   double fitMin, double fitMax, double amplitude,
                   double mean, double sigma)
{
  const double safeAmplitude = std::max(amplitude, 1.e-9);
  const double safeSigma = std::max(sigma, 1.e-9);
  TF1 *fit = nullptr;
  if (UseGaussianOnlyFit(path)) {
    fit = new TF1(name.c_str(), "gaus", fitMin, fitMax);
    fit->SetParameters(safeAmplitude, mean, safeSigma);
    fit->SetParNames("Amplitude", "Mean", "Sigma");
    fit->SetParLimits(0, 0, 10.0 * safeAmplitude);
    fit->SetParLimits(1, fitMin, fitMax);
    fit->SetParLimits(2, 0.2 * safeSigma, 5.0 * safeSigma);
  }
  else {
    fit = new TF1(name.c_str(), LeftTailCrystalBall, fitMin, fitMax, 5);
    fit->SetParameters(safeAmplitude, mean, safeSigma, 0.5, 3.0);
    fit->SetParNames("Amplitude", "Mean", "Sigma", "Alpha", "N");
    fit->SetParLimits(0, 0, 10.0 * safeAmplitude);
    fit->SetParLimits(1, fitMin, fitMax);
    fit->SetParLimits(2, 0.2 * safeSigma, 5.0 * safeSigma);
    fit->SetParLimits(3, 0.05, 8.0);
    fit->SetParLimits(4, 1.01, 30.0);
  }

  fit->SetNpx(500);
  fit->SetLineWidth(2);
  return fit;
}

TF1 *CreateInitializedPeakFit(const std::string &name, const std::string &path,
                              TH1D *hist, double fitMin, double fitMax,
                              double amplitude, double mean, double sigma)
{
  if (UseGaussianOnlyFit(path))
    return CreatePeakFit(name, path, fitMin, fitMax, amplitude, mean, sigma);

  TF1 prefit((name + "_prefit").c_str(), "gaus", fitMin, fitMax);
  prefit.SetParameters(amplitude, mean, sigma);
  prefit.SetParLimits(0, 0, 10.0 * std::max(amplitude, 1.e-9));
  prefit.SetParLimits(1, fitMin, fitMax);
  prefit.SetParLimits(2, 0.2 * std::max(sigma, 1.e-9),
                      5.0 * std::max(sigma, 1.e-9));
  hist->Fit(&prefit, "RQN");

  const double prefitAmplitude = prefit.GetParameter(0);
  const double prefitMean = prefit.GetParameter(1);
  const double prefitSigma = std::min(std::abs(prefit.GetParameter(2)),
                                      2.0 * std::max(sigma, 1.e-9));

  return CreatePeakFit(name, path, fitMin, fitMax,
                       prefitAmplitude, prefitMean, prefitSigma);
}

void GetPeakFitRange(TH1D *hist, const std::string &path,
                     double mean, double sigma,
                     double &fitMin, double &fitMax)
{
  const double axisMin = hist->GetXaxis()->GetXmin();
  const double axisMax = hist->GetXaxis()->GetXmax();

  if (UseGaussianOnlyFit(path)) {
    fitMin = std::max(axisMin, mean - 5.0 * sigma);
    fitMax = std::min(axisMax, mean + 5.0 * sigma);
    return;
  }

  fitMin = std::max(axisMin, mean - 15.0 * sigma);
  fitMax = std::min(axisMax, mean + 4.0 * sigma);
}

struct PulserPoint {
  double vamp = 0;
  double mean = 0;
  double stddev = 0;
};

struct BlankCalibration {
  std::string path;
  TH1D *hist = nullptr;
  TLegend *legend = nullptr;
  double mean = 0;
  double meanError = 0;
  double conversion = 0;
  double conversionError = 0;
};

TLegend *CreatePeakFitLegend(TF1 *fit)
{
  const bool isTailed = fit->GetNpar() > 3;
  auto legend = new TLegend(0.12, isTailed ? 0.42 : 0.54, 0.58, 0.88);
  legend->SetFillColor(0);
  legend->SetFillStyle(0);
  legend->SetBorderSize(0);
  legend->SetTextSize(0.035);

  char line[128];
  if (isTailed) {
    legend->AddEntry(static_cast<TObject *>(nullptr), "Left-tail Crystal Ball", "");
  }
  else {
    legend->AddEntry(static_cast<TObject *>(nullptr), "Gaussian", "");
  }

  std::snprintf(line, sizeof(line), "A = %.3g", fit->GetParameter(0));
  legend->AddEntry(static_cast<TObject *>(nullptr), line, "");
  std::snprintf(line, sizeof(line), "#mu = %.3f", fit->GetParameter(1));
  legend->AddEntry(static_cast<TObject *>(nullptr), line, "");
  std::snprintf(line, sizeof(line), "#sigma = %.3f", fit->GetParameter(2));
  legend->AddEntry(static_cast<TObject *>(nullptr), line, "");
  if (isTailed) {
    std::snprintf(line, sizeof(line), "#alpha = %.3f", fit->GetParameter(3));
    legend->AddEntry(static_cast<TObject *>(nullptr), line, "");
  }
  const double reducedChi2 =
      fit->GetNDF() > 0 ? fit->GetChisquare() / fit->GetNDF() : 0;
  std::snprintf(line, sizeof(line), "#chi^{2}/NDF = %.3g/%d = %.3g",
                fit->GetChisquare(), fit->GetNDF(), reducedChi2);
  legend->AddEntry(static_cast<TObject *>(nullptr), line, "");
  return legend;
}

TLegend *CreatePol1ParameterLegend(TF1 *fit, const char *title)
{
  auto legend = new TLegend(0.12, 0.66, 0.50, 0.88);
  legend->SetFillColor(0);
  legend->SetFillStyle(0);
  legend->SetBorderSize(0);
  legend->SetTextSize(0.035);

  char line[128];
  legend->AddEntry(static_cast<TObject *>(nullptr), title, "");
  std::snprintf(line, sizeof(line), "p0 = %.4g", fit->GetParameter(0));
  legend->AddEntry(static_cast<TObject *>(nullptr), line, "");
  std::snprintf(line, sizeof(line), "p1 = %.4g", fit->GetParameter(1));
  legend->AddEntry(static_cast<TObject *>(nullptr), line, "");
  const double reducedChi2 =
      fit->GetNDF() > 0 ? fit->GetChisquare() / fit->GetNDF() : 0;
  std::snprintf(line, sizeof(line), "#chi^{2}/NDF = %.3g/%d = %.3g",
                fit->GetChisquare(), fit->GetNDF(), reducedChi2);
  legend->AddEntry(static_cast<TObject *>(nullptr), line, "");
  return legend;
}

void SetFiveSigmaRange(TH1D *hist)
{
  const double mean = hist->GetMean();
  const double sigma = hist->GetStdDev();
  if (sigma <= 0)
    return;

  const double axisMin = hist->GetXaxis()->GetXmin();
  const double axisMax = hist->GetXaxis()->GetXmax();
  const double xMin = std::max(axisMin, mean - 5.0 * sigma);
  const double xMax = std::min(axisMax, mean + 5.0 * sigma);
  hist->GetXaxis()->SetRangeUser(xMin, xMax);
}

std::vector<std::string> ReadSpeList(const char *listFile)
{
  std::vector<std::string> paths;
  std::ifstream input(listFile);
  if (!input.is_open()) {
    std::cerr << "Cannot open list file: " << listFile << std::endl;
    return paths;
  }

  std::string line;
  while (std::getline(input, line)) {
    line = Trim(line);
    if (line.empty() || line[0] == '#')
      continue;

    std::string path;
    const auto tab = line.find('\t');
    if (tab != std::string::npos) {
      // Current list.txt format: "<type>\t<relative_path>".
      path = Trim(line.substr(tab + 1));
    }
    else {
      // Backward-compatible with old one-column path lists, and with
      // whitespace-separated "<type> <path>" lists.
      std::istringstream ss(line);
      std::string first;
      std::string second;
      ss >> first;
      if (ss >> second)
        path = second;
      else
        path = first;
    }

    if (!path.empty() && !IsOldPath(path))
      paths.push_back(path);
  }
  return paths;
}

TH1D *ReadSpeHistogram(const std::string &path)
{
  std::ifstream input(path);
  if (!input.is_open()) {
    std::cerr << "Cannot open SPE file: " << path << std::endl;
    return nullptr;
  }

  std::string line;
  for (int iLine = 0; iLine < kHeaderLines; ++iLine) {
    if (!std::getline(input, line)) {
      std::cerr << "Header ended early in: " << path << std::endl;
      return nullptr;
    }
  }

  auto hist = new TH1D(HistNameFromPath(path).c_str(), path.c_str(),
                       kNChannels, 0, kNChannels);
  hist->SetDirectory(nullptr);
  hist->GetXaxis()->SetTitle("Channel");
  hist->GetYaxis()->SetTitle("Counts");

  double counts = 0;
  for (int channel = 0; channel < kNChannels; ++channel) {
    if (!(input >> counts)) {
      std::cerr << "Only read " << channel << " channels from: " << path
                << std::endl;
      delete hist;
      return nullptr;
    }
    if (ShouldSkipLowPulserChannel(path, channel))
      continue;
    if (ShouldSkipLowTargetChannel(path, channel))
      continue;
    hist->SetBinContent(channel + 1, counts);
  }

  return hist;
}

} // namespace

void draw_all_spe(const char *listFile = "list.txt",
                  const char *outRootFile = "spe_histograms.root",
                  const char *outPdfFile = "all_spe.pdf",
                  int energyBinFactor = 3)
{
  // energyBinFactor: divide kNChannels by this value for energy histogram bins.
  //   1 → 4096 bins (full resolution)
  //   4 → 1024 bins (default, faster rendering)
  //   8 →  512 bins
  const int kEnergyBins = kNChannels / std::max(energyBinFactor, 1);

  //gStyle->SetOptStat(0);

  auto top = new LKDrawingGroup("top");
  auto group0 = top -> CreateGroup("blank");
  auto group1 = top -> CreateGroup("pulser");
  auto group2 = top -> CreateGroup("target");
  auto group3 = top -> CreateGroup("target_energy");
  auto group4 = top -> CreateGroup("target_thickness");

  // Analytic approximation: E(x) = kAlphaEnergy * (1 - x/kAlphaRange)^(1/kPowerN)
  // Inverse:                x(E) = kAlphaRange * (1 - (E/kAlphaEnergy)^kPowerN)
  // Fitted to NIST ASTAR numerical integration for 3.271 MeV alpha in CD2 (rho=1.06 g/cm3)
  constexpr double kAlphaEnergy = 3.271;  // MeV, initial alpha energy (Gd-148)
  constexpr double kAlphaRange  = 16.57;  // um, range in CD2
  constexpr double kThicknessRange = 20;
  constexpr double kPowerN      = 1.753;  // power-law exponent

  auto EnergyToThickness = [](double eMeV) -> double {
    if (eMeV <= 0) return kAlphaRange;
    if (eMeV >= kAlphaEnergy) return 0.0;
    return kAlphaRange * (1.0 - std::pow(eMeV / kAlphaEnergy, kPowerN));
  };

  // Jacobian |dx/dE| for sigma conversion: dx = |dx/dE| * dE
  // dx/dE = -kAlphaRange * kPowerN / kAlphaEnergy * (E/kAlphaEnergy)^(kPowerN-1)
  auto dxdE = [](double eMeV) -> double {
    if (eMeV <= 0) return 0.0;
    return kAlphaRange * kPowerN / kAlphaEnergy
           * std::pow(eMeV / kAlphaEnergy, kPowerN - 1.0);
  };

  const auto paths = ReadSpeList(listFile);
  if (paths.empty()) {
    std::cerr << "No SPE files found in " << listFile << std::endl;
    return;
  }

  std::vector<TH1D *> hists;
  hists.reserve(paths.size());

  for (const auto &path : paths) {
    auto hist = ReadSpeHistogram(path);
    if (hist)
      hists.push_back(hist);
  }

  if (hists.empty()) {
    std::cerr << "No histograms were created." << std::endl;
    return;
  }

  const int colors[] = {kBlack,       kRed + 1,    kBlue + 1,   kGreen + 2,
                        kMagenta + 1, kCyan + 2,   kOrange + 7, kViolet + 1,
                        kAzure + 7,   kSpring + 5, kPink + 7,   kTeal + 3};
  const int nColors = sizeof(colors) / sizeof(colors[0]);

  double maxCounts = 0;
  for (auto hist : hists)
    maxCounts = std::max(maxCounts, hist->GetMaximum());

  std::vector<TF1 *> fits;
  fits.reserve(hists.size());
  std::vector<TLegend *> fitLegends;
  std::vector<PulserPoint> pulserPoints;
  std::vector<BlankCalibration> blankCalibrations;
  std::vector<TH1D *> targetEnergyHists;
  std::vector<TH1D *> targetThicknessHists;
  std::vector<TF1 *> targetEnergyFits;
  std::vector<TLegend *> targetEnergyFitLegends;

  for (size_t iHist = 0; iHist < hists.size(); ++iHist) {
    auto hist = hists[iHist];
    //hist->SetLineColor(colors[iHist % nColors]);
    hists[iHist]->SetTitle(paths[iHist].c_str());
    hist->SetLineWidth(1);
    //hist->SetMaximum(maxCounts * 2.0);
    hist->SetMinimum(0.5);
    //hist->Draw(iHist == 0 ? "hist" : "hist same");

    int data_type = IsBlankPath(paths[iHist]) ? 0 : (IsPulserPath(paths[iHist]) ? 1 : 2);

    LKDrawingGroup *group = group2;
    if (data_type == 0)      group = group0;
    else if (data_type == 1) group = group1;

    auto draw = group -> CreateDrawing();
    if (data_type == 1) {
      double vamp = 0;
      if (GetPulserVamp(paths[iHist], vamp))
        pulserPoints.push_back({vamp, hist->GetMean(), hist->GetStdDev()});
	      else
	        std::cerr << "Cannot parse pulser amplitude from: " << paths[iHist]
	                  << std::endl;
	      draw -> Add(hist);
	      SetFiveSigmaRange(hist);
	      continue;
    }

    const auto initial = GetPeakInitialParameters(hist);

    if (initial.sigma > 0) {
      double fitMin = 0;
      double fitMax = 0;
      GetPeakFitRange(hist, paths[iHist], initial.mean, initial.sigma,
                      fitMin, fitMax);
      auto fit = CreateInitializedPeakFit(FitNameFromPath(paths[iHist]),
                                          paths[iHist], hist, fitMin, fitMax,
                                          initial.amplitude, initial.mean,
                                          initial.sigma);
      //fit->SetLineColor(colors[iHist % nColors]);
      hist->Fit(fit, "RQ");
      //auto draw = top -> CreateDrawing();
      auto fitLegend = CreatePeakFitLegend(fit);
      draw -> Add(hist);
      draw -> Add(fit,"samel");
      draw -> Add(fitLegend,"same");
      draw -> SetLegendCorner(1);
      fits.push_back(fit);
      fitLegends.push_back(fitLegend);
      if (data_type == 0) {
        BlankCalibration calibration;
        calibration.path = paths[iHist];
        calibration.hist = hist;
        calibration.legend = fitLegend;
        calibration.mean = fit->GetParameter(1);
        calibration.meanError = fit->GetParError(1);
        blankCalibrations.push_back(calibration);
      }
    }
    else {
      std::cerr << "Skipping fit for zero-width spectrum: " << paths[iHist]
                << std::endl;
    }

	    SetFiveSigmaRange(hist);
	  }

  TGraphErrors *pulserGraph = nullptr;
  TF1 *pulserFit = nullptr;
  TGraphErrors *pulserGraphNo1V = nullptr;
  TF1 *pulserFitNo1V = nullptr;
  TH2D *pulserFrame = nullptr;
  TH2D *pulserFrameNo1V = nullptr;
  TLegend *pulserFitLegend = nullptr;
  TLegend *pulserFitLegendNo1V = nullptr;
  if (!pulserPoints.empty()) {
    std::sort(pulserPoints.begin(), pulserPoints.end(),
              [](const PulserPoint &a, const PulserPoint &b) {
                return a.vamp < b.vamp;
              });

    pulserGraph = new TGraphErrors(static_cast<int>(pulserPoints.size()));
    pulserGraph->SetName("g_pulser_mean_vs_vamp");
    pulserGraph->SetTitle("Pulser;V_{amp} (V);Mean channel");
    pulserGraph->SetMarkerStyle(20);
    pulserGraph->SetMarkerSize(1.0);
    pulserGraph->SetLineWidth(2);

    std::vector<PulserPoint> pulserPointsNo1V;
    for (size_t iPoint = 0; iPoint < pulserPoints.size(); ++iPoint) {
      pulserGraph->SetPoint(static_cast<int>(iPoint), pulserPoints[iPoint].vamp,
                            pulserPoints[iPoint].mean);
      pulserGraph->SetPointError(static_cast<int>(iPoint), 0,
                                 pulserPoints[iPoint].stddev);
      if (std::abs(pulserPoints[iPoint].vamp - 1.0) > 1.e-9)
        pulserPointsNo1V.push_back(pulserPoints[iPoint]);
    }

    const double fitMin = 0.0;
    const double fitMax = 2.1;
    pulserFit = new TF1("f_pulser_mean_vs_vamp", "pol1", fitMin, fitMax);
    pulserFit->SetNpx(500);
    pulserFit->SetLineWidth(2);
    pulserGraph->Fit(pulserFit, "RQ");
    pulserFitLegend = CreatePol1ParameterLegend(pulserFit, "All points");

    if (pulserPointsNo1V.size() >= 2) {
      pulserGraphNo1V = new TGraphErrors(static_cast<int>(pulserPointsNo1V.size()));
      pulserGraphNo1V->SetName("g_pulser_mean_vs_vamp_no1V");
      pulserGraphNo1V->SetTitle("Pulser without 1 V;V_{amp} (V);Mean channel");
      pulserGraphNo1V->SetMarkerStyle(24);
      pulserGraphNo1V->SetMarkerSize(1.0);

      for (size_t iPoint = 0; iPoint < pulserPointsNo1V.size(); ++iPoint) {
        pulserGraphNo1V->SetPoint(static_cast<int>(iPoint),
                                  pulserPointsNo1V[iPoint].vamp,
                                  pulserPointsNo1V[iPoint].mean);
        pulserGraphNo1V->SetPointError(static_cast<int>(iPoint), 0,
                                       pulserPointsNo1V[iPoint].stddev);
      }

      pulserFitNo1V = new TF1("f_pulser_mean_vs_vamp_no1V", "pol1",
                              fitMin, fitMax);
      pulserFitNo1V->SetNpx(500);
      pulserFitNo1V->SetLineColor(kRed + 1);
      pulserFitNo1V->SetLineStyle(2);
      pulserFitNo1V->SetLineWidth(2);
      pulserGraphNo1V->Fit(pulserFitNo1V, "RQ");
      pulserFitLegendNo1V = CreatePol1ParameterLegend(pulserFitNo1V, "Without 1 V");
    }

    pulserFrame = new TH2D("h_pulser_frame",
                           "Pulser;V_{amp} (V);Mean channel",
                           100, 0, 2.1, 100, 0, 900);
    pulserFrame->SetDirectory(nullptr);
    pulserFrame->SetStats(0);

    auto draw = group1 -> CreateDrawing();
    draw -> Add(pulserFrame, "");
    draw -> Add(pulserGraph, "P same");
    draw -> Add(pulserFit, "samel");
    draw -> Add(pulserFitLegend, "same");
    draw -> SetLegendCorner(1);

	  if (pulserGraphNo1V && pulserFitNo1V) {
      pulserFrameNo1V = new TH2D("h_pulser_frame_no1V",
                                 "Pulser without 1 V;V_{amp} (V);Mean channel",
                                 100, 0, 2.1, 100, 0, 900);
      pulserFrameNo1V->SetDirectory(nullptr);
      pulserFrameNo1V->SetStats(0);

      auto drawNo1V = group1 -> CreateDrawing();
      drawNo1V -> Add(pulserFrameNo1V, "");
      drawNo1V -> Add(pulserGraphNo1V, "P same");
      drawNo1V -> Add(pulserFitNo1V, "samel");
      drawNo1V -> Add(pulserFitLegendNo1V, "same");
	      drawNo1V -> SetLegendCorner(1);
	    }
	  }

  if (pulserFitNo1V && !blankCalibrations.empty()) {
    const double gd148AlphaMeV = 3.271;
    const double p0 = pulserFitNo1V->GetParameter(0);
    const double p0Error = pulserFitNo1V->GetParError(0);
    std::map<std::string, double> conversionByDate;

    std::ofstream calibrationOutput("blank_energy_calibration.txt");
    calibrationOutput << "# energy = (Mean - [p0]) * [conversion]\n";
    calibrationOutput << "# p0(no 1V pulser fit) = " << p0
                      << " +/- " << p0Error << " channel\n";
    calibrationOutput << "# Gd-148 alpha = " << gd148AlphaMeV << " MeV\n";
    calibrationOutput << "# path\tmean_channel\tmean_error\tconversion_MeV_per_channel\tconversion_error\n";

    std::cout << "Blank energy calibration using no-1V pulser p0 = " << p0
              << std::endl;
    for (auto &calibration : blankCalibrations) {
      const double correctedChannel = calibration.mean - p0;
      const double correctedChannelError =
          std::sqrt(calibration.meanError * calibration.meanError +
                    p0Error * p0Error);
      calibration.conversion = gd148AlphaMeV / correctedChannel;
      calibration.conversionError =
          gd148AlphaMeV * correctedChannelError /
          (correctedChannel * correctedChannel);

      if (calibration.legend) {
        char line[128];
        std::snprintf(line, sizeof(line), "E = (#mu - p0) #times conv");
        calibration.legend->AddEntry(static_cast<TObject *>(nullptr), line, "");
        std::snprintf(line, sizeof(line), "conv = %.6f MeV/ch",
                      calibration.conversion);
        calibration.legend->AddEntry(static_cast<TObject *>(nullptr), line, "");
        std::snprintf(line, sizeof(line), "p0 = %.4f ch", p0);
        calibration.legend->AddEntry(static_cast<TObject *>(nullptr), line, "");
      }

      calibrationOutput << calibration.path << "\t"
                        << calibration.mean << "\t"
                        << calibration.meanError << "\t"
                        << calibration.conversion << "\t"
                        << calibration.conversionError << "\n";
      const auto date = DateFromPath(calibration.path);
      if (!date.empty())
        conversionByDate[date] = calibration.conversion;
      std::cout << "  " << calibration.path
                << ": conversion = " << calibration.conversion
                << " MeV/channel, energy = (Mean - [" << p0
                << "]) * [" << calibration.conversion << "]" << std::endl;

      // ── group3: blank channel → energy ───────────────────────────────────
      if (calibration.hist) {
        const double eMin = (0.0 - p0) * calibration.conversion;
        const double eMax = (static_cast<double>(kNChannels) - p0) * calibration.conversion;
        auto blankEnergyHist = new TH1D(
            ("blank_" + EnergyHistNameFromPath(calibration.path)).c_str(),
            (calibration.path + ";Energy (MeV);Counts").c_str(),
            kEnergyBins, eMin, eMax);
        blankEnergyHist->SetDirectory(nullptr);
        blankEnergyHist->SetLineWidth(1);

        for (int iBin = 1; iBin <= kNChannels; ++iBin) {
          const double counts = calibration.hist->GetBinContent(iBin);
          if (counts <= 0)
            continue;
          const double channel = calibration.hist->GetBinCenter(iBin);
          const double energy = (channel - p0) * calibration.conversion;
          const int eBin = blankEnergyHist->FindBin(energy);
          if (eBin >= 1 && eBin <= blankEnergyHist->GetNbinsX())
            blankEnergyHist->Fill(energy, counts);
        }

        auto drawBlankE = group3 -> CreateDrawing();
        drawBlankE -> Add(blankEnergyHist, "hist");

        const auto initialBlank = GetPeakInitialParameters(blankEnergyHist);
        TF1 *blankEFit = nullptr;
        if (initialBlank.sigma > 0) {
          double fitMin = 0, fitMax = 0;
          GetPeakFitRange(blankEnergyHist, calibration.path,
                          initialBlank.mean, initialBlank.sigma, fitMin, fitMax);
          blankEFit = CreateInitializedPeakFit(
              "blank_" + EnergyFitNameFromPath(calibration.path),
              calibration.path, blankEnergyHist, fitMin, fitMax,
              initialBlank.amplitude, initialBlank.mean, initialBlank.sigma);
          blankEnergyHist->Fit(blankEFit, "RQ");
          auto blankELegend = CreatePeakFitLegend(blankEFit);
          drawBlankE -> Add(blankEFit, "samel");
          drawBlankE -> Add(blankELegend, "same");
          drawBlankE -> SetLegendCorner(1);
          targetEnergyFits.push_back(blankEFit);
          targetEnergyFitLegends.push_back(blankELegend);
        }
        SetFiveSigmaRange(blankEnergyHist);
        targetEnergyHists.push_back(blankEnergyHist);

        // ── group4: blank energy → thickness (should peak near 0 µm) ────────
        {
          constexpr int kThicknessBins = 100;
          auto blankThickHist = new TH1D(
              ("blank_ht_" + EnergyHistNameFromPath(calibration.path)).c_str(),
              (calibration.path + ";Target Thickness (#mum);Counts").c_str(),
              kThicknessBins, 0.0, kThicknessRange);
          blankThickHist->SetDirectory(nullptr);
          blankThickHist->SetLineWidth(1);

          // Blank peak sits at exactly kAlphaEnergy by construction.
          // Do NOT apply the eBinCenter >= kAlphaEnergy cut here — that would
          // clip the right half of the Gaussian and bias the peak to > 0 um.
          // Instead, clamp E to [0, kAlphaEnergy] so the full peak maps to
          // thickness = 0 (right edge is clamped, not discarded).
          for (int iBin = 1; iBin <= blankEnergyHist->GetNbinsX(); ++iBin) {
            const double eBinCenter = blankEnergyHist->GetBinCenter(iBin);
            const double counts     = blankEnergyHist->GetBinContent(iBin);
            if (counts <= 0 || eBinCenter <= 0)
              continue;
            const double eClamp    = std::min(eBinCenter, kAlphaEnergy);
            const double thickness = EnergyToThickness(eClamp);
            const int thickBin = blankThickHist->FindBin(thickness);
            if (thickBin >= 1 && thickBin <= kThicknessBins)
              blankThickHist->Fill(thickness, counts);
          }

          auto drawBlankT = group4 -> CreateDrawing();
          drawBlankT -> Add(blankThickHist, "hist");

          if (blankEFit) {
            const double ePeak  = blankEFit->GetParameter(1);
            const double eSigma = std::abs(blankEFit->GetParameter(2));
            if (ePeak > 0 && ePeak < kAlphaEnergy) {
              const double xPeak  = EnergyToThickness(ePeak);
              const double xSigma = dxdE(ePeak) * eSigma;

              const double fitMin = std::max(0.0, xPeak - 5.0 * xSigma);
              const double fitMax = std::min(kThicknessRange, xPeak + 5.0 * xSigma);
              const double amp    = blankThickHist->GetMaximum();

              auto blankThickFit = new TF1(
                  ("blank_fthick_" + EnergyHistNameFromPath(calibration.path)).c_str(),
                  "gaus", fitMin, fitMax);
              blankThickFit->SetParameters(amp, xPeak, std::max(xSigma, 1.e-4));
              blankThickFit->SetParNames("Amplitude", "Mean", "Sigma");
              blankThickFit->SetParLimits(0, 0, 10.0 * std::max(amp, 1.e-9));
              blankThickFit->SetParLimits(1, fitMin, fitMax);
              blankThickFit->SetParLimits(2, 0.2 * std::max(xSigma, 1.e-4),
                                             5.0 * std::max(xSigma, 1.e-4));
              blankThickFit->SetNpx(500);
              blankThickFit->SetLineWidth(2);
              blankThickHist->Fit(blankThickFit, "RQ");

              auto blankThickLegend = new TLegend(0.12, 0.54, 0.58, 0.88);
              blankThickLegend->SetFillColor(0);
              blankThickLegend->SetFillStyle(0);
              blankThickLegend->SetBorderSize(0);
              blankThickLegend->SetTextSize(0.035);
              char line[128];
              blankThickLegend->AddEntry(static_cast<TObject *>(nullptr), "Gaussian (Blank)", "");
              std::snprintf(line, sizeof(line), "#mu = %.4f #mum", blankThickFit->GetParameter(1));
              blankThickLegend->AddEntry(static_cast<TObject *>(nullptr), line, "");
              std::snprintf(line, sizeof(line), "#sigma = %.4f #mum", blankThickFit->GetParameter(2));
              blankThickLegend->AddEntry(static_cast<TObject *>(nullptr), line, "");
              const double reducedChi2 = blankThickFit->GetNDF() > 0
                  ? blankThickFit->GetChisquare() / blankThickFit->GetNDF() : 0;
              std::snprintf(line, sizeof(line), "#chi^{2}/NDF = %.3g/%d = %.3g",
                            blankThickFit->GetChisquare(), blankThickFit->GetNDF(), reducedChi2);
              blankThickLegend->AddEntry(static_cast<TObject *>(nullptr), line, "");

              drawBlankT -> Add(blankThickFit, "samel");
              drawBlankT -> Add(blankThickLegend, "same");
              drawBlankT -> SetLegendCorner(1);
            }
          }

          blankThickHist->GetXaxis()->SetRangeUser(0, 5);
          targetThicknessHists.push_back(blankThickHist);
        }
      }
    }

    for (size_t iHist = 0; iHist < hists.size(); ++iHist) {
      const auto &path = paths[iHist];
      if (IsBlankPath(path) || IsPulserPath(path))
        continue;

      const auto date = DateFromPath(path);
      auto conversionDate = date;
      if (date == "20260529")
        conversionDate = "20260528";

      const auto conversionIt = conversionByDate.find(conversionDate);
      if (conversionIt == conversionByDate.end()) {
        std::cerr << "Skipping energy histogram, no blank conversion for date "
                  << conversionDate << ": " << path << std::endl;
        continue;
      }

      const double conversion = conversionIt->second;
      const double eMin = (0.0 - p0) * conversion;
      const double eMax = (static_cast<double>(kNChannels) - p0) * conversion;
      auto energyHist = new TH1D(EnergyHistNameFromPath(path).c_str(),
                                 (path + ";Energy (MeV);Counts").c_str(),
                                 kEnergyBins, eMin, eMax);
      energyHist->SetDirectory(nullptr);
      energyHist->SetLineWidth(1);

      for (int iBin = 1; iBin <= kNChannels; ++iBin) {
        const double counts = hists[iHist]->GetBinContent(iBin);
        if (counts <= 0)
          continue;
        const double channel = hists[iHist]->GetBinCenter(iBin);
        const double energy = (channel - p0) * conversion;
        const int eBin = energyHist->FindBin(energy);
        if (eBin >= 1 && eBin <= energyHist->GetNbinsX())
          energyHist->Fill(energy, counts);
      }

      auto draw = group3 -> CreateDrawing();
      draw -> Add(energyHist, "hist");

      const auto initial = GetPeakInitialParameters(energyHist);
      if (initial.sigma > 0) {
        double fitMin = 0;
        double fitMax = 0;
        GetPeakFitRange(energyHist, path, initial.mean, initial.sigma,
                        fitMin, fitMax);
        auto fit = CreateInitializedPeakFit(EnergyFitNameFromPath(path), path,
                                            energyHist, fitMin, fitMax,
                                            initial.amplitude, initial.mean,
                                            initial.sigma);
        energyHist->Fit(fit, "RQ");

        auto fitLegend = CreatePeakFitLegend(fit);
        draw -> Add(fit, "samel");
        draw -> Add(fitLegend, "same");
        draw -> SetLegendCorner(1);

        targetEnergyFits.push_back(fit);
        targetEnergyFitLegends.push_back(fitLegend);
      }
      else {
        std::cerr << "Skipping energy fit for zero-width spectrum: "
                  << path << std::endl;
      }

      SetFiveSigmaRange(energyHist);
      targetEnergyHists.push_back(energyHist);

      // ── group4: remaining energy → target thickness (analytic approximation) ──
      {
        // x(E) = kAlphaRange * (1 - (E/kAlphaEnergy)^kPowerN)
        // Histogram x-axis: thickness in um (0 ~ kAlphaRange)
        //constexpr int kThicknessBins = 100;//kNChannels/4;
        //constexpr int kThicknessBins = kNChannels/4;
        constexpr int kThicknessBins = 100;
        auto thickHist = new TH1D(
            ("ht_" + EnergyHistNameFromPath(path)).c_str(),
            (path + ";Target Thickness (#mum);Counts").c_str(),
            //kThicknessBins, 0.0, kAlphaRange);
            kThicknessBins, 0.0, kThicknessRange);
        thickHist->SetDirectory(nullptr);
        thickHist->SetLineWidth(1);

        // Fill: transform each energy bin to a thickness bin
        for (int iBin = 1; iBin <= energyHist->GetNbinsX(); ++iBin) {
          const double eBinCenter = energyHist->GetBinCenter(iBin);
          const double counts     = energyHist->GetBinContent(iBin);
          if (counts <= 0 || eBinCenter <= 0 || eBinCenter >= kAlphaEnergy)
            continue;
          const double thickness = EnergyToThickness(eBinCenter);
          const int thickBin = thickHist->FindBin(thickness);
          if (thickBin >= 1 && thickBin <= kThicknessBins)
            thickHist->Fill(thickness, counts);
        }

        auto drawThick = group4 -> CreateDrawing();
        drawThick -> Add(thickHist, "hist");

        // Fit the thickness peak if the energy fit succeeded
        if (!targetEnergyFits.empty()) {
          auto eFit = targetEnergyFits.back();
          const double ePeak  = eFit->GetParameter(1);
          const double eSigma = std::abs(eFit->GetParameter(2));

          if (ePeak > 0 && ePeak < kAlphaEnergy) {
            const double xPeak  = EnergyToThickness(ePeak);
            const double xSigma = dxdE(ePeak) * eSigma;

            // Gaussian fit around the converted peak position
            const double fitMin = std::max(0.0, xPeak - 5.0 * xSigma);
            const double fitMax = std::min(kAlphaRange, xPeak + 5.0 * xSigma);
            const double amp    = thickHist->GetMaximum();

            auto thickFit = new TF1(
                ("fthick_" + EnergyHistNameFromPath(path)).c_str(),
                "gaus", fitMin, fitMax);
            thickFit->SetParameters(amp, xPeak, xSigma);
            thickFit->SetParNames("Amplitude", "Mean", "Sigma");
            thickFit->SetParLimits(0, 0, 10.0 * std::max(amp, 1.e-9));
            thickFit->SetParLimits(1, fitMin, fitMax);
            thickFit->SetParLimits(2, 0.2 * std::max(xSigma, 1.e-9),
                                      5.0 * std::max(xSigma, 1.e-9));
            thickFit->SetNpx(500);
            thickFit->SetLineWidth(2);
            thickHist->Fit(thickFit, "RQ");

            // Legend
            auto thickLegend = new TLegend(0.12, 0.54, 0.58, 0.88);
            thickLegend->SetFillColor(0);
            thickLegend->SetFillStyle(0);
            thickLegend->SetBorderSize(0);
            thickLegend->SetTextSize(0.035);
            char line[128];
            thickLegend->AddEntry(static_cast<TObject *>(nullptr), "Gaussian", "");
            std::snprintf(line, sizeof(line), "#mu = %.4f #mum",
                          thickFit->GetParameter(1));
            thickLegend->AddEntry(static_cast<TObject *>(nullptr), line, "");
            std::snprintf(line, sizeof(line), "#sigma = %.4f #mum",
                          thickFit->GetParameter(2));
            thickLegend->AddEntry(static_cast<TObject *>(nullptr), line, "");
            const double reducedChi2 = thickFit->GetNDF() > 0
                ? thickFit->GetChisquare() / thickFit->GetNDF() : 0;
            std::snprintf(line, sizeof(line), "#chi^{2}/NDF = %.3g/%d = %.3g",
                          thickFit->GetChisquare(), thickFit->GetNDF(), reducedChi2);
            thickLegend->AddEntry(static_cast<TObject *>(nullptr), line, "");
            std::snprintf(line, sizeof(line),
                          "x(E) = R(1-(E/E_{0})^{n}), n=%.3f", kPowerN);
            thickLegend->AddEntry(static_cast<TObject *>(nullptr), line, "");

            drawThick -> Add(thickFit, "samel");
            drawThick -> Add(thickLegend, "same");
            drawThick -> SetLegendCorner(1);
          }
        }

        // Zoom x-axis to peak region
        thickHist->GetXaxis()->SetRangeUser(0, kAlphaRange * 0.5);
        targetThicknessHists.push_back(thickHist);
      }
    }
  }

  TFile output(outRootFile, "RECREATE");
  for (auto hist : hists)
    hist->Write();
  for (auto hist : targetEnergyHists)
    hist->Write();
  for (auto hist : targetThicknessHists)
    hist->Write();
  for (auto fit : fits)
    fit->Write();
  for (auto fit : targetEnergyFits)
    fit->Write();
  for (auto legend : fitLegends)
    legend->Write();
  for (auto legend : targetEnergyFitLegends)
    legend->Write();
  if (pulserFrame)
    pulserFrame->Write();
  if (pulserFrameNo1V)
    pulserFrameNo1V->Write();
  if (pulserGraph)
    pulserGraph->Write();
  if (pulserGraphNo1V)
    pulserGraphNo1V->Write();
  if (pulserFit)
    pulserFit->Write();
  if (pulserFitNo1V)
    pulserFitNo1V->Write();
  if (pulserFitLegend)
    pulserFitLegend->Write();
  if (pulserFitLegendNo1V)
    pulserFitLegendNo1V->Write();
  top->Write("top");
  output.Close();

  const int nColumns = 5;
  const int nRows = 5;
  const int perPage = nColumns * nRows;

  const std::string pdfOpen = std::string(outPdfFile) + "(";
  const std::string pdfClose = std::string(outPdfFile) + ")";

  for (size_t first = 0; first < hists.size(); first += perPage) {
    //pages.Clear();
    //pages.Divide(nColumns, nRows, 0.001, 0.001);

    const size_t last = std::min(first + perPage, hists.size());
    for (size_t iHist = first; iHist < last; ++iHist) {
      //pages.cd(static_cast<int>(iHist - first + 1));
      //gPad->SetLogy();
      //hists[iHist]->Draw("hist");
      //lk_debug << hists[iHist]->GetEntries() << endl;
      //top -> AddHist(hists[iHist],"hist");
        //break;
    }

    top -> Draw();

    //if (first == 0 && last == hists.size())
    //  pages.SaveAs(outPdfFile);
    //else if (first == 0)
    //  pages.SaveAs(pdfOpen.c_str());
    //else if (last == hists.size())
    //  pages.SaveAs(pdfClose.c_str());
    //else
    //  pages.SaveAs(outPdfFile);
  }

  std::cout << "Created " << hists.size() << " histograms from " << listFile
            << std::endl;
  std::cout << "Wrote " << outRootFile << " and " << outPdfFile << std::endl;
}
