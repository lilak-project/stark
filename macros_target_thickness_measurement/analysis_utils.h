#ifndef ANALYSIS_UTILS_H
#define ANALYSIS_UTILS_H

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "TH1D.h"

// ── Constants ─────────────────────────────────────────────────────────────────
constexpr int    kNChannels   = 4096;
constexpr int    kHeaderLines = 12;
constexpr double kAlphaEnergy = 3.271;  // MeV, Gd-148 alpha

// ── Data type ─────────────────────────────────────────────────────────────────
enum class SpeType { kBlank, kPulser, kCh2, kCarbon, kCd2, kUnknown };

inline std::string SpeTypeName(SpeType t)
{
  switch (t) {
    case SpeType::kBlank:  return "blank";
    case SpeType::kPulser: return "pulser";
    case SpeType::kCh2:    return "ch2";
    case SpeType::kCarbon: return "carbon";
    case SpeType::kCd2:    return "cd2";
    default:               return "unknown";
  }
}

inline SpeType SpeTypeFromString(const std::string &s)
{
  if (s == "blank")  return SpeType::kBlank;
  if (s == "pulser") return SpeType::kPulser;
  if (s == "ch2")    return SpeType::kCh2;
  if (s == "carbon") return SpeType::kCarbon;
  if (s == "cd2")    return SpeType::kCd2;
  return SpeType::kUnknown;
}

// ── List entry ────────────────────────────────────────────────────────────────
struct SpeEntry {
  SpeType     type = SpeType::kUnknown;
  std::string path;
};

// ── Path utilities ────────────────────────────────────────────────────────────
inline std::string Trim(const std::string &s)
{
  const auto a = s.find_first_not_of(" \t\r\n");
  if (a == std::string::npos) return "";
  return s.substr(a, s.find_last_not_of(" \t\r\n") - a + 1);
}

inline std::string DateFromPath(const std::string &path)
{
  const std::string prefix = "data/";
  const auto begin = path.find(prefix);
  if (begin == std::string::npos) return "";
  const auto dateBegin = begin + prefix.size();
  const auto dateEnd   = path.find('/', dateBegin);
  if (dateEnd == std::string::npos) return "";
  return path.substr(dateBegin, dateEnd - dateBegin);
}

inline std::string HistNameFromPath(std::string path)
{
  std::replace(path.begin(), path.end(), '/', '_');
  std::replace(path.begin(), path.end(), '.', '_');
  std::replace(path.begin(), path.end(), '-', '_');
  return "h_" + path;
}

inline bool GetPulserVamp(const std::string &path, double &vamp)
{
  const auto slash = path.find_last_of('/');
  const auto name  = slash == std::string::npos ? path : path.substr(slash + 1);
  const auto dash  = name.find('-');
  if (dash == std::string::npos) return false;

  const char *beg = name.c_str() + dash + 1;
  char *end = nullptr;
  const double value = std::strtod(beg, &end);
  if (end == beg) return false;

  const std::string unit(end, std::min<size_t>(2, std::strlen(end)));
  if (unit == "mV") { vamp = value * 1.e-3; return true; }
  if (!unit.empty() && unit[0] == 'V') { vamp = value; return true; }
  return false;
}

// ── Channel cut helpers (mirrors draw_all_spe.C) ──────────────────────────────

// Skip low channels for 1V and 1.6V pulser (noisy region near zero)
inline bool ShouldSkipLowPulserChannel(const SpeEntry &e, int channel)
{
  if (e.type != SpeType::kPulser) return false;
  double vamp = 0;
  if (!GetPulserVamp(e.path, vamp)) return false;
  const bool is1V  = vamp > 0.9  && vamp < 1.1;
  const bool is1V6 = vamp > 1.5  && vamp < 1.7;
  return channel < 50 && (is1V || is1V6);
}

// Skip low channels for target data (below detector threshold)
inline bool ShouldSkipLowTargetChannel(const SpeEntry &e, int channel)
{
  return channel < 300 &&
         e.type != SpeType::kBlank &&
         e.type != SpeType::kPulser;
}

// ── Peak utilities ────────────────────────────────────────────────────────────
struct PeakParams {
  double amplitude = 0;
  double mean      = 0;
  double sigma     = 0;
};

inline double EstimatePeakSigma(TH1D *hist)
{
  const int    maxBin = hist->GetMaximumBin();
  const double halfMax = 0.5 * hist->GetBinContent(maxBin);
  if (halfMax <= 0) return hist->GetStdDev();

  int lb = maxBin, rb = maxBin;
  while (lb > 1                 && hist->GetBinContent(lb) > halfMax) --lb;
  while (rb < hist->GetNbinsX() && hist->GetBinContent(rb) > halfMax) ++rb;

  const double fwhm = hist->GetBinCenter(rb) - hist->GetBinCenter(lb);
  return fwhm > 0 ? fwhm / 2.35482 : hist->GetStdDev();
}

inline PeakParams GetPeakParams(TH1D *hist)
{
  const int maxBin = hist->GetMaximumBin();
  return { hist->GetBinContent(maxBin),
           hist->GetBinCenter(maxBin),
           EstimatePeakSigma(hist) };
}

// Set x-axis range to mean ± 5σ (from draw_all_spe.C)
inline void SetFiveSigmaRange(TH1D *hist)
{
  const double mean  = hist->GetMean();
  const double sigma = hist->GetStdDev();
  if (sigma <= 0) return;
  hist->GetXaxis()->SetRangeUser(
      std::max(hist->GetXaxis()->GetXmin(), mean - 5.0 * sigma),
      std::min(hist->GetXaxis()->GetXmax(), mean + 5.0 * sigma));
}

// ── List file reader ──────────────────────────────────────────────────────────
// Format: "<type>\t<relative_path>"
inline std::vector<SpeEntry> ReadSpeList(const char *listFile)
{
  std::vector<SpeEntry> entries;
  std::ifstream f(listFile);
  if (!f.is_open()) {
    std::cerr << "Cannot open list file: " << listFile << std::endl;
    return entries;
  }
  std::string line;
  while (std::getline(f, line)) {
    line = Trim(line);
    if (line.empty() || line[0] == '#') continue;
    const auto tab = line.find('\t');
    if (tab == std::string::npos) continue;
    SpeEntry e;
    e.type = SpeTypeFromString(Trim(line.substr(0, tab)));
    e.path = Trim(line.substr(tab + 1));
    if (e.type != SpeType::kUnknown && !e.path.empty())
      entries.push_back(e);
  }
  return entries;
}

// ── SPE histogram reader ──────────────────────────────────────────────────────
inline TH1D *ReadSpeHistogram(const SpeEntry &entry)
{
  std::ifstream f(entry.path);
  if (!f.is_open()) {
    std::cerr << "Cannot open SPE file: " << entry.path << std::endl;
    return nullptr;
  }

  std::string line;
  for (int i = 0; i < kHeaderLines; ++i) {
    if (!std::getline(f, line)) {
      std::cerr << "Header ended early: " << entry.path << std::endl;
      return nullptr;
    }
  }

  const std::string name  = HistNameFromPath(entry.path);
  const std::string title = entry.path;
  auto hist = new TH1D(name.c_str(), title.c_str(),
                       kNChannels, 0, kNChannels);
  hist->SetDirectory(nullptr);
  hist->GetXaxis()->SetTitle("Channel");
  hist->GetYaxis()->SetTitle("Counts");

  double counts = 0;
  for (int ch = 0; ch < kNChannels; ++ch) {
    if (!(f >> counts)) {
      std::cerr << "Only read " << ch << " channels from: " << entry.path << std::endl;
      delete hist;
      return nullptr;
    }
    // Apply channel cuts from original draw_all_spe.C
    if (ShouldSkipLowPulserChannel(entry, ch)) continue;
    if (ShouldSkipLowTargetChannel(entry, ch)) continue;
    hist->SetBinContent(ch + 1, counts);
  }
  return hist;
}

#endif // ANALYSIS_UTILS_H
