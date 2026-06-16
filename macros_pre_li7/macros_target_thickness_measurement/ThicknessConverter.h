#ifndef THICKNESS_CONVERTER_H
#define THICKNESS_CONVERTER_H

/**
 * ThicknessConverter.h
 *
 * Remaining alpha energy -> target thickness using tabulated NIST ASTAR
 * mass stopping powers.  Compound targets use Bragg additivity:
 *
 *   S_mix(E) = sum_i w_i S_i(E)
 *
 * with S in MeV cm2/mg.  Thickness is obtained by numerical integration:
 *
 *   x(E) = integral_E^E0 dE' / (S_mix(E') * rho * 0.1)
 *
 * where rho is in g/cm3 and x is in um.  The factor 0.1 converts
 * (MeV cm2/mg)*(g/cm3) to MeV/um.
 *
 * CD2 uses the same approximation as alpha_cd2_table.py:
 * carbon 75% + deuterium 25% by mass, with D stopping approximated by
 * the ASTAR hydrogen mass stopping power.
 */

#include <algorithm>
#include <cmath>
#include <string>
#include <vector>

#include "TH1D.h"
#include "TGraph.h"

namespace ThicknessConv {
  constexpr double kE0 = 3.271;  // MeV, Gd-148 alpha energy

  // Densities in g/cm3.  Adjust here if the actual target density differs.
  constexpr double kDensityCarbon = 2.267;  // graphite
  constexpr double kDensityCH2    = 0.94;   // polyethylene
  constexpr double kDensityCD2    = 1.06;   // deuterated polyethylene

  constexpr double kMassC = 12.0107;
  constexpr double kMassH = 1.00794;
  constexpr double kMassD = 2.01410;

  constexpr double kCH2WeightC = kMassC / (kMassC + 2.0 * kMassH);
  constexpr double kCH2WeightH = 1.0 - kCH2WeightC;
  constexpr double kCD2WeightC = kMassC / (kMassC + 2.0 * kMassD);
  constexpr double kCD2WeightD = 1.0 - kCD2WeightC;

  constexpr double kThicknessRange = 30.0;  // um, generic plotting upper limit

  enum class Material { kCarbon, kCH2, kCD2 };

  // NIST ASTAR alpha mass stopping power tables [MeV, MeV cm2/mg].
  constexpr int kNTable = 45;
  constexpr double kEnergy[kNTable] = {
      0.10, 0.125, 0.15, 0.175, 0.20, 0.225, 0.25, 0.275, 0.30,
      0.35, 0.40, 0.45, 0.50, 0.55, 0.60, 0.65, 0.70, 0.75,
      0.80, 0.85, 0.90, 0.95, 1.00, 1.25, 1.50, 1.75, 2.00,
      2.25, 2.50, 2.75, 3.00, 3.50, 4.00, 4.50, 5.00, 5.50,
      6.00, 6.50, 7.00, 7.50, 8.00, 8.50, 9.00, 9.50, 10.00,
  };
  constexpr double kStopCarbon[kNTable] = {
      4.918, 5.051, 5.097, 5.094, 4.999, 4.870, 4.727, 4.577, 4.426,
      4.130, 3.858, 3.617, 3.405, 3.218, 3.052, 2.904, 2.771, 2.650,
      2.541, 2.442, 2.352, 2.269, 2.194, 1.898, 1.682, 1.517, 1.389,
      1.286, 1.202, 1.131, 1.072, 0.975, 0.899, 0.838, 0.789, 0.747,
      0.712, 0.681, 0.654, 0.629, 0.608, 0.588, 0.570, 0.554, 0.540,
  };
  constexpr double kStopHydrogen[kNTable] = {
      15.91, 16.18, 15.96, 15.51, 14.93, 14.29, 13.63, 12.97, 12.34,
      11.16, 10.10, 9.181, 8.386, 7.699, 7.104, 6.586, 6.132, 5.732,
      5.378, 5.063, 4.781, 4.528, 4.298, 3.373, 2.749, 2.309, 1.990,
      1.748, 1.559, 1.408, 1.284, 1.086, 0.942, 0.832, 0.746, 0.678,
      0.622, 0.576, 0.537, 0.503, 0.473, 0.448, 0.425, 0.405, 0.387,
  };
}

inline const char *MaterialName(ThicknessConv::Material material)
{
  using Material = ThicknessConv::Material;
  switch (material) {
    case Material::kCarbon: return "carbon";
    case Material::kCH2:    return "ch2";
    case Material::kCD2:    return "cd2";
  }
  return "unknown";
}

inline double LogLogInterpolate(double energy, const double *table)
{
  using namespace ThicknessConv;
  const double e = std::clamp(energy, 0.001, 10.0);

  int idx = 0;
  if (e <= kEnergy[0]) {
    idx = 0;
  }
  else if (e >= kEnergy[kNTable - 1]) {
    idx = kNTable - 2;
  }
  else {
    const double *upper = std::upper_bound(kEnergy, kEnergy + kNTable, e);
    idx = std::max(0, int(upper - kEnergy) - 1);
  }

  const double e1 = kEnergy[idx];
  const double e2 = kEnergy[idx + 1];
  const double s1 = table[idx];
  const double s2 = table[idx + 1];
  const double t = (std::log(e) - std::log(e1)) / (std::log(e2) - std::log(e1));
  return std::exp(std::log(s1) + t * (std::log(s2) - std::log(s1)));
}

inline double AStarMassStoppingPower(ThicknessConv::Material material,
                                     double energyMeV)
{
  using namespace ThicknessConv;
  const double sC = LogLogInterpolate(energyMeV, kStopCarbon);
  const double sH = LogLogInterpolate(energyMeV, kStopHydrogen);
  switch (material) {
    case Material::kCarbon: return sC;
    case Material::kCH2:    return kCH2WeightC * sC + kCH2WeightH * sH;
    case Material::kCD2:    return kCD2WeightC * sC + kCD2WeightD * sH;
  }
  return sC;
}

inline double MaterialDensity(ThicknessConv::Material material)
{
  using namespace ThicknessConv;
  switch (material) {
    case Material::kCarbon: return kDensityCarbon;
    case Material::kCH2:    return kDensityCH2;
    case Material::kCD2:    return kDensityCD2;
  }
  return kDensityCD2;
}

inline double EnergyLossPerUm(ThicknessConv::Material material, double energyMeV)
{
  return AStarMassStoppingPower(material, energyMeV)
       * MaterialDensity(material) * 0.1;
}

inline double EnergyToThickness(ThicknessConv::Material material, double eMeV)
{
  using namespace ThicknessConv;
  if (eMeV >= kE0) return 0.0;
  const double eMin = std::max(eMeV, 0.001);
  const int nSteps = 800;
  const double h = (kE0 - eMin) / nSteps;
  if (h <= 0) return 0.0;

  double sum = 0.0;
  for (int i = 0; i <= nSteps; ++i) {
    const double e = eMin + i * h;
    const double weight = (i == 0 || i == nSteps) ? 1.0 : (i % 2 == 0 ? 2.0 : 4.0);
    sum += weight / std::max(EnergyLossPerUm(material, e), 1.e-12);
  }
  return h * sum / 3.0;
}

// Backward-compatible default: CD2.
inline double EnergyToThickness(double eMeV)
{
  return EnergyToThickness(ThicknessConv::Material::kCD2, eMeV);
}

inline double ThicknessJacobian(ThicknessConv::Material material, double eMeV)
{
  if (eMeV <= 0) return 0.0;
  return 1.0 / std::max(EnergyLossPerUm(material, eMeV), 1.e-12);
}

inline double ThicknessJacobian(double eMeV)
{
  return ThicknessJacobian(ThicknessConv::Material::kCD2, eMeV);
}

inline TGraph *MakeAStarThicknessGraph(ThicknessConv::Material material,
                                       const char *name,
                                       int nPoints = 300,
                                       double eMin = 0.1,
                                       double eMax = ThicknessConv::kE0)
{
  auto graph = new TGraph(nPoints);
  graph->SetName(name);
  graph->SetTitle(name);
  for (int i = 0; i < nPoints; ++i) {
    const double e = eMin + (eMax - eMin) * i / std::max(nPoints - 1, 1);
    graph->SetPoint(i, e, EnergyToThickness(material, e));
  }
  return graph;
}

inline double MaterialThicknessRange(ThicknessConv::Material material)
{
  return std::min(ThicknessConv::kThicknessRange,
                  1.05 * EnergyToThickness(material, 0.1));
}

inline TH1D *MakeThicknessHist(const TH1D *hE, const std::string &name,
                               int nBins = 100,
                               double xMax = ThicknessConv::kThicknessRange,
                               ThicknessConv::Material material = ThicknessConv::Material::kCD2)
{
  using namespace ThicknessConv;
  auto hT = new TH1D(name.c_str(),
                     (std::string(hE->GetTitle()) +
                      ";Target Thickness (#mum);Counts").c_str(),
                     nBins, 0.0, xMax);
  hT->SetDirectory(nullptr);
  hT->SetLineWidth(2);

  const int nE = hE->GetNbinsX();
  for (int i = 1; i <= nE; ++i) {
    const double e = hE->GetBinCenter(i);
    const double c = hE->GetBinContent(i);
    if (c <= 0 || e <= 0) continue;
    const double x = EnergyToThickness(material, std::min(e, kE0));
    const int tb = hT->FindBin(x);
    if (tb >= 1 && tb <= nBins)
      hT->Fill(x, c);
  }
  return hT;
}

#endif // THICKNESS_CONVERTER_H
