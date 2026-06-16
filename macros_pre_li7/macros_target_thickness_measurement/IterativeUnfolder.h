#ifndef ITERATIVE_UNFOLDER_H
#define ITERATIVE_UNFOLDER_H

/**
 * IterativeUnfolder
 *
 * Iterative unfolding following the method of:
 *   "4.8 Spectrum Corrections" (iterative MC correction technique)
 *
 * Algorithm per iteration k:
 *   1. Smear current estimate:  h_sim = h_true^k ⊗ IRF
 *   2. Correction factor:       C(x)  = h_meas(x) / h_sim(x)
 *   3. Update:                  h_true^(k+1) = h_true^k * C
 *
 * The IRF is modeled as a Gaussian with width sigma_irf (from blank measurement).
 * This is equivalent to Richardson-Lucy deconvolution when the IRF is symmetric.
 *
 * Usage:
 *   IterativeUnfolder unfolder(hMeasured, sigmaIRF);
 *   TH1D* hUnfolded = unfolder.Unfold(10);  // 10 iterations
 */

#include <cmath>
#include <string>
#include <vector>

#include "TH1D.h"

class IterativeUnfolder {
public:
  // ── Constructor ─────────────────────────────────────────────────────────────
  // hMeas   : measured histogram (will not be modified)
  // sigmaIRF: Gaussian sigma of the instrument response (from blank)
  IterativeUnfolder(const TH1D *hMeas, double sigmaIRF)
      : fSigmaIRF(sigmaIRF)
  {
    // Deep-copy the measured histogram
    fHMeas = static_cast<TH1D *>(hMeas->Clone("_unfolder_meas"));
    fHMeas->SetDirectory(nullptr);

    // Initialise estimate with the measured histogram
    fHCurrent = static_cast<TH1D *>(hMeas->Clone("_unfolder_current"));
    fHCurrent->SetDirectory(nullptr);

    fNBins = hMeas->GetNbinsX();
  }

  ~IterativeUnfolder()
  {
    delete fHMeas;
    delete fHCurrent;
  }

  // ── Main interface ───────────────────────────────────────────────────────────
  // Performs nIter iterations and returns a NEW histogram (caller owns it).
  // histName / histTitle are applied to the returned histogram.
  TH1D *Unfold(int nIter = 5,
               const char *histName  = "h_unfolded",
               const char *histTitle = "Unfolded")
  {
    // Reset estimate to measured distribution before starting
    for (int i = 1; i <= fNBins; ++i)
      fHCurrent->SetBinContent(i, fHMeas->GetBinContent(i));

    for (int iter = 0; iter < nIter; ++iter)
      DoOneIteration();

    TH1D *result = static_cast<TH1D *>(fHCurrent->Clone(histName));
    result->SetDirectory(nullptr);
    result->SetTitle(histTitle);
    return result;
  }

  // ── Auto stopping: iterate until chi2 stops improving ───────────────────────
  // Stops when chi2 increases by more than chi2Tolerance compared to best so far.
  // maxIter: safety limit.
  // chi2History: if non-null, filled with chi2 per iteration (caller owns).
  // Returns the unfolded histogram at the best iteration (caller owns it).
  TH1D *UnfoldAuto(int maxIter = 30,
                   double chi2Tolerance = 0.0,
                   std::vector<double> *chi2History = nullptr,
                   const char *histName  = "h_unfolded_auto",
                   const char *histTitle = "Unfolded (auto)")
  {
    // Reset
    for (int i = 1; i <= fNBins; ++i)
      fHCurrent->SetBinContent(i, fHMeas->GetBinContent(i));

    double bestChi2  = 1.e30;
    int    bestIter  = 0;
    TH1D  *bestSnap  = nullptr;

    for (int iter = 0; iter < maxIter; ++iter) {
      DoOneIteration();

      // chi2 = sum_i [ (h_sim_i - h_meas_i)^2 / h_meas_i ]
      // where h_sim = h_current ⊗ IRF
      TH1D *hSim = Convolve(fHCurrent, "_unfolder_chi2_sim");
      double chi2 = 0.0;
      int    ndf  = 0;
      for (int i = 1; i <= fNBins; ++i) {
        const double meas = fHMeas->GetBinContent(i);
        if (meas <= 0) continue;
        const double sim  = hSim->GetBinContent(i);
        chi2 += (sim - meas) * (sim - meas) / meas;
        ++ndf;
      }
      delete hSim;

      if (chi2History) chi2History->push_back(chi2);

      if (chi2 < bestChi2 - chi2Tolerance) {
        bestChi2 = chi2;
        bestIter = iter + 1;
        delete bestSnap;
        bestSnap = static_cast<TH1D *>(fHCurrent->Clone("_unfolder_best"));
        bestSnap->SetDirectory(nullptr);
      } else {
        // chi2 stopped improving — stop here
        break;
      }
    }

    fBestIter  = bestIter;
    fBestChi2  = bestChi2;

    TH1D *result = static_cast<TH1D *>(bestSnap->Clone(histName));
    result->SetDirectory(nullptr);
    result->SetTitle(histTitle);
    delete bestSnap;
    return result;
  }

  int    GetBestIter() const { return fBestIter; }
  double GetBestChi2() const { return fBestChi2; }

  // Retrieve the smeared version of the current estimate (for diagnostics)
  TH1D *GetSmeared(const char *name = "h_smeared_estimate") const
  {
    return Convolve(fHCurrent, name);
  }

  // Change the IRF sigma after construction
  void SetSigmaIRF(double sigma) { fSigmaIRF = sigma; }
  double GetSigmaIRF()     const { return fSigmaIRF; }

private:
  // ── One iteration ────────────────────────────────────────────────────────────
  void DoOneIteration()
  {
    // Step 1: smear the current estimate
    TH1D *hSim = Convolve(fHCurrent, "_unfolder_sim");

    // Step 2 & 3: correction factor, then update bin by bin
    for (int i = 1; i <= fNBins; ++i) {
      const double sim = hSim->GetBinContent(i);
      if (sim <= 0) continue;

      const double meas    = fHMeas->GetBinContent(i);
      const double current = fHCurrent->GetBinContent(i);
      const double corrected = current * (meas / sim);
      fHCurrent->SetBinContent(i, std::max(corrected, 0.0));

      // Propagate errors (approximate: scale by same correction factor)
      const double errCurrent = fHCurrent->GetBinError(i);
      fHCurrent->SetBinError(i, errCurrent * (meas / sim));
    }

    delete hSim;
  }

  // ── Discrete Gaussian convolution ───────────────────────────────────────────
  // For each output bin j, accumulates contributions from all input bins i
  // weighted by a Gaussian kernel G(x_j - x_i ; 0, sigmaIRF).
  TH1D *Convolve(const TH1D *hIn, const char *name) const
  {
    TH1D *hOut = static_cast<TH1D *>(hIn->Clone(name));
    hOut->SetDirectory(nullptr);
    hOut->Reset();

    const double sigma = fSigmaIRF;
    if (sigma <= 0) {
      // No smearing: just copy
      for (int i = 1; i <= fNBins; ++i)
        hOut->SetBinContent(i, hIn->GetBinContent(i));
      return hOut;
    }

    // Truncate kernel at ±5σ for efficiency
    const double xMin = hIn->GetXaxis()->GetXmin();
    const double xMax = hIn->GetXaxis()->GetXmax();
    const double binW = (xMax - xMin) / fNBins;
    const int    kHalf = static_cast<int>(std::ceil(5.0 * sigma / binW));

    for (int j = 1; j <= fNBins; ++j) {
      const double xj = hIn->GetBinCenter(j);
      double sum    = 0.0;
      double weight = 0.0;

      const int iLo = std::max(1,        j - kHalf);
      const int iHi = std::min(fNBins,   j + kHalf);

      for (int i = iLo; i <= iHi; ++i) {
        const double xi = hIn->GetBinCenter(i);
        const double dx = xj - xi;
        const double w  = std::exp(-0.5 * dx * dx / (sigma * sigma));
        sum    += hIn->GetBinContent(i) * w;
        weight += w;
      }

      hOut->SetBinContent(j, weight > 0 ? sum / weight : 0.0);
    }

    // Normalise to preserve total counts
    const double inIntegral  = hIn->Integral();
    const double outIntegral = hOut->Integral();
    if (outIntegral > 0 && inIntegral > 0)
      hOut->Scale(inIntegral / outIntegral);

    return hOut;
  }

  // ── Members ──────────────────────────────────────────────────────────────────
  TH1D  *fHMeas    = nullptr;
  TH1D  *fHCurrent = nullptr;
  double fSigmaIRF = 0.0;
  int    fNBins    = 0;
  int    fBestIter = 0;
  double fBestChi2 = 0.0;
};

#endif // ITERATIVE_UNFOLDER_H
