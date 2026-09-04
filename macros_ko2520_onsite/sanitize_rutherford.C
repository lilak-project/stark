#include "TFile.h"
#include "TH1.h"
#include "TH1D.h"
#include "TString.h"
#include "TSystem.h"

static TH1D* cleanHistogram(TH1* source)
{
  if (!source) return nullptr;
  TH1D* clean = new TH1D(source->GetName(),source->GetTitle(),
                         source->GetNbinsX(),
                         source->GetXaxis()->GetXmin(),
                         source->GetXaxis()->GetXmax());
  clean->SetDirectory(nullptr);
  clean->SetLineColor(source->GetLineColor());
  clean->SetLineWidth(source->GetLineWidth());
  clean->SetMarkerColor(source->GetMarkerColor());
  clean->SetMarkerStyle(source->GetMarkerStyle());
  clean->SetMarkerSize(source->GetMarkerSize());
  for (int bin=0;bin<=source->GetNbinsX()+1;bin++) {
    clean->SetBinContent(bin,source->GetBinContent(bin));
    clean->SetBinError(bin,source->GetBinError(bin));
  }
  return clean;
}

void sanitize_rutherford(const char* inputFileName, const char* outputFileName)
{
  const Bool_t oldAddDirectory = TH1::AddDirectoryStatus();
  TH1::AddDirectory(kFALSE);
  TFile input(inputFileName,"READ");
  if (input.IsZombie()) {
    Error("sanitize_rutherford","cannot open %s",inputFileName);
    gSystem->Exit(2);
  }

  TH1D* data = cleanHistogram((TH1*)input.Get("w1_2_junction_pattern_data"));
  TH1D* fit = cleanHistogram((TH1*)input.Get("w1_2_rutherford_fit"));
  TH1D* ratio = cleanHistogram((TH1*)input.Get("w1_2_data_over_fit"));
  if (!data || !fit || !ratio) {
    Error("sanitize_rutherford","required histograms not found in %s",inputFileName);
    gSystem->Exit(3);
  }
  input.Close();

  TFile output(outputFileName,"RECREATE");
  if (output.IsZombie()) {
    Error("sanitize_rutherford","cannot create %s",outputFileName);
    gSystem->Exit(4);
  }
  data->Write();
  fit->Write();
  ratio->Write();
  output.Close();
  TH1::AddDirectory(oldAddDirectory);
  Printf("sanitised Rutherford histograms: %s",outputFileName);
}
