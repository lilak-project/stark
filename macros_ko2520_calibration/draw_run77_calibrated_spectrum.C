#include "TCanvas.h"
#include "TFile.h"
#include "TF1.h"
#include "TGraph.h"
#include "TH1.h"
#include "TH1D.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLine.h"
#include "TObjArray.h"
#include "TStyle.h"
#include "TSystem.h"

#include <fstream>
#include <vector>

namespace {
void collect_energy_histograms(TObject *object, std::vector<TH1 *> &histograms)
{
    if (object == nullptr)
        return;

    auto histogram = dynamic_cast<TH1 *>(object);
    if (histogram != nullptr) {
        TString name = histogram->GetName();
        if (histogram->GetDimension() == 1 && name.BeginsWith("stage6_energy_"))
            histograms.push_back(histogram);
        return;
    }

    auto array = dynamic_cast<TObjArray *>(object);
    if (array == nullptr)
        return;
    for (auto item : *array)
        collect_energy_histograms(item, histograms);
}

double peak_position(TH1 *histogram, double low, double high)
{
    auto first = histogram->GetXaxis()->FindBin(low);
    auto last = histogram->GetXaxis()->FindBin(high);
    auto maximumBin = first;
    for (auto bin = first + 1; bin <= last; ++bin)
        if (histogram->GetBinContent(bin) > histogram->GetBinContent(maximumBin))
            maximumBin = bin;
    return histogram->GetBinContent(maximumBin) > 0
        ? histogram->GetBinCenter(maximumBin) : -1;
}

void draw_reference(double energy, double maximum, Color_t color)
{
    auto line = new TLine(energy, 0, energy, maximum);
    line->SetLineColor(color);
    line->SetLineStyle(2);
    line->SetLineWidth(2);
    line->Draw();
}

void draw_horizontal_reference(double energy, double maximum, Color_t color)
{
    auto line = new TLine(0, energy, maximum, energy);
    line->SetLineColor(color);
    line->SetLineStyle(2);
    line->SetLineWidth(2);
    line->Draw();
}
}

void draw_run77_calibrated_spectrum()
{
    const double gdEnergy = 3.1822;
    const double amEnergy = 5.486;
    const char *inputName = "ko2520_0077.si_calibration_stage_6.root";
    const char *outputPath = "data_calibration/run_0077/validation";

    TFile input(inputName, "READ");
    if (input.IsZombie()) {
        Error("draw_run77_calibrated_spectrum", "Cannot open %s", inputName);
        return;
    }

    auto top = dynamic_cast<TObjArray *>(input.Get("top"));
    if (top == nullptr) {
        Error("draw_run77_calibrated_spectrum", "LKDrawingGroup 'top' was not found");
        return;
    }

    std::vector<TH1 *> channelHistograms;
    collect_energy_histograms(top, channelHistograms);
    if (channelHistograms.empty()) {
        Error("draw_run77_calibrated_spectrum", "No stage-6 energy histograms were found");
        return;
    }

    auto sum = new TH1D("run77_calibrated_sum",
        "Run 77 calibrated silicon spectrum;Energy (MeV);Counts / 0.04 MeV",
        500, 0, 20);
    sum->SetDirectory(nullptr);
    sum->Sumw2();

    std::vector<double> channelIndex;
    std::vector<double> gdPeaks;
    std::vector<double> amPeaks;
    std::vector<TString> channelNames;
    int usedChannels = 0;
    for (auto histogram : channelHistograms) {
        if (histogram->GetEntries() < 100)
            continue;
        sum->Add(histogram);
        auto gdPeak = peak_position(histogram, 2.8, 3.55);
        auto amPeak = peak_position(histogram, 5.1, 5.85);
        if (gdPeak > 0 && amPeak > 0) {
            channelIndex.push_back(usedChannels);
            gdPeaks.push_back(gdPeak);
            amPeaks.push_back(amPeak);
            channelNames.push_back(histogram->GetName());
        }
        ++usedChannels;
    }

    TF1 gdFit("gd_fit", "gaus(0)+pol1(3)", 2.9, 3.42);
    gdFit.SetParameters(sum->GetMaximum(), gdEnergy, 0.05, 0, 0);
    sum->Fit(&gdFit, "RQ0");
    TF1 amFit("am_fit", "gaus(0)+pol1(3)", 5.2, 5.72);
    amFit.SetParameters(sum->GetMaximum(), amEnergy, 0.05, 0, 0);
    sum->Fit(&amFit, "RQ0+");

    gSystem->mkdir(outputPath, true);
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);

    TCanvas canvas("validation", "Run 77 calibration validation", 1500, 1050);
    canvas.Divide(2, 2, 0.002, 0.002);

    canvas.cd(1);
    gPad->SetLogy();
    gPad->SetGridx();
    sum->SetLineColor(kBlack);
    sum->SetLineWidth(2);
    sum->GetXaxis()->SetRangeUser(0, 7);
    sum->Draw("hist");
    draw_reference(gdEnergy, sum->GetMaximum(), kBlue + 1);
    draw_reference(amEnergy, sum->GetMaximum(), kRed + 1);
    TLegend legend(0.55, 0.67, 0.88, 0.88);
    legend.SetBorderSize(0);
    legend.SetFillStyle(0);
    legend.AddEntry(sum, Form("sum of %d channels", usedChannels), "l");
    legend.AddEntry((TObject *) nullptr, "^{148}Gd: 3.1822 MeV", "");
    legend.AddEntry((TObject *) nullptr, "^{241}Am: 5.486 MeV", "");
    legend.Draw();

    canvas.cd(2);
    auto gdZoom = dynamic_cast<TH1D *>(sum->Clone("gd_zoom"));
    gdZoom->SetTitle("^{148}Gd peak;Energy (MeV);Counts / 0.04 MeV");
    gdZoom->GetXaxis()->SetRangeUser(2.75, 3.6);
    gdZoom->Draw("hist");
    gdFit.SetLineColor(kBlue + 1);
    gdFit.Draw("same");
    draw_reference(gdEnergy, gdZoom->GetMaximum() * 1.05, kBlue + 1);
    TLatex text;
    text.SetNDC();
    text.SetTextSize(0.04);
    text.DrawLatex(0.14, 0.84, Form("fit mean = %.4f MeV", gdFit.GetParameter(1)));
    text.DrawLatex(0.14, 0.78, Form("sigma = %.4f MeV", TMath::Abs(gdFit.GetParameter(2))));

    canvas.cd(3);
    auto amZoom = dynamic_cast<TH1D *>(sum->Clone("am_zoom"));
    amZoom->SetTitle("^{241}Am peak;Energy (MeV);Counts / 0.04 MeV");
    amZoom->GetXaxis()->SetRangeUser(5.0, 5.95);
    amZoom->Draw("hist");
    amFit.SetLineColor(kRed + 1);
    amFit.Draw("same");
    draw_reference(amEnergy, amZoom->GetMaximum() * 1.05, kRed + 1);
    text.DrawLatex(0.14, 0.84, Form("fit mean = %.4f MeV", amFit.GetParameter(1)));
    text.DrawLatex(0.14, 0.78, Form("sigma = %.4f MeV", TMath::Abs(amFit.GetParameter(2))));

    canvas.cd(4);
    gPad->SetGridy();
    TGraph gdGraph(channelIndex.size(), channelIndex.data(), gdPeaks.data());
    TGraph amGraph(channelIndex.size(), channelIndex.data(), amPeaks.data());
    gdGraph.SetTitle("Channel-by-channel peak positions;Calibrated channel index;Peak energy (MeV)");
    gdGraph.SetMarkerStyle(20);
    gdGraph.SetMarkerSize(0.55);
    gdGraph.SetMarkerColor(kBlue + 1);
    gdGraph.SetMinimum(2.7);
    gdGraph.SetMaximum(6.0);
    gdGraph.Draw("AP");
    amGraph.SetMarkerStyle(20);
    amGraph.SetMarkerSize(0.55);
    amGraph.SetMarkerColor(kRed + 1);
    amGraph.Draw("P same");
    auto lastChannel = channelIndex.empty() ? 1. : channelIndex.back();
    draw_horizontal_reference(gdEnergy, lastChannel, kBlue + 1);
    draw_horizontal_reference(amEnergy, lastChannel, kRed + 1);
    TLegend graphLegend(0.72, 0.72, 0.89, 0.88);
    graphLegend.SetBorderSize(0);
    graphLegend.SetFillStyle(0);
    graphLegend.AddEntry(&gdGraph, "^{148}Gd", "p");
    graphLegend.AddEntry(&amGraph, "^{241}Am", "p");
    graphLegend.Draw();

    canvas.SaveAs(Form("%s/run77_calibrated_spectrum.svg", outputPath));
    canvas.SaveAs(Form("%s/run77_calibrated_spectrum.pdf", outputPath));

    TFile output(Form("%s/run77_calibrated_spectrum.root", outputPath), "RECREATE");
    sum->Write();
    gdGraph.Write("gd_channel_peak_positions");
    amGraph.Write("am_channel_peak_positions");
    canvas.Write();
    output.Close();

    std::ofstream report(Form("%s/run77_calibration_validation.txt", outputPath));
    report << "input " << inputName << "\n";
    report << "histograms_found " << channelHistograms.size() << "\n";
    report << "channels_summed " << usedChannels << "\n";
    report << "channels_with_both_peaks " << channelIndex.size() << "\n";
    report << "gd_reference_MeV " << gdEnergy << "\n";
    report << "gd_fit_mean_MeV " << gdFit.GetParameter(1) << "\n";
    report << "gd_fit_sigma_MeV " << TMath::Abs(gdFit.GetParameter(2)) << "\n";
    report << "am_reference_MeV " << amEnergy << "\n";
    report << "am_fit_mean_MeV " << amFit.GetParameter(1) << "\n";
    report << "am_fit_sigma_MeV " << TMath::Abs(amFit.GetParameter(2)) << "\n";
    report.close();

    Info("draw_run77_calibrated_spectrum",
        "summed %d channels; Gd mean %.4f MeV, Am mean %.4f MeV",
        usedChannels, gdFit.GetParameter(1), amFit.GetParameter(1));
}
