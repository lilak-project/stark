#include "LKBinning.h"
#include "LKDrawingGroup.h"
#include "LKParameterContainer.h"
#include "TCanvas.h"
#include "TH2D.h"
#include "TGraphErrors.h"
#include "TLegend.h"
#include "TBox.h"
#include "TError.h"
#include "TMultiGraph.h"
#include "TString.h"
#include "TSystem.h"
#include "TAxis.h"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <map>
#include <set>
#include <sstream>
#include <tuple>
#include <vector>

struct DepVPoint {
    double voltage = 0;
    double energy = 0;
    double error = 0;
    double rawMean = 0;
};

struct DepVGraph {
    TGraphErrors *graph = nullptr;
    TString option = "L SAME";
};

using DepVKey = std::tuple<TString, int, int, int>; // detector type, idx, side, channel

LKBinning bnnADC(500, 0, 2800);

std::map<std::tuple<TString, int>, double> MakeBestBiasMap()
{
    // Initial values are selected from the minimum average junction resolution.
    // Edit these values directly when the preferred operating bias changes.
    return {
        {std::make_tuple("BB10", 1), 25.0},
        {std::make_tuple("BB10", 2), 40.0},
        {std::make_tuple("BB10", 3), 25.0},
        {std::make_tuple("BB10", 4), 25.0},
        {std::make_tuple("BB10", 5), 25.0},
        {std::make_tuple("BB10", 6), 35.0},
        {std::make_tuple("BB10", 7), 20.0},
        {std::make_tuple("BB10", 8), 20.0},
        {std::make_tuple("QQQ5", 1), 112.0},
        {std::make_tuple("QQQ5", 2), 112.0},
        {std::make_tuple("QQQ5", 3), 112.0},
        {std::make_tuple("QQQ5", 4), 112.0},
    };
}

TString DetectorTitle(TString detType, int detectorNumber)
{
    return Form("%s #%d", detType.Data(), detectorNumber);
}

TString SideSuffix(int side)
{
    return side == 0 ? "J" : "O";
}

int DepVColor(int i)
{
    static const int colors[] = {
        kBlack, kRed + 1, kBlue + 1, kGreen + 2, kMagenta + 1, kCyan + 2,
        kOrange + 7, kViolet + 1, kAzure + 7, kSpring + 5, kPink + 7, kTeal + 3
    };
    return colors[i % (sizeof(colors) / sizeof(colors[0]))];
}

int DepVLineStyle(int i)
{
    static const int styles[] = {
        1, 2, 3, 4, 5, 6, 7, 8, 9, 10
    };
    return styles[(i / 5) % (sizeof(styles) / sizeof(styles[0]))];
}

bool ParseDepVFileName(TString fileName, TString &detType, double &voltage)
{
    if (!fileName.EndsWith(".dat"))
        return false;
    if (fileName.EndsWith("_lc.dat"))
        return false;
    if (!fileName.EndsWith("_er.dat"))
        return false;

    fileName.ReplaceAll("_er.dat", "");
    const auto underscore = fileName.Last('_');
    if (underscore == kNPOS)
        return false;

    detType = fileName(0, underscore);
    TString voltageText = fileName(underscore + 1, fileName.Length() - underscore - 1);
    voltage = voltageText.Atof();
    return !detType.IsNull();
}

bool ParseLeakageCurrentFileName(TString fileName, TString &detType, int &idx)
{
    if (!fileName.EndsWith("_lc.dat"))
        return false;

    fileName.ReplaceAll("_lc.dat", "");
    auto idxPos = fileName.Index("_N");
    int offset = 2;
    if (idxPos == kNPOS) {
        idxPos = fileName.Index("_idx");
        offset = 4;
    }
    if (idxPos == kNPOS)
        return false;

    detType = fileName(0, idxPos);
    TString idxText = fileName(idxPos + offset, fileName.Length() - idxPos - offset);
    idx = idxText.Atoi();
    return !detType.IsNull() && idx >= 0;
}

void DrawMultiGraph(
LKDrawingGroup* group,
    TString title,
    const std::vector<DepVGraph> &graphs,
    TString outputName,
    const LKBinning &bnnVoltage,
    const LKBinning &bnnY,
    TString yTitle,
    bool saveImages,
    double bestBias = -1,
    double bestBiasBandHalfWidth = 2.5)
{
    if (graphs.empty())
        return;

    //e_debug << outputName << endl;
    auto drawing = group -> CreateDrawing(outputName);

    //auto canvas = new TCanvas(Form("c_%s", outputName.Data()), title, 1200, 800);
    auto frameBnnX = bnnVoltage;
    auto frameBnnY = bnnY;
    auto frame = (frameBnnX * frameBnnY).NewH2(Form("frame_%s", outputName.Data()), title + yTitle);
    frame->SetStats(false);
    //frame->Draw("AXIS");
    drawing -> Add(frame, "hist");

    if (bestBias >= 0 && std::isfinite(bestBias) && bestBiasBandHalfWidth > 0) {
        auto bestBiasBand = new TBox(
            bestBias - bestBiasBandHalfWidth,
            bnnY.GetX1(),
            bestBias + bestBiasBandHalfWidth,
            bnnY.GetX2()
        );
        bestBiasBand->SetFillColorAlpha(kYellow, 0.35);
        bestBiasBand->SetLineColor(kYellow + 1);
        bestBiasBand->SetLineStyle(1);
        bestBiasBand->SetLineWidth(1);
        drawing -> Add(bestBiasBand, "f");
    }

    auto legend = new TLegend(0.15, 0.12, 0.85, 0.33);
    legend->SetBorderSize(0);
    legend->SetFillStyle(0);
    if (0) {}
    //if (outputName.Contains("QQQ5_side0")) legend->SetNColumns(6);
    //else if (graphs.size() > 80) legend->SetNColumns(10);
    //else if (graphs.size() > 60) legend->SetNColumns(8);
    //else if (graphs.size() > 35) legend->SetNColumns(6);
    else if (graphs.size() > 18) legend->SetNColumns(8);
    else legend->SetNColumns(4);
    legend->SetTextSize(0.045);

    int iLine = 0;
    for (auto const &graphSpec : graphs) {
        auto graph = graphSpec.graph;
        if (graph == nullptr)
            continue;
        if (graphSpec.option.Contains("L")) {
            graph->SetLineStyle(DepVLineStyle(iLine));
            ++iLine;
        }
        //graph->Draw("PL SAME");
        drawing -> Add(graph, graphSpec.option);
        legend->AddEntry(graph, graph->GetTitle(), "lp");
    }

    frame->GetXaxis()->SetTitleOffset(1.2);
    frame->GetYaxis()->SetTitleOffset(1.7);
    //legend->Draw();
    drawing -> Add(legend);
    //drawing -> SetCreateLegend();
    //canvas->Modified();
    //canvas->Update();

    drawing -> SetLeftMargin(0.12);
    drawing -> SetGridx();
    //drawing -> SetGridy();

    if (saveImages) {
        //gSystem->mkdir("figures_depv_energy", true);
        //canvas->SaveAs(Form("figures_depv_energy/%s.png", outputName.Data()));
        //canvas->SaveAs(Form("figures_depv_energy/%s.pdf", outputName.Data()));
    }
}

TGraphErrors *CreateVoltageAverageGraph(
    TString name,
    TString title,
    const std::vector<TGraphErrors *> &graphs)
{
    std::map<double, std::vector<double>> valuesByVoltage;
    for (auto graph : graphs) {
        if (graph == nullptr)
            continue;
        for (int i = 0; i < graph->GetN(); ++i) {
            double x = 0;
            double y = 0;
            graph->GetPoint(i, x, y);
            valuesByVoltage[x].push_back(y);
        }
    }

    auto averageGraph = new TGraphErrors(valuesByVoltage.size());
    averageGraph->SetName(name);
    averageGraph->SetTitle(title);
    averageGraph->SetLineColor(kBlack);
    averageGraph->SetMarkerColor(kBlack);
    averageGraph->SetMarkerStyle(20);
    averageGraph->SetMarkerSize(0.8);

    int iPoint = 0;
    for (auto const &entry : valuesByVoltage) {
        const auto &values = entry.second;
        if (values.empty())
            continue;

        double sum = 0;
        for (auto value : values)
            sum += value;
        const double mean = sum / values.size();

        double variance = 0;
        for (auto value : values)
            variance += (value - mean) * (value - mean);
        const double sigma = values.size() > 1 ? std::sqrt(variance / (values.size() - 1)) : 0;
        const double error = values.size() > 0 ? sigma / std::sqrt(double(values.size())) : 0;

        averageGraph->SetPoint(iPoint, entry.first, mean);
        averageGraph->SetPointError(iPoint, 0, error);
        ++iPoint;
    }
    return averageGraph;
}

std::map<std::tuple<TString, int>, TGraphErrors *> LoadLeakageCurrentGraphs(
    TString dataDir,
    TString selectedDetType,
    int selectedDetIdx)
{
    std::map<std::tuple<TString, int>, TGraphErrors *> leakageGraphs;

    auto dir = gSystem->OpenDirectory(dataDir);
    if (dir == nullptr)
        return leakageGraphs;

    while (auto entry = gSystem->GetDirEntry(dir)) {
        TString fileName(entry);
        TString detType;
        int idx = -1;
        if (!ParseLeakageCurrentFileName(fileName, detType, idx))
            continue;
        if (!selectedDetType.IsNull() && detType != selectedDetType)
            continue;
        if (selectedDetIdx >= 0 && idx != selectedDetIdx)
            continue;

        std::ifstream input(Form("%s/%s", dataDir.Data(), fileName.Data()));
        if (!input.is_open())
            continue;

        std::vector<DepVPoint> points;
        TString line;
        while (line.ReadLine(input)) {
            line = line.Strip(TString::kBoth);
            if (line.IsNull() || line.BeginsWith("#"))
                continue;

            std::istringstream parser(line.Data());
            double voltage = 0;
            double current = 0;
            if (!(parser >> voltage >> current))
                continue;
            points.push_back({voltage, current, 0});
        }
        std::sort(points.begin(), points.end(), [](const auto &a, const auto &b) {
            return a.voltage < b.voltage;
        });

        auto graph = new TGraphErrors(points.size());
        graph->SetName(Form("g_lc_%s_idx%d", detType.Data(), idx));
        graph->SetTitle("leakage current");
        graph->SetLineColor(kBlack);
        graph->SetMarkerColor(kBlack);
        graph->SetMarkerStyle(20);
        graph->SetMarkerSize(0.8);
        for (int i = 0; i < int(points.size()); ++i) {
            graph->SetPoint(i, points[i].voltage, points[i].energy);
            graph->SetPointError(i, 0, 0);
        }
        leakageGraphs[std::make_tuple(detType, idx)] = graph;
    }
    gSystem->FreeDirectory(dir);
    return leakageGraphs;
}

void draw_depv_energy(
    TString dataDir = "data_depv",
    TString parameterFileName = "ana_depv_calibration.conf",
    bool saveImages = true,
    TString selectedDetType = "",
    //TString selectedDetType = "BB10",
    int selectedDetIdx = -1,
    bool draw_mean_graph = false,
    bool draw_resolution_graph = true,
    bool draw_raw_mean_graph = true
    //int selectedDetIdx = -1
    )
{
    auto par = LKParameterContainer(parameterFileName);
    std::map<TString, LKBinning> bnnVoltageByType = {
        //{"BB10", LKBinning(200, 0, 70)},
        //{"QQQ5", LKBinning(200, 0, 220)},
        {"BB10", LKBinning(200, 0, 100)},
        {"QQQ5", LKBinning(200, 0, 250)},
    };
    std::map<TString, LKBinning> bnnEnergyByType = {
        {"BB10", LKBinning(500, 0, 7)},
        {"QQQ5", LKBinning(500, 0, 7)},
    };
    std::map<TString, LKBinning> bnnRawMeanByType = {
        {"BB10", bnnADC},
        {"QQQ5", bnnADC},
    };
    std::map<TString, LKBinning> bnnResolutionByType = {
        {"BB10", LKBinning(500, 0, 0.15)},
        {"QQQ5", LKBinning(500, 0, 0.15)},
    };
    std::map<TString, LKBinning> bnnLeakageCurrentByType = {
        {"BB10", LKBinning(500, 0, 1)},
        {"QQQ5", LKBinning(500, 0, 120)},
    };
    double bestBiasBandHalfWidth = 2.5;
    auto bestBiasByDetector = MakeBestBiasMap();
    std::map<std::tuple<TString, int>, LKBinning> bnnResolutionByTypeSide = {
        {std::make_tuple("BB10", 0), LKBinning(500, 0, 0.15)},
        {std::make_tuple("BB10", 1), LKBinning(500, 0, 0.15)},
        {std::make_tuple("QQQ5", 0), LKBinning(500, 0, 0.15)},
        {std::make_tuple("QQQ5", 1), LKBinning(500, 0, 0.15)},
    };
    for (auto &entry : bnnVoltageByType)
        par.UpdatePar(entry.second, Form("bnnVoltage%s", entry.first.Data()));
    for (auto &entry : bnnEnergyByType)
        par.UpdatePar(entry.second, Form("bnnEnergy%s", entry.first.Data()));
    for (auto &entry : bnnRawMeanByType)
        par.UpdatePar(entry.second, Form("bnnRawMean%s", entry.first.Data()));
    for (auto &entry : bnnResolutionByType) {
        par.UpdatePar(entry.second, Form("bnnResolution%s", entry.first.Data()));
        bnnResolutionByTypeSide[std::make_tuple(entry.first, 0)] = entry.second;
        bnnResolutionByTypeSide[std::make_tuple(entry.first, 1)] = entry.second;
    }
    for (auto &entry : bnnLeakageCurrentByType)
        par.UpdatePar(entry.second, Form("bnnLeakageCurrent%s", entry.first.Data()));
    par.UpdatePar(bestBiasBandHalfWidth, "bestBiasBandHalfWidth");
    for (auto &entry : bnnResolutionByTypeSide) {
        auto detType = std::get<0>(entry.first);
        auto side = std::get<1>(entry.first);
        par.UpdatePar(entry.second, Form("bnnResolution%s%s", detType.Data(), SideSuffix(side).Data()));
    }

    selectedDetType = selectedDetType.Strip(TString::kBoth);
    selectedDetType.ToUpper();
    auto leakageCurrentGraphs = LoadLeakageCurrentGraphs(dataDir, selectedDetType, selectedDetIdx);

    std::map<DepVKey, std::vector<DepVPoint>> data;
    std::set<TString> detTypes;
    std::set<std::tuple<TString, int>> detectorSet;

    auto dir = gSystem->OpenDirectory(dataDir);
    if (dir == nullptr) {
        ::Error("draw_depv_energy", "Cannot open data directory %s", dataDir.Data());
        return;
    }

    while (auto entry = gSystem->GetDirEntry(dir)) {
        TString fileName(entry);
        TString detType;
        double voltage = 0;
        if (!ParseDepVFileName(fileName, detType, voltage))
            continue;
        if (!selectedDetType.IsNull() && detType != selectedDetType)
            continue;

        std::ifstream input(Form("%s/%s", dataDir.Data(), fileName.Data()));
        if (!input.is_open())
            continue;

        detTypes.insert(detType);
        TString line;
        while (line.ReadLine(input)) {
            line = line.Strip(TString::kBoth);
            if (line.IsNull() || line.BeginsWith("#"))
                continue;

            std::istringstream parser(line.Data());
            std::vector<double> values;
            double value = 0;
            while (parser >> value)
                values.push_back(value);

            int idx = -1;
            int side = -1;
            int channel = -1;
            double amp = 0;
            double mean = 0;
            double error = 0;
            double rawMean = 0;
            if (values.size() >= 11) {
                // Temporary CAAC format:
                // cobo asad aget chan det_number side channel corrected_fit_amp corrected_fit_energy corrected_fit_fwhm raw_fit_mean
                idx = int(values[4]);
                side = int(values[5]);
                channel = int(values[6]);
                amp = values[7];
                mean = values[8];
                error = values[9];
                rawMean = values[10];
            }
            else if (values.size() >= 10) {
                // Temporary CAAC format without raw_fit_mean.
                idx = int(values[4]);
                side = int(values[5]);
                channel = int(values[6]);
                amp = values[7];
                mean = values[8];
                error = values[9];
                rawMean = mean;
            }
            else if (values.size() >= 7) {
                // Current format:
                // det_number side channel corrected_fit_amp corrected_fit_energy corrected_fit_fwhm raw_fit_mean
                idx = int(values[0]);
                side = int(values[1]);
                channel = int(values[2]);
                amp = values[3];
                mean = values[4];
                error = values[5];
                rawMean = values[6];
            }
            else if (values.size() >= 6) {
                // Old format:
                // det_number side channel corrected_fit_amp corrected_fit_energy corrected_fit_fwhm
                idx = int(values[0]);
                side = int(values[1]);
                channel = int(values[2]);
                amp = values[3];
                mean = values[4];
                error = values[5];
                rawMean = mean;
            }
            else
                continue;
            if (selectedDetIdx >= 0 && idx != selectedDetIdx)
                continue;

            data[DepVKey(detType, idx, side, channel)].push_back({voltage, mean, error, rawMean});
            detectorSet.insert({detType, idx});
        }
    }
    gSystem->FreeDirectory(dir);

    for (auto &entry : data) {
        auto &points = entry.second;
        std::sort(points.begin(), points.end(), [](const auto &a, const auto &b) {
            return a.voltage < b.voltage;
        });
    }

    auto top = new LKDrawingGroup("top");

    int count = 0;
    /*
    for (auto const &detType : detTypes) {
        auto group = top -> CreateGroup(Form("%s_%d",detType.Data(),count++));
        group -> SetCanvasDivision(2,3);
        auto bnnVoltage = bnnVoltageByType.count(detType) ? bnnVoltageByType[detType] : LKBinning(200, 0, 220);
        auto bnnEnergy = bnnEnergyByType.count(detType) ? bnnEnergyByType[detType] : LKBinning(500, 0, 7);
        for (int side : {0, 1}) {
            auto bnnResolution = bnnResolutionByTypeSide.count(std::make_tuple(detType, side)) ? bnnResolutionByTypeSide[std::make_tuple(detType, side)] : LKBinning(500, 0, 0.15);
            std::vector<TGraphErrors *> meanChannelGraphs;
            std::vector<TGraphErrors *> resolutionChannelGraphs;
            int iGraph = 0;
            for (auto const &entry : data) {
                if (std::get<0>(entry.first) != detType || std::get<2>(entry.first) != side)
                    continue;

                const int idx = std::get<1>(entry.first);
                const int channel = std::get<3>(entry.first);
                const auto &points = entry.second;
                auto meanGraph = new TGraphErrors(points.size());
                meanGraph->SetName(Form("g_mean_%s_idx%d_side%d_ch%d", detType.Data(), idx, side, channel));
                meanGraph->SetTitle(Form("#%d s%d ch%d", idx, side, channel));
                meanGraph->SetLineColor(DepVColor(iGraph));
                auto resolutionGraph = new TGraphErrors(points.size());
                resolutionGraph->SetName(Form("g_resolution_%s_idx%d_side%d_ch%d", detType.Data(), idx, side, channel));
                resolutionGraph->SetTitle(Form("#%d s%d ch%d", idx, side, channel));
                resolutionGraph->SetLineColor(DepVColor(iGraph));
                for (int i = 0; i < int(points.size()); ++i) {
                    meanGraph->SetPoint(i, points[i].voltage, points[i].energy);
                    meanGraph->SetPointError(i, 0, points[i].error);
                    const double resolution = points[i].energy > 0 ? points[i].error / points[i].energy : 0;
                    resolutionGraph->SetPoint(i, points[i].voltage, resolution);
                    resolutionGraph->SetPointError(i, 0, 0);
                }
                meanChannelGraphs.push_back(meanGraph);
                resolutionChannelGraphs.push_back(resolutionGraph);
                ++iGraph;
            }
            std::vector<DepVGraph> meanGraphs;
            std::vector<DepVGraph> resolutionGraphs;
            for (auto graph : meanChannelGraphs)
                meanGraphs.push_back({graph, "L SAME"});
            for (auto graph : resolutionChannelGraphs)
                resolutionGraphs.push_back({graph, "L SAME"});
            meanGraphs.push_back({
                CreateVoltageAverageGraph(Form("g_mean_avg_%s_side%d", detType.Data(), side), "avrg.", meanChannelGraphs),
                "P SAME"
            });
            resolutionGraphs.push_back({
                CreateVoltageAverageGraph(Form("g_resolution_avg_%s_side%d", detType.Data(), side), "avrg.", resolutionChannelGraphs),
                "P SAME"
            });
            if (draw_mean_graph)
                DrawMultiGraph(
                    group,
                    Form("%s %s all channels mean", detType.Data(), (side==0?"junction":"ohmic")),
                    meanGraphs,
                    Form("%s_side%d_all_mean", detType.Data(), side),
                    bnnVoltage,
                    bnnEnergy,
                    "Corrected fit energy",
                    saveImages
                );
            if (draw_resolution_graph)
                DrawMultiGraph(
                    group,
                    Form("%s %s all channels resolution", detType.Data(), (side==0?"junction":"ohmic")),
                    resolutionGraphs,
                    Form("%s_side%d_all_resolution", detType.Data(), side),
                    bnnVoltage,
                    bnnResolution,
                    ";Bias [V];Energy resolution (FWHM / mean) [%]",
                    saveImages
                );
        }
    }
    */

    for (auto const &detector : detectorSet) {
        const auto detType = std::get<0>(detector);
        const int idx = std::get<1>(detector);
        const auto detectorTitle = DetectorTitle(detType, idx);
        auto group = top -> CreateGroup(Form("all_%s_%d",detType.Data(),count++));
        group -> SetCanvasDivision(2,3);
        //group -> SetCanvasSize(800,800);
        auto bnnVoltage = bnnVoltageByType.count(detType) ? bnnVoltageByType[detType] : LKBinning(200, 0, 220);
            auto bnnEnergy = bnnEnergyByType.count(detType) ? bnnEnergyByType[detType] : LKBinning(500, 0, 7);
            auto bnnRawMean = bnnRawMeanByType.count(detType) ? bnnRawMeanByType[detType] : bnnADC;
            auto bnnLeakageCurrent = bnnLeakageCurrentByType.count(detType) ? bnnLeakageCurrentByType[detType] : LKBinning(500, 0, 120);
            const auto bestBiasIt = bestBiasByDetector.find(std::make_tuple(detType, idx));
            const double bestBias = bestBiasIt != bestBiasByDetector.end() ? bestBiasIt->second : -1;
        for (int side : {0, 1}) {
            auto bnnResolution = bnnResolutionByTypeSide.count(std::make_tuple(detType, side)) ? bnnResolutionByTypeSide[std::make_tuple(detType, side)] : LKBinning(500, 0, 0.15);
            std::vector<TGraphErrors *> rawMeanChannelGraphs;
            std::vector<TGraphErrors *> meanChannelGraphs;
            std::vector<TGraphErrors *> resolutionChannelGraphs;
            int iGraph = 0;
            for (auto const &entry : data) {
                if (std::get<0>(entry.first) != detType || std::get<1>(entry.first) != idx || std::get<2>(entry.first) != side)
                    continue;

                const int channel = std::get<3>(entry.first);
                const auto &points = entry.second;
                auto rawMeanGraph = new TGraphErrors(points.size());
                rawMeanGraph->SetName(Form("g_raw_mean_%s_idx%d_side%d_ch%d_detail", detType.Data(), idx, side, channel));
                rawMeanGraph->SetTitle(Form("ch%d", channel));
                rawMeanGraph->SetLineColor(DepVColor(iGraph));
                auto meanGraph = new TGraphErrors(points.size());
                meanGraph->SetName(Form("g_mean_%s_idx%d_side%d_ch%d_detail", detType.Data(), idx, side, channel));
                meanGraph->SetTitle(Form("ch%d", channel));
                meanGraph->SetLineColor(DepVColor(iGraph));
                auto resolutionGraph = new TGraphErrors(points.size());
                resolutionGraph->SetName(Form("g_resolution_%s_idx%d_side%d_ch%d_detail", detType.Data(), idx, side, channel));
                resolutionGraph->SetTitle(Form("ch%d", channel));
                resolutionGraph->SetLineColor(DepVColor(iGraph));
                for (int i = 0; i < int(points.size()); ++i) {
                    rawMeanGraph->SetPoint(i, points[i].voltage, points[i].rawMean);
                    rawMeanGraph->SetPointError(i, 0, 0);
                    meanGraph->SetPoint(i, points[i].voltage, points[i].energy);
                    meanGraph->SetPointError(i, 0, points[i].error);
                    const double resolution = points[i].energy > 0 ? points[i].error / points[i].energy : 0;
                    resolutionGraph->SetPoint(i, points[i].voltage, resolution);
                    resolutionGraph->SetPointError(i, 0, 0);
                }
                rawMeanChannelGraphs.push_back(rawMeanGraph);
                meanChannelGraphs.push_back(meanGraph);
                resolutionChannelGraphs.push_back(resolutionGraph);
                ++iGraph;
            }
            std::vector<DepVGraph> rawMeanGraphs;
            std::vector<DepVGraph> meanGraphs;
            std::vector<DepVGraph> resolutionGraphs;
            for (auto graph : rawMeanChannelGraphs)
                rawMeanGraphs.push_back({graph, "L SAME"});
            for (auto graph : meanChannelGraphs)
                meanGraphs.push_back({graph, "L SAME"});
            for (auto graph : resolutionChannelGraphs)
                resolutionGraphs.push_back({graph, "L SAME"});
            rawMeanGraphs.push_back({
                CreateVoltageAverageGraph(Form("g_raw_mean_avg_%s_idx%d_side%d", detType.Data(), idx, side), "avrg.", rawMeanChannelGraphs),
                "P SAME"
            });
            meanGraphs.push_back({
                CreateVoltageAverageGraph(Form("g_mean_avg_%s_idx%d_side%d", detType.Data(), idx, side), "avrg.", meanChannelGraphs),
                "P SAME"
            });
            resolutionGraphs.push_back({
                CreateVoltageAverageGraph(Form("g_resolution_avg_%s_idx%d_side%d", detType.Data(), idx, side), "avrg.", resolutionChannelGraphs),
                "P SAME"
            });
            if (draw_raw_mean_graph)
                DrawMultiGraph(
                    group,
                    Form("%s %s raw mean", detectorTitle.Data(), (side==0?"junction":"ohmic")),
                    rawMeanGraphs,
                    Form("%s_idx%d_side%d_raw_mean", detType.Data(), idx, side),
                    bnnVoltage,
                    bnnRawMean,
                    ";Bias [V];Raw amplitude",
                    saveImages,
                    bestBias,
                    bestBiasBandHalfWidth
                );
            if (draw_mean_graph)
                DrawMultiGraph(
                    group,
                    Form("%s %s mean", detectorTitle.Data(), (side==0?"junction":"ohmic")),
                    meanGraphs,
                    Form("%s_idx%d_side%d_mean", detType.Data(), idx, side),
                    bnnVoltage,
                    bnnEnergy,
                    "Corrected fit energy",
                    saveImages,
                    bestBias,
                    bestBiasBandHalfWidth
                );
            if (draw_resolution_graph)
                DrawMultiGraph(
                    group,
                    Form("%s %s resolution", detectorTitle.Data(), (side==0?"junction":"ohmic")),
                    resolutionGraphs,
                    Form("%s_idx%d_side%d_resolution", detType.Data(), idx, side),
                    bnnVoltage,
                    bnnResolution,
                    ";Bias [V];Energy resolution (FWHM / mean) [%]",
                    saveImages,
                    bestBias,
                    bestBiasBandHalfWidth
                );
        }
        auto leakageCurrentGraph = leakageCurrentGraphs[std::make_tuple(detType, idx)];
        if (leakageCurrentGraph != nullptr) {
            std::vector<DepVGraph> leakageGraphs = {
                {leakageCurrentGraph, "PL SAME"}
            };
            DrawMultiGraph(
                group,
                Form("%s leakage current", detectorTitle.Data()),
                leakageGraphs,
                Form("%s_idx%d_leakage_current", detType.Data(), idx),
                bnnVoltage,
                bnnLeakageCurrent,
                ";Bias [V];Leakage current",
                saveImages,
                bestBias,
                bestBiasBandHalfWidth
            );
        }
    }

    //top -> Draw("viewer");
    top -> Draw("");
    top -> Save("png");
}
