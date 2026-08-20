#include <TString.h>
#include <TSystem.h>
#include <TSystemDirectory.h>
#include <TSystemFile.h>
#include <TList.h>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>

namespace {

using ResKey = std::tuple<std::string, int, double, int>; // detector set, detector number, voltage, side

bool ParseDepVFileName(const TString &fileName, TString &detectorSet, double &voltage)
{
    if (!fileName.EndsWith("_er.dat"))
        return false;

    TString stem(fileName);
    stem.ReplaceAll("_er.dat", "");

    auto lastUnderscore = stem.Last('_');
    if (lastUnderscore == kNPOS)
        return false;

    detectorSet = stem(0, lastUnderscore);
    TString voltageText = stem(lastUnderscore + 1, stem.Length() - lastUnderscore - 1);
    if (!voltageText.IsFloat())
        return false;

    voltage = voltageText.Atof();
    return true;
}

double Mean(const std::vector<double> &values)
{
    if (values.empty())
        return 0;
    double sum = 0;
    for (auto value : values)
        sum += value;
    return sum / values.size();
}

} // namespace

void make_x6_depv_latex_table(TString dataDir = "data_depv",
                              TString outputFile = "depv_x6_best_bias_table.tex")
{
    std::set<std::string> selectedSets = {"X6_1", "X6_2"};
    std::map<ResKey, std::vector<double>> resolutions;

    TSystemDirectory directory(dataDir, dataDir);
    auto fileList = directory.GetListOfFiles();
    if (fileList == nullptr) {
        std::cerr << "Cannot open data directory: " << dataDir << std::endl;
        return;
    }

    TIter next(fileList);
    while (auto object = next()) {
        auto file = dynamic_cast<TSystemFile *>(object);
        if (file == nullptr || file->IsDirectory())
            continue;

        TString fileName = file->GetName();
        TString detectorSet;
        double voltage = 0;
        if (!ParseDepVFileName(fileName, detectorSet, voltage))
            continue;
        if (selectedSets.find(detectorSet.Data()) == selectedSets.end())
            continue;

        TString path = dataDir + "/" + fileName;
        std::ifstream input(path.Data());
        if (!input.is_open()) {
            std::cerr << "Cannot open input file: " << path << std::endl;
            continue;
        }

        std::string line;
        while (std::getline(input, line)) {
            if (line.empty() || line[0] == '#')
                continue;

            int detNumber = -1;
            int side = -1;
            int channel = -1;
            double amp = 0;
            double energy = 0;
            double fwhm = 0;
            double rawMean = 0;

            std::istringstream lineStream(line);
            lineStream >> detNumber >> side >> channel >> amp >> energy >> fwhm;
            if (!lineStream)
                continue;
            lineStream >> rawMean; // optional column; not needed for this table

            if (energy <= 0 || fwhm <= 0 || !std::isfinite(energy) || !std::isfinite(fwhm))
                continue;

            double resolution = 100. * fwhm / energy;
            if (resolution <= 0 || resolution >= 200.)
                continue;

            resolutions[ResKey(detectorSet.Data(), detNumber, voltage, side)].push_back(resolution);
        }
    }

    struct Row {
        std::string detectorSet;
        int detNumber = -1;
        double bestBias = 0;
        double junctionResolution = 0;
        int junctionChannels = 0;
        double ohmicResolution = 0;
        int ohmicChannels = 0;
    };

    std::vector<Row> rows;
    std::set<std::tuple<std::string, int>> detectors;
    for (auto const &entry : resolutions)
        detectors.insert(std::make_tuple(std::get<0>(entry.first), std::get<1>(entry.first)));

    for (auto const &detector : detectors) {
        auto detectorSet = std::get<0>(detector);
        auto detNumber = std::get<1>(detector);

        double bestBias = 0;
        double bestResolution = 0;
        int bestChannels = 0;
        bool hasBest = false;

        for (auto const &entry : resolutions) {
            auto const &key = entry.first;
            if (std::get<0>(key) != detectorSet || std::get<1>(key) != detNumber || std::get<3>(key) != 0)
                continue;

            double meanResolution = Mean(entry.second);
            if (!hasBest || meanResolution < bestResolution) {
                hasBest = true;
                bestBias = std::get<2>(key);
                bestResolution = meanResolution;
                bestChannels = entry.second.size();
            }
        }

        if (!hasBest)
            continue;

        auto ohmicKey = ResKey(detectorSet, detNumber, bestBias, 1);
        auto ohmicIt = resolutions.find(ohmicKey);

        Row row;
        row.detectorSet = detectorSet;
        row.detNumber = detNumber;
        row.bestBias = bestBias;
        row.junctionResolution = bestResolution;
        row.junctionChannels = bestChannels;
        if (ohmicIt != resolutions.end()) {
            row.ohmicResolution = Mean(ohmicIt->second);
            row.ohmicChannels = ohmicIt->second.size();
        }
        rows.push_back(row);
    }

    std::ofstream output(outputFile.Data());
    if (!output.is_open()) {
        std::cerr << "Cannot open output file: " << outputFile << std::endl;
        return;
    }

    output << "% X6 depletion-voltage summary from " << dataDir << "/X6_1_*_er.dat and "
           << dataDir << "/X6_2_*_er.dat\n";
    output << "% Best operation bias is selected from the minimum average junction resolution:\n";
    output << "% resolution = corrected_fit_fwhm / corrected_fit_energy * 100.\n";
    output << "\\begin{table}[htbp]\n";
    output << "\\centering\n";
    output << "\\caption{Best operation bias summary for the re-analysed X6 detectors.}\n";
    output << "\\label{tab:x6-best-bias}\n";
    output << "\\begin{tabular}{c c c c c c c}\n";
    output << "\\hline\n";
    output << "Set & Det.\\# & Best bias [V] & Junction res. [\\%] & Junction ch. & Ohmic res. [\\%] & Ohmic ch. \\\\\n";
    output << "\\hline\n";

    std::string previousSet;
    output << std::fixed << std::setprecision(2);
    for (auto const &row : rows) {
        if (!previousSet.empty() && previousSet != row.detectorSet)
            output << "\\hline\n";
        previousSet = row.detectorSet;

        output << row.detectorSet << " & "
               << row.detNumber << " & "
               << std::setprecision(0) << row.bestBias << " & "
               << std::setprecision(2) << row.junctionResolution << " & "
               << row.junctionChannels << " & ";
        if (row.ohmicChannels > 0)
            output << row.ohmicResolution << " & " << row.ohmicChannels;
        else
            output << "-- & 0";
        output << " \\\\\n";
    }

    output << "\\hline\n";
    output << "\\end{tabular}\n";
    output << "\\end{table}\n";

    std::cout << "Wrote " << outputFile << " with " << rows.size() << " rows." << std::endl;
}
