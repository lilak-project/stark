#include "LKSiCalibrationHistogramBuilder.h"

#include "LKSiChannel.h"

#include "TClonesArray.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TTree.h"
#include "TSpectrum.h"

#include <algorithm>
#include <fstream>

ClassImp(LKSiCalibrationHistogramBuilder)

LKSiCalibrationHistogramBuilder::LKSiCalibrationHistogramBuilder()
    : TNamed("LKSiCalibrationHistogramBuilder", "")
{
    fGateEnergies = {3182.69, 5485.56};
}

void LKSiCalibrationHistogramBuilder::Clear(Option_t *)
{
    fRawEnergyHistMap.clear();
    fRawEnergySumHistMap.clear();
    fRawLeftRightHistMap.clear();
    fRawPositionEnergyHistMap.clear();
    fGateLeftRightHistMap.clear();
    fGatePositionEnergyHistMap.clear();
    fDefaultGateRanges.clear();
    fC0Parameters.clear();
    fC1Parameters.clear();
    fC2Parameters.clear();
}

bool LKSiCalibrationHistogramBuilder::Configure(TString detectorMapFileName, TString channelMapFileName)
{
    return fMapping.Load(detectorMapFileName, channelMapFileName);
}

void LKSiCalibrationHistogramBuilder::SetEnergyBinning(int nbin, double x1, double x2)
{
    fNBinEnergy = nbin;
    fBinEnergy1 = x1;
    fBinEnergy2 = x2;
}

void LKSiCalibrationHistogramBuilder::SetPositionBinning(int nbin, double x1, double x2)
{
    fNBinPosition = nbin;
    fBinPosition1 = x1;
    fBinPosition2 = x2;
}

void LKSiCalibrationHistogramBuilder::SetGateConfig(const std::vector<double> &gateEnergies, double expectedResolution, double entriesCut)
{
    fGateEnergies = gateEnergies;
    fExpectedResolution = expectedResolution;
    fEntriesCut = entriesCut;
}

void LKSiCalibrationHistogramBuilder::SetCalibrationStage(int stage)
{
    fApplyC0 = (stage >= 1);
    fApplyC1 = (stage >= 1);
    fApplyC2 = (stage >= 2);
    if (stage <= 0)
        fHistogramPrefix = "raw";
    else if (stage == 1)
        fHistogramPrefix = "c1";
    else
        fHistogramPrefix = "c2";
}

bool LKSiCalibrationHistogramBuilder::LoadC0Parameters(TString fileName)
{
    fC0Parameters.clear();
    std::ifstream input(fileName.Data());
    if (!input.is_open())
        return false;

    std::string line;
    std::getline(input, line);
    int det = -1;
    int side = -1;
    int strip = -1;
    double entries = 0;
    double itcpt = 0;
    double slope = 0;
    while (input >> det >> side >> strip >> entries >> itcpt >> slope)
        fC0Parameters[C0Key(det, side, strip)] = {itcpt, slope};
    return true;
}

bool LKSiCalibrationHistogramBuilder::LoadC1Parameters(TString fileName)
{
    fC1Parameters.clear();
    std::ifstream input(fileName.Data());
    if (!input.is_open())
        return false;

    std::string line;
    std::getline(input, line);
    int det = -1;
    int side = -1;
    int strip = -1;
    int nPoints = 0;
    double itcptL = 0;
    double slopeL = 0;
    double itcptR = 0;
    double slopeR = 0;
    while (input >> det >> side >> strip >> nPoints >> itcptL >> slopeL >> itcptR >> slopeR) {
        fC1Parameters[C1Key(det, side, strip, 0)] = {itcptL, slopeL};
        fC1Parameters[C1Key(det, side, strip, 1)] = {itcptR, slopeR};
    }
    return true;
}

bool LKSiCalibrationHistogramBuilder::LoadC2Parameters(TString fileName)
{
    fC2Parameters.clear();
    std::ifstream input(fileName.Data());
    if (!input.is_open())
        return false;

    std::string line;
    std::getline(input, line);
    int det = -1;
    int side = -1;
    int strip = -1;
    double entries = 0;
    double b0 = 0;
    double b1 = 0;
    double b2 = 0;
    while (input >> det >> side >> strip >> entries >> b0 >> b1 >> b2)
        fC2Parameters[C2Key(det, side, strip)] = {b0, b1, b2};
    return true;
}

TString LKSiCalibrationHistogramBuilder::MakeBaseName(const char *prefix, int detNumber, int side, int strip, int gate) const
{
    auto fullPrefix = TString(prefix);
    if (!fHistogramPrefix.IsNull())
        fullPrefix = Form("%s_%s", fHistogramPrefix.Data(), prefix);
    if (gate >= 0)
        return Form("%s_det%03d_side%d_strip%02d_gate%d", fullPrefix.Data(), detNumber, side, strip, gate);
    return Form("%s_det%03d_side%d_strip%02d", fullPrefix.Data(), detNumber, side, strip);
}

TH1D *LKSiCalibrationHistogramBuilder::MakeHist1(const char *prefix, int detNumber, int side, int strip, int gate)
{
    auto name = MakeBaseName(prefix, detNumber, side, strip, gate);
    auto hist = new TH1D(name, name, fNBinEnergy, fBinEnergy1, fBinEnergy2);
    hist->SetDirectory(nullptr);
    return hist;
}

TH2D *LKSiCalibrationHistogramBuilder::MakeHist2(const char *prefixX, const char *prefixY, int detNumber, int side, int strip, int gate)
{
    auto name = MakeBaseName(Form("%s_%s", prefixX, prefixY), detNumber, side, strip, gate);
    auto hist = new TH2D(name, name, fNBinPosition, fBinPosition1, fBinPosition2, fNBinEnergy, fBinEnergy1, fBinEnergy2);
    if (TString(prefixX).Contains("left") || TString(prefixY).Contains("right")) {
        delete hist;
        hist = new TH2D(name, name, fNBinEnergy, fBinEnergy1, fBinEnergy2, fNBinEnergy, fBinEnergy1, fBinEnergy2);
    }
    hist->SetDirectory(nullptr);
    return hist;
}

bool LKSiCalibrationHistogramBuilder::IsResistiveChannel(const LKSiliconMapping::ChannelInfo *channelInfo, const LKSiliconMapping::DetectorInfo *detectorInfo) const
{
    if (channelInfo == nullptr || detectorInfo == nullptr)
        return false;
    return detectorInfo->detType.EqualTo("X6", TString::kIgnoreCase) && channelInfo->side == 1;
}

void LKSiCalibrationHistogramBuilder::AllocateHistograms()
{
    fRawEnergyHistMap.clear();
    fRawEnergySumHistMap.clear();
    fRawLeftRightHistMap.clear();
    fRawPositionEnergyHistMap.clear();
    fGateLeftRightHistMap.clear();
    fGatePositionEnergyHistMap.clear();
    fDefaultGateRanges.clear();
    auto numChannels = fMapping.GetNumChannels();
    for (auto iChannel = 0; iChannel < numChannels; ++iChannel) {
        auto channel = fMapping.GetChannelByVectorIndex(iChannel);
        if (channel == nullptr)
            continue;
        auto detector = fMapping.FindDetectorByIndex(channel->detIndex);
        if (detector == nullptr)
            continue;

        StripKey key = {channel->detNumber, channel->side, channel->strip};
        if (IsResistiveChannel(channel, detector)) {
            if (fRawEnergySumHistMap.count(key) == 0) {
                fRawEnergySumHistMap[key] = MakeHist1("raw_esum", channel->detNumber, channel->side, channel->strip);
                fRawLeftRightHistMap[key] = MakeHist2("raw_left", "raw_right", channel->detNumber, channel->side, channel->strip);
                fRawPositionEnergyHistMap[key] = MakeHist2("raw_rpos", "raw_esum", channel->detNumber, channel->side, channel->strip);
                for (int gate = 0; gate < (int) fGateEnergies.size(); ++gate) {
                    auto gateKey = std::make_tuple(channel->detNumber, channel->side, channel->strip, gate);
                    fGateLeftRightHistMap[gateKey] = MakeHist2("gate_left", "gate_right", channel->detNumber, channel->side, channel->strip, gate);
                    fGatePositionEnergyHistMap[gateKey] = MakeHist2("gate_rpos", "gate_esum", channel->detNumber, channel->side, channel->strip, gate);
                }
            }
        }
        else {
            if (fRawEnergyHistMap.count(key) == 0)
                fRawEnergyHistMap[key] = MakeHist1("raw_energy", channel->detNumber, channel->side, channel->strip);
        }
    }
}

void LKSiCalibrationHistogramBuilder::FillMappedChannel(int cobo, int asad, int aget, int chan, double energy1, double energy2)
{
    auto channelInfo = fMapping.FindChannel(cobo, asad, aget, chan);
    if (channelInfo == nullptr)
        return;
    auto detectorInfo = fMapping.FindDetectorByIndex(channelInfo->detIndex);
    if (detectorInfo == nullptr)
        return;

    StripKey key = {channelInfo->detNumber, channelInfo->side, channelInfo->strip};
    if (IsResistiveChannel(channelInfo, detectorInfo) && energy2 >= 0) {
        if (fApplyC1) {
            CalibrateC1(channelInfo->detNumber, channelInfo->side, channelInfo->strip, 0, energy2);
            CalibrateC1(channelInfo->detNumber, channelInfo->side, channelInfo->strip, 1, energy1);
        }
        auto sum = energy1 + energy2;
        if (sum <= 0)
            return;
        auto pos = (energy1 - energy2) / sum;
        if (fApplyC2)
            CalibrateC2(channelInfo->detNumber, channelInfo->side, channelInfo->strip, sum, pos);
        if (fRawEnergySumHistMap.count(key)) fRawEnergySumHistMap[key]->Fill(sum);
        if (fRawLeftRightHistMap.count(key)) fRawLeftRightHistMap[key]->Fill(energy1, energy2);
        if (fRawPositionEnergyHistMap.count(key)) fRawPositionEnergyHistMap[key]->Fill(pos, sum);
    }
    else {
        if (fApplyC0)
            CalibrateC0(channelInfo->detNumber, channelInfo->side, channelInfo->strip, energy1);
        if (fRawEnergyHistMap.count(key))
            fRawEnergyHistMap[key]->Fill(energy1);
    }
}

void LKSiCalibrationHistogramBuilder::FillFromSiChannel(const LKSiChannel *channel)
{
    if (channel == nullptr)
        return;
    FillMappedChannel(channel->GetCobo(), channel->GetAsad(), channel->GetAget(), channel->GetChan(), channel->GetEnergy(), channel->GetEnergy2());
}

void LKSiCalibrationHistogramBuilder::ComputeDefaultGateRanges()
{
    fDefaultGateRanges.assign(fGateEnergies.size(), std::make_pair(0.0, 0.0));
    TSpectrum spectrum((int) fGateEnergies.size() + 2);
    for (auto &entry : fRawEnergySumHistMap) {
        auto hist = entry.second;
        if (hist == nullptr || hist->GetEntries() < fEntriesCut)
            continue;
        auto numPeaks = spectrum.Search(hist, 5, "goff nodraw");
        if (numPeaks < (int) fGateEnergies.size())
            continue;
        auto peaks = spectrum.GetPositionX();
        std::vector<double> sortedPeaks(peaks, peaks + numPeaks);
        std::sort(sortedPeaks.begin(), sortedPeaks.end());
        for (int iGate = 0; iGate < (int) fGateEnergies.size() && iGate < (int) sortedPeaks.size(); ++iGate) {
            auto center = sortedPeaks[iGate];
            auto width = center * fExpectedResolution * 8.0;
            if (fDefaultGateRanges[iGate].second <= fDefaultGateRanges[iGate].first) {
                fDefaultGateRanges[iGate] = std::make_pair(center - width, center + width);
            }
        }
    }
}

void LKSiCalibrationHistogramBuilder::ResetGateHistograms()
{
    for (auto &entry : fGateLeftRightHistMap)
        entry.second->Reset();
    for (auto &entry : fGatePositionEnergyHistMap)
        entry.second->Reset();
}

void LKSiCalibrationHistogramBuilder::FillGateHistogramsFromSiChannel(const LKSiChannel *channel)
{
    if (channel == nullptr)
        return;
    auto channelInfo = fMapping.FindChannel(channel->GetCobo(), channel->GetAsad(), channel->GetAget(), channel->GetChan());
    if (channelInfo == nullptr)
        return;
    auto detectorInfo = fMapping.FindDetectorByIndex(channelInfo->detIndex);
    if (!IsResistiveChannel(channelInfo, detectorInfo))
        return;

    auto energyR = channel->GetEnergy();
    auto energyL = channel->GetEnergy2();
    auto sum = energyR + energyL;
    if (sum <= 0)
        return;
    auto pos = (energyR - energyL) / sum;

    for (int gate = 0; gate < (int) fDefaultGateRanges.size(); ++gate) {
        auto range = fDefaultGateRanges[gate];
        if (sum < range.first || sum > range.second)
            continue;
        auto gateKey = std::make_tuple(channelInfo->detNumber, channelInfo->side, channelInfo->strip, gate);
        if (fGateLeftRightHistMap.count(gateKey)) fGateLeftRightHistMap[gateKey]->Fill(energyR, energyL);
        if (fGatePositionEnergyHistMap.count(gateKey)) fGatePositionEnergyHistMap[gateKey]->Fill(pos, sum);
    }
}

bool LKSiCalibrationHistogramBuilder::BuildFromTree(TTree *tree, TString branchName)
{
    if (tree == nullptr)
        return false;

    AllocateHistograms();

    TClonesArray *array = nullptr;
    tree->SetBranchAddress(branchName.Data(), &array);
    auto numEntries = tree->GetEntries();
    for (Long64_t iEntry = 0; iEntry < numEntries; ++iEntry) {
        tree->GetEntry(iEntry);
        if (array == nullptr)
            continue;
        auto numChannels = array->GetEntriesFast();
        for (int iChannel = 0; iChannel < numChannels; ++iChannel)
            FillFromSiChannel((LKSiChannel *) array->At(iChannel));
    }

    ComputeDefaultGateRanges();
    ResetGateHistograms();

    for (Long64_t iEntry = 0; iEntry < numEntries; ++iEntry) {
        tree->GetEntry(iEntry);
        if (array == nullptr)
            continue;
        auto numChannels = array->GetEntriesFast();
        for (int iChannel = 0; iChannel < numChannels; ++iChannel)
            FillGateHistogramsFromSiChannel((LKSiChannel *) array->At(iChannel));
    }

    return true;
}

bool LKSiCalibrationHistogramBuilder::BuildFromTreeFile(TString fileName, TString treeName, TString branchName)
{
    TFile file(fileName, "read");
    if (file.IsZombie())
        return false;
    auto tree = (TTree *) file.Get(treeName);
    if (tree == nullptr)
        return false;
    return BuildFromTree(tree, branchName);
}

bool LKSiCalibrationHistogramBuilder::Write(TString fileName)
{
    TFile file(fileName, "recreate");
    if (file.IsZombie())
        return false;
    for (auto &entry : fRawEnergyHistMap) entry.second->Write();
    for (auto &entry : fRawEnergySumHistMap) entry.second->Write();
    for (auto &entry : fRawLeftRightHistMap) entry.second->Write();
    for (auto &entry : fRawPositionEnergyHistMap) entry.second->Write();
    for (auto &entry : fGateLeftRightHistMap) entry.second->Write();
    for (auto &entry : fGatePositionEnergyHistMap) entry.second->Write();
    file.Close();
    return true;
}

TH1D *LKSiCalibrationHistogramBuilder::GetRawEnergyHist(int detNumber, int side, int strip) const
{
    auto key = std::make_tuple(detNumber, side, strip);
    auto it = fRawEnergyHistMap.find(key);
    return it == fRawEnergyHistMap.end() ? nullptr : it->second;
}

TH1D *LKSiCalibrationHistogramBuilder::GetRawEnergySumHist(int detNumber, int side, int strip) const
{
    auto key = std::make_tuple(detNumber, side, strip);
    auto it = fRawEnergySumHistMap.find(key);
    return it == fRawEnergySumHistMap.end() ? nullptr : it->second;
}

TH2D *LKSiCalibrationHistogramBuilder::GetRawLeftRightHist(int detNumber, int side, int strip) const
{
    auto key = std::make_tuple(detNumber, side, strip);
    auto it = fRawLeftRightHistMap.find(key);
    return it == fRawLeftRightHistMap.end() ? nullptr : it->second;
}

TH2D *LKSiCalibrationHistogramBuilder::GetRawPositionEnergyHist(int detNumber, int side, int strip) const
{
    auto key = std::make_tuple(detNumber, side, strip);
    auto it = fRawPositionEnergyHistMap.find(key);
    return it == fRawPositionEnergyHistMap.end() ? nullptr : it->second;
}

TH2D *LKSiCalibrationHistogramBuilder::GetGateLeftRightHist(int detNumber, int side, int strip, int gate) const
{
    auto key = std::make_tuple(detNumber, side, strip, gate);
    auto it = fGateLeftRightHistMap.find(key);
    return it == fGateLeftRightHistMap.end() ? nullptr : it->second;
}

TH2D *LKSiCalibrationHistogramBuilder::GetGatePositionEnergyHist(int detNumber, int side, int strip, int gate) const
{
    auto key = std::make_tuple(detNumber, side, strip, gate);
    auto it = fGatePositionEnergyHistMap.find(key);
    return it == fGatePositionEnergyHistMap.end() ? nullptr : it->second;
}

bool LKSiCalibrationHistogramBuilder::CalibrateC0(int detNumber, int side, int strip, double &energy) const
{
    auto iterator = fC0Parameters.find(C0Key(detNumber, side, strip));
    if (iterator == fC0Parameters.end())
        return false;
    energy = iterator->second[1] * energy + iterator->second[0];
    return true;
}

bool LKSiCalibrationHistogramBuilder::CalibrateC1(int detNumber, int side, int strip, int lr, double &energy) const
{
    auto iterator = fC1Parameters.find(C1Key(detNumber, side, strip, lr));
    if (iterator == fC1Parameters.end())
        return false;
    energy = iterator->second[1] * energy + iterator->second[0];
    return true;
}

bool LKSiCalibrationHistogramBuilder::CalibrateC2(int detNumber, int side, int strip, double &energySum, double position) const
{
    auto iterator = fC2Parameters.find(C2Key(detNumber, side, strip));
    if (iterator == fC2Parameters.end())
        return false;
    auto scale = iterator->second[0] + iterator->second[1] * position + iterator->second[2] * position * position;
    if (scale == 0)
        return false;
    energySum = energySum / scale * 5.48556;
    return true;
}
