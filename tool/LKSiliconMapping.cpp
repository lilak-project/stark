#include "LKSiliconMapping.h"

#include "LKLogger.h"
#include "TSystem.h"

#include <fstream>
#include <iostream>
#include <sstream>

ClassImp(LKSiliconMapping)

namespace
{
std::vector<TString> TokenizeColumns(const TString &line)
{
    std::vector<TString> columns;
    std::istringstream stream(line.Data());
    std::string token;
    while (stream >> token)
        columns.push_back(TString(token));
    return columns;
}
}

void LKSiliconMapping::DetectorInfo::Print() const
{
    std::cout
        << "type=" << detType
        << " det_idx=" << detIndex
        << " det_number=" << detNumber
        << " ring_type=" << ringType
        << " phi_number=" << phiNumber
        << " phi=" << phi
        << std::endl;
}

void LKSiliconMapping::ChannelInfo::Print() const
{
    std::cout
        << "caac=(" << cobo << "," << asad << "," << aget << "," << chan << ")"
        << " det_idx=" << detIndex
        << " det_number=" << detNumber
        << " side=" << side
        << " strip=" << strip
        << std::endl;
}

LKSiliconMapping::LKSiliconMapping()
    : TNamed("LKSiliconMapping", "")
{
    Clear();
}

TString LKSiliconMapping::FindFirstMappingFile(TString mappingPath, TString suffix) const
{
    auto resolvedPath = mappingPath;
    gSystem->ExpandPathName(resolvedPath);
    void *dirHandle = gSystem->OpenDirectory(resolvedPath);
    if (dirHandle == nullptr)
        return "";

    TString found = "";
    while (auto entry = gSystem->GetDirEntry(dirHandle)) {
        TString name(entry);
        if (name == "." || name == "..")
            continue;
        if (name.EndsWith(suffix)) {
            found = TString::Format("%s/%s", resolvedPath.Data(), name.Data());
            break;
        }
    }
    gSystem->FreeDirectory(dirHandle);
    return found;
}

bool LKSiliconMapping::Load(TString mappingPath)
{
    auto resolvedPath = mappingPath;
    gSystem->ExpandPathName(resolvedPath);

    auto detectorFile = TString::Format("%s/detector_mapping.txt", resolvedPath.Data());
    auto channelFile = TString::Format("%s/channel_mapping.txt", resolvedPath.Data());

    if (gSystem->AccessPathName(detectorFile))
        detectorFile = FindFirstMappingFile(mappingPath, "_detector_mapping.txt");
    if (gSystem->AccessPathName(channelFile))
        channelFile = FindFirstMappingFile(mappingPath, "_channel_mapping.txt");

    if (detectorFile.IsNull() || channelFile.IsNull()) {
        lk_error << "Cannot find mapping files in " << mappingPath << std::endl;
        lk_error << "  detector mapping: " << detectorFile << std::endl;
        lk_error << "  channel mapping : " << channelFile << std::endl;
        return false;
    }
    return Load(detectorFile, channelFile);
}

void LKSiliconMapping::Clear(Option_t *)
{
    fDetectorMappingFileName = "";
    fChannelMappingFileName = "";
    fDetectorStyleIsNew = false;
    fChannelStyleIsNew = false;
    fMaxDetectorNumber = -1;
    fDetectors.clear();
    fChannels.clear();
    fDetectorVectorIndexByDetIndex.clear();
    fDetectorVectorIndexByDetNumber.clear();
    fChannelVectorIndexByCAAC.clear();
}

bool LKSiliconMapping::Load(TString detectorFileName, TString channelFileName)
{
    Clear();
    if (!LoadDetectorMapping(detectorFileName))
        return false;
    if (!LoadChannelMapping(channelFileName))
        return false;
    return true;
}

bool LKSiliconMapping::LoadDetectorMapping(TString fileName)
{
    fDetectorMappingFileName = fileName;

    std::ifstream input(fileName.Data());
    if (!input.is_open()) {
        lk_error << "Cannot open detector mapping file " << fileName << std::endl;
        return false;
    }

    TString firstLine = "";
    std::string rawLine;
    while (std::getline(input, rawLine)) {
        firstLine = TString(rawLine).Strip(TString::kBoth);
        if (!firstLine.IsNull())
            break;
    }

    input.clear();
    input.seekg(0);

    fDetectorStyleIsNew = firstLine.BeginsWith("detector", TString::kIgnoreCase);
    int lineIndex = 0;
    while (std::getline(input, rawLine)) {
        ++lineIndex;
        TString line(rawLine);
        line = line.Strip(TString::kBoth);
        if (line.IsNull() || line.BeginsWith("#"))
            continue;
        if (fDetectorStyleIsNew && lineIndex <= 2)
            continue;

        auto columns = TokenizeColumns(line);
        DetectorInfo info;
        bool ok = fDetectorStyleIsNew ? ParseNewDetectorRow(columns, info) : ParseLegacyDetectorRow(columns, info);
        if (!ok)
            continue;
        fDetectors.push_back(info);
    }

    RebuildLookupTables();
    return true;
}

bool LKSiliconMapping::LoadChannelMapping(TString fileName)
{
    fChannelMappingFileName = fileName;

    std::ifstream input(fileName.Data());
    if (!input.is_open()) {
        lk_error << "Cannot open channel mapping file " << fileName << std::endl;
        return false;
    }

    TString firstLine = "";
    std::string rawLine;
    while (std::getline(input, rawLine)) {
        firstLine = TString(rawLine).Strip(TString::kBoth);
        if (!firstLine.IsNull())
            break;
    }

    input.clear();
    input.seekg(0);

    fChannelStyleIsNew = firstLine.BeginsWith("channel", TString::kIgnoreCase);
    int lineIndex = 0;
    int dataIndex = 0;
    while (std::getline(input, rawLine)) {
        ++lineIndex;
        TString line(rawLine);
        line = line.Strip(TString::kBoth);
        if (line.IsNull() || line.BeginsWith("#"))
            continue;
        if (fChannelStyleIsNew && lineIndex <= 2)
            continue;

        auto columns = TokenizeColumns(line);
        ChannelInfo info;
        bool ok = fChannelStyleIsNew ? ParseNewChannelRow(columns, info) : ParseLegacyChannelRow(columns, info, dataIndex);
        if (!ok)
            continue;

        auto detector = FindDetectorByIndex(info.detIndex);
        if (detector != nullptr) {
            info.detNumber = detector->detNumber;
            if (info.detType.IsNull())
                info.detType = detector->detType;
        }

        fChannelVectorIndexByCAAC[MakeChannelKey(info.cobo, info.asad, info.aget, info.chan)] = fChannels.size();
        fChannels.push_back(info);
        ++dataIndex;
    }

    return true;
}

bool LKSiliconMapping::ParseNewDetectorRow(const std::vector<TString> &columns, DetectorInfo &info) const
{
    if (columns.size() < 15)
        return false;

    info.detType = columns[0];
    info.detIndex = columns[1].Atoi();
    info.detNumber = columns[2].Atoi();
    info.detThickness = columns[3].Atof();
    info.detWidth = columns[4].Atof();
    info.detHeight = columns[5].Atof();
    info.ringNumber = columns[6].Atoi();
    info.ringType = columns[7];
    info.ringRadius = columns[8].Atof();
    info.ringZ = columns[9].Atof();
    info.dEE = columns[10];
    info.phiNumber = columns[11].Atoi();
    info.phi = columns[12].Atof();
    info.phi1 = columns[13].Atof();
    info.phi2 = columns[14].Atof();
    info.isLegacy = false;
    return true;
}

bool LKSiliconMapping::ParseLegacyDetectorRow(const std::vector<TString> &columns, DetectorInfo &info) const
{
    if (columns.size() < 16)
        return false;

    info.detType = columns[0];
    info.detIndex = columns[1].Atoi();
    info.detNumber = info.detIndex;
    info.cobo = columns[2].Atoi();
    info.asad = columns[3].Atoi();
    info.zapJNo = columns[4];
    info.zapONo = columns[5];
    info.markNo = columns[6];
    info.markID = columns[7].Atoi();
    info.fb = columns[8].Atoi();
    info.ringNumber = columns[9].Atoi();
    info.dEE = columns[10] == "1" ? "E" : "dE";
    info.polarID = columns[11].Atoi();
    info.ringRadius = 5. * columns[12].Atof();
    info.ringZ = 10. * columns[13].Atof();
    info.phi = columns[14].Atof();
    info.phiNumber = columns[15].Atoi();
    info.ringType = Form("%d%s", info.ringNumber, info.dEE.Data());
    info.isLegacy = true;
    return true;
}

bool LKSiliconMapping::ParseNewChannelRow(const std::vector<TString> &columns, ChannelInfo &info) const
{
    if (columns.size() < 9)
        return false;

    info.channelIndex = columns[0].Atoi();
    info.cobo = columns[1].Atoi();
    info.asad = columns[2].Atoi();
    info.aget = columns[3].Atoi();
    info.chan = columns[4].Atoi();
    info.detIndex = columns[5].Atoi();
    info.phiNumber = columns[6].Atoi();
    info.side = columns[7].Atoi();
    info.strip = columns[8].Atoi();
    info.isLegacy = false;
    return true;
}

bool LKSiliconMapping::ParseLegacyChannelRow(const std::vector<TString> &columns, ChannelInfo &info, int rowIndex) const
{
    if (columns.size() < 12)
        return false;

    info.channelIndex = rowIndex;
    info.cobo = columns[0].Atoi();
    info.asad = columns[1].Atoi();
    info.aget = columns[2].Atoi();
    info.chan = columns[3].Atoi();
    info.chan2 = columns[4].Atoi();
    info.detType = columns[5];
    info.detIndex = columns[6].Atoi();
    info.detNumber = info.detIndex;
    info.side = columns[7].Atoi();
    info.strip = columns[8].Atoi();
    info.lr = columns[9].Atoi();
    info.detRadius = columns[10].Atof();
    info.detDistance = columns[11].Atof();
    info.phiNumber = -1;
    info.isLegacy = true;
    return true;
}

void LKSiliconMapping::RebuildLookupTables()
{
    int maxDetIndex = -1;
    fMaxDetectorNumber = -1;
    for (const auto &detector : fDetectors) {
        if (maxDetIndex < detector.detIndex)
            maxDetIndex = detector.detIndex;
        if (fMaxDetectorNumber < detector.detNumber)
            fMaxDetectorNumber = detector.detNumber;
    }

    fDetectorVectorIndexByDetIndex.assign(maxDetIndex + 1, -1);
    fDetectorVectorIndexByDetNumber.assign(fMaxDetectorNumber + 1, -1);
    for (size_t iDetector = 0; iDetector < fDetectors.size(); ++iDetector) {
        auto detIndex = fDetectors[iDetector].detIndex;
        auto detNumber = fDetectors[iDetector].detNumber;
        if (detIndex >= 0)
            fDetectorVectorIndexByDetIndex[detIndex] = iDetector;
        if (detNumber >= 0)
            fDetectorVectorIndexByDetNumber[detNumber] = iDetector;
    }
}

Long64_t LKSiliconMapping::MakeChannelKey(int cobo, int asad, int aget, int chan) const
{
    return (((Long64_t)cobo) << 24) | (((Long64_t)asad) << 16) | (((Long64_t)aget) << 8) | (Long64_t)chan;
}

const LKSiliconMapping::DetectorInfo *LKSiliconMapping::GetDetectorByVectorIndex(int index) const
{
    if (index < 0 || index >= (int) fDetectors.size())
        return nullptr;
    return &fDetectors[index];
}

const LKSiliconMapping::DetectorInfo *LKSiliconMapping::FindDetectorByIndex(int detIndex) const
{
    if (detIndex < 0 || detIndex >= (int) fDetectorVectorIndexByDetIndex.size())
        return nullptr;
    auto vectorIndex = fDetectorVectorIndexByDetIndex[detIndex];
    return GetDetectorByVectorIndex(vectorIndex);
}

const LKSiliconMapping::DetectorInfo *LKSiliconMapping::FindDetectorByNumber(int detNumber) const
{
    if (detNumber < 0 || detNumber >= (int) fDetectorVectorIndexByDetNumber.size())
        return nullptr;
    auto vectorIndex = fDetectorVectorIndexByDetNumber[detNumber];
    return GetDetectorByVectorIndex(vectorIndex);
}

const LKSiliconMapping::ChannelInfo *LKSiliconMapping::GetChannelByVectorIndex(int index) const
{
    if (index < 0 || index >= (int) fChannels.size())
        return nullptr;
    return &fChannels[index];
}

int LKSiliconMapping::FindChannelIndex(int cobo, int asad, int aget, int chan) const
{
    auto iterator = fChannelVectorIndexByCAAC.find(MakeChannelKey(cobo, asad, aget, chan));
    if (iterator == fChannelVectorIndexByCAAC.end())
        return -1;
    return iterator->second;
}

const LKSiliconMapping::ChannelInfo *LKSiliconMapping::FindChannel(int cobo, int asad, int aget, int chan) const
{
    return GetChannelByVectorIndex(FindChannelIndex(cobo, asad, aget, chan));
}

int LKSiliconMapping::FindDetectorIndex(int cobo, int asad, int aget, int chan) const
{
    auto channel = FindChannel(cobo, asad, aget, chan);
    if (channel == nullptr)
        return -1;
    return channel->detIndex;
}

int LKSiliconMapping::FindDetectorNumber(int cobo, int asad, int aget, int chan) const
{
    auto channel = FindChannel(cobo, asad, aget, chan);
    if (channel == nullptr)
        return -1;
    return channel->detNumber;
}

void LKSiliconMapping::Print(Option_t *) const
{
    lk_info << "LKSiliconMapping" << std::endl;
    lk_info << "  Detector file: " << fDetectorMappingFileName << std::endl;
    lk_info << "  Channel file : " << fChannelMappingFileName << std::endl;
    lk_info << "  Detector style: " << (fDetectorStyleIsNew ? "new" : "legacy") << std::endl;
    lk_info << "  Channel style : " << (fChannelStyleIsNew ? "new" : "legacy") << std::endl;
    lk_info << "  Number of detectors: " << fDetectors.size() << std::endl;
    lk_info << "  Number of channels : " << fChannels.size() << std::endl;
    lk_info << "  Max detector number: " << fMaxDetectorNumber << std::endl;
}

void LKSiliconMapping::PrintDetectors(int numPrint) const
{
    if (numPrint < 0)
        numPrint = fDetectors.size();
    if (numPrint > (int) fDetectors.size())
        numPrint = fDetectors.size();

    std::cout << "# First detectors" << std::endl;
    for (auto iDetector = 0; iDetector < numPrint; ++iDetector) {
        auto detector = GetDetectorByVectorIndex(iDetector);
        if (detector == nullptr)
            continue;
        std::cout << "det[" << iDetector << "] ";
        detector->Print();
    }
}

void LKSiliconMapping::PrintChannels(int numPrint) const
{
    if (numPrint < 0)
        numPrint = fChannels.size();
    if (numPrint > (int) fChannels.size())
        numPrint = fChannels.size();

    std::cout << "# First channels" << std::endl;
    for (auto iChannel = 0; iChannel < numPrint; ++iChannel) {
        auto channel = GetChannelByVectorIndex(iChannel);
        if (channel == nullptr)
            continue;
        std::cout << "ch[" << iChannel << "] ";
        channel->Print();
    }
}

void LKSiliconMapping::PrintChannelLookup(int cobo, int asad, int aget, int chan) const
{
    std::cout << "# Lookup by CAAC" << std::endl;
    std::cout << "query: (" << cobo << "," << asad << "," << aget << "," << chan << ")" << std::endl;

    auto channel = FindChannel(cobo, asad, aget, chan);
    if (channel == nullptr) {
        std::cout << "channel not found" << std::endl;
        return;
    }

    std::cout << "query result ";
    channel->Print();

    auto detector = FindDetectorByIndex(channel->detIndex);
    if (detector != nullptr) {
        std::cout << "matched detector ";
        detector->Print();
    }
}

void LKSiliconMapping::PrintTest(int queryCobo, int queryAsad, int queryAget, int queryChan, int numPrintDetectors, int numPrintChannels) const
{
    Print();
    std::cout << std::endl;
    PrintDetectors(numPrintDetectors);
    std::cout << std::endl;
    PrintChannels(numPrintChannels);
    std::cout << std::endl;
    PrintChannelLookup(queryCobo, queryAsad, queryAget, queryChan);
}
