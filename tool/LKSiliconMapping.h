#ifndef LKSILICONMAPPING_HH
#define LKSILICONMAPPING_HH

#include "TNamed.h"
#include "TString.h"

#include <unordered_map>
#include <vector>

class LKSiliconMapping : public TNamed
{
    public:
        struct DetectorInfo
        {
            TString detType = "";
            int detIndex = -1;
            int detNumber = -1;
            double detThickness = 0;
            double detWidth = 0;
            double detHeight = 0;
            int ringNumber = -1;
            TString ringType = "";
            double ringRadius = 0;
            double ringZ = 0;
            TString dEE = "";
            int phiNumber = -1;
            double phi = 0;
            double phi1 = 0;
            double phi2 = 0;

            int cobo = -1;
            int asad = -1;
            TString zapJNo = "";
            TString zapONo = "";
            TString markNo = "";
            int markID = -1;
            int fb = -1;
            int polarID = -1;
            bool isLegacy = false;

            void Print() const;
        };

        struct ChannelInfo
        {
            int channelIndex = -1;
            int cobo = -1;
            int asad = -1;
            int aget = -1;
            int chan = -1;
            int detIndex = -1;
            int detNumber = -1;
            int phiNumber = -1;
            int side = -1;
            int strip = -1;

            int chan2 = -1;
            int lr = -1;
            TString detType = "";
            double detRadius = 0;
            double detDistance = 0;
            bool isLegacy = false;

            void Print() const;
        };

    public:
        LKSiliconMapping();
        virtual ~LKSiliconMapping() {}

        bool Load(TString mappingPath);
        bool Load(TString detectorFileName, TString channelFileName);
        bool LoadDetectorMapping(TString fileName);
        bool LoadChannelMapping(TString fileName);
        void Clear(Option_t *option="");
        void Print(Option_t *option="") const;
        void PrintDetectors(int numPrint = 20) const;
        void PrintChannels(int numPrint = 20) const;
        void PrintChannelLookup(int cobo, int asad, int aget, int chan) const;
        void PrintTest(int queryCobo = 0, int queryAsad = 0, int queryAget = 3, int queryChan = 0, int numPrintDetectors = 20, int numPrintChannels = 20) const;

        int GetNumDetectors() const { return fDetectors.size(); }
        int GetNumChannels() const { return fChannels.size(); }
        int GetMaxDetectorNumber() const { return fMaxDetectorNumber; }
        bool IsNewDetectorStyle() const { return fDetectorStyleIsNew; }
        bool IsNewChannelStyle() const { return fChannelStyleIsNew; }

        const DetectorInfo *GetDetectorByVectorIndex(int index) const;
        const DetectorInfo *FindDetectorByIndex(int detIndex) const;
        const DetectorInfo *FindDetectorByNumber(int detNumber) const;
        const ChannelInfo *GetChannelByVectorIndex(int index) const;
        const ChannelInfo *FindChannel(int cobo, int asad, int aget, int chan) const;

        int FindChannelIndex(int cobo, int asad, int aget, int chan) const;
        int FindDetectorIndex(int cobo, int asad, int aget, int chan) const;
        int FindDetectorNumber(int cobo, int asad, int aget, int chan) const;

    private:
        TString FindFirstMappingFile(TString mappingPath, TString suffix) const;
        bool ParseNewDetectorRow(const std::vector<TString> &columns, DetectorInfo &info) const;
        bool ParseLegacyDetectorRow(const std::vector<TString> &columns, DetectorInfo &info) const;
        bool ParseNewChannelRow(const std::vector<TString> &columns, ChannelInfo &info) const;
        bool ParseLegacyChannelRow(const std::vector<TString> &columns, ChannelInfo &info, int rowIndex) const;
        void RebuildLookupTables();
        Long64_t MakeChannelKey(int cobo, int asad, int aget, int chan) const;

    private:
        TString fDetectorMappingFileName = "";
        TString fChannelMappingFileName = "";
        bool fDetectorStyleIsNew = false;
        bool fChannelStyleIsNew = false;
        int fMaxDetectorNumber = -1;

        std::vector<DetectorInfo> fDetectors;
        std::vector<ChannelInfo> fChannels;
        std::vector<int> fDetectorVectorIndexByDetIndex;
        std::vector<int> fDetectorVectorIndexByDetNumber;
        std::unordered_map<Long64_t, int> fChannelVectorIndexByCAAC;

    ClassDef(LKSiliconMapping, 1);
};

#endif
