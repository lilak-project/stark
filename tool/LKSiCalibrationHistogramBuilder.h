#ifndef LKSICALIBRATIONHISTOGRAMBUILDER_HH
#define LKSICALIBRATIONHISTOGRAMBUILDER_HH

#include "TNamed.h"
#include "TString.h"

#include "LKSiliconMapping.h"

#include <array>
#include <map>
#include <tuple>
#include <vector>

class TH1D;
class TH2D;
class TTree;
class TClonesArray;
class TDirectory;
class LKSiChannel;

class LKSiCalibrationHistogramBuilder : public TNamed
{
    public:
        LKSiCalibrationHistogramBuilder();
        virtual ~LKSiCalibrationHistogramBuilder() {}

        bool Configure(TString detectorMapFileName, TString channelMapFileName);
        void Clear(Option_t *option="");
        void SetEnergyBinning(int nbin, double x1, double x2);
        void SetPositionBinning(int nbin, double x1, double x2);
        void SetGateConfig(const std::vector<double> &gateEnergies, double expectedResolution=0.01, double entriesCut=100.0);
        void SetHistogramPrefix(TString prefix) { fHistogramPrefix = prefix; }
        void SetApplyC0(bool value=true) { fApplyC0 = value; }
        void SetApplyC1(bool value=true) { fApplyC1 = value; }
        void SetApplyC2(bool value=true) { fApplyC2 = value; }
        void SetCalibrationStage(int stage);
        bool LoadC0Parameters(TString fileName);
        bool LoadC1Parameters(TString fileName);
        bool LoadC2Parameters(TString fileName);

        bool BuildFromTree(TTree *tree, TString branchName="SiChannel");
        bool BuildFromTreeFile(TString fileName, TString treeName="event", TString branchName="SiChannel");
        bool Write(TString fileName);

        TH1D *GetRawEnergyHist(int detNumber, int side, int strip) const;
        TH1D *GetRawEnergySumHist(int detNumber, int side, int strip) const;
        TH2D *GetRawLeftRightHist(int detNumber, int side, int strip) const;
        TH2D *GetRawPositionEnergyHist(int detNumber, int side, int strip) const;
        TH2D *GetGateLeftRightHist(int detNumber, int side, int strip, int gate) const;
        TH2D *GetGatePositionEnergyHist(int detNumber, int side, int strip, int gate) const;

        const std::vector<double> &GetGateEnergies() const { return fGateEnergies; }
        const std::vector<std::pair<double,double>> &GetDefaultGateRanges() const { return fDefaultGateRanges; }
        LKSiliconMapping *GetMapping() { return &fMapping; }

    private:
        using StripKey = std::tuple<int,int,int>;
        using C0Key = std::tuple<int,int,int>;
        using C1Key = std::tuple<int,int,int,int>;
        using C2Key = std::tuple<int,int,int>;

        TString MakeBaseName(const char *prefix, int detNumber, int side, int strip, int gate=-1) const;
        TH1D *MakeHist1(const char *prefix, int detNumber, int side, int strip, int gate=-1);
        TH2D *MakeHist2(const char *prefixX, const char *prefixY, int detNumber, int side, int strip, int gate=-1);
        bool IsResistiveChannel(const LKSiliconMapping::ChannelInfo *channelInfo, const LKSiliconMapping::DetectorInfo *detectorInfo) const;
        void FillMappedChannel(int cobo, int asad, int aget, int chan, double energy1, double energy2);
        void FillFromSiChannel(const LKSiChannel *channel);
        void ComputeDefaultGateRanges();
        void AllocateHistograms();
        void ResetGateHistograms();
        void FillGateHistogramsFromSiChannel(const LKSiChannel *channel);
        bool CalibrateC0(int detNumber, int side, int strip, double &energy) const;
        bool CalibrateC1(int detNumber, int side, int strip, int lr, double &energy) const;
        bool CalibrateC2(int detNumber, int side, int strip, double &energySum, double position) const;

    private:
        LKSiliconMapping fMapping;

        int fNBinEnergy = 400;
        double fBinEnergy1 = 0.0;
        double fBinEnergy2 = 20000.0;
        int fNBinPosition = 200;
        double fBinPosition1 = -1.0;
        double fBinPosition2 = 1.0;
        double fExpectedResolution = 0.01;
        double fEntriesCut = 100.0;
        TString fHistogramPrefix = "raw";
        bool fApplyC0 = false;
        bool fApplyC1 = false;
        bool fApplyC2 = false;

        std::vector<double> fGateEnergies;
        std::vector<std::pair<double,double>> fDefaultGateRanges;
        std::map<C0Key, std::array<double,2>> fC0Parameters;
        std::map<C1Key, std::array<double,2>> fC1Parameters;
        std::map<C2Key, std::array<double,3>> fC2Parameters;

        std::map<StripKey, TH1D*> fRawEnergyHistMap;
        std::map<StripKey, TH1D*> fRawEnergySumHistMap;
        std::map<StripKey, TH2D*> fRawLeftRightHistMap;
        std::map<StripKey, TH2D*> fRawPositionEnergyHistMap;
        std::map<std::tuple<int,int,int,int>, TH2D*> fGateLeftRightHistMap;
        std::map<std::tuple<int,int,int,int>, TH2D*> fGatePositionEnergyHistMap;

    ClassDef(LKSiCalibrationHistogramBuilder, 1);
};

#endif
