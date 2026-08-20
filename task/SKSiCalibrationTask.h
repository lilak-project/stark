#ifndef SKSICALIBRATIONTASK_HH
#define SKSICALIBRATIONTASK_HH

#include "LKTask.h"

#include "TClonesArray.h"
#include "TString.h"

#include <array>
#include <map>
#include <utility>
#include <tuple>
#include <vector>

class TH1D;
class TH2D;
class TPaveText;
class LKSiliconArray;

class SKSiCalibrationTask : public LKTask
{
    public:
        SKSiCalibrationTask();
        virtual ~SKSiCalibrationTask();

        bool Init();
        void Exec(Option_t *);
        bool EndOfRun();

    private:
        using StripKey = std::tuple<int,int,int>;
        using LRKey = std::tuple<int,int,int,int>;

        TH1D *GetEnergyHist(int det, int side, int strip);
        TH1D *GetEnergySumHist(int det, int side, int strip);
        TH2D *GetLeftRightHist(int det, int side, int strip);
        TH2D *GetPositionEnergyHist(int det, int side, int strip);

        bool LoadC0Parameters(TString fileName);
        bool LoadC1Parameters(TString fileName);
        bool LoadC2Parameters(TString fileName);
        bool LoadC3Parameters(TString fileName);
        bool LoadParametersForStage();

        double ApplyLinear(const std::array<double,2> &par, double value) const;
        bool ApplyC2(int det, int side, int strip, double position, double &energy) const;
        bool ApplyC3(int det, int side, int strip, double &energy) const;
        bool ParseSourceEnergies(TString sourceEnergies);
        double FitPeak(TH1D *hist, TString fitName="fit") const;
        void ApplyPaveTextStyle(TPaveText *text) const;
        void AddC0FitAnnotations(TH1D *hist, const std::vector<double> &amps, const std::vector<double> &means, const std::vector<double> &sigmas) const;
        void AddC2FitAnnotation(TH2D *hist, const std::array<double,3> &fitPar, const std::array<double,3> &c2Par, double entries) const;
        bool GetGateDisplayRange(double axisMin, double axisMax, double &rangeMin, double &rangeMax) const;
        void SetEnergyGateDisplayRange(TH1D *hist) const;
        void SetEnergyGateDisplayRange(TH2D *hist) const;
        void SetEnergyGateDisplayRangeXY(TH2D *hist) const;
        void AddEnergyPositionGateLines(TH2D *hist) const;
        void AddEnergyPositionGateLines();
        void ApplyStageEnergyZoom();

        void FitAndWriteC0();
        void FitAndWriteC1();
        void FitAndWriteC2();
        void FitAndWriteC3();
        void WriteHistograms();
        void WriteC0ParameterContainer() const;
        void WriteC1ParameterContainer() const;
        void WriteC2ParameterContainer() const;
        void WriteC3ParameterContainer() const;
        void WriteCombinedParameterContainer() const;
        void WriteCombinedCalibrationFile();
        void InitializeExpectedChannels();
        void WriteChannelSummary() const;
        double GetCalibratedEnergyAxisMax() const;
        TString GetSourceEnergyListText() const;
        TString MakeFileName(TString suffix) const;
        TString MakeHistName(TString quantity, int det, int side, int strip) const;

    private:
        TClonesArray *fSiChannelArray = nullptr;
        LKSiliconArray *fSiArray = nullptr;

        int fStage = 0;
        double fSourceEnergy = 0;
        TString fSourceEnergyText = "";
        std::vector<double> fSourceEnergies;
        double fC2ReferenceEnergy = 5.486;
        double fExpectedResolution = 0.01;
        double fGateDisplayRangeScale = 3.;
        double fEntriesCut = 100;
        TString fOutputPath = "data_calibration";
        TString fOutputPrefix = "";
        TString fInputC0FileName = "";
        TString fInputC1FileName = "";
        TString fInputC2FileName = "";
        TString fInputC3FileName = "";
        TString fListFileName = "";

        std::map<StripKey, TH1D *> fEnergyHistMap;
        std::map<StripKey, TH1D *> fEnergySumHistMap;
        std::map<StripKey, TH2D *> fLeftRightHistMap;
        std::map<StripKey, TH2D *> fPositionEnergyHistMap;
        std::vector<std::pair<StripKey, TH1D *>> fFitHistArray;

        std::map<StripKey, bool> fExpectedStandaloneMap;
        std::map<StripKey, bool> fExpectedPairedMap;

        std::map<StripKey, std::array<double,2>> fC0Parameters;
        std::map<LRKey, std::array<double,2>> fC1Parameters;
        std::map<StripKey, std::array<double,3>> fC2Parameters;
        std::map<StripKey, std::array<double,2>> fC3Parameters;

    ClassDef(SKSiCalibrationTask, 1)
};

#endif
