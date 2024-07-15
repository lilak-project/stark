#ifndef LKSIARRAYPLANE_HH
#define LKSIARRAYPLANE_HH

#include "LKEvePlane.h"
#include "LKPad.h"
#include "TPad.h"
#include "LKParameterContainer.h"
#include <vector>
#include "LKSiDetector.h"
#include "LKSiChannel.h"

class SKSiArrayPlane : public LKEvePlane
{
    public:
        SKSiArrayPlane();
        SKSiArrayPlane(const char *name, const char *title);
        virtual ~SKSiArrayPlane() {};

        //virtual void Clear(Option_t *option = "");
        //virtual void Print(Option_t *option = "") const;
        virtual bool Init();
        virtual void Print(Option_t *option="") const;

        //virtual TVector3 GetPositionError(int padID);

    public:
        virtual void FillDataToHistEventDisplay1(Option_t *option="");
        //virtual void FillDataToHistEventDisplay2(Option_t *option="");

        virtual bool SetDataFromBranch();

        virtual TPad* Get3DEventPad();

        virtual TH2* GetHistUserDrawing(Option_t *option="-1");
        virtual TH2* GetHistControlDataDP(Option_t *option="-1");
        virtual TH2* GetHistEventDisplay1(Option_t *option="-1");
        virtual TH2* GetHistEventDisplay2(Option_t *option="-1");
        virtual TH1D* GetHistChannelBuffer();

        virtual int FindPadID(int cobo, int asad, int aget, int chan);
        //virtual LKPad* FindPad(int cobo, int asad, int aget, int chan);
        virtual int FindPadIDFromHistEventDisplay1Bin(int hbin);
        //virtual int FindZFromHistEventDisplay2Bin(int hbin);

        int GetNumSiDetectors() const { return fDetectorArray -> GetEntries(); }
        LKSiDetector *GetSiDetector(int idx);
        LKSiDetector *GetSiDetector(int cobo, int asad, int aget, int chan);
        int FindSiDetectorID(int cobo, int asad, int aget, int chan);

        LKSiDetector *GetSiDetectorData(int idx);

        LKSiChannel *GetSiChannel(int idx);
        LKSiChannel *GetSiChannel(int cobo, int asad, int aget, int chan);
        int FindSiChannelID(int cobo, int asad, int aget, int chan) { return FindPadID(cobo, asad, aget, chan); }
        bool SetSiChannelData(LKSiChannel* siChannel, GETChannel* channel);

        virtual TCanvas *GetCanvas(Option_t *option="");
        virtual TH2* GetHist(Option_t *option="");
        void ExecMouseClickEventOnPad(TVirtualPad *pad, double xOnClick, double yOnClick);
        virtual void ClickedEventDisplay1(double xOnClick, double yOnClick);
        virtual void ClickedEventDisplay2(double xOnClick, double yOnClick);
        virtual void ClickedJOSideDisplay(int jo);
        virtual void ClickedUserDrawing(double xOnClick, double yOnClick);
        virtual void ClickedControlDataDP(double xOnClick, double yOnClick);

        virtual void UpdateAll();
        virtual void UpdateEventDisplay1();
        virtual void UpdateEventDisplay2();
        virtual void UpdateJunctionOhmic();
        //virtual void UpdateChannelBuffer();
        virtual void UpdateDataDisplays();
        virtual void UpdateUserDrawing();

        bool AddUserDrawings(TString label, int detID, int joID, TObjArray* userDrawingArray, int leastNDraw=1);

    protected:
        TString fMappingFileName;
        TString fDetectorParName;
        int fNumCobo = 6;
        int fNumAsad = 4;
        int fNumAget = 4;
        int fNumChan = 68;

        int ****fMapCAACToChannelIndex; ///< [cobo][asad][aget][chan] to channel-index mapping
        int ****fMapCAACToDetectorIndex; ///< [cobo][asad][aget][chan] to detector-index mapping
        int *fMapBinToDetector; ///< histogram bin to pad-id mapping

        int fCurrentView = 1;
        double fThreshold = 300;

        int fSelDetID = 0;
        int fSelJOID = -1;

        int fMaxLayerIndex = -INT_MAX;
        double fMinPhi = DBL_MAX;
        double fMaxPhi = -DBL_MAX;
        /// If histogram bin size is smaller than 1, mouse click position might not be pointing 
        /// inside the bin since position is returned as integer.
        /// This parameter will increase resolution of histogram bin.
        double fLayerYScale = 100;
        double fLayerYOffset = 0.05;

    protected:
        int fPrevEventID = -1;
        int fPrevEventID2 = -1;

        TObjArray *fDetectorArray = nullptr;
        TObjArray *fChannelArray = nullptr;
        LKParameterContainer fDetectorTypeArray;

        TH2* fFrameEventDisplay1 = nullptr;
        int fMaxDetectors = 40;

        TPad* fPadJOSideDisplay[2];
        //TPad* fPadJSideDisplay = nullptr;
        //TPad* fPadOSideDisplay = nullptr;

        TPad* fPadCtrlUserDrawing = nullptr;
        TH2* fHistCtrlUserDrawing = nullptr;
        int fNumUserDrawing = 0;
        LKParameterContainer fParUserDrawing;
        TObjArray* fUserDrawingArrayCollection = nullptr;
        int ****fUserDrawingArrayIndex; // [4][7][40][2] menu-tab, menu-number, detector-number, junction/ohmic
        int ****fUserDrawingLeastNDraw; // [4][7][40][2]
        TString fUserDrawingName[4][7];

        TPad* fPadControlDataDP = nullptr;
        int fBinCtrlPrevPage;
        int fBinCtrlCurrPage;
        int fBinCtrlNextPage;
        TPad* fPadDataDisplayFull;
        TPad* fPadDataDisplayTwo[2];
        TPad* fPadDataDisplayThree[3];
        int fNumDDX = 2;
        int fNumDDY = 3;
        int fNumDataDisplays = fNumDDX*fNumDDY;
        TPad* fPadDataDisplaySmall[24];
        TH2D* fHistControlDataDP = nullptr;
        int fNumUDPage = 0;

        int fSelUDTab = 0;
        int fSelUDBin = -1;
        int fNumUDLabels = 0;
        int fUDGoodIndex[10][16];
        int fNumDrawingsInArray = 0;
        int fNumGoodDrawings = 0;
        int fSelUDArrayID = 0;
        int fSelUDLeastNDraw = 1;
        int fSelUDPage = 0;
        int fNumSelUDGroup = 0;
        TObjArray* fSelUDArray = nullptr;
        bool fSignalUDTabChange = true;
        bool fSignalUDBinChange = true;

        //double fCtrlLabelSize = 0.18;
        double fCtrlLabelOffset = 0.10;

        TGraph* fGSelJOSideDisplay = nullptr;

    protected:
        int kUpdateUserDrawing   = 6;
        int kUpdateDataDisplays  = 7;
        int kUpdateJunctionOhmic = 8;

    private:
        int fCountHistArray1 = 0;
        int fCountHistArray2 = 0;
        TClonesArray* fHistArray1 = nullptr;
        TClonesArray* fHistArray2 = nullptr;

        int fNumArrayChannelDrawing  = 0;
        int fNumArrayDetectorDrawing = 0;
        int fNumArrayPairDrawing     = 0;
        int fNumArrayEventDrawing    = 0;
        LKParameterContainer* fArrayChannelDrawing  = nullptr;
        LKParameterContainer* fArrayPairDrawing     = nullptr;
        LKParameterContainer* fArrayDetectorDrawing = nullptr;
        LKParameterContainer* fArrayEventDrawing    = nullptr;

    ClassDef(SKSiArrayPlane, 1)
};

#endif
