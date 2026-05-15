#ifndef ELARK_HH
#define ELARK_HH

#include "LKDetectorPlane.h"
#include "LKPad.h"
#include "TPad.h"
#include "TH2Poly.h"
#include "LKParameterContainer.h"
#include <vector>
#include "LKSiDetector.h"
#include "LKSiChannel.h"
#include "LKDrawing.h"

class ELARK : public LKDetectorPlane
{
    public:
        ELARK();
        ELARK(const char *name, const char *title);
        virtual ~ELARK() {};

        virtual bool Init();
        virtual bool Init2();
        virtual void Print(Option_t *option="") const;

    public:
        virtual bool SetDataFromBranch() { return true; }
        virtual void FillDataToHist() {}

        TH2Poly* NewHistPlane(TString name);
        LKDrawing* NewHistPlaneDrawing(TString name); 

        virtual int FindPadID(int cobo, int asad, int aget, int chan);

        int GetNumSiDetectors() const { return fDetectorArray -> GetEntries(); }
        LKSiDetector *GetSiDetector(int idx);
        LKSiDetector *GetSiDetector(int cobo, int asad, int aget, int chan);
        int FindSiDetectorID(int cobo, int asad, int aget, int chan);
        int FindEPairDetectorID(int det, int i=-1);

        LKSiDetector *GetSiDetectorData(int idx);

        LKSiChannel *GetSiChannel(int idx);
        LKSiChannel *GetSiChannel(int cobo, int asad, int aget, int chan);
        int FindSiChannelID(int cobo, int asad, int aget, int chan) { return FindPadID(cobo, asad, aget, chan); }
        bool SetSiChannelData(LKSiChannel* siChannel, GETChannel* channel);

        void FireStrip(int det, int side, int strip, double energy);
        void ClearFiredFlags();
        int GetNumFiredDetectors();
        LKSiDetector* GetFiredDetector(int iFired);

        int GetEPairDetID(int pair, int i) { return fdEEPairMapping[pair][i]; }

        double CalculatePairSolidAngle(double theta1, double theta2, int numDetectors=1);
        double Calculate16RingSolidAngle(double theta1, double theta2, int numDetectors=1);

    protected:
        TString fMappingFileName;
        TString fDetectorParName;
        int fNumCobo = 6;
        int fNumAsad = 4;
        int fNumAget = 4;
        int fNumChan = 68;

        int ****fMapCAACToChannelIndex; ///< [cobo][asad][aget][chan] to channel-index mapping
        int ****fMapCAACToDetectorIndex; ///< [cobo][asad][aget][chan] to detector-index mapping

        int fMaxLayerIndex = -INT_MAX;
        double fMinPhi = DBL_MAX;
        double fMaxPhi = -DBL_MAX;

        int fdEEPairMapping[40][2];
        /// If histogram bin size is smaller than 1, mouse click position might not be pointing 
        /// inside the bin since position is returned as integer.
        /// This parameter will increase resolution of histogram bin.
        double fLayerYScale = 100;
        double fUserPhiScale = 1.;
        double fLayerYOffset = 0.05;

    protected:
        TObjArray *fDetectorArray = nullptr;
        TObjArray *fChannelArray = nullptr;
        LKParameterContainer fDetectorTypeArray;

    ClassDef(ELARK, 1)
};

#endif
