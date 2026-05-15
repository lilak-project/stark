#ifndef LKMTEMERGER_HH
#define LKMTEMERGER_HH

#include "LKLogger.h"
#include "TFile.h"
#include "TFile.h"
#include "TTree.h"
#include <vector>

#define NUM_MTE_SOURCES 10

#ifndef LILAK_VERSION
//#define e_info cout<<"\033[0;32m==\033[0m "
//#define e_error cout<<"\033[0;31mER\033[0m "
//#define e_test cout<<"\033[0;36mtest\033[0m "
#endif

/// Master trigger electronics
class LKMTEMerger
{
    public:
        LKMTEMerger(TFile* outputFile);
        LKMTEMerger(TString outputFileName);
        virtual ~LKMTEMerger() {};

        void SetKeySource  (int no, TString name); ///< Set Key channel bit number and name
        void SetInputSource(int no, TString name); ///< Set channel bit number and name (not key channel)
        void SetMTEKobraTimeWindowFactor(double factor) { fMTEKobraTimeWindowFactor = factor; }
        void SetTestKobraDAQEntry(int testEntry) { fTestKobraDAQEntry = testEntry; }

        void GetTimeOffset  (TString fileName="data/mte_time_offset.txt"); ///< Read and store time offset data from fileName
        void FindTimeOffset (TString fileName="data/mte_time_offset.txt"); ///< Read single source trees and find time offset between them, and save out to fileName
        bool ReadMTE        (TString inputFileName); ///< Read MTE to create separate trees for different sources
        bool MapKobra       (TString fileName, TString kobraName="Kobra"); ///< Read Kobra DAQ data and match time-stamp with MTE Kobra tree by matching time-diffs between events
        void WriteSummary   (bool writeToFile=true); ///< Write summary tree with entry matching data
        bool TestKobraEntry (int maxEntry=0); ///< Read Kobra DAQ to find matching entry of the Key source entry. Print up to maxEntry of cobo entries
        int  GetKobraEntry  (int entryCobo); ///< XXX Get koobra daq entry regarding to cobo daq entry

    public:
        TTree *GetKobraTree() { return fKobraTree; }
        bool ConfigureKobraFile(TString fileName, TString kobraName="Kobra");

    private:
        TString fOutputFileName = "";

        int fKeySourceNumber;
        TString fInputSourceName[NUM_MTE_SOURCES];
        std::vector<int> fInputSourceArray;
        int fEntryOffset[NUM_MTE_SOURCES] = {0};
        double fMTETimeOffset[NUM_MTE_SOURCES] = {0};

        // MTE tree
        int    bTrigNum = 0;
        int    bTrigType = 0;
        double bTrigTime = 0;

        // Summary tree
        int    bKeyNum  = 0;
        double bKeyTime = 0;
        int    bKeyType = 0;
        int    bMTEEntry[NUM_MTE_SOURCES] = {0};
        int    bDAQEntry[NUM_MTE_SOURCES] = {0};

        TTree* fMTETree[NUM_MTE_SOURCES];
        int fMTEEntries[NUM_MTE_SOURCES];
        int fMaxEntries = 0;

        double fMTEKobraTimeWindowFactor = 1./50000;

        TFile *fOutputFile = nullptr;
        TTree *fTreeSummary = nullptr;

        TFile *fKobraFile = nullptr;
        TTree *fKobraTree = nullptr;
        int fIKobraTree = -1;
        int fTestKobraDAQEntry = 5;

#ifndef LILAK_VERSION
    ClassDef(LKMTEMerger, 1)
#endif
};

#endif
