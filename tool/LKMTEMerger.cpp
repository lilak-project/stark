#include "LKMTEMerger.h"
#include <stdio.h>
#include <iostream>
#include <limits.h>
#include <float.h>
using namespace std;

#ifndef LILAK_VERSION
ClassImp(LKMTEMerger)
#endif

LKMTEMerger::LKMTEMerger(TFile* outputFile)
{
    fOutputFile = outputFile;
}

LKMTEMerger::LKMTEMerger(TString outputFileName)
{
    fOutputFileName = outputFileName;
}

void LKMTEMerger::SetKeySource(int no, TString name)
{
    fInputSourceArray.push_back(no);
    fInputSourceName[no] = name;
    fKeySourceNumber = no;
}

void LKMTEMerger::SetInputSource(int no, TString name)
{
    fInputSourceArray.push_back(no);
    fInputSourceName[no] = name;
}

bool LKMTEMerger::ReadMTE(TString inputFileName)
{
    FILE *fileMTE;
    unsigned int fileSize;
    char data[16];
    int extPattern[NUM_MTE_SOURCES];
    int iTmp;
    double fTmp;
    int numEvents;

    if (fOutputFile==nullptr)
        fOutputFile = new TFile(fOutputFileName,"recreate");

    for (auto iSource : fInputSourceArray)
    {
        fMTETree[iSource] = new TTree(Form("mte_%s",fInputSourceName[iSource].Data()),"");
        fMTETree[iSource] -> Branch("num",&bTrigNum);
        fMTETree[iSource] -> Branch("time",&bTrigTime);
        fMTETree[iSource] -> Branch("type",&bTrigType);
    }

    if (inputFileName.IsNull()) {
        e_error << "Cannot find " << inputFileName << endl;
        return false;
    }
    e_info << "MTE file: " << inputFileName << endl;
    fileMTE = fopen(inputFileName, "rb");
    if (fileMTE==nullptr) {
        e_error << "Cannot open " << inputFileName << endl;
        return false;
    }

    // get file size to know # of events, 1 event = 32 byte
    fseek(fileMTE, 0L, SEEK_END);
    fileSize = ftell(fileMTE);
    fclose(fileMTE);
    numEvents = fileSize / 16;
    fileMTE = fopen(inputFileName, "rb");

    for (int evt=0; evt<numEvents; evt++)
    {
        fread(data, 1, 16, fileMTE);

        // trigger logic #
        memcpy(&bTrigNum, data, 4);

        // trigger time
        iTmp = data[4] & 0xFF;
        fTmp = iTmp;
        bTrigTime = fTmp * 0.008;        // trig_ftime = 8 ns unit
        iTmp = data[5] & 0xFF;
        fTmp = iTmp;
        bTrigTime = bTrigTime + fTmp;
        iTmp = data[6] & 0xFF;
        fTmp = iTmp;
        fTmp = fTmp * 256.0;
        bTrigTime = bTrigTime + fTmp;
        iTmp = data[7] & 0xFF;
        fTmp = iTmp;
        fTmp = fTmp * 256.0 * 256.0;
        bTrigTime = bTrigTime + fTmp;
        iTmp = data[8] & 0xFF;
        fTmp = iTmp;
        fTmp = fTmp * 256.0 * 256.0 * 256.0;
        bTrigTime = bTrigTime + fTmp;
        iTmp = data[9] & 0xFF;
        fTmp = iTmp;
        fTmp = fTmp * 256.0 * 256.0 * 256.0 * 256.0;
        bTrigTime = bTrigTime + fTmp;
        iTmp = data[NUM_MTE_SOURCES] & 0xFF;
        fTmp = iTmp;
        fTmp = fTmp * 256.0 * 256.0 * 256.0 * 256.0 * 256.0;
        bTrigTime = bTrigTime + fTmp;

        // trigger type
        bTrigType = data[11] & 0xFF;

        // external trigger pattern
        extPattern[0] = data[12] & 0x1;
        extPattern[1] = (data[12] >> 1) & 0x1;
        extPattern[2] = (data[12] >> 2) & 0x1;
        extPattern[3] = (data[12] >> 3) & 0x1;
        extPattern[4] = (data[12] >> 4) & 0x1;
        extPattern[5] = (data[12] >> 5) & 0x1;
        extPattern[6] = (data[12] >> 6) & 0x1;
        extPattern[7] = (data[12] >> 7) & 0x1;
        extPattern[8] = data[13] & 0x1;
        extPattern[9] = (data[13] >> 1) & 0x1;

        for (auto iSource : fInputSourceArray)
        {
            if (extPattern[iSource]>0)
                fMTETree[iSource] -> Fill();
        }
    }

    fclose(fileMTE);

    fOutputFile -> cd();
    for (auto iSource : fInputSourceArray) {
        fMTEEntries[iSource] = fMTETree[iSource] -> GetEntries();
        if (fMaxEntries<fMTEEntries[iSource])
            fMaxEntries=fMTEEntries[iSource];
        fMTETree[iSource] -> Write();
    }
    e_info << "Output file: " << fOutputFile->GetName() << endl;

    return true;
}

void LKMTEMerger::WriteSummary(bool writeToFile)
{
    fOutputFile -> cd();
    fTreeSummary = new TTree("mte_summary","");
    fTreeSummary -> Branch("num" ,&bKeyNum);
    fTreeSummary -> Branch("type",&bKeyType);
    fTreeSummary -> Branch("time",&bKeyTime);
    for (auto iSource : fInputSourceArray) {
        if (fKeySourceNumber==iSource)
            continue;
        fTreeSummary -> Branch(Form("mte_id%d",iSource), &bMTEEntry[iSource]);
        fTreeSummary -> Branch(Form("daq_id%d",iSource), &bDAQEntry[iSource]);
    }

    double fTimeWindowCut = 5;

    for (auto key_entry=0; key_entry<fMTEEntries[fKeySourceNumber]; ++key_entry)
    {
        fMTETree[fKeySourceNumber] -> GetEntry(key_entry);
        bKeyNum  = bTrigNum;
        bKeyTime = bTrigTime;
        bKeyType = bTrigType;
        for (auto iSource : fInputSourceArray)
        {
            if (fKeySourceNumber==iSource)
                continue;

            bMTEEntry[iSource] = -1;
            bDAQEntry[iSource] = -1;
            double time_diff = DBL_MAX;
            auto n = fMTEEntries[iSource];
            for (auto entry=0; entry<n; ++entry)
            {
                fMTETree[iSource] -> GetEntry(entry);
                time_diff = bKeyTime - bTrigTime - fMTETimeOffset[iSource];
                if (abs(time_diff)<fTimeWindowCut)
                {
                    bMTEEntry[iSource] = entry + 1;
                    bDAQEntry[iSource] = entry + fEntryOffset[iSource];
                }
                else if (time_diff<0)
                {
                    break;
                }
            }
        }

        fTreeSummary -> Fill();
    }

    if (writeToFile) {
        fOutputFile -> cd();
        fTreeSummary -> Write();
    }
}

void LKMTEMerger::GetTimeOffset(TString fileName)
{
    int iSource;
    double offset;
    e_info << "Reading " << fileName << endl;
    ifstream offsetFile(fileName);
    while (offsetFile >> iSource >> offset) {
        fMTETimeOffset[iSource] = offset;
        e_info << fInputSourceName[iSource] << " time-offset = " << fMTETimeOffset[iSource] << endl;
    }
}

void LKMTEMerger::FindTimeOffset(TString fileName)
{
    ///////////////////////////////////////////////////////////////////////
    // Find tree with maximum bTrigNum at entry 0
    int    num [NUM_MTE_SOURCES] = {0};
    double time[NUM_MTE_SOURCES] = {0};
    int    type[NUM_MTE_SOURCES] = {0};
    for (auto iSource : fInputSourceArray) {
        num [iSource] = 0;
        time[iSource] = 0;
        type[iSource] = 0;
    }
    for (auto iSource : fInputSourceArray)
    {
        fMTETree[iSource] -> GetEntry(0);
        num [iSource] = bTrigNum;
        time[iSource] = bTrigTime;
        type[iSource] = bTrigType;
    }
    int numAtMax = 0;
    double timeAtMax = 0;
    int iSourceAtMax = 0;
    for (auto iSource : fInputSourceArray)
    {
        if (numAtMax<num[iSource])  {
            numAtMax = num[iSource];
            timeAtMax = time[iSource];
            iSourceAtMax = iSource;
        }
    }

    ///////////////////////////////////////////////////////////////////////
    // For each tree, find closest time with timeAtMax
    double timeAtMin[NUM_MTE_SOURCES] = {0};
    timeAtMin[iSourceAtMax] = timeAtMax;
    for (auto iSource : fInputSourceArray)
    {
        if (iSourceAtMax==iSource)
            continue;
        double time_diff_min = DBL_MAX;
        for (auto entry=0; entry<fMaxEntries; ++entry)
        {
            fMTETree[iSource] -> GetEntry(entry);
            time[iSource] = bTrigTime;
            double time_diff = abs(time[iSource]-timeAtMax);
            if (time_diff_min>time_diff)
            {
                time_diff_min = time_diff;
                timeAtMin[iSource] = bTrigTime;
            }
            else if (time_diff_min<time_diff)
                break;
        }
    }

    ofstream offsetFile(fileName);
    for (auto iSource : fInputSourceArray) {
        fMTETimeOffset[iSource] = timeAtMin[fKeySourceNumber] - timeAtMin[iSource];
        e_info     << fInputSourceName[iSource] << " time-offset = " << fMTETimeOffset[iSource] << endl;
        offsetFile << iSource << " " << fMTETimeOffset[iSource]  << endl;
    }
}

bool LKMTEMerger::ConfigureKobraFile(TString fileName, TString kobraName)
{
    if (fKobraTree!=nullptr)
        return true;

    for (auto iSource : fInputSourceArray) {
        if (fInputSourceName[iSource]==kobraName) {
            fIKobraTree = iSource;
            break;
        }
    }
    if (fIKobraTree<0) {
        e_error << "MTE tree for Kobra do not exist!" << endl;
        return false;
    }

    fKobraFile = new TFile(fileName,"read");
    if (fKobraFile->IsOpen()==false) {
        e_error << "Kobra file do not exist! " << fileName << endl;
        return false;
    }
    else
        e_info << "Kobra file: " << fileName << endl;

    fKobraTree = (TTree*) fKobraFile -> Get("midas_data");

    return true;
}

bool LKMTEMerger::MapKobra(TString fileName, TString kobraName)
{
    ConfigureKobraFile(fileName, kobraName);

    e_info << "Mapping Kobra entry to Cobo entry ..." << endl;

    int eventidT1, ref_pulse, scaler, scaler_prev;
    fKobraTree -> SetBranchAddress("scaler",    &scaler);

    scaler_prev = INT_MAX;
    fKobraTree -> GetEntry(fTestKobraDAQEntry);
    scaler_prev = scaler;
    fKobraTree -> GetEntry(fTestKobraDAQEntry);
    fKobraTree -> GetEntry(fTestKobraDAQEntry+1);
    int dscaler = scaler - scaler_prev;

    //cout.precision(9);

    double time_prev = DBL_MAX;
    double time_diff = 0;
    int entry;
    for (entry=0; entry<fMTEEntries[fIKobraTree]; ++entry)
    {
        time_prev = bTrigTime;
        fMTETree[fIKobraTree] -> GetEntry(entry);

        if (entry>0) {
            time_diff = (bTrigTime - time_prev)*100;

            double time_window = dscaler*fMTEKobraTimeWindowFactor; // XXX
            if (abs(dscaler-time_diff)<time_window) {
                e_info << "Matching event! cobo-entry=" << entry << ", dsacaler-time_diff=" << dscaler - time_diff << ", time-window=" << time_window << endl;
                break;
            }
        }
    }

    fEntryOffset[fIKobraTree] = fTestKobraDAQEntry+2 - entry;
    e_info << fInputSourceName[fIKobraTree] << " " << fEntryOffset[fIKobraTree] << endl;

    return true;
}

bool LKMTEMerger::TestKobraEntry(int maxEntry)
{
    int eventidT1, ref_pulse, scaler;
    fKobraTree -> SetBranchAddress("eventidT1", &eventidT1);
    fKobraTree -> SetBranchAddress("ref_pulse", &ref_pulse);
    fKobraTree -> SetBranchAddress("scaler",    &scaler);

    int entryKobra;
    int numSummary = fTreeSummary -> GetEntries();
    if (maxEntry>0 && numSummary>maxEntry)
        numSummary = maxEntry;
    for (int entryCobo=0; entryCobo<numSummary; ++entryCobo)
    {
        eventidT1=0;
        ref_pulse=0;
        scaler=0;

        fTreeSummary -> GetEntry(entryCobo);
        entryKobra = bDAQEntry[fIKobraTree];
        if (entryKobra<0)
            continue;

        fKobraTree -> GetEntry(entryKobra);
        e_test << "cobo=" << entryCobo << " kobra=" << entryKobra <<  " eventidT1=" << eventidT1 << " ref_pulse=" << ref_pulse << " scaler=" << scaler << endl;
    }
    return true;
}


int LKMTEMerger::GetKobraEntry(int entryCobo)
{
    fTreeSummary -> GetEntry(entryCobo);
    auto entryKobra = bDAQEntry[fIKobraTree];
    if (entryKobra>=0)
        fKobraTree -> GetEntry(entryKobra);
    return entryKobra;
}
