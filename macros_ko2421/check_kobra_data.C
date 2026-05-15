void check_kobra_data()
{
    TString fileNameDAQ = "data_kobra_daq/output00661.root";
    TString fileNameAna = "data_kobra_ana/output00000661.root";

    if (1)
    {
        auto fileDAQ = new TFile(fileNameDAQ,"read");
        auto treeDAQ = (TTree*) fileDAQ -> Get("midas_data");
        //cout << fileNameDAQ << " " << treeDAQ -> GetEntries() << endl;
        //treeDAQ -> Scan("HW_trigger_after_latch:HW_trigger_before_latch");
        int scaler, HW_trigger_after_latch, HW_trigger_before_latch;
        treeDAQ -> SetBranchAddress("scaler",&scaler);
        treeDAQ -> SetBranchAddress("HW_trigger_after_latch",&HW_trigger_after_latch);
        treeDAQ -> SetBranchAddress("HW_trigger_before_latch",&HW_trigger_before_latch);
        auto numAna = 30;
        uint64_t t1;
        uint64_t t2;
        for (auto i=0; i<numAna; ++i)
        {
            treeDAQ -> GetEntry(i);
            t2 = t1;
            t1 = scaler;
            cout << HW_trigger_after_latch << " " << HW_trigger_before_latch << " " <<  t1 << " dt=" << t1 - t2 << endl;
        }
    }

    if (1)
    {
        auto fileAna = new TFile(fileNameAna,"read");
        auto treeAna = (TTree*) fileAna -> Get("kobra");
        //cout << fileNameAna << " " << treeAna -> GetEntries() << endl;
        //treeAna -> Scan("scaler.hwtriga:scaler.hwtrigb");
        TClonesArray *array = nullptr;
        treeAna -> SetBranchAddress("scaler",&array);
        //auto numAna = treeAna -> GetEntries();
        auto numAna = 30;
        uint64_t t1;
        uint64_t t2;
        for (auto i=0; i<numAna; ++i)
        {
            treeAna -> GetEntry(i);
            auto n = array -> GetEntries();
            if (n!=1) continue;
            auto scaler_data = (TScalerData*) array -> At(0);
            t2 = t1;
            t1 = scaler_data->ts;
            cout << scaler_data->hwtriga << " " << scaler_data->hwtrigb << " " << t1 << " dt=" << t1 - t2 << endl;
        }
    }
}
