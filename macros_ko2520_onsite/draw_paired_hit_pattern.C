void draw_paired_hit_pattern()
{
    auto par = new LKParameterContainer("common_run.mac");
    auto runID = par -> GetParInt("LKRun/RunID");
    auto file = new TFile(Form("data_reco/ko2520_%04d.reco.root",runID));
    auto tree = (TTree*) file -> Get("event");


    TString nameTop = Form("ko2520_%04d_hit_pattern",runID);
    auto top = new LKDrawingGroup(nameTop);
    TH2D *hist[4][4] = {0};
    TString asadTitle;
    TString agetTitle;
    for (auto asad=0; asad<4; ++asad)
    {
        if (asad==0) asadTitle = "[QQQ5]"; 
        else asadTitle = "[X6]"; 
        for (auto aget=0; aget<4; ++aget)
        {
            if (aget==0) agetTitle = "(Ohmic)"; 
            else agetTitle = "(Junction)"; 
            TString cut;
            cut = cut + Form("RawData.fAsad==%d",asad);
            cut = cut + Form("&&RawData.fAget==%d",aget);
            TString histName = Form("hist_%d_%d",asad,aget);
            TString histTitle = asadTitle + " " + agetTitle + " " + Form("AsAd=%d Aget=%d;channel (strip);energy",asad,aget);
            cout << histName << " " << histTitle << endl;
            if (asad!=0 && aget!=0)
                hist[asad][aget] = new TH2D(histName,histTitle,32,0,32,200,0,2000);
            else
                hist[asad][aget] = new TH2D(histName,histTitle,64,0,64,200,0,2000);
            //if (asad==1&&aget==2) top -> CreateDrawing() -> Add(hist[asad][aget]);
            //if (asad==1) top -> CreateDrawing() -> Add(hist[asad][aget]);
            top -> CreateDrawing() -> Add(hist[asad][aget]);
        }
    }

    double energySum = 0;
    int asadp = -1;
    int agetp = -1;
    int chanp = 0;

    TClonesArray *bRawDataArray = nullptr;
    tree -> SetBranchAddress("RawData",&bRawDataArray);
    auto numEvents = tree -> GetEntries();
    cout << numEvents << endl;
    for (auto event=0; event<numEvents; ++event)
    {
        tree -> GetEntry(event);
        auto numRawData = bRawDataArray -> GetEntries();
        if (event%5000==0) cout << "event=" << event << " #channels=" << numRawData << endl;

        energySum = 0;
        asadp = -1;
        agetp = -1;
        chanp = -1;

        if (numRawData>3) continue;

        for (auto iRawData=0; iRawData<numRawData; ++iRawData)
        {
            auto rawData = (GETChannel*) bRawDataArray -> At(iRawData);
            auto asad0 = rawData -> GetAsad();
            auto aget0 = rawData -> GetAget();

            auto chan = rawData -> GetChan();
            if (chan==56) continue;
            if (chan==45) continue;
            if (chan==22) continue;
            if (chan==11) continue;
            auto chan2 = chan;
            if (chan2>56) chan2 = chan2-1;
            if (chan2>45) chan2 = chan2-1;
            if (chan2>22) chan2 = chan2-1;
            if (chan2>11) chan2 = chan2-1;
            auto chanf = chan2;

            auto energy = rawData -> GetEnergy();
            auto energyf = energy;

            if (asad0!=0 && aget0!=0)
            {
                if (chanp<0) {
                    energySum = energy;
                    asadp = asad0;
                    agetp = aget0;
                    chanp = chan2;
                    continue;
                }
                else if (asadp==asad0 && agetp==aget0 && abs(chanp-chan2)==1){
                    energySum += energy;
                    energyf = energySum;
                    if (chanf%2==1) chanf = chanf-1;
                    chanf = chanf/2;
                }
            }

            if (event%5000==0) cout << "  asad=" << asad0 << " aget=" << aget0 << endl;
            for (auto asad=0; asad<4; ++asad)
            {
                for (auto aget=0; aget<4; ++aget)
                {
                    if (asad0==asad && aget0==aget) {
                        hist[asad][aget] -> Fill(chanf,energyf);
                    }
                }
            }
        }
    }
    top -> Draw("viewer");
    top -> GetCanvas() -> SaveAs(TString("figures/")+nameTop+".png");
    top -> GetCanvas() -> SaveAs(TString("figures/")+nameTop+".pdf");
}
