bool is_ohmic(int asad, int aget) { if (aget==0) return true; return false; }
bool is_qqq5(int asad, int aget) { if (asad==0 && (aget==2 || aget==3)) return true; return false; }
bool is_x6(int asad, int aget) { if (is_ohmic(asad,aget)) return false; if (is_qqq5(asad,aget)) return false; return true; }

void draw_update_hit_pattern()
{
    gStyle -> SetPalette(kRainBow);

    int multRange1 = 2;
    int multRange2 = 5;
    int printAfter = 50000;
    double emax = 4200;

    auto par = new LKParameterContainer("get_common.mac");
    auto runID = par -> GetParInt("LKRun/RunID");
    cout << runID << endl;

    TChain *tree = nullptr;

    if (runID==8000) {
        tree = new TChain("event");
        tree -> Add(Form("data_reco/ko2520_%04d.reco.root",83));
        tree -> Add(Form("data_reco/ko2520_%04d.reco.root",84));
        tree -> Add(Form("data_reco/ko2520_%04d.reco.root",85));
        tree -> Add(Form("data_reco/ko2520_%04d.reco.root",86));
        tree -> Add(Form("data_reco/ko2520_%04d.reco.root",87));
        tree -> Add(Form("data_reco/ko2520_%04d.reco.root",88));
        tree -> Add(Form("data_reco/ko2520_%04d.reco.root",89));
        tree -> Add(Form("data_reco/ko2520_%04d.reco.root",90));
    }
    else if (runID == 8003) {
        // Production runs acquired on 2026-09-04.
        tree = new TChain("event");
        for (int inputRun:{134,135,136,139,141,142,146,147,148})
            tree->Add(Form("data_reco/ko2520_%04d.reco.root",inputRun));
    }
    else if (runID == 8004)
    {
        tree = new TChain("event");
        for (int inputRun:{134, 135, 136, 139, 141, 142}) {
            cout << Form("data_reco/ko2520_%04d.reco.root",inputRun) << endl;
            TString fileName = Form("data_reco/ko2520_%04d.reco.root",inputRun);
            tree->Add(fileName);
        }
    }
    else {
        TString fileName = Form("data_reco/ko2520_%04d.reco.root",runID);
        cout << fileName << endl;
        auto file = new TFile(fileName);
        tree = (TChain*) file -> Get("event");
    }

    TString nameTop = Form("ko2520_%04d",runID);
    TString nameTop1 = Form("ko2520_%04d_paired",runID);
    TString nameTop2 = Form("ko2520_%04d_energy",runID);
    TString nameTop3 = Form("ko2520_%04d_et",runID);
    TString nameTop4 = Form("ko2520_%04d_np",runID);
    TString nameTop5 = Form("ko2520_%04d_q",runID);
    auto top1 = new LKDrawingGroup(nameTop1);
    top1 -> CreateGroup("top1_1");
    top1 -> CreateGroup("top1_2");
    top1 -> CreateGroup("top1_3");
    top1 -> CreateGroup("top1_4");
    top1 -> CreateGroup("top1_5");
    auto top2 = new LKDrawingGroup(nameTop2);
    auto top3 = new LKDrawingGroup(nameTop3);
    auto top4 = new LKDrawingGroup(nameTop4);
    auto top5 = new LKDrawingGroup(nameTop5);

    //TH2D *hist1[4][4] = {0};
    TH2D *hist1[2][50] = {0};
    TH1D *hist2[4][4] = {0};
    TH2D *hist3[4][4] = {0};
    TH2D *hist4[4][4] = {0};
    TString asadTitle;
    TString agetTitle;

    LKSiliconMapping mapping;
    mapping.Load("mapping_ko2520_0904");

    for (auto id=0; id<mapping.GetNumDetectors(); ++id)
    {
        //auto id = group*10 + lcid;
        auto det = mapping.FindDetectorByIndex(id);
        if (det==nullptr) continue;
        int detNum = det->detNumber;
        TString detType = det->detType;  // "X6", "QQQ5", "BB10"
        TString histName1 = Form("hist1_%d",id);
        TString headerTitle = Form("[%d|%d,%s] M[%d,%d] ",id,detNum,detType.Data(),multRange1,multRange2);
        if (detType=="QQQ5") hist1[0][id] = new TH2D(histName1+"_J",headerTitle + "_J;strip;energy",32,0,32,200,0,emax);
        else hist1[0][id] = new TH2D(histName1+"_J",headerTitle + "_J;strip;energy",8,0,8,200,0,emax);
        hist1[1][id] = new TH2D(histName1+"_O",headerTitle + "_O;strip;energy",4,0,4,200,0,emax);
        LKDrawing *draw1;
             if (id>= 0 && id<12) { draw1 = top1 -> GetGroup(0) -> CreateDrawing(); draw1 -> Add(hist1[0][id]); }
        else if (id>=12 && id<24) { draw1 = top1 -> GetGroup(1) -> CreateDrawing(); draw1 -> Add(hist1[0][id]); }
        else if (id>=24 && id<36) { draw1 = top1 -> GetGroup(2) -> CreateDrawing(); draw1 -> Add(hist1[0][id]); }
        else if (id>=36 && id<48) { draw1 = top1 -> GetGroup(3) -> CreateDrawing(); draw1 -> Add(hist1[0][id]); }
             if (id>= 0 && id<12) { draw1 = top1 -> GetGroup(0) -> CreateDrawing(); draw1 -> Add(hist1[1][id]); }
        else if (id>=12 && id<24) { draw1 = top1 -> GetGroup(1) -> CreateDrawing(); draw1 -> Add(hist1[1][id]); }
        else if (id>=24 && id<36) { draw1 = top1 -> GetGroup(2) -> CreateDrawing(); draw1 -> Add(hist1[1][id]); }
        else if (id>=36 && id<48) { draw1 = top1 -> GetGroup(3) -> CreateDrawing(); draw1 -> Add(hist1[1][id]); }
        else if (detType=="X6")   for (auto i=0; i<8 ; ++i) { auto line = new TLine(8 *i,0,8 *i,4200); line -> SetLineColor(kBlue); draw1 -> Add(line); }
        else if (detType=="QQQ5") for (auto i=0; i<2 ; ++i) { auto line = new TLine(32*i,0,32*i,4200); line -> SetLineColor(kRed);  draw1 -> Add(line); }
    }

    for (auto asad=0; asad<4; ++asad)
    {
        if (asad==0) asadTitle = "QQQ5"; 
        else asadTitle = "X6"; 
        for (auto aget=0; aget<4; ++aget)
        {
            if (aget==0) agetTitle = "(Ohmic)"; 
            else agetTitle = "(Junction)"; 
            TString headerTitle = Form("[%04d] M[%d,%d] ",runID,multRange1,multRange2) + asadTitle + " " + agetTitle + " ";
            TString histName1 = Form("hist_%d_%d",asad,aget);
            TString histName2 = Form("histe_%d_%d",asad,aget);
            TString histName3 = Form("histet_%d_%d",asad,aget);
            TString histName4 = Form("histnp_%d_%d",asad,aget);
            TString histTitle1 = headerTitle + Form("AsAd=%d Aget=%d;strip;energy",asad,aget);
            TString histTitle2 = headerTitle + Form("AsAd=%d Aget=%d;energy",asad,aget);
            TString histTitle3 = headerTitle + Form("AsAd=%d Aget=%d;time;energy",asad,aget);
            TString histTitle4 = headerTitle + Form("AsAd=%d Aget=%d;channel;energy",asad,aget);

            //if      (is_ohmic(asad,aget)) hist1[asad][aget] = new TH2D(histName1,histTitle1,48,0,48,200,0,emax);
            //else if (is_qqq5 (asad,aget)) hist1[asad][aget] = new TH2D(histName1,histTitle1,64,0,64,200,0,emax);
            //else if (is_x6   (asad,aget)) hist1[asad][aget] = new TH2D(histName1,histTitle1,32,0,32,200,0,emax);
            hist2[asad][aget] = new TH1D(histName2,histTitle2,200,0,emax);           
            hist3[asad][aget] = new TH2D(histName3,histTitle3,100,0,500,200,0,4200); 
            hist4[asad][aget] = new TH2D(histName4,histTitle4,68,0,68,200,0,4200);   


            //if (asad==1&&aget==2)
            {
                //auto draw1 = top1 -> CreateDrawing(); draw1 -> Add(hist1[asad][aget]);
                auto draw2 = top2 -> CreateDrawing(); draw2 -> Add(hist2[asad][aget]);
                auto draw3 = top3 -> CreateDrawing(); draw3 -> Add(hist3[asad][aget]);
                auto draw4 = top4 -> CreateDrawing(); draw4 -> Add(hist4[asad][aget]);
            }
        }
    }

    auto draw5 = top5 -> CreateDrawing();
    auto hist5 = new TH2D("hist","",200,0,emax,200,0,emax);
    draw5 -> Add(hist5);

    double energySum = 0;
    int asadp = -1;
    int agetp = -1;
    int chanp = 0;

    TClonesArray *bRawDataArray = nullptr;
    TClonesArray *bSiChannelArray = nullptr;
    tree -> SetBranchAddress("RawData",&bRawDataArray);
    tree -> SetBranchAddress("SiChannel",&bSiChannelArray);

    auto numEvents = tree -> GetEntries();
    //numEvents = 10000;
    cout << numEvents << endl;

    for (auto event=0; event<numEvents; ++event)
    {
        tree -> GetEntry(event);
        auto numRawData = bRawDataArray -> GetEntries();
        auto numSiChannel = bSiChannelArray -> GetEntries();
        if (event%printAfter==0) cout << "event=" << event << " #channels=" << numRawData << endl;

        energySum = 0;
        asadp = -1;
        agetp = -1;
        chanp = -1;

        //if (!(numRawData>=multRange1 && numRawData<=multRange2)) continue;

        for (auto iRawData=0; iRawData<numRawData; ++iRawData)
        {
            auto rawData = (GETChannel*) bRawDataArray -> At(iRawData);
            auto asad0 = rawData -> GetAsad();
            auto aget0 = rawData -> GetAget();
            auto chan = rawData -> GetChan();
            auto energy = rawData -> GetEnergy();

            for (auto asad=0; asad<4; ++asad) {
                for (auto aget=0; aget<4; ++aget) {
                    if (asad0==asad && aget0==aget) {
                        hist4[asad][aget] -> Fill(chan,energy);
                        hist2[asad][aget] -> Fill(energy);
                    }
                }
            }
        }

        double q_energy1 = -1;
        double q_energy2 = -1;

        for (auto iSiChannel=0; iSiChannel<numSiChannel; ++iSiChannel)
        {
            auto siChannel = (LKSiChannel*) bSiChannelArray -> At(iSiChannel);
            auto id0 = siChannel -> GetDetID();
            if (event%printAfter==0)
                cout << iSiChannel << ": "
                    << " id=" << siChannel -> GetDetID()
                    << " side=" << siChannel -> GetSide()
                    << " ttl=" << siChannel -> GetTheta1()
                    << " phi=" << siChannel -> GetPhi0()
                    << " r=" << siChannel -> GetRadius()
                    << " asad=" << siChannel -> GetAsad()
                    << " aget=" << siChannel -> GetAget() << endl;
            auto side = siChannel -> GetSide();
            auto strip = siChannel -> GetStrip();
            auto energy1 = siChannel -> GetEnergy1();
            auto energySum = siChannel -> GetEnergySum();

            if (id0>=0 && id0<50) {
                if (side==0 && hist1[0][id0]!=nullptr) {
                    if (energy1>=0 && id0>=32 && id0<=35) hist1[0][id0] -> Fill(strip,energy1);
                    else if (energySum>=0) hist1[0][id0] -> Fill(strip,energySum);
                }
                if (side==1 && energy1>=0   && hist1[1][id0]!=nullptr) hist1[1][id0] -> Fill(strip,energy1);
            }

            auto detID = siChannel -> GetDetID();
            if (event%1000==0) {
                //if (detID==34 || detID==35) cout << ">> " << detID << endl;
                //else cout << detID << endl;
            }

            if (detID==34) { q_energy1 = energy1; } //cout << detID << " " << energySum << " " << energy1 << endl; }//cout << q_energy1 << endl; }
            if (detID==35) { q_energy2 = energy1; } //cout << detID << " " << energySum << " " << energy1 << endl; }//cout << q_energy2 << endl; }
        }
        if (q_energy1>0 && q_energy2>0) {
            //cout << q_energy1 << " " << q_energy2 << endl;
            hist5 -> Fill(q_energy1,q_energy2);
        }
    }

    TString outName = TString("data_reco/")+nameTop+".root";
    cout << outName << endl;
    auto out = new TFile(outName,"recreate");
    top1 -> Draw(); //top1 -> GetCanvas() -> Write();
    top1 -> GetGroup(0) -> GetCanvas() -> Write();
    top1 -> GetGroup(1) -> GetCanvas() -> Write();
    top1 -> GetGroup(2) -> GetCanvas() -> Write();
    top1 -> GetGroup(3) -> GetCanvas() -> Write();
    //top1 -> GetGroup(4) -> GetCanvas() -> Write();
    top2 -> Draw(); top2 -> GetCanvas() -> Write();
    //top3 -> Draw(); top3 -> GetCanvas() -> Write();
    top4 -> Draw(); top4 -> GetCanvas() -> Write();
    top5 -> Draw(); top5 -> GetCanvas() -> Write();
    out -> Close();
}
