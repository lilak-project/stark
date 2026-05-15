void check_fit()
{
    double c0 = -1.90397;
    double c1 =  9.39922;

    auto anaMax = new LKChannelAnalyzer();
    anaMax -> Print();

    auto anaPulse = new LKChannelAnalyzer();
    anaPulse -> SetPulse("../../../common/pulseReference.root");
    anaPulse -> Print();

    auto file = new TFile("../data/stark_0013.conv.root");
    auto tree = (TTree*) file -> Get("event");
    TClonesArray* array = nullptr;
    tree -> SetBranchAddress("RawData",&array);
    auto numEvents = tree -> GetEntries();
    cout << "Number of events: " << numEvents << endl;

    if (1)
    {
        auto file_cal = new TFile("amplitude_slope.root","recreate");
        auto tree_cal = new TTree("hit","");
        double bSlope;
        double bAmplitude;
        tree_cal -> Branch("slope",&bSlope);
        tree_cal -> Branch("amplitude",&bAmplitude);

        auto cvsNormal = new TCanvas("cvsNormal","",1600,900);
        //int nCvs = 12; cvsNormal -> Divide(4,3);
        int nCvsNormal = 6;
        cvsNormal -> Divide(3,2);
        int iCvsNormal = 1; 

        numEvents = 1000;
        for (auto iEvent=0; iEvent<numEvents; ++iEvent)
        {
            bool breakFlag = false;
            tree -> GetEntry(iEvent);
            auto numChannels = array -> GetEntries();
            if (iEvent%100==0)
                cout << "Event-" << iEvent << ", Number of hits: " << numChannels << endl;

            for (auto iChannel=0; iChannel<numChannels; ++iChannel)
            {
                auto channel = (GETChannel *) array -> At(iChannel);
                if (channel->GetAget()==0) {
                    anaPulse -> SetDataIsInverted(false);
                    anaMax -> SetDataIsInverted(false);
                }
                else {
                    anaPulse -> SetDataIsInverted(true);
                    anaMax -> SetDataIsInverted(true);
                }
                auto wf = channel->GetWaveformY();

                bool continueFlag = false;
                for (auto t=0; t<512; ++t) {
                    if (wf[t]==4096) {
                        lk_debug << wf[t] << endl;
                        continueFlag = true;
                        break;
                    }
                }
                if (continueFlag)
                    continue;

                anaPulse -> Analyze(wf);
                auto hist = anaPulse -> NewHistBuffer();
                auto numRecoHits = anaPulse -> GetNumHits();
                if (numRecoHits!=1)
                    continue;

                double tbHit = anaPulse -> GetTbHit(0);
                double amplitude = anaPulse -> GetAmplitude(0);
                double pedestal = anaPulse -> GetPedestal();

                if (amplitude<500 || amplitude>3500)
                    continue;

                double t1 = tbHit-10; if (t1<0) t1 = 0;
                double t2 = tbHit+10; if (t2>512) t2 = 512;
                TF1 *fit = new TF1(Form("fit2_%d_%d",iEvent,iChannel),"[2]+(x>[0]&&x<([0]+6))*([1]/(([0]+6)-[0]))*(x-[0])",t1,t2);
                auto v1 = anaPulse -> Eval(tbHit, tbHit, amplitude);
                auto v2 = anaPulse -> Eval(tbHit+6, tbHit, amplitude);
                bSlope = (v2-v1)/6.;
                bAmplitude = amplitude;
                tree_cal -> Fill();
                fit -> SetParameters(tbHit,v2,pedestal);
                hist -> Fit(fit,"QN0");
                fit -> SetLineColor(kBlue);

                if (iCvsNormal<nCvsNormal) {
                    cvsNormal -> cd(iCvsNormal++);
                    anaPulse -> Draw("new");
                    anaPulse -> GetHistBuffer() -> GetXaxis() -> SetRangeUser(tbHit-50,tbHit+50);
                    fit -> Draw("samel");
                }
                //else breakFlag = true;
            }
            if (breakFlag)
                break;
        }

        file_cal -> cd();
        tree_cal -> Write();
        new TCanvas();
        auto hist = new TH2D("hist","",150,0,400,150,0,4000);
        tree_cal -> Draw("amplitude:slope>>hist","","colz");
        auto f1 = new TF1("f1","pol1",0,400);
        hist -> Fit("f1");
        c0 = f1 -> GetParameter(0);
        c1 = f1 -> GetParameter(1);
        cout << c0 << " " << c1 << endl;
        cout << file_cal -> GetName() << endl;
    }

    auto cvs = new TCanvas("cvs","",1600,900);
    int nCvs = 6;
    cvs -> Divide(3,2);
    int iCvs = 1; 

    for (auto iEvent=0; iEvent<numEvents; ++iEvent)
    {
        tree -> GetEntry(iEvent);
        auto numChannels = array -> GetEntries();
        if (iEvent%10000==0) cout << "Event-" << iEvent << ", Number of hits: " << numChannels << endl;
        for (auto iChannel=0; iChannel<numChannels; ++iChannel)
        {
            if (iCvs==nCvs)
                return;

            auto channel = (GETChannel *) array -> At(iChannel);
            if (channel->GetAget()==0) {
                anaPulse -> SetDataIsInverted(false);
                anaMax -> SetDataIsInverted(false);
            }
            else {
                anaPulse -> SetDataIsInverted(true);
                anaMax -> SetDataIsInverted(true);
            }

            auto wf = channel->GetWaveformY();
            bool continueFlag = true;
            for (auto t=0; t<512; ++t) {
                if (wf[t]==4095) {
                    continueFlag = false;
                    break;
                }
            }
            if (continueFlag)
                continue;

            cout << iEvent << " " << iChannel << endl;
            anaMax -> Analyze(wf);
            cvs -> cd(iCvs++);
            anaMax -> Draw("new");

            auto hist = anaMax -> GetHistBuffer();
            auto numRecoHits = anaMax -> GetNumHits();
            if (numRecoHits>=1)
            {
                double tbHit = anaMax -> GetTbHit(0);
                double amplitude = anaMax -> GetAmplitude(0);
                double pedestal = anaMax -> GetPedestal();

                hist -> GetXaxis() -> SetRangeUser(tbHit-50,tbHit+50);

                double t1 = tbHit-10; if (t1<0) t1 = 0;
                double t2 = tbHit+5; if (t2>512) t2 = 512;
                TF1 *fit = new TF1(Form("fit_%d_%d",iEvent,iChannel),"[3]+(x>[0]&&x<[1])*([2]/([1]-[0]))*(x-[0])+(x>[1]&&x<[1]+4)*[2]",t1,t2);
                fit -> SetParameters(tbHit-5,tbHit,amplitude,pedestal);
                hist -> Fit(fit,"QN0");
                double a0 = fit -> GetParameter(0);
                double a1 = fit -> GetParameter(1);
                double a2 = fit -> GetParameter(2);
                double slope = a2 / (a1-a0);
                double reco_amp = c0 + c1*slope;
                cout << a0 << " " << a1 << " " << a2 << " a=" << amplitude << " s=" << slope << " r=" << reco_amp << endl;
                fit -> SetLineColor(kBlue);
                fit -> Draw("samel");
            }
        }
    }
}
