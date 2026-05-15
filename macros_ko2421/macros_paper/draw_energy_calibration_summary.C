TCanvas *cvs = nullptr;

void draw_energy_calibration_summary()
{
    cout << "execute using --web option to save as accurate pdf file" << endl;
    cout << "use SaveCanvas() to save eps, png, pdf" << endl;

    //gStyle -> SetPalette(kRainBow);

    bool pull_channels = true;
    bool set_custom_label = true;

    int numJ = 8;
    int numO = 4;
    double min = 1.05;
    double max = 6.50;
    const int oooohmic = 0;
    const int junction = 1;

    auto file = new TFile("energy_calibration_summary.root");
    auto c1 = (TCanvas*) file -> Get("c1"); c1 -> Draw();
    auto list = c1 -> GetListOfPrimitives();

    TH2D* hist_old[2];
    hist_old[junction] = (TH2D*) (((TPad*) list -> At(0)) -> GetListOfPrimitives() -> FindObject("h0"));
    hist_old[oooohmic] = (TH2D*) (((TPad*) list -> At(1)) -> GetListOfPrimitives() -> FindObject("h1"));

    TObjArray *lines[2];
    lines[oooohmic] = new TObjArray();
    lines[junction] = new TObjArray();

    for (auto jo : {oooohmic,junction})
    {
        auto hist = hist_old[jo];
        int numStrips = numO;
        if (jo==junction) numStrips = numJ;

        hist -> GetYaxis() -> SetNdivisions(505);
        hist -> GetXaxis() -> SetTitleSize(0.055);
        hist -> GetYaxis() -> SetTitleSize(0.055);
        //hist -> GetXaxis() -> SetLabelSize(0.070);
        hist -> GetXaxis() -> SetLabelOffset(0.01);
        hist -> GetYaxis() -> SetLabelOffset(0.005);
        hist -> GetXaxis() -> SetLabelSize(0.065);
        hist -> GetYaxis() -> SetLabelSize(0.055);
        hist -> GetXaxis() -> SetTitleOffset(1.2);
        hist -> GetYaxis() -> SetTitleOffset(1.3);
        hist -> GetXaxis() -> SetTickSize(0.00001);
        hist -> GetYaxis() -> SetTickSize(.01);
        hist -> GetYaxis() -> SetRangeUser(min,max);
        if (pull_channels) hist -> GetXaxis() -> SetRangeUser(0,37*numStrips);
        hist -> SetStats(0);
        if (jo==junction) hist -> SetTitle(";;Junction Energy (Mev)");
        if (jo==oooohmic) hist -> SetTitle(";;Ohmic Energy (Mev)");
        hist -> GetYaxis() -> SetTitleOffset(0.4);
        //hist -> GetYaxis() -> CenterTitle();

        if (set_custom_label) for (auto bin=1; bin<=40*numStrips; ++bin) hist -> GetXaxis() -> SetBinLabel(bin,""); 
    }

    TH2D* hist_new[2];
    hist_new[oooohmic] = (TH2D*) hist_old[oooohmic] -> Clone("hist_o");
    hist_new[junction] = (TH2D*) hist_old[junction] -> Clone("hist_j");
    hist_new[oooohmic] -> Reset("ICES");
    hist_new[junction] -> Reset("ICES");

    for (auto jo : {oooohmic,junction})
    {
        int numStrips = numO;
        if (jo==junction) numStrips = numJ;
        for (auto det=0; det<12; ++det)
        {
            //int det0 = det+16;
            int det0 = det;
            int binx1 = det*numStrips+1;
            int binx2 = (det+1)*numStrips;
            int binx3 = det0*numStrips+1;
            if (pull_channels && det0>=21) { binx3 = binx3 - numStrips; }
            if (pull_channels && det0>=36) { binx3 = binx3 - 2*numStrips; }
            int numBinsY = hist_old[jo] -> GetYaxis() -> GetNbins();
            for (auto biny=0; biny<numBinsY; ++biny) {
                for (auto strip=0; strip<numStrips; ++strip) {
                    auto content = hist_old[jo] -> GetBinContent(binx1+strip,biny);
                    hist_new[jo] -> SetBinContent(binx3+strip,biny,content);
                }
            }
        }
        for (auto det=12; det<28; ++det)
        {
            //int det0 = det-12;
            int det0 = det+12;
            int binx1 = det*numStrips+1;
            int binx2 = (det+1)*numStrips;
            int binx3 = det0*numStrips+1;
            if (pull_channels && det0>=21) { binx3 = binx3 - numStrips; }
            if (pull_channels && det0>=36) { binx3 = binx3 - 2*numStrips; }
            int numBinsY = hist_old[jo] -> GetYaxis() -> GetNbins();
            for (auto biny=0; biny<numBinsY; ++biny) {
                for (auto strip=0; strip<numStrips; ++strip) {
                    auto content = hist_old[jo] -> GetBinContent(binx1+strip,biny);
                    hist_new[jo] -> SetBinContent(binx3+strip,biny,content);
                    //cout << binx3+strip << " " << biny << endl;
                }
            }
        }
        for (auto det=28; det<40; ++det)
        {
            int det0 = det-16;
            int binx1 = det*numStrips+1;
            int binx2 = (det+1)*numStrips;
            int binx3 = det0*numStrips+1;
            if (pull_channels && det0>=21) { binx3 = binx3 - numStrips; }
            if (pull_channels && det0>=36) { binx3 = binx3 - 2*numStrips; }
            int numBinsY = hist_old[jo] -> GetYaxis() -> GetNbins();
            for (auto biny=0; biny<numBinsY; ++biny) {
                for (auto strip=0; strip<numStrips; ++strip) {
                    auto content = hist_old[jo] -> GetBinContent(binx1+strip,biny);
                    hist_new[jo] -> SetBinContent(binx3+strip,biny,content);
                    //cout << binx3+strip << " " << biny << endl;
                }
            }
        }
    }

    for (auto jo : {0,1})
    {
        auto hist = hist_new[jo];
        int numStrips = numO;
        if (jo==junction) numStrips = numJ;

        for (auto det=0; det<40; ++det)
        {
            int binx1 = det*numStrips+1;
            int binx2 = (det+1)*numStrips;
            double x1 = hist -> GetXaxis() -> GetBinLowEdge(binx1);
            double x2 = hist -> GetXaxis() -> GetBinUpEdge (binx2);
            if (det!=40) {
                auto line = new TLine(x2,min,x2,max);
                line -> SetLineStyle(2);
                //line -> SetLineColor(kGray+2);
                lines[jo] -> Add(line);
            }
            int ringIdx, detIdx;
            TString comment;
            int det0 = det;
            if (pull_channels && det0>=20) { det0 = det0 + 1; }
            if (pull_channels && det0>=34) { det0 = det0 + 2; }
            if (det0>=0  && det0<12) { ringIdx = 2; detIdx = det0+1; }
            if (det0>=12 && det0<24) { ringIdx = 3; detIdx = det0-12+1; }
            if (det0>=24 && det0<40) { ringIdx = 1; detIdx = det0-24+1; }
            if (det0>32) comment = "(CSD)";
            TString label = Form("%d - %d",ringIdx,detIdx);
            if (set_custom_label) hist -> GetXaxis() -> SetBinLabel(floor(0.5*(binx1+binx2))+1,label); 
            //if (set_custom_label) hist -> GetXaxis() -> SetBinLabel(binx1,label); 
            //if (set_custom_label) hist -> GetXaxis() -> SetLabelSize(0.5);
            if (set_custom_label) hist -> GetXaxis() -> LabelsOption("v");

            if (jo==oooohmic)
            {
                if ((!pull_channels&&(det>=16&&det<24)) || (pull_channels&&(det>=16&&det<23)))
                {
                    int binx = binx1;
                    int numBinsY = hist -> GetYaxis() -> GetNbins();
                    for (auto biny=0; biny<numBinsY; ++biny)
                    {
                        auto content = hist -> GetBinContent(binx,biny);
                        for (auto strip=0; strip<numStrips; ++strip)
                        {
                            hist -> SetBinContent(binx+strip,biny,content);
                        }
                    }
                }
            }
            if (!pull_channels && det==20)
            {
                int binx = binx1;
                int numBinsY = hist -> GetYaxis() -> GetNbins();
                for (auto biny=0; biny<numBinsY; ++biny)
                {
                    for (auto strip=0; strip<numStrips; ++strip)
                    {
                        hist -> SetBinContent(binx+strip,biny,0);
                    }
                }
            }
        }
    }

    for (auto oldnew : {1})
    {
        cvs = new TCanvas(Form("cvs_%d",oldnew),Form("cvs_%d",oldnew),1400,800);
        cvs -> Divide(1,2,0,0);
        cvs -> cd(1) -> SetMargin(0.05, 0.01, 0.0, 0.15);
        cvs -> cd(2) -> SetMargin(0.05, 0.01, 0.15, 0.0);
        for (auto jo : {oooohmic,junction})
        {
            if (jo==oooohmic) cvs -> cd(1) -> SetLogz();
            if (jo==junction) cvs -> cd(2) -> SetLogz();
            if (oldnew==0) hist_old[jo] -> Draw("col");
            if (oldnew==1) hist_new[jo] -> Draw("col");
            lines[jo] -> Draw("samel");
        }
    }
}

void SaveCanvas()
{
    cvs -> SaveAs("elark_energy_calibration_summary.eps");
    cvs -> SaveAs("elark_energy_calibration_summary.png");
    cvs -> SaveAs("elark_energy_calibration_summary.pdf");
}
