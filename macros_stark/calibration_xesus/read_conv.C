#include "read_conv.h"

void read_conv(int RunN=253){
    ofstream txtfile("dump_data.txt");

    auto fChannelAnalyzer = new LKChannelAnalyzer();
    auto fStarkPlane = new SKSiArrayPlane();
    fStarkPlane -> AddPar("../config_stark.mac");
    fStarkPlane -> Init();
    auto siChannel1 = new LKSiChannel();

    auto file = new TFile(Form("../data/stark_0%d.conv.root", RunN));
    auto tree = (TTree *) file->Get("event");

    TFile* ofile = new TFile(Form("../data/stark_0%d.canv.root",RunN),"recreate");
    X6tree = new TTree("X6tree","X6tree");
    treedef();

    TClonesArray *channelArray = nullptr;
    tree->SetBranchAddress("RawData", &channelArray);

    //Histogram definition in read_conv.h
    def();

    //reading starkmap
    readmap();

    auto numEvents = tree->GetEntries();
    cout << "Nentries: "<< tree->GetEntries() << endl;
    //numEvents = 160000;
    for (iEvent=0; iEvent<numEvents; ++iEvent){
        for(int i=0;i<36;i++){
            for(int j=0;j<4;j++){
                //X6_strip_ohmic[i][j]=0;
                E_ohmic[i][j]=0;
            }
            for(int j=0;j<8;j++){
                //X6_det[i][j]=0;
                //X6_strip[i][j]=0;
                //X6_position[i][j]=0;
                E_junction[i][j][0]=0;
                E_junction[i][j][1]=0;
            }
        }
        for(int i=0;i<12;i++){
            CSD_ohmic[i]=0;
            for(int j=0;j<8;j++){
                CSD_junction[i][j]=0;
            }
        }
        tree->GetEntry(iEvent); // channelArray update for event number = iEvent

        auto numChannels = channelArray->GetEntries();
        X6Mult = numChannels;
        //cout << iEvent << "\tmultiplicity" << numChannels << endl;
        for (auto iChannel=0; iChannel<numChannels; ++iChannel)
        {
            auto channel = (GETChannel*) channelArray->At(iChannel);
            //channel -> Print(); return;
            //channel -> Draw(); return;
            auto cobo   = channel->GetCobo();
            auto asad   = channel->GetAsad();
            auto aget   = channel->GetAget();
            auto chan   = channel->GetChan();
            auto energy = channel->GetEnergy();
            auto time   = channel->GetTime();
            auto pedestal = channel->GetPedestal();
            int* buffer   = channel->GetWaveformY();

            {
                bool good = fStarkPlane -> SetSiChannelData(siChannel1, channel);
                if (good) {
                    if (siChannel1->GetSide()==0) fChannelAnalyzer -> SetDataIsInverted(true);
                    else fChannelAnalyzer -> SetDataIsInverted(false);
                    auto data = channel -> GetWaveformY();
                    fChannelAnalyzer -> Analyze(data);
                    auto numRecoHits = fChannelAnalyzer -> GetNumHits();
                    if (numRecoHits>=1)
                        energy = fChannelAnalyzer -> GetAmplitude(0);
                    else
                        energy = 0;
                }
                else {
                    energy = 0;
                }
            }

            X6cobo[iChannel] = cobo;
            X6asad[iChannel] = asad;
            X6aget[iChannel] = aget;
            X6chan[iChannel] = chan;
            X6nrj[iChannel]  = energy;

            if(stark_map[asad][aget][chan]<0) continue;
            else{
                //cout << i <<"\t"<< Cobo[j] <<"\t"<< Asad[j] <<"\t"<< Aget[j] <<"\t"<< channel[j]<<"\t\t"<< energy[j] << endl;
                X6_det[iChannel]   = stark_map[asad][aget][chan]/20;
                X6_strip[iChannel] = (stark_map[asad][aget][chan]%20);
            }

            //cout << iEvent <<" "<< iChannel <<" "<< X6Mult <<" "<< cobo <<" "<< asad <<" "<< aget <<" "<< chan <<" | "<< time <<" "<< energy << endl;

            hitpattern->Fill(asad*500+aget*100+chan,abs(energy));
            //if(asad==2) cout << iEvent << " " << iChannel << " " << numChannels << " " << cobo << " " << asad << " " << aget << " " << chan << " " << time << " " << energy << " " << pedestal << endl;
            if(chan==11 || chan==22 || chan==45 || chan==56) continue;
            if(chan>11 && chan<22) chan-=1;
            if(chan>22 && chan<45) chan-=2;
            if(chan>45 && chan<56) chan-=3;
            if(chan>56 && chan<68) chan-=4;

            if(aget==0 && chan<48){//ohmic
                det=0;
                det_back = 0;
                strip_back = 0;
                if(chan>=0&&chan<16) det = asad*12 + int(chan/4) + 8;
                if(chan>=16&&chan<32) det = asad*12 + int((chan-16)/4) + 4;
                if(chan>=32&&chan<48) det = asad*12 + int((chan-32)/4);
                strip = chan%4;
                det_back = det;
                strip_back = strip;
                //cout << iEvent <<" "<< asad << " " << chan << " " << det << " " << strip << endl;
                E_ohmic[det][strip]=abs(energy);
                //cout << iEvent <<" "<< iChannel <<" "<< X6Mult <<" | "<< cobo <<" "<< asad <<" "<< aget <<" "<< chan <<" "<< energy <<" | "<< det <<" "<< strip <<" "<< CSD_junction[det][strip] <<" "<< CSD_ohmic[det] << endl;
            }else{//junction
                det = asad*12+(aget-1)*4+int(chan/16);
                strip = int(chan/2) - int(chan/16)*8;
                if(chan%2==0) E_junction[det][strip][0]=abs(energy);
                if(chan%2==1) E_junction[det][strip][1]=abs(energy);
                //cout << iEvent <<" "<< iChannel <<" "<< X6Mult <<" | "<< cobo <<" "<< asad <<" "<< aget <<" "<< chan <<" "<< energy <<" | "<< det <<" "<< strip <<" "<< CSD_junction[det][strip] <<" "<< CSD_ohmic[det] << endl;
            }
            if(asad==2){
                det = 0;
                strip = 0;
                if(aget==0){//ohmic 
                    if(chan>=0&&chan<16) det = int(chan/4) + 8;
                    if(chan>=16&&chan<32) det = int((chan-16)/4) + 4;
                    if(chan>=32&&chan<48) det = int((chan-32)/4);
                    CSD_ohmic[det]=abs(energy);
                    if(chan%4==2 && CSD_ohmic[det]>0){
                        hCSD_ohmic[det]->Fill(CSD_ohmic[det]);
                    }
                    //cout << iEvent <<" "<< iChannel <<" "<< X6Mult <<" | "<< cobo <<" "<< asad <<" "<< aget <<" "<< chan <<" "<< energy <<" | "<< det <<" "<< strip <<" "<< CSD_junction[det][strip] <<" "<< CSD_ohmic[det] << endl;
                }
                if(aget==1){//CSD junction 
                    det = int(chan/8);
                    strip = chan - det*8;
                    CSD_junction[det][strip]=abs(energy);
                    if(CSD_junction[det][strip]>0){
                        hCSD_junction[det][strip]->Fill(CSD_junction[det][strip]);
                    }
                    //cout << iEvent <<" "<< iChannel <<" "<< X6Mult <<" | "<< cobo <<" "<< asad <<" "<< aget <<" "<< chan <<" "<< energy <<" | "<< det <<" "<< strip <<" "<< CSD_junction[det][strip] <<" "<< CSD_ohmic[det] << endl;
                }
                //cout << iEvent <<" "<< iChannel <<" "<< X6Mult <<" | "<< cobo <<" "<< asad <<" "<< aget <<" "<< chan <<" "<< energy <<" | "<< det <<" "<< strip <<" "<< CSD_junction[det][strip] <<" "<< CSD_ohmic[det] << endl;
                //cout << iEvent << " " << X6Mult << " " << det << " " << strip << " " << CSD_junction[det][strip] << " " << CSD_ohmic[det] << endl;
            }
            //cout << iEvent <<" "<< iChannel <<" "<< X6Mult <<" | "<< cobo <<" "<< asad <<" "<< aget <<" "<< chan <<" "<< energy <<" | "<< det <<" "<< strip <<" "<< CSD_junction[det][strip] <<" "<< CSD_ohmic[det] << endl;
            //txtfile << iEvent << " " << iChannel << " " << numChannels << " " << cobo << " " << asad << " " << aget << " " << chan << " " << time << " " << energy << " " << pedestal << endl;
        }

        for(int i=0;i<36;i++){
            for(int j=0;j<8;j++){
                //if(E_ohmic[i][j]>0 && E_junction[i][j][0]>0 && E_junction[i][j][1]>0)
                if(E_junction[i][j][0]>0 && E_junction[i][j][1]>0)
                {
                    double pos = ((E_junction[i][j][0]-E_junction[i][j][1])*1000)/(E_junction[i][j][0]+E_junction[i][j][1]);
                    auto esum = E_junction[i][j][0]+E_junction[i][j][1];
                    LvsR[i][j]->Fill(E_junction[i][j][0],E_junction[i][j][1]);
                    EvsPos[i][j]->Fill(pos/1000,esum);
                    StripvsPos[i][j]->Fill(pos/1000,j+1);
                }
            }
            for(int j=0;j<4;j++){
                if(E_ohmic[i][j]>0){
                    E_ohmic_raw[i][j]->Fill(E_ohmic[i][j]);
                }
            }
        }

        X6tree->Fill();
        if(iEvent%10000==0) cout << iEvent << endl;
    }

    gStyle->SetPalette(kCool);
    TCanvas* chitpattern=new TCanvas("chitpattern","chitpattern",5,5,1000,600);
    //TCanvas* chitpattern = LKPainter::GetPainter() -> CanvasFull("chitpattern");
    chitpattern->SetLogz();
    hitpattern->Draw("colz");

    for(int i=0;i<36;i++){
        cX6[i]=new TCanvas(Form("cX6_%d",i),Form("cX6_%d",i),5,5,1000,600);
        //cX6[i]= LKPainter::GetPainter() -> CanvasFull(Form("cX6_%d",i));
        cX6[i]->Divide(4,5);
        for(int j=0;j<4;j++){
            cX6[i]->cd(j+1);
            E_ohmic_raw[i][j]->GetXaxis()->SetRangeUser(200,700);
            E_ohmic_raw[i][j]->Draw();
        }
        for(int j=0;j<8;j++){
            cX6[i]->cd(j*2+5);
            LvsR[i][j]->Draw("colz");
            cX6[i]->cd(j*2+6);
            EvsPos[i][j]->Draw("colz");
        }
    }

    for(int i=0;i<12;i++){
        cCSD[i]=new TCanvas(Form("cCSD_%d",i),Form("cCSD_%d",i),5,5,1000,600);
        //cCSD[i]= LKPainter::GetPainter() -> CanvasFull(Form("cCSD_%d",i));
        cCSD[i]->Divide(3,3);
        cCSD[i]->cd(1);
        hCSD_ohmic[i]->Draw();
        for(int j=0;j<8;j++){
            cCSD[i]->cd(j+2);
            hCSD_junction[i][j]->Draw();
        }
    }

    X6tree->Write();
    hitpattern->Write();
    for(int i=0;i<36;i++) cX6[i]->Write();
    for(int i=0;i<36;i++) for(int j=0;j<8;j++) EvsPos[i][j]->Write();
    for(int i=0;i<36;i++) for(int j=0;j<4;j++) E_ohmic_raw[i][j]->Write();
    for(int i=0;i<12;i++) for(int j=0;j<8;j++) hCSD_junction[i][j]->Write();
    ofile->Close();
    ofile->ls();
}
