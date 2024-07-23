#include "X6_stark.h"

void susoX6_stark(string filename="253", int FirstEntry=0,int LastEntry=0){

    cout << "Running susoX6_stark("<< FirstEntry <<","<< LastEntry <<")"<<endl;
    int i,j,k,l;
    int det, strip;
    double nrj;

    //TFile* ifile = new TFile("run_0260_X6_et.root","read");
    TFile* ifile = new TFile(Form("../data/stark_0%s.canv.root",filename.c_str()),"read");
    tr = (TTree*)ifile->Get("X6tree");
    treeaddress();

    TFile *newfile;//  TFile* ifile = new TFile(Form("%s.root",filename.c_str()),"read");string filename="run_0260_X6_et"
    newfile = new TFile(Form("stark_0%s.root",filename.c_str()),"recreate");
    starktree = new TTree("starktree","starktree");
    starktreedef();

    //Setting process limits ----------------------
    int N = tr->GetEntries();
    cout << "Entries: "<< N << endl;
    int stop = ::LastEntry(LastEntry);
    int init = ::FirstEntry(FirstEntry, stop);  
    cout << "Processing data....  from entry " << init << " to entry "<< stop << endl;

    //Initializing STARKMap, calibrations and histograms
    readmap();
    readcal();
    histos();

    /*gup[11][3]	= 0.0015614806;	//for E>3MeV
      gdwn[11][3]	= 0.0016073192; //for E>3MeV
      p0[11][3]	= 0.6598079439;
      p1[11][3]	=-0.6855971902;

      b0[11][3]	= 1.05563137234289273e+00;
      b1[11][3]	=-5.49029762003765132e-03;
      b2[11][3]	= 7.26722109436991925e-05;*/

    //cout << "i \tCobo[j]\tAsad[j]\tAget[j]\tchannel[j]\tenergy[j]" << endl;
    for(i=FirstEntry; i<stop; i++){ //(1)
        if(filename=="253" && i==277103) continue;	//file 253, entry with too big multiplicity
        tr->GetEntry(i);

        STARK_Mult	= 0;
        X6Mult_up	= 0;
        X6Mult_dw	= 0;
        X6Mult_back	= 0;
        for(j=0;j<X6Mult;j++){//(2)
            X6det[j]			= -1;
            X6strip[j]		= -1;

            if(stark_map[X6asad[j]][X6aget[j]][X6chan[j]]<0) continue;   
            else{//(3)
                X6det[j]	= stark_map[X6asad[j]][X6aget[j]][X6chan[j]]/20;
                X6strip[j]	= stark_map[X6asad[j]][X6aget[j]][X6chan[j]]%20;
            }//(3)
        }//(2)

        for(j=0;j<X6Mult;j++){//(2)  
                              //tr->Show(j);
            X6strip_back[j]		= 0;
            X6energy_up[j]		= 0;
            X6energy_up_cal[j]	= 0;
            X6energy_dw[j]		= 0;
            X6energy_dw_cal[j]	= 0;
            X6energy_back[j]		= 0;
            X6energy_back_cal[j]	= 0;
            X6flag[j]			= 0;

            //cout << i <<"\t"<< j <<"\t"<< X6det[j] <<"\t"<< X6strip[j];
            //cout << i <<"\t"<< j <<"\t\t"<< X6_det[j] <<"\t"<< X6det[j] <<"\t\t"<< X6_strip[j] <<"\t"<< X6strip[j] << endl;

            if(X6strip[j]<16 && (X6strip[j]%2)==0){ //(3)		//UPSTREAM
                                                    //cout << "\tFound upstream....   ";
                for(k=0;k<X6Mult;k++){ //(4)
                    if(X6flag[k]!=0) continue;
                    if(X6strip[k]<16 && (X6strip[k]%2)!=0){ //(5)	//DOWNSTREAM
                        if(X6det[j]==X6det[k] && X6strip[j]==(X6strip[k]-1)){ //(6)
                                                                              //cout << "Found matching downstream....   ";
                            for(l=0;l<X6Mult;l++){//(7)
                                if(X6det[j]==X6det[l] && X6strip[l]>15){//(8)	//BACK signal
                                                                        //cout << "Found matching back signal....   ";

                                    X6det2[STARK_Mult]		 = X6det[j];
                                    X6strip_front[STARK_Mult]	 = X6strip[j]/2;
                                    X6strip_back[STARK_Mult]	 = X6strip[l]-16;

                                    X6energy_up[STARK_Mult]	 = abs(X6nrj[j]);
                                    X6energy_dw[STARK_Mult]	 = abs(X6nrj[k]);
                                    X6energy_back[STARK_Mult]	 = X6nrj[l];

                                    X6energy_up_cal[STARK_Mult]	 = abs(X6nrj[j])*gup[X6det[j]][X6strip[j]/2];
                                    X6energy_dw_cal[STARK_Mult]	 = abs(X6nrj[k])*gdwn[X6det[k]][X6strip[k]/2];
                                    X6energy_back_cal[STARK_Mult]	 = aback[X6det[l]][X6strip_back[l]] + bback[X6det[l]][X6strip_back[l]]*X6nrj[l];

                                    //- - ENERGY - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                                    X6energy_raw[STARK_Mult]	= X6energy_up[STARK_Mult] + X6energy_dw[STARK_Mult];
                                    X6energy_cal[STARK_Mult]	= (X6energy_up_cal[STARK_Mult] + X6energy_dw_cal[STARK_Mult]);

                                    //- - POSITION - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                                    //X6position[STARK_Mult]	= (1+(X6energy_up[STARK_Mult] - X6energy_dw[STARK_Mult])/X6energy_raw[STARK_Mult])*75/2.;
                                    X6position[STARK_Mult]	= (1+(X6energy_up_cal[STARK_Mult] - X6energy_dw_cal[STARK_Mult])/X6energy_cal[STARK_Mult])*75/2.;
                                    X6position_cal[STARK_Mult]	= X6position[STARK_Mult];

                                    //- - BALLISTIC DEFICIT CORRECTION - - - - - - - - - - - - - - - - - - - - - - -
                                    X6energy[STARK_Mult]		= X6energy_cal[STARK_Mult]/(b0[X6det[j]][X6strip[j]/2]+b1[X6det[j]][X6strip[j]/2]*X6position[STARK_Mult]+b2[X6det[j]][X6strip[j]/2]*X6position[STARK_Mult]*X6position[STARK_Mult]);

                                    //cout << setprecision(8);
                                    //cout << i <<"\t"<< X6det[j] <<" "<< X6strip[j] <<"\t"<< b0[X6det[j]][X6strip[j]/2] <<"\t"<< b1[X6det[j]][X6strip[j]/2] <<"\t"<< b2[X6det[j]][X6strip[j]/2] <<"\t"<< X6energy_cal[j] <<"\t"<< X6energy[STARK_Mult] <<"\t"<< X6position[STARK_Mult] << endl;

                                    //- - HISTOGRAMS - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                                    Signal[X6det[j]][X6strip[j]]->Fill(X6energy_up[STARK_Mult]);
                                    Signal[X6det[k]][X6strip[k]]->Fill(X6energy_dw[STARK_Mult]);
                                    Signal[X6det[l]][X6strip[l]]->Fill(X6energy_back[STARK_Mult]);
                                    UpvsDown[X6det[j]][X6strip[j]/2]->Fill(X6energy_up[STARK_Mult],X6energy_dw[STARK_Mult]);

                                    Energy[X6det[j]][X6strip[l]-8]->Fill(X6energy_back_cal[STARK_Mult]);
                                    Eraw_backs->Fill(X6det[j]*4+X6strip[l]-16,X6energy_back[STARK_Mult]);
                                    Ecal_backs->Fill(X6det[j]*4+X6strip[l]-16,X6energy_back_cal[STARK_Mult]);

                                    Eraw_strip->Fill(X6det[j]*8+X6strip[j]/2,X6energy_raw[STARK_Mult]);
                                    Ecal_strip->Fill(X6det[j]*8+X6strip[j]/2,X6energy_cal[STARK_Mult]);
                                    Ebal_strip->Fill(X6det[j]*8+X6strip[j]/2,X6energy[STARK_Mult]);

                                    Energy[X6det[j]][X6strip[j]/2]->Fill(X6energy[STARK_Mult]);

                                    Evspos_cal[X6det[j]][X6strip[j]/2]->Fill(X6position[STARK_Mult],X6energy_cal[STARK_Mult]);
                                    Evspos_bal[X6det[j]][X6strip[j]/2]->Fill(X6position[STARK_Mult],X6energy[STARK_Mult]);
                                    Evspos_tot->Fill(X6position_cal[STARK_Mult],X6energy[STARK_Mult]);

                                    X6flag[j]++;
                                    X6flag[k]++;
                                    X6flag[l]++;
                                    STARK_Mult++;
                                    break;
                                } //(8)
                            } //(7)
                        } //(6)
                    } //(5)
                    if(X6flag[k]!=0) break;
                } //(4)
            } //(3)
        } //(2)

        //cout << endl;
        starktree->Fill();

        if((i+1)%1 == 0 || i%(N-1) == 0){ 
            cout << fixed <<"\rProcessing: Events: " << N << "\tDone: " << (i+1) << "  (%): " << setprecision(2) << (((float) (i+1))/N)*100.;
            cout.flush();
        }
    } //(1)    
    cout << endl;

    //ifile->Close();
    //TFile *newfile;//  TFile* ifile = new TFile(Form("%s.root",filename.c_str()),"read");string filename="run_0260_X6_et"
    newfile->cd();
    starktree->Write();

    Evspos_tot->Write();
    Eraw_strip->Write();
    Ecal_strip->Write();
    Ebal_strip->Write();
    Eraw_backs->Write();
    Ecal_backs->Write();

    newfile->mkdir("Evspos_bal");
    newfile->cd("Evspos_bal");
    for (int i=0; i<32; i++){
        for (int j=0; j<8; j++)	Evspos_bal[i][j]->Write();
    }
    newfile->mkdir("Evspos_cal");
    newfile->cd("Evspos_cal");
    for (int i=0; i<32; i++){
        for (int j=0; j<8; j++)	Evspos_cal[i][j]->Write();
    }
    newfile->mkdir("UpvsDown");
    newfile->cd("UpvsDown");
    for (int i=0; i<32; i++){
        for (int j=0; j<8; j++)	UpvsDown[i][j]->Write();
    }
    newfile->mkdir("Energy");
    newfile->cd("Energy");
    for (int i=0; i<32; i++){
        for (int j=0; j<12; j++)	Energy[i][j]->Write();
    }
    newfile->mkdir("Signals");
    newfile->cd("Signals");
    for (int i=0; i<32; i++){
        for (int j=0; j<20; j++)	Signal[i][j]->Write();
    }
    //for(int i=0;i<640;i++)    Test[i]->Write();
    //gROOT->ProcessLine(".q");
    //cout << "Output:   "<< Form("/home/cens-alpha-00/cydaq/data/X6_tamudata_run0260_%d.root",flag) <<endl;
    cout << newfile -> GetName() << endl;
}
