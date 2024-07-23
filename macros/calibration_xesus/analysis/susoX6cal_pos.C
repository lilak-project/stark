#include "TSystem.h"
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TH2.h"
#include "TF1.h"
#include "TGraph.h"
#include "TString.h"
#include "TSpectrum.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <algorithm>


Double_t fpeaks(Double_t *x, Double_t *par) {
    for (Int_t p=0;p<4;p++) {
        Double_t norm  = par[3*p]; // "height" or "area"
        Double_t mean  = par[3*p+1];
        Double_t sigma = par[3*p+2];
#if defined(__PEAKS_C_FIT_AREAS__)
        norm /= sigma * (TMath::Sqrt(TMath::TwoPi())); // "area"
#endif /* defined(__PEAKS_C_FIT_AREAS__) */
        Double_t result = norm*TMath::Gaus(x[0],mean,sigma);
    }
    return result;
}


double interstrip[3]={0.25*75,0.50*75,0.75*75};

void susoX6cal_pos(){

    int i,j,k,l,n;

    time_t rawtime;
    struct tm * timeinfo;
    char buffer[80];

    time (&rawtime);
    timeinfo = localtime(&rawtime);
    strftime(buffer,sizeof(buffer),"%Y-%m-%d",timeinfo);
    std::string str(buffer);

    std::string calfile = "X6_position_" + str + ".cal";
    cout << calfile << endl;
    ofstream outfile1(calfile.c_str());

    //TFile* ifile = new TFile("../X6_cribdata_bal-det12ss4.root","read");
    TFile* ifile1 = new TFile("../../stark_0199.root","read");
    TFile* ifile2 = new TFile("../../stark_0253.root","read");

    TCanvas* c1[8];
    TCanvas* c2[8];
    for (int j=0; j<8; j++){
        c1[j] = new TCanvas(Form("X6 Calibration %d",j),Form("X6 Calibration %d",j),100,10*j,2200,1100);
        c1[j]->Divide(8,4);  

        c2[j] = new TCanvas(Form("X6 Calibration Fit %d",j),Form("X6 Calibration Fit %d",j),100+j*10,10,2200,1100);
        c2[j]->Divide(8,4);
    }

    TH1D* Posbackstrip[32][8][4];
    TH1D* Posback[32][8];
    TGraph* gr[32][8];
    //Parameters variables
    double pos[3];
    double offset[32][8];
    double slope[32][8];

    for (i=0; i<32; i++){
        for (j=0; j<8; j++){
            c1[i/4]->cd(8*(i%4)+j+1);
            gPad->SetGridy();
            gPad->SetGridx();

            //Posback[i][j]= new TH1D();
            for (k=0; k<4; k++){       
                cout << Form("Pos vs back det%d strip%d back%d",i+1,j+1,k+1) << endl;
                if(i<12) Posbackstrip[i][j][k] = (TH1D*)ifile1->Get(Form("Posvsbackstrip/Pos vs back det%d strip%d back%d",i+1,j+1,k+1));
                else	 Posbackstrip[i][j][k] = (TH1D*)ifile2->Get(Form("Posvsbackstrip/Pos vs back det%d strip%d back%d",i+1,j+1,k+1));

                if(k==0) Posback[i][j] = (TH1D*)Posbackstrip[i][j][k]->Clone();
                else	 Posback[i][j]->Add(Posbackstrip[i][j][k],pow(-1,k));
            }
            //Posbackstrip[i][j][0]->Draw("");
            if((i==20 && j==1) || (i==22 && j==7))Posback[i][j]->Rebin(2);
            Posback[i][j]->Draw("");
            Posback[i][j]->SetTitle(Form("Pos vs back det%d strip%d",i+1,j+1));
            Posback[i][j]->SetStats(0);

            //Find interstrip positions
            n=0;
            for (k=0; k<Posback[i][j]->GetNbinsX(); k++){
                if(i==0 || i==8){ //BACK SIGNAL MISSING - different conditions
                    if(Posback[i][j]->GetBinContent(k)*Posback[i][j]->GetBinContent(k+1)<0){
                        pos[n] = Posback[i][j]->GetBinLowEdge(k+1);
                        //cout << i <<"\t"<< j <<"\t"<< pos[n] << endl;
                        n++;
                    } else if(k<111 && Posback[i][j]->GetBinContent(k)>0 && Posback[i][j]->GetBinContent(k-1)>0 && Posback[i][j]->GetBinContent(k+1)==0 && Posback[i][j]->GetBinContent(k+2)==0){
                        pos[n] = Posback[i][j]->GetBinCenter(k);
                        //cout << i <<"\t"<< j <<"\t"<< k <<"\t"<< pos[n] << endl;
                        n++;
                    } else if(Posback[i][j]->GetBinContent(k-1)==0 && Posback[i][j]->GetBinContent(k)>0 && Posback[i][j]->GetBinContent(k+1)>0){
                        pos[n] = Posback[i][j]->GetBinCenter(k);
                        //cout << i <<"\t"<< j <<"\t"<< pos[n] << endl;
                        n++;
                    }
                    if(n>2) break;
                } else {
                    if(Posback[i][j]->GetBinContent(k)*Posback[i][j]->GetBinContent(k+1)<0){
                        pos[n] = Posback[i][j]->GetBinLowEdge(k+1);
                        cout << i <<"\t"<< j <<"\t"<< pos[n] << endl;
                        n++;
                    } else if(Posback[i][j]->GetBinContent(k)==0 && Posback[i][j]->GetBinContent(k-1)*Posback[i][j]->GetBinContent(k+1)<0){
                        pos[n] = Posback[i][j]->GetBinCenter(k);
                        cout << i <<"\t"<< j <<"\t"<< pos[n] << endl;
                        n++;
                    }
                    if(n>2) break;
                }       
            }	

            c2[i/4]->cd(8*(i%4)+j+1);
            gPad->SetGridy();
            gPad->SetGridx();
            gr[i][j]=new TGraph(3,pos,interstrip);
            gr[i][j]->SetMarkerStyle(20);   
            gr[i][j]->SetMarkerSize(0.9);
            gr[i][j]->Fit("pol1","Q");

            //Get fit function
            TF1 *fgr=(TF1*)gr[i][j]->GetFunction("pol1");
            fgr->SetLineWidth(1);   
            fgr->SetLineColor(4);  
            fgr->SetLineStyle(1); 
            gr[i][j]->Draw("AP");
            gr[i][j]->SetTitle(Form("Position calibration X6 det%d strip%d",i+1,j+1));
            gr[i][j]->GetYaxis()->SetRangeUser(0.,75.);
            gr[i][j]->GetXaxis()->SetLimits(0.,75.);
            gr[i][j]->GetXaxis()->SetNdivisions(504,kFALSE);
            gr[i][j]->GetYaxis()->SetNdivisions(504,kFALSE);
            gr[i][j]->GetXaxis()->SetTitle("Position (au)");  
            gr[i][j]->GetYaxis()->SetTitle("Position (mm)");
            gr[i][j]->GetYaxis()->SetTitleOffset(1.11);

            //c2->Modified();

            offset[i][j] = fgr->GetParameter(0);
            slope[i][j]  = fgr->GetParameter(1);

            cout <<"det:" << i <<"\tstrip: "<< j <<"\t"<< offset[i][j] << "\t" << slope[i][j] << endl;
            outfile1 << offset[i][j] << "     \t" << slope[i][j] << endl;
        }
    }

    for (int j=0; j<8; j++){
        c1[j]->SaveAs(Form("X6 position Calibration %d.png",j+1));
        c2[j]->SaveAs(Form("X6 position Calibration Fit %d.png",j+1));
    }
}
