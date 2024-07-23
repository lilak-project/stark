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


//Double_t fpeaks(Double_t *x, Double_t *par) {
//   for (Int_t p=0;p<4;p++) {
//      Double_t norm  = par[3*p]; // "height" or "area"
//      Double_t mean  = par[3*p+1];
//      Double_t sigma = par[3*p+2];
//#if defined(__PEAKS_C_FIT_AREAS__)
//      norm /= sigma * (TMath::Sqrt(TMath::TwoPi())); // "area"
//#endif /* defined(__PEAKS_C_FIT_AREAS__) */
//      Double_t result = norm*TMath::Gaus(x[0],mean,sigma);
//   }
//   return result;
//}


double en[2]={3.182,5.486};
//double en[3]={3.0744,5.400496,5.7194901};

void susoX6cal_back(){

    time_t rawtime;
    struct tm * timeinfo;
    char buffer[80];

    time (&rawtime);
    timeinfo = localtime(&rawtime);
    strftime(buffer,sizeof(buffer),"%Y-%m-%d",timeinfo);
    std::string str(buffer);

    std::string calfile = "X6_back_" + str + ".cal";
    cout << calfile << endl;
    ofstream outfile1(calfile.c_str());

    //TFile* ifile = new TFile("../X6_cribdata_bal-det12ss4.root","read");
    TFile* ifile1 = new TFile("../stark_0199.root","read");
    TFile* ifile2 = new TFile("../stark_0253.root","read");

    TCanvas* c1[8];
    TCanvas* c2[8];
    for (int j=0; j<8; j++){
        c1[j] = new TCanvas(Form("X6 Calibration %d",j),Form("X6 Calibration %d",j),300,10*j,1600,1100);
        c1[j]->Divide(4,4);  

        c2[j] = new TCanvas(Form("X6 Calibration Fit %d",j),Form("X6 Calibration Fit %d",j),300+j*10,10,1600,1100);
        c2[j]->Divide(4,4);
    }

    TH1D* Eback[32][4];
    TGraph* gr[32][4];
    //Parameters variables
    double chback[4];
    double offset[32][4];
    double slope[32][4];

    //for (int i=0; i<32; i++){
    for (int i=12; i<32; i++){
        for (int j=0; j<4; j++){
            c1[i/4]->cd(4*(i%4)+j+1);
            gPad->SetGridy();
            gPad->SetGridx();

            cout << Form("Back signal det%d strip%d",i+1,j+9);// << endl;
            if(i<12)	Eback[i][j]=(TH1D*)ifile1->Get(Form("Signals/Back signal det%d strip%d",i+1,j+9));
            else	Eback[i][j]=(TH1D*)ifile2->Get(Form("Signals/Back signal det%d strip%d",i+1,j+9));
            Eback[i][j]->Rebin(2);
            Eback[i][j]->Draw("hist");

            if((i==0 && j==3)||(i==8 && j==3)){
                offset[i][j] = 0.000;
                slope[i][j]  = 0.000;

                outfile1 << setprecision(6) << offset[i][j] << "\t" << slope[i][j] << endl;
                continue;
            }

            //Find peaks
            Int_t npeaks=2;
            Float_t fit_range;
            TSpectrum *s = new TSpectrum(npeaks);
            Int_t nfound = s->Search(Eback[i][j],3,"new",0.05);
            //Eback[i][j]->Draw("hist");
            cout <<"\t\tnpeaks:  "<< npeaks <<"\t" << "nfound:  "<< nfound << endl;

            //Order peaks 
            Double_t* xpeakr = s->GetPositionX();
            sort(xpeakr,xpeakr+nfound);      
            //for (Int_t p=0;p<nfound;p++)	cout << p <<"\t"<< s->GetPositionX()[p] <<"\t"<< xpeakr[p] << endl;		

            //fit 4-gaussians
            TF1 *g[3];
            Double_t param[12];
            for (Int_t p=0;p<nfound;p++){
                g[p] = new TF1(Form("g%d",p),"gaus",xpeakr[p]*0.96,xpeakr[p]*1.05);
                /*if(p==0)  g[p] = new TF1(Form("g%d",p),"gaus",xpeakr[p]*0.96,xpeakr[p]*1.05);
                  if(p==1)  g[p] = new TF1(Form("g%d",p),"gaus",xpeakr[p]*0.95,xpeakr[p]*1.02);
                  if(p==2)  g[p] = new TF1(Form("g%d",p),"gaus",xpeakr[p]*0.98,xpeakr[p]*1.05);*/
                g[p]->SetParLimits(0,5,Eback[i][j]->GetMaximum());		//fix amplitude to positive values
                g[p]->SetParLimits(1,xpeakr[p]*0.99,xpeakr[p]*1.01);		//fix mean around peak
                g[p]->SetParLimits(2,2,7);					//fix sigma to positive values
                Eback[i][j]->Fit(g[p],"QR+");  

                chback[npeaks-nfound+p] = g[p]->GetParameter(1);
                cout << p <<"\t"<< s->GetPositionX()[p] <<"\t"<< xpeakr[p] <<"\t"<< chback[npeaks-nfound+p] << endl;
            }

            c2[i/4]->cd(4*(i%4)+j+1);
            gPad->SetGridy();
            gPad->SetGridx();
            gr[i][j]=new TGraph(npeaks,chback,en);
            //if(nfound<4)	gr[i][j]->RemovePoint(0);
            //if(nfound<3)	gr[i][j]->RemovePoint(0);
            gr[i][j]->SetMarkerStyle(20);   
            gr[i][j]->SetMarkerSize(0.9);
            gr[i][j]->Fit("pol1","Q");

            //Get fit function
            TF1 *fgr=(TF1*)gr[i][j]->GetFunction("pol1");
            fgr->SetLineWidth(1);   
            fgr->SetLineColor(4);  
            fgr->SetLineStyle(1); 
            gr[i][j]->Draw("AP");
            gr[i][j]->SetTitle(Form("Calibration X6 back det%d strip%d",i+1,j+1));
            gr[i][j]->GetYaxis()->SetRangeUser(3.0,6.);
            gr[i][j]->GetXaxis()->SetLimits(0,1250);
            //gr[i][j]->GetXaxis()->SetNdivisions(503,kFALSE);
            gr[i][j]->GetXaxis()->SetTitle("Energy (ch)");  
            gr[i][j]->GetYaxis()->SetTitle("Energy (MeV)");
            //gr[i][j]->GetYaxis()->SetTitleOffset(1.44);

            //c2->Modified();

            offset[i][j] = fgr->GetParameter(0);
            slope[i][j]  = fgr->GetParameter(1);

            cout <<"det:" << i <<"\tstrip: "<< j <<"\t"<< offset[i][j] << "\t" << slope[i][j] << endl;
            outfile1 << offset[i][j] << "\t" << slope[i][j] << endl;
        }
    }

    for (int j=0; j<8; j++){
        c1[j]->SaveAs(Form("X6 back Calibration %d.png",j+1));
        c2[j]->SaveAs(Form("X6 back Calibration Fit %d.png",j+1));
    }
}
