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

void X6_Ecal_front(){

  time_t rawtime;
  struct tm * timeinfo;
  char buffer[80];

  time (&rawtime);
  timeinfo = localtime(&rawtime);
  strftime(buffer,sizeof(buffer),"%Y-%m-%d",timeinfo);
  std::string str(buffer);

  std::string calfile = "X6_Ecal_front_" + str + ".cal";
  cout << calfile << endl;
  ofstream outfile1(calfile.c_str());
  
  std::string mkdir = "mkdir Ecal_front_" + str;
  gSystem->Exec(mkdir.c_str());

  //TFile* ifile = new TFile("../X6_cribdata_bal-det12ss4.root","read");
  TFile* ifile1 = new TFile("../stark_0199.root","read");
  TFile* ifile2 = new TFile("../stark_0253.root","read");
  
  TCanvas* c1[8];
  TCanvas* c2[8];
  TCanvas* c3[8];
  for (int j=0; j<8; j++){
    c1[j] = new TCanvas(Form("X6 Calibration %d",j),Form("X6 Calibration %d",j),300,10*j,2200,1100);
    c1[j]->Divide(8,4);  
  
    c2[j] = new TCanvas(Form("X6 Calibration Projection %d",j),Form("X6 Calibration Projection %d",j),300+j*10,10,2200,1100);
    c2[j]->Divide(8,4);
  
    c3[j] = new TCanvas(Form("X6 Calibration Fit %d",j),Form("X6 Calibration Fit %d",j),300+j*10,10,2200,1100);
    c3[j]->Divide(8,4);
  }
  
  TGraph* gr[32][8];
  //Parameters variables
  double xpeakr[4];
  double offset[32][8];
  double slope[32][8];
  
  TH2D* epc;
  TH1D* proj; 
   
  //for (int i=0; i<32; i++){
  for (int i=12; i<32; i++){
      cout    << "= Detector "<< i+1 <<" ========================================="<<endl;
      for (int j=0; j<8; j++){
          cout    << "- Strip "<< j+1 <<" --------------------------------------------"<<endl;
          c1[i/4]->cd(8*(i%4)+j+1);
          gPad->SetGridy();
          gPad->SetGridx();

          cout << Form("Evspos_cal/E vs pos det%d strip%d",i+1,j+1);// << endl;
          if(i<12)	epc = (TH2D*)ifile1->Get(Form("Evspos_cal/E vs pos det%d strip%d",i+1,j+1));
          else	epc = (TH2D*)ifile2->Get(Form("Evspos_cal/E vs pos det%d strip%d",i+1,j+1));
          epc->Draw("colz");
          //Eback[i][j]->Rebin(2);
          //Eback[i][j]->Draw("hist");

          c2[i/4]->cd(8*(i%4)+j+1);
          gPad->SetGridx();
          proj=(TH1D*)epc->ProjectionY(Form("Y-proj%d det%d strip%d",15,i+1,j+1),0,epc->GetXaxis()->FindBin(15));
          if((i==19 && j==1)||(i==28 && j==0)||(i==30 && j==1)||(i==30 && j==7))	proj->Rebin(6);
          else if(i==0 && j==4)	proj->Rebin(8);
          else			proj->Rebin(5);
          //if(i==31 && j==5) proj->Rebin(2);
          proj->Draw();

          //Find peaks
          Int_t npeaks=2;
          Float_t fit_range;
          TSpectrum *s = new TSpectrum(npeaks);
          Int_t nfound = s->Search(proj,2,"new",0.1);
          //Eback[i][j]->Draw("hist");
          cout <<"\t\tnpeaks:  "<< npeaks <<"\t" << "nfound:  "<< nfound << endl;

          //Order peaks 
          xpeakr[0] = s->GetPositionX()[0];
          xpeakr[1] = s->GetPositionX()[1];
          sort(xpeakr,xpeakr+nfound);      
          //for (Int_t p=0;p<nfound;p++)	cout << p <<"\t"<< s->GetPositionX()[p] <<"\t"<< xpeakr[p] << endl;		

          //fit 4-gaussians
          /*TF1 *g[3];
            Double_t param[12];
            for (Int_t p=0;p<nfound;p++){
            g[p] = new TF1(Form("g%d",p),"gaus",xpeakr[p]*0.96,xpeakr[p]*1.05);
            g[p]->SetParLimits(0,5,Eback[i][j]->GetMaximum());		//fix amplitude to positive values
            g[p]->SetParLimits(1,xpeakr[p]*0.99,xpeakr[p]*1.01);		//fix mean around peak
            g[p]->SetParLimits(2,2,7);					//fix sigma to positive values
            Eback[i][j]->Fit(g[p],"QR+");  

            chback[npeaks-nfound+p] = g[p]->GetParameter(1);
            cout << p <<"\t"<< s->GetPositionX()[p] <<"\t"<< xpeakr[p] <<"\t"<< chback[npeaks-nfound+p] << endl;
            }*/

          c3[i/4]->cd(8*(i%4)+j+1);
          gPad->SetGridy();
          gPad->SetGridx();
          gr[i][j]=new TGraph(npeaks,xpeakr,en);
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
          outfile1 << offset[i][j] << "    \t" << slope[i][j] << endl;
      }
  }

  for (int j=0; j<8; j++){
      {
          std::string dir = "Ecal_front_"+ str + Form("/X6 front Evspos det%d-%d.png",4*j+1,4*j+4);
          c1[j]->SaveAs(dir.c_str());
      }
      {
          std::string dir = "Ecal_front_"+ str + Form("/X6 front Projection 0-15 det%d-%d.png",4*j+1,4*j+4);
          c2[j]->SaveAs(dir.c_str());
      }
      {
          std::string dir = "Ecal_front_"+ str + Form("/X6 front Calibration Fit det%d-%d.png",4*j+1,4*j+4);
          c3[j]->SaveAs(dir.c_str());
      }
  }
}
