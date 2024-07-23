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

using namespace std;

Int_t ci = 314; 	// color index
//TColor *color = new TColor(ci, 169/256., 196/256., 256/256.);

Int_t npeaks=2;
//Alpha source energies in MeV
Float_t en[2]={3.182,5.486};
//Float_t en[3]={4.89619,5.23426,5.56150};


int X6_ballistic(){

  int i,j,k,l,p,n;
  
  time_t rawtime;
  struct tm * timeinfo;
  char buffer[80];

  time (&rawtime);
  timeinfo = localtime(&rawtime);
  strftime(buffer,sizeof(buffer),"%Y-%m-%d",timeinfo);
  std::string str(buffer);

  std::string calfile = "X6_ballistic_" + str + ".cal";
  cout << calfile << endl;
  ofstream outfile(calfile.c_str());
  
  std::string mkdir = "mkdir ballistic_" + str;
  gSystem->Exec(mkdir.c_str());
  
  //Graph variables
  float x[2][10];
  float mean[2][10];
  float FWHM[2][10];
  
  float fitpar[2][3];
  float calpar[32][8][3];
    
  TFile* fin1 = new TFile("../stark_0199.root","read");
  TFile* fin2 = new TFile("../stark_0253.root","read");

  TCanvas *ballistic;
  ballistic = new TCanvas("X6 Ballistic","X6 Ballistic",300,100,1600,1100);
  ballistic->Divide(2,1);
  
  TCanvas *projection;
  projection = new TCanvas("Projection","Projection",100,10,2100,1100);
  projection->Divide(5,2);

  TMultiGraph *Alphas;
  TGraphErrors *Alpha[2];

  TH2D* epc;
  TH1D* proj; 

  TF1 *g[3];
  TF1* fitf[3];
  TF1* myfit = new TF1("myfit","[0]+[1]*x+[2]*x*x",0.0,1.2);  
  
  int det=32;
  int strip=8;
  //for (i=0; i<det; i++){
  for (i=12; i<det; i++){
  //for (i=0; i<det; i++){
     cout    << "= Detector "<< i+1 <<" ========================================="<<endl;
     for (j=0; j<strip; j++){
     //for (j=0; j<strip; j++){
	cout    << "- Strip "<< j+1 <<" --------------------------------------------"<<endl;
	
	//ballistic = new TCanvas(Form("X6 Ballistic det%d strip%d",i+1,j+1),Form("X6 Ballistic det%d strip%d",i+1,j+1),300,10*j,1600,1100);
	//ballistic->Divide(2,1);
	ballistic->cd(1);
	
	if(i<12) epc = (TH2D*)fin1->Get(Form("Evspos_cal/E vs pos det%d strip%d",i+1,j+1));
	else	 epc = (TH2D*)fin2->Get(Form("Evspos_cal/E vs pos det%d strip%d",i+1,j+1));
	epc->Draw("colz");

	n=0;
	for (k=15; k<64; k+=5){
	   //cout    << "\nProjecting bin X="<< k <<" ----------------------------------"<<endl;
	   projection->cd((k-15)/5+1);
	   gPad->SetGridx();
	   proj=(TH1D*)epc->ProjectionY(Form("Y-proj%d det%d strip%d",k,i+1,j+1),k-2,k+2);
	   proj->Rebin(20);
	   if(i>27) proj->Rebin(3);
	   //if(i==28 && j==6) proj->Rebin(2);
	   //if(i==31 && j==5) proj->Rebin(2);
	   proj->Draw();
	   proj->GetXaxis()->SetRangeUser(2.,6.);

	   //Find peaks
	   Int_t nfound;
	   Float_t fit_range;
	   TSpectrum *s = new TSpectrum(npeaks);
	   nfound = s->Search(proj,3,"new",0.05);
	   if((i==7&&j==4)||(i==20&&j==3)||(i==21&&j==5)||(i==31&&j==6)) nfound = s->Search(proj,1,"new",0.05);
	   
	   //Order peaks 
	   Double_t* xpeaks = s->GetPositionX();
	   sort(xpeaks,xpeaks+nfound);  
	   //for (p=0;p<nfound;p++) cout << xpeaks[p] << endl;
	   
	   if(nfound!=2) continue;
	   for (p=0;p<nfound;p++){
		g[p] = new TF1(Form("g%d",p),"gaus",xpeaks[p]*0.96,xpeaks[p]*1.05);
		g[p]->SetParLimits(0,5,proj->GetMaximum());		//fix amplitude to positive values
		g[p]->SetParLimits(1,xpeaks[p]*0.99,xpeaks[p]*1.01);		//fix mean around peak
		g[p]->SetParLimits(2,2,7);					//fix sigma to positive values
		proj->Fit(g[p],"QR+");  
    			
		x[p][n]    = k;
		mean[p][n] = g[p]->GetParameter(1); 
		FWHM[p][n] = 2.35*g[p]->GetParameter(2)/sqrt(g[p]->Integral(mean[p][n]-2*g[p]->GetParameter(2),mean[p][n]+2*g[p]->GetParameter(2)));
	   }
	   n++;
	}
	   
	ballistic->cd(2);
	gPad->SetGridx();
	gPad->SetGridy();
	
	Alphas = new TMultiGraph();
	Alphas->SetTitle(Form("Ballistic deficit correction det%d strip%d",i+1,j+1));
	for(Int_t p=0;p<npeaks;p++){   
	   Alpha[p] = new TGraphErrors(n,x[p],mean[p],0,FWHM[p]);
	   Alpha[p]->SetMarkerColor(2+p*2);
	   Alpha[p]->SetMarkerSize(1.2);
	   Alpha[p]->SetMarkerStyle(20);
	   Alphas->Add(Alpha[p]);
	   
	   if(i==7  && j==7 && p==0) Alpha[p]->RemovePoint(n-1);
	   if(i==7  && j==7 && p==0) Alpha[p]->RemovePoint(n-2);
	   if(i==7  && j==7 && p==1) Alpha[p]->RemovePoint(n-2);
	   if(i==20 && j==7 && p==1) Alpha[p]->RemovePoint(n-2);
	   if(i==28 && j==1 && p==1) Alpha[p]->RemovePoint(n-1);
	   if(i==28 && j==3 && p==1) Alpha[p]->RemovePoint(n-3);
	   if(i==28 && j==6 && p==1) Alpha[p]->RemovePoint(n-3);
	   if(i==29 && j==1 && p==0) Alpha[p]->RemovePoint(n-1);
	   if(i==29 && j==3 && p==1) Alpha[p]->RemovePoint(n-1);
	   if(i==29 && j==4 && p==1) Alpha[p]->RemovePoint(4);
	   if(i==29 && j==4 && p==1) Alpha[p]->RemovePoint(3);
	   if(i==30 && j==2 && p==0) Alpha[p]->RemovePoint(n-3);
	   if(i==30 && j==2 && p==1) Alpha[p]->RemovePoint(n-3);
	   if(i==31 && j==7 && p==0) Alpha[p]->RemovePoint(n-1);
	   
	   Alpha[p]->Fit("pol2","Q+");
	   fitf[p] = Alpha[p]->GetFunction("pol2");
       if (fitf[p]==nullptr) {
           cout << "!!!!!!!!!!!!!!" << endl;
           continue;
       }
	   fitf[p]->SetLineColor(2+p*2);
	   
	   for(l=0;l<3;l++){
	      fitpar[p][l] = fitf[p]->GetParameter(l)/en[p];
	      //cout << fitpar[p][l] << "\t";
	   }
	   //cout << endl;
	}
	Alphas->Draw("AP");
	Alphas->GetYaxis()->SetRangeUser(0.,7.);
	
	for(l=0;l<3;l++){
	   calpar[i][j][l] = (fitpar[0][l] + fitpar[1][l])/2.;
	   if((i==7&&j==7)||(i==28&&j==0)) calpar[i][j][l] = fitpar[0][l];
	}
	
	//Printing results
	outfile << fixed << setprecision(10);
	outfile << calpar[i][j][0] << "    " << calpar[i][j][1] << "    " << calpar[i][j][2] << endl;

	//cout    << "\n    b0    " << "    " << "    b1    " << "    " << "    b2  " << endl;
	//cout    << fixed << setprecision(10);
	//cout    << calpar[i][j][0] << "    " << calpar[i][j][1] << "    " << calpar[i][j][2] << endl;

	//Saving canvas
	std::string dir = "ballistic_"+ str + Form("/X6 Ballistic det%d strip%d.png",i+1,j+1);
	ballistic->SaveAs(dir.c_str());
	//ballistic->SaveAs(Form("ballistic/X6 Ballistic det%d strip%d",i+1,j+1));
      }
   }
       
   return 0;
}
