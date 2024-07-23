#include "TFile.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TH2.h"
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
double en[2]={3.182,5.486};
//double en[3]={3.0744,5.400496,5.7194901};
//Float_t en[4]={3.0744,5.05698,5.400496,5.7194901};

double x, y;
TCanvas* c1;

void Point() {
   int event = gPad->GetEvent();
   int px    = gPad->GetEventX();
   int py    = gPad->GetEventY();
   x = gPad->AbsPixeltoX(px);
   y = gPad->AbsPixeltoY(py); 
   
   if(event==1) { // left mouse button click
      c1->WaitPrimitive();
      //cout <<"X:  "<< x << endl;
      //cout <<"Y:  "<< y << endl;
      return;
   }
}

int X6_Eflat(int det=32, int strip=8){

   TString name;
   int i,j,k, bin;

   //Parameters variables
   double x0,x1,y0,y1,gx,gy;
   //double x0[4];    double x1[4];    double y0[4];    double y1[4];
   //double gx[4];    double gy[4];    double p0[4];    double p1[4];	//The 4th index is for the mean values

   //char *w =new char[1];
   
   time_t rawtime;
   struct tm * timeinfo;
   char buffer[80];

   time (&rawtime);
   timeinfo = localtime(&rawtime);
   strftime(buffer,sizeof(buffer),"%Y-%m-%d",timeinfo);
   std::string str(buffer);

   std::string calfile = "X6_Eflat_" + str + ".cal";
   cout << calfile << endl;
   ofstream outfile1;
   outfile1.open(calfile.c_str());
   
   /*std::string inputfile = "X6_Eimproved_inputs" + str + ".cal";
   cout << calfile << endl;
   ofstream outfile2;
   outfile2.open(inputfile.c_str(), std::ios_base::app);*/
        
   std::string mkdir = "mkdir Eflat_" + str;
   gSystem->Exec(mkdir.c_str());

   c1 = new TCanvas("X6 E-calibration","X6 E-calibration",300,50,700,700);
   gPad->SetBorderMode(0);
   gPad->SetTopMargin(0.03);	//percentage
   gPad->SetBottomMargin(0.077);
   gPad->SetRightMargin(0.03);
   gPad->SetLeftMargin(0.092);
   gPad->SetGridx();
   gPad->SetGridy();
   
   TCanvas *projection;
   projection = new TCanvas("Projection","Projection",1000,50,700,700);
   projection->Divide(1,2);
   
   TFile* fin1 = new TFile("../stark_0199.root","read");
   TFile* fin2 = new TFile("../stark_0253.root","read");

   TH2D* lvh[32][8];
   gStyle->SetPalette(1);
   gStyle->SetPalette(1);
   
   //for (i=det; i<det+1; i++){
   for (i=0; i<32; i++){
     cout    << "= Detector "<< i+1 <<" ========================================="<<endl;
     //for (j=strip; j<strip+1; j++){
     for (j=0; j<8; j++){
        c1->cd();
	cout    << "- Strip "<< j+1 <<" --------------------------------------------"<<endl;
	if(i<12) lvh[i][j] = (TH2D*)fin1->Get(Form("UpvsDown/Up vs Down det%d strip%d",i+1,j+1));
	else	 lvh[i][j] = (TH2D*)fin2->Get(Form("UpvsDown/Up vs Down det%d strip%d",i+1,j+1));
	//lvh[i][j]->SetStats(0);
	lvh[i][j]->SetTitleOffset(0.5);
	lvh[i][j]->GetYaxis()->SetTitleOffset(1.);
	lvh[i][j]->GetYaxis()->SetTitleOffset(1.33);
	lvh[i][j]->GetXaxis()->SetRangeUser(0,4000);
	lvh[i][j]->GetYaxis()->SetRangeUser(0,4000);
	lvh[i][j]->SetMarkerStyle(20);
	lvh[i][j]->SetMarkerSize(0.64);
	lvh[i][j]->Draw("");
	c1->SetLogz();
	
	//Saving canvas
	std::string dir = "Eflat_"+ str + Form("/X6 Eflat det%d strip%d.png",i+1,j+1);
	c1->SaveAs(dir.c_str());
	
	//cout    << "\nProjecting bin X="<< k <<" ----------------------------------"<<endl;
	for (k=0; k<2; k++){
	   projection->cd(k+1);
	   y0 = 500.;	if((i==0 && j==4)||(i==8 && j==7)||(i==12 && j==2)||(i==19 && j==2)) y0 = 400.;
	   x1 = 500.;	if((i==12 && j==2)) x1 = 400.;
	   gPad->SetGridx();
	   if(k==0){
	      proj=(TH1D*)lvh[i][j]->ProjectionX(Form("%c-proj%d det%d strip%d",(char)(k+88),(int)y0,i+1,j+1),lvh[i][j]->GetYaxis()->FindBin(y0)-2,lvh[i][j]->GetYaxis()->FindBin(y0)+2);
	      proj->SetTitle(Form("%c-proj%d det%d strip%d",(char)(k+88),(int)y0,i+1,j+1));
	   } else{
	      proj=(TH1D*)lvh[i][j]->ProjectionY(Form("%c-proj%d det%d strip%d",(char)(k+88),(int)x1,i+1,j+1),lvh[i][j]->GetXaxis()->FindBin(x1)-2,lvh[i][j]->GetXaxis()->FindBin(x1)+2);
	      proj->SetTitle(Form("%c-proj%d det%d strip%d",(char)(k+88),(int)x1,i+1,j+1));
	   }
	   if(i==26 && j==0 && k==0)	proj->Rebin(3);
	   else				proj->Rebin(2);
	   //if(i==28 && j==6) proj->Rebin(2);
	   //if(i==31 && j==5) proj->Rebin(2);
	   proj->Draw("hist");
	   //proj->GetXaxis()->SetRangeUser(2.,6.);
	   
	   //Find peaks
	   Int_t nfound;
	   Float_t fit_range;
	   TSpectrum *s = new TSpectrum(1);
	   nfound = s->Search(proj,1,"new",0.4);
	   //if((i==7&&j==4)||(i==20&&j==3)||(i==21&&j==5)||(i==31&&j==6)) nfound = s->Search(proj,1,"new",0.05);
	   //proj->Draw("hist");
	   if(k==0)	x0 = s->GetPositionX()[0];
	   else		y1 = s->GetPositionX()[0];
	   cout << "Peak: "<< s->GetPositionX()[0] << endl;
	}
	  
	gx = y0*(y1-y0)/(x0*y1-x1*y0);
	gy = y0*(x0-x1)/(x0*y1-x1*y0);
	   
	outfile1 << fixed << setprecision(8);
	outfile1 << gx << "\t" << gy << endl;
	
	std::string dir = "Eflat_"+ str + Form("/X6 Eflat det%d strip%d-proj.png",i+1,j+1);
	projection->SaveAs(dir.c_str());
     }
   }
       
   return 0;
}
