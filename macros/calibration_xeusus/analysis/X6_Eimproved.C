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

int X6_Eimproved(int det=32, int strip=8){

   TString name;
   int i,j,k;

   //Parameters variables
   double x0[4];    double x1[4];    double y0[4];    double y1[4];
   double gx[4];    double gy[4];    double p0[4];    double p1[4];	//The 4th index is for the mean values

   //char *w =new char[1];
   
   time_t rawtime;
   struct tm * timeinfo;
   char buffer[80];

   time (&rawtime);
   timeinfo = localtime(&rawtime);
   strftime(buffer,sizeof(buffer),"%Y-%m-%d",timeinfo);
   std::string str(buffer);

   std::string calfile = "X6_Eimproved" + str + ".cal";
   cout << calfile << endl;
   ofstream outfile1;
   outfile1.open(calfile.c_str(), std::ios_base::app);
   
   std::string inputfile = "X6_Eimproved_inputs" + str + ".cal";
   cout << calfile << endl;
   ofstream outfile2;
   outfile2.open(inputfile.c_str(), std::ios_base::app);
        
   std::string mkdir = "mkdir Ecal_" + str;
   gSystem->Exec(mkdir.c_str());

   c1 = new TCanvas("X6 E-calibration","X6 E-calibration",300,10,700,700);
   gPad->SetBorderMode(0);
   gPad->SetTopMargin(0.03);	//percentage
   gPad->SetBottomMargin(0.077);
   gPad->SetRightMargin(0.03);
   gPad->SetLeftMargin(0.092);
   gPad->SetGridx();
   gPad->SetGridy();
   
   TFile* fin1 = new TFile("../stark_0199.root","read");
   TFile* fin2 = new TFile("../stark_0253.root","read");

   TH2D* lvh[32][8];
   gStyle->SetPalette(1);
   gStyle->SetPalette(1);
   
   //for (i=det; i<det+1; i++){
   for (i=0; i<det; i++){
     cout    << "= Detector "<< i+1 <<" ========================================="<<endl;
     //for (j=strip; j<strip+1; j++){
     for (j=0; j<strip; j++){
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
	std::string dir = "Ecal_"+ str + Form("/X6 Ecal det%d strip%d.png",i+1,j+1);
	c1->SaveAs(dir.c_str());
	
	//NO DEAD CHANNELS
	/*if((i==2 && j==7)){	//DEAD CHANNELS
	  gx[4] = 0.;
	  gy[4] = 0.;
	  p0[4] = 0.;
	  p1[4] = 0.;
	  
	  outfile1 << fixed << setprecision(8);
	  outfile1 << gx[4] << "\t" << gy[4] << "\t" << p0[4] << "\t" << p1[4] << endl;
	  
	  outfile2 << fixed << setprecision(8);
	  outfile2 << 0.000 << "\t" << 0.000 << "\t" << 0.000 << "\t" << 0.000 << endl;
	  continue;
	}*/
	//c1->AddExec("Point()","Point()");c1->WaitPrimitive();
	
	for(k=-1; k<4; k++){
	  c1->AddExec("Point()","Point()");
	  if(k==-1) c1->Update();
	  else{
	    if(k%2==0){
	      x0[k/2] = x;
	      y0[k/2] = y;
	      cout    << "\n- Alpha ="<< en[k/2] <<" MeV ----------------------------------\n"<<endl;
	      cout << "   + Upstream x and y coords for "  << en[k/2] <<" MeV:\t" << x0[k/2] <<"\t"<< y0[k/2] << endl;
	    } else {
	      x1[k/2] = x;
	      y1[k/2] = y;
	      cout << "   + Downstream x and y coords for "<< en[k/2] <<" MeV:\t" << x1[k/2] <<"\t"<< y1[k/2] << endl;
	      
	      outfile2 << fixed << setprecision(8);
	      outfile2 << x0[k/2] << "\t" << y0[k/2] << "\t" << x1[k/2] << "\t" << y1[k/2] << endl;
	      
	      gx[k/2] = en[k/2]*(y1[k/2]-y0[k/2])/(x0[k/2]*y1[k/2]-x1[k/2]*y0[k/2]);
	      gy[k/2] = en[k/2]*(x0[k/2]-x1[k/2])/(x0[k/2]*y1[k/2]-x1[k/2]*y0[k/2]);
	      p0[k/2] = (gy[k/2]*y0[k/2]-gx[k/2]*x0[k/2])/(gx[k/2]*x0[k/2]+gy[k/2]*y0[k/2]);
	      p1[k/2] = (gy[k/2]*y1[k/2]-gx[k/2]*x1[k/2])/(gx[k/2]*x1[k/2]+gy[k/2]*y1[k/2]);

	      cout << fixed << setprecision(5);
	      cout << "   " << endl;
	      cout << "   +----------+----------+----------+----------+" << endl;
	      cout << "   |    gx    |    gy    |    p0    |    p1    |" << endl;
	      cout << "   +----------+----------+----------+----------+" << endl;
	      cout << "   |  " << gx[k/2] << " |  " << gy[k/2] << " |  " << p0[k/2] << " | " << p1[k/2] << " |" << endl;
	      cout << "   +----------+----------+----------+----------+" << endl;
	    }
	    //cout << k <<"  X:  "<< x << endl;
	    //cout << k <<"  Y:  "<< y << endl;
	  }
	  if(k<4) c1->WaitPrimitive();
	  else break;
	}
	outfile2 << endl;

	cout    << "\n- Mean values --------------------------------------------------"<<endl;
	gx[4] = (gx[0]+gx[1])/2.;
	gy[4] = (gy[0]+gy[1])/2.;

	p0[4] = (p0[0]+p0[1])/2.;
	p1[4] = (p1[0]+p1[1])/2.;
	//3-point calibration
	/*gx[4] = (gx[0]+gx[1]+gx[2])/3.;
	gy[4] = (gy[0]+gy[1]+gy[2])/3.;

	p0[4] = (p0[0]+p0[1]+p0[2])/3.;
	p1[4] = (p1[0]+p1[1]+p1[2])/3.;*/
	//4-point calibration
	/*gx[4] = (gx[0]+gx[1]+gx[2]+gx[3])/4.;
	gy[4] = (gy[0]+gy[1]+gy[2]+gy[3])/4.;

	p0[4] = (p0[0]+p0[1]+p0[2]+p0[3])/4.;
	p1[4] = (p1[0]+p1[1]+p1[2]+p1[3])/4.;*/
       
	cout << fixed << setprecision(5);
	cout << "   " << endl;
	cout << " +----------+----------+----------+----------+" << endl;
	cout << " |    gx    |    gy    |    p0    |    p1    |" << endl;
	cout << " +----------+----------+----------+----------+" << endl;
	cout << " |  " << gx[4] << " |  " << gy[4] << " |  " << p0[4] << " | " << p1[4] << " |" << endl;
	cout << " +----------+----------+----------+----------+" << endl;
	outfile1 << fixed << setprecision(8);
	outfile1 << gx[4] << "\t" << gy[4] << "\t" << p0[4] << "\t" << p1[4] << endl;
     }
   }
       
   return 0;
}
