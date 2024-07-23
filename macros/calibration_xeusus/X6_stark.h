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
#include "TMath.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <algorithm>


TTree* tr;
TTree *starktree;

//[Asad][Aget][Chan]
int stark_map[4][4][68];

//** Intput leaf variables *************************************************************
int iEvent;
int X6Mult;
int X6cobo[20];
int X6asad[20];
int X6aget[20];
int X6chan[20];
double X6nrj[20];
   
int STARK_Mult;
int X6_det[16];
int X6_strip[16];

//Set input tree address
void treeaddress(){
  tr->SetBranchAddress("Event",	&iEvent);
  tr->SetBranchAddress("Mult",	&X6Mult);
  tr->SetBranchAddress("cobo",	X6cobo);
  tr->SetBranchAddress("asad",	X6asad);
  tr->SetBranchAddress("aget",	X6aget);
  tr->SetBranchAddress("chan",	X6chan);
  tr->SetBranchAddress("energy",X6nrj);
  tr->SetBranchAddress("X6det",	X6_det);
  tr->SetBranchAddress("X6strip",X6_strip);
}

//** Output leaf variables ************************************************************
int		X6Evt;
//Int_t		X6Mult;
Int_t		X6Mult_up;
Int_t		X6Mult_dw;
Int_t		X6Mult_back;
int		X6det[16];		//[X6Mult]
int		X6det2[16];		//[X6Mult]
int		X6strip[16];		//[X6Mult]
Int_t		X6strip_front[642];	//[X6Mult]
Int_t		X6strip_back[642];	//[X6Mult]
Float_t		X6energy_up[642];	//[X6Mult]
Float_t		X6energy_up_cal[642];	//[X6Mult]
Float_t		X6energy_dw[642];	//[X6Mult]
Float_t		X6energy_dw_cal[642];	//[X6Mult]
Float_t		X6energy_back[642];	//[X6Mult] 
Float_t		X6energy_back_cal[642];	//[X6Mult] 
Float_t		X6energy_raw[642];	//[X6Mult]
Float_t		X6energy_cal[642];	//[X6Mult]
Float_t		X6energy[642];		//[X6Mult]
Float_t		X6position[642];	//[X6Mult]
Float_t		X6position_cal[642];	//[X6Mult]
Int_t		X6flag[642];		//[X6Mult]

void starktreedef(){
  starktree->Branch("Event",	&iEvent,	"iEvent/I");
  starktree->Branch("Mult",	&X6Mult,	"X6Mult/I");
  starktree->Branch("cobo",	X6cobo,		"X6cobo[X6Mult]/I");
  starktree->Branch("asad",	X6asad,		"X6asad[X6Mult]/I");
  starktree->Branch("aget",	X6aget,		"X6aget[X6Mult]/I");
  starktree->Branch("chan",	X6chan,		"X6chan[X6Mult]/I");
  starktree->Branch("energy",	X6nrj,		"X6nrj[X6Mult]/D");

  //starktree->Branch("X6Mult",	&X6_Mult,	"X6_Mult/I");
  //starktree->Branch("X6det",		X6_det,			"X6_det[X6Mult]/I");
  //starktree->Branch("X6strip",	X6_strip,	"X6_strip[X6Mult]/I");
  //starktree->Branch("X6strip_b",	X6_stripback,	"X6_strip_bX6_strip_b[X6Mult]/I");
  //starktree->Branch("X6pos",	X6_position,	"X6_position[X6Mult]/I");
  //starktree->Branch("X6Eup",	X6_Eup,		"X6_Eup[X6Mult]/I");
  //starktree->Branch("X6Edwn",	X6_Edwn,	"X6_Edwn[X6Mult]/I");
  //starktree->Branch("X6E",	X6_E,		"X6_E[X6Mult]/I");
  //starktree->Branch("X6E_back",	X6_Eback,	"X6_Eback[X6Mult]/D");
    
  
  starktree->Branch("det",		X6det,			"X6det[X6Mult]/I"); 
  starktree->Branch("strip",		X6strip,		"X6strip[X6Mult]/I"); 
  starktree->Branch("X6Mult",		&STARK_Mult,		"STARK_Mult/I");
  starktree->Branch("X6det",		X6det2,			"X6det2[STARK_Mult]/I");
  starktree->Branch("X6strip_front",	X6strip_front,		"X6strip_front[STARK_Mult]/I");
  starktree->Branch("X6strip_back",	X6strip_back,		"X6strip_back[STARK_Mult]/I"); 
  starktree->Branch("X6energy_up",	X6energy_up,		"X6energy_up[STARK_Mult]/F");
  starktree->Branch("X6energy_up_cal",	X6energy_up_cal,	"X6energy_up_cal[STARK_Mult]/F");
  starktree->Branch("X6energy_dw",	X6energy_dw,		"X6energy_dw[STARK_Mult]/F");
  starktree->Branch("X6energy_dw_cal",	X6energy_dw_cal,	"X6energy_dw_cal[STARK_Mult]/F");
  starktree->Branch("X6energy_back",	X6energy_back, 		"X6energy_back[STARK_Mult]/F");
  starktree->Branch("X6energy_back_cal",X6energy_back_cal,	"X6energy_back_cal[STARK_Mult]/F"); 
  starktree->Branch("X6energy_raw",	X6energy_raw,		"X6energy_raw[STARK_Mult]/F");
  starktree->Branch("X6energy_cal",	X6energy_cal,		"X6energy_cal[STARK_Mult]/F");
  starktree->Branch("X6energy_bal",	X6energy,		"X6energy[STARK_Mult]/F");
  starktree->Branch("X6position",	X6position,		"X6position[STARK_Mult]/F");
  starktree->Branch("X6position_cal",	X6position_cal,		"X6position_cal[STARK_Mult]/F");
  starktree->Branch("X6flag",		X6flag,			"X6flag[X6Mult]/I");
}

//** Setting first & last entries to be processed *************************************
int FirstEntry(int FirstEntry, int LastEntry){
  if(FirstEntry<0){
    cout << "ERROR: First entry < 0. Starting at 0 instead" << endl;
    return 0;
  } else if(FirstEntry>LastEntry){
    cout << "ERROR: First entry > LastEntry. Processing from " << LastEntry-1 << " instead" << endl;
    return LastEntry-1;
  } else {
    return FirstEntry;
  }
}

int LastEntry(int LastEntry){
  int N = tr->GetEntries();
  if(LastEntry==0)	return N;		//Til the last entry by default
  else if(LastEntry>N){
    cout << "ERROR: Last entry > NEntries!!! End at "<< N << endl;
    return N;
  } else if(LastEntry<0){
    cout << "ERROR: Last entry < 0!!! Please set LastEntry between 0 and "<< N <<". Processing only entry 0" << endl;
    return 1;
  } else 		return LastEntry;
}
   
//** Calibration coefficients *********************************************************
//int thr_up[642],thr_dwn[642];					//Inner barrel thresholds
float aback[32][4],bback[32][4];				//inner barrel backs: offset & slope
float gdwn[32][8],gup[32][8],p0[32][8],p1[32][8];		//inner barrel e628_analysis pos coefs
float b0[32][8],b1[32][8],b2[32][8];				//ballisitic def

//HISTOGRAMS ==============================================
TH1D* Energy[32][12];
TH1D* Signal[32][20];
//TH1D* Test[600];

TH2D* Evspos_cal[32][8];
TH2D* Evspos_bal[32][8];
TH2D* UpvsDown[32][8];
TH2D* hitpattern = new TH2D("hitpattern","hitpattern",8,1,9,76,1,76);
//TH2D* Evspos_tot = new TH2D("E vs pos","E vs pos",75,1,76,6000,0,6000);  
TH2D* Evspos_tot = new TH2D("E vs pos","E vs pos",75,1,76,5000,0,6.);  

TH2D* Ebal_strip = new TH2D("E vs strip","E vs strip",256,0,256,1024,0,6.);
TH2D* Ecal_strip = new TH2D("Ecal vs strip","Ecal vs strip",256,0,256,1024,0,6.);
TH2D* Eraw_strip = new TH2D("Eraw vs strip","Eraw vs strip",256,0,256,1024,0,4096.);
TH2D* Ecal_backs = new TH2D("GM Backs vs strip","GM Backs vs strip",128,0,128,1024,0,6.);
TH2D* Eraw_backs = new TH2D("Raw Backs vs strip","Raw Backs vs strip",128,0,128,1250,0,1250);

void histos(){
  int i,j;
  for(i=0;i<32;i++){
    for(j=0;j<8;j++){
      //Signal[2*j]	= new TH1D(Form("Downstream det%d strip%d",i+1,j+1),	Form("Downstream det%d strip%d",i+1,j+1),	1000,0,7.);
      //Signal[2*j+1]	= new TH1D(Form("Upstream det%d strip%d",i+1,j+1),	Form("Upstream det%d strip%d",i+1,j+1),		1000,0,7.);
      Signal[i][2*j]	= new TH1D(Form("Upstream det%d strip%d",i+1,j+1),	Form("Upstream det%d strip%d",i+1,j+1),		6000,0,6000);
      Signal[i][2*j+1]	= new TH1D(Form("Downstream det%d strip%d",i+1,j+1),	Form("Downstream det%d strip%d",i+1,j+1),	6000,0,6000);
    
      Energy[i][j]	= new TH1D(Form("E det%d strip%d-front",i+1,j+1),	Form("E det%d strip%d-front",i+1,j+1),		2000,0,7.);
      //Evspos_cal[j]	= new TH2D(Form("E vs pos strip%d",i+1,j+1),		Form("E vs pos det%d strip%d",i+1,j+1),		75,1,76,1024,0,7.);
      Evspos_cal[i][j]	= new TH2D(Form("E vs pos det%d strip%d",i+1,j+1),	Form("E vs pos det%d strip%d",i+1,j+1),		75,1,76,6000,0,7);
      //Evspos_bal[j]	= new TH2D(Form("Ebal vs pos strip%d",i+1,j+1),		Form("Ebal vs pos det%d strip%d",i+1,j+1),	75,1,76,1024,0,7.);
      Evspos_bal[i][j]	= new TH2D(Form("Ebal vs pos det%d strip%d",i+1,j+1),	Form("Ebal vs pos det%d strip%d",i+1,j+1),	75,1,76,6000,0,7);
      //UpvsDown[j]	= new TH2D(Form("Up vs Down det%d strip%d",i+1,j+1),	Form("Up vs Down det%d strip%d",i+1,j+1),	500,0,7.,500,0,7.);
      UpvsDown[i][j]	= new TH2D(Form("Up vs Down det%d strip%d",i+1,j+1),	Form("Up vs Down det%d strip%d",i+1,j+1),	400,0,4000,400,0,4000);
    }
    for(j=8;j<12;j++){
      Signal[i][j+8]	= new TH1D(Form("Back signal det%d strip%d",i+1,j+1), 	Form("Back signal det%d strip%d",i+1,j+1),	1250,0,1250);
      //Energy[i][j]	= new TH1D(Form("E det%d strip%d-back",i+1,j-7),	Form("E det%d strip%d-back",i+1,j-7),		6000,0,6000);
      Energy[i][j]	= new TH1D(Form("E det%d strip%d-back",i+1,j-7),	Form("E det%d strip%d-back",i+1,j-7),		2000,0,7.);
    }
  }
  //for(i=0;i<600;i++)
  //    Test[i]		= new TH1D(Form("E ch%d (det%d strip%d)",i,i/20,i%20),	Form("E ch%d (det%d strip%d)",i,i/20,i%20),	6000,0,6000);
}


void readmap(string calfile="susostarkmap.txt"){
  ifstream inputfile;
  int asad, aget, chan, index;
  for(int i=0; i<4; i++){
    for(int j=0; j<4; j++){
      for(int k=0; k<68; k++)
	stark_map[i][j][k] = -1; 
    }  
  }
  
  inputfile.open(calfile.c_str());
  for(int j=0; j<640; j++){
    if (inputfile.is_open()){
      inputfile >> asad >> aget >> chan >> index;	//Reading x and w from file
      stark_map[asad][aget][chan] = index;		//Reading x and w from file
    } else {
      cout << "ERROR opening mapping file " << calfile.c_str() << endl;
      return;
    }
  }
  inputfile.close();
}


void readcal(){
   int i;   
   //** CALIBRATION PARAMETERS ***********************************************************
   ifstream inf0, inf1, inf2;
   inf0.open("parameters/X6_back.cal");
   for (i=0; i<128; i++){
     if (inf0.is_open()){
       inf0 >> aback[i/4][i%4] >> bback[i/4][i%4];
     } else {
       aback[i/4][i%4] = 0.;
       bback[i/4][i%4] = 0.;
     }
   }
   //inf0.close();
   /*cout << "X6 BACKS" << endl;
   for (i=0; i<128; i++){
     cout << i << "  " << aback[i/4][i%4] << "  " << bback[i/4][i%4] << endl;
   }*/
   
   inf1.open("parameters/X6_Eimproved.cal");
   for (i=0; i<256; i++){
     if (inf1.is_open()){
       inf1 >> gup[i/8][i%8]  >> gdwn[i/8][i%8] >> p0[i/8][i%8] >> p1[i/8][i%8];
     } else {
       gup[i/8][i%8]	= 1.0;
       gdwn[i/8][i%8]	= 1.0;
       p0[i/8][i%8]	= 0.0;
       p1[i/8][i%8]	= 1.0;
     }
   }
   //inf1.close();
   cout << fixed << setprecision(5);
   /*cout << "X6 CALIBRATIONS" << endl;
   for (i=0; i<256; i++){
     cout << i << "   " << gup[i/8][i%8] << "  " << gdwn[i/8][i%8] << "  " << p0[i/8][i%8] << "  " << p1[i/8][i%8] << endl;
   }*/

   inf2.open("parameters/X6_ballistic.cal");
   for (i=0; i<256; i++){
     if (inf2.is_open()){
       inf2 >> b0[i/8][i%8] >> b1[i/8][i%8] >> b2[i/8][i%8];
     } else {
       b0[i/8][i%8] = 1.0;
       b1[i/8][i%8] = 0.0;
       b2[i/8][i%8] = 0.0;
     }
   }
   inf2.close();
   /*cout << "X6 BALLISTIC" << endl;
   for (i=0; i<256; i++){
     cout << i << " " << b0[i/8][i%8] << " " << b1[i/8][i%8] << " " << b2[i/8][i%8] << endl;
   }*/
}
