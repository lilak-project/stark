
    int det,strip;
    int det_back,strip_back;

    //[Asad][Aget][Chan]
    int stark_map[3][4][71];

    int iEvent;
    int X6Mult;
    int X6cobo[20];
    int X6asad[20];
    int X6aget[20];
    int X6chan[20];
    double X6nrj[20];
   
    int X6_Mult;
    int X6_det[16];
    int X6_strip[16];
    int X6_stripback[16];
    double X6_position[16];
    double X6_E[16];
    double X6_Eup[16];
    double X6_Edwn[16];
    double X6_Eback[16];


    int E_ohmic[36][4];
    int E_junction[36][8][2];
    int CSD_ohmic[12];
    int CSD_junction[12][8];
    TH1D* hCSD_ohmic[12];
    TH1D* hCSD_junction[12][8];
    TCanvas* cCSD[12];
    TH2D* hitpattern = new TH2D("hitpattern","hitpattern;ch;energy",2000,0,2000,300,0,3000);
    TH2D* LvsR[36][8];
    TH2D* EvsPos[36][8];
    TH2D* StripvsPos[36][8];
    TH1D* E_ohmic_raw[36][4];
    TTree* X6tree;
    TCanvas* cX6[36];


void treedef(){
    X6tree->Branch("Event",	&iEvent,	"iEvent/I");
    X6tree->Branch("Mult",	&X6Mult,	"X6Mult/I");
    X6tree->Branch("cobo",	X6cobo,		"X6cobo[X6Mult]/I");
    X6tree->Branch("asad",	X6asad,		"X6asad[X6Mult]/I");
    X6tree->Branch("aget",	X6aget,		"X6aget[X6Mult]/I");
    X6tree->Branch("chan",	X6chan,		"X6chan[X6Mult]/I");
    X6tree->Branch("energy",	X6nrj,		"X6nrj[X6Mult]/D");

    //X6tree->Branch("X6Mult",	&X6_Mult,	"X6_Mult/I");
    X6tree->Branch("X6det",	X6_det,		"X6_det[X6Mult]/I");
    X6tree->Branch("X6strip",	X6_strip,	"X6_strip[X6Mult]/I");
    //X6tree->Branch("X6strip_b",	X6_stripback,	"X6_strip_bX6_strip_b[X6Mult]/I");
    //X6tree->Branch("X6pos",	X6_position,	"X6_position[X6Mult]/I");
    //X6tree->Branch("X6Eup",	X6_Eup,		"X6_Eup[X6Mult]/I");
    //X6tree->Branch("X6Edwn",	X6_Edwn,	"X6_Edwn[X6Mult]/I");
    //X6tree->Branch("X6E",	X6_E,		"X6_E[X6Mult]/I");
    //X6tree->Branch("X6E_back",	X6_Eback,	"X6_Eback[X6Mult]/D");
}
 
void def(){
    for(int i=0;i<36;i++){
    	for(int j=0;j<4;j++){
		E_ohmic_raw[i][j] = new TH1D(Form("E_ohmic_raw_%d_%d",i,j),Form("E_ohmic_raw_%d_%d;channel;counts",i,j),3000,0,3000);
    	}
    	for(int j=0;j<8;j++){
		LvsR[i][j] = new TH2D(Form("LvsR_%d_%d",i,j),Form("LvsR_%d_%d;Right;Left",i,j),500,0,3000,500,0,3000);
		EvsPos[i][j] = new TH2D(Form("EvsPos_%d_%d",i,j),Form("EvsPos_%d_%d;Pos(ratio);Energy",i,j),400,-1,1,300,0,3000);
		StripvsPos[i][j] = new TH2D(Form("StripvsPos_%d_%d",i,j),Form("StripvsPos_%d_%d;Pos(ratio);Energy",i,j),400,-1,1,8,1,9);
    	}
    }
    for(int i=0;i<12;i++){
    	CSD_ohmic[i]=0;
	hCSD_ohmic[i] = new TH1D(Form("CSD_ohmic_%d",i),Form("CSD_ohmic_%d;channel;counts",i),3000,0,3000);
    	for(int j=0;j<8;j++){
    		CSD_junction[i][j]=0;
		hCSD_junction[i][j] = new TH1D(Form("CSD_junction%d_%d",i,j),Form("CSD_junction_%d_%d;channel;counts",i,j),3000,0,3000);
    	}
    }
}


void readmap(string calfile="susostarkmap.txt"){
  ifstream inputfile;
  int asad, aget, chan, index;
  for(int i=0; i<4; i++){
    for(int j=0; j<4; j++){
      for(int k=0; k<71; k++)
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



/*
void susostarkmap(string calfile="starkmap_output.txt"){
  ifstream inputfile;
  ofstream outputfile;
  int cobo, asad, aget, chanfpn, chan, det, junc_ohm, strip, updwn, index;
  inputfile.open(calfile.c_str());
  outputfile.open("susostarkmap.txt");
  for(int j=0; j<640; j++){
    if (inputfile.is_open()){
      inputfile >> cobo >> asad >> aget >> chanfpn >> chan >> det >> junc_ohm >> strip >> updwn;	//Reading from file
      //cout << cobo <<"\t"<< asad <<"\t"<< aget <<"\t"<< chanfpn <<"\t"<< chan <<"\t"<< det <<"\t"<< junc_ohm <<"\t"<< strip <<"\t"<< updwn << endl;
      if(junc_ohm==1) index = 20*(det-1)+2*strip+updwn;
      else            index = 20*(det-1)+strip+16;	
      cout << asad <<"\t"<< aget <<"\t"<< chanfpn <<"\t"<< index << endl;
      outputfile << asad <<"\t"<< aget <<"\t"<< chanfpn <<"\t"<< index << endl;
    } else {
      cout << "ERROR opening mapping file " << calfile.c_str() << endl;
      return;
    }
  }
  inputfile.close();
}
*/
