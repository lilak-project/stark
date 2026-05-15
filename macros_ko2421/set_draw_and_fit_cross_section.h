bool fShowAll = true;
TString fMainName, fMainTitle, fAnaName;
TTree* fTree = nullptr;
ELARK* fELARK = nullptr;
TGraph* fGraphQE = nullptr;
TGraph* fGraphRF = nullptr;
TGraph* fGraphFC[100];
int fNumFresco = 0;
LKDrawingGroup* fTop = nullptr;
LKDrawingGroup* fSave = nullptr;
TFile* fSaveFit = nullptr;
LKParameterContainer* fPar = nullptr;
TGraph* fGraphFCDivRF = nullptr;
TGraph* fGraphThetaError[8];
LKBinning bnnE, bnnQ, bnnL, bnnC, bnnT, bnnT0, bnnSplitX12, bnnSplitX16, bnnSplitX12Com, bnnSplitX16Com;
ofstream fFileFitParameterOut;
ifstream fFileFitParameterIn;
int ***fHistParameters;
double ***fFitRange;
double ****fFitParameters;
vector<int> fProjectionBins;
vector<int> fIDArray12, fIDArray16;
TString fDataXPath = "data_x";
bool fPickup = false;
int fBackgroundType16 = -1;
int fBackgroundType12 = -1;
bool fAllDraw1 = false;
double fThetaCoMScaleRange1 = 44;
double fThetaCoMScaleRange2 = 48;
double fSolidAngleRatio = 0.8;
int fNumPairs = 28;
int fMaxPoints = 50;
int fNumVariables = 10;
int fNumVE = 6;
double**** fDataValues;
double*** fDataRange;
double fQProjTestSigmaFactor, fQProjSigmaFactor, fBeamEnergy, fSelectBeamEnergy, fHistQEntryCut;
double fQProjSigmaFactor2, fQProjSigmaFactor3;
double fBinWidthForFit = 0;
double fNumAtomsPerUnitArea = 1.72e+19;
double fBeamCount;
double fmbarn = 1e-27; // >> mbarn
TString fCrossSectionTitle = "#theta_{CoM.} (deg);#it{d}#sigma/#it{d}#Omega (mb/sr)";
bool fUseSavedLinearBG = true;
int fEnergyIndex;
double bIndex[3];
double bValue[7];
double fFinalThetaRange[2][2];
double fGlobalBeamNormalizationFactor = 0;
TGraph* fQTEdge;
double fSignalRatioCut = 0.2;
TString fNameFileX;
TString fBeamEnergyString;
int fHistEntriesCut = 800;
TGraph* fGraphK[2] = {0};
bool fUseScalingWithPrevData = false;
TString fTagBnn;
TGraph *fGraphParameterGGAR = nullptr;
TF1* fFitGGAR = nullptr;
bool fUseSecondPeak = false;

TFile *fParFile;
TTree *fParTree;
TGraph* fGraphFix[2];
double ScaleFunction0(double *x, double *par) { return (par[0] * fGraphFix[0]->Eval(x[0])); }
double ScaleFunction1(double *x, double *par) { return (par[0] * fGraphFix[1]->Eval(x[0])); }

///////////////////////////////////////////////////////////////////////////////////////////////////////////
double          GetFirstBoundary(TH1D* histQ, double thresholdFactor=0.01);
LKDrawingGroup* AnalyzeWithBeamEnergyAndPairID(double beamEnergy, int pairID=-1, bool addToTop=true);
bool            Initialize(int energyIndex, int optionIndex);
TGraphErrors*   FitAmplitude1To2(int idx, TGraphErrors* graph1, TGraph* graph2, TF1* f1, double &scale);
void            DrawWriteDrawings(bool justDraw = true);
void            DrawWriteEach();
void            WriteParameters();
LKDrawingGroup* DrawCrossSection(int selPair);
void            DrawXPair(int pairID);
bool            ShouldIncludeThetaBin(double theta_com);
TString MakeName(int pairID, TString header, TString tag);
TString FitET(double beamEnergy, int pairID);
double ThetaComToLab(double theta_com) { return (180. - theta_com)/2.; }
double ThetaLabToCom(double theta_lab) { return (180. - 2*theta_lab); }
double GetMaxValueAround(TH1D *histE, double mean1K, int binRange);
void InitFitParameters2(int energyIndex, int pairID);
///////////////////////////////////////////////////////////////////////////////////////////////////////////

double CalculatePairSolidAngle(double theta1, double theta2, int numDetectors);
double Calculate16RingSolidAngle(double theta1, double theta2, int numDetectors);
