//////////////////////////////////////////////////////////////////////////////////
void InitParameter(int run, TString dataDir="", bool isTest=false);
void StartWriteC0Parameters();
void StartWriteC1Parameters();
void StartWriteCEParameters();
void StartWriteCPParameters();
void StartWriteC2Parameters();
void StartWriteC3Parameters(TString fileName="");
void FillC0Parameters(int det, int side, int strip, double entries, double itcpt, double slope);
void FillC1Parameters(int det, int side, int strip, double entries, double itcptL, double slopeL, double itcptR, double slopeR);
void FillCEParameters(int det, int side, int strip, int    nPoints, double itcptL, double slopeL, double itcptR, double slopeR);
void FillCPParameters(int det, int side, int strip, double x[2][5]);
void FillC2Parameters(int det, int side, int strip, double entries, double b0, double b1, double b2);
void FillC3Parameters(int det, int side, int strip, double entries, double itcpt, double slope);
TString EndWriteParameters();
void GetC0Parameters(int calRun=-1);
void GetC1Parameters(int calRun=-1);
void GetC2Parameters(int calRun=-1);
void GetC3Parameters(int calRun=-1);
void GetC3Parameters(TString calDataName);
void GetCEParameters(int calRun=-1);
void GetCPParameters(int calRun=-1);
double CalibrateC0(int det, int side, int strip, double &energy);
double CalibrateC1(int det, int side, int strip, int lr, double &energy);
double CalibrateCE(int det, int side, int strip, int lr, double &energy);
double CalibrateC2(int det, int side, int strip, double &sum, double &pos);
double CalibrateC3(int det, int side, int strip, double &energy);

//////////////////////////////////////////////////////////////////////////////////
TString fC0ParDataName;
TString fC1ParDataName;
TString fCEParDataName;
TString fCPParDataName;
TString fC2ParDataName;
TString fC3ParDataName;
TString fGlobalParDataName;
double fC0Parameters[40][2][8][2];
double fC1Parameters[40][2][8][2][2];
double fCEParameters[40][2][8][2][2];
double fCPParameters[40][2][8][2][5];
double fC2Parameters[40][2][8][3];
double fC3Parameters[40][2][8][2];
TString fDataParDir;
ofstream fOStreamCalParameters;
ifstream fIStreamCalParameters;

//////////////////////////////////////////////////////////////////////////////////
void InitParameter(int run, TString dataDir, bool isTest)
{
    if (!dataDir.IsNull()) fDataParDir = dataDir;
    fC0ParDataName  = Form("%s/stark_%04d.c0.dat",fDataParDir.Data(),run);
    fC1ParDataName  = Form("%s/stark_%04d.c1.dat",fDataParDir.Data(),run);
    fCEParDataName  = Form("%s/stark_%04d.ce.dat",fDataParDir.Data(),run);
    fCPParDataName  = Form("%s/stark_%04d.cp.dat",fDataParDir.Data(),run);
    fC2ParDataName  = Form("%s/stark_%04d.c2.dat",fDataParDir.Data(),run);
    fC3ParDataName  = Form("%s/stark_%04d.c3.dat",fDataParDir.Data(),run);
    GetC0Parameters(run);
    GetC1Parameters(run);
    GetC2Parameters(run);
    GetCEParameters(run);
}

void StartWriteC0Parameters()
{
    fGlobalParDataName = fC0ParDataName;
    fOStreamCalParameters.open(fC0ParDataName);
    fOStreamCalParameters
        << "#det"   << "\t"
        << "side"   << "\t"
        << "strip"  << "\t"
        << "entries"<< "\t"
        << "itcpt"  << "\t"
        << "slope"  << endl;
}

void FillC0Parameters(int det, int side, int strip, double entries, double itcpt, double slope)
{
    fOStreamCalParameters
        << det     << "\t"
        << side    << "\t"
        << strip   << "\t"
        << entries << "\t"
        << itcpt   << "\t"
        << slope   << endl;
    fC0Parameters[det][side][strip][0] = itcpt;
    fC0Parameters[det][side][strip][1] = slope;
}

void StartWriteC1Parameters()
{
    fGlobalParDataName = fC1ParDataName;
    fOStreamCalParameters.open(fC1ParDataName);
    fOStreamCalParameters
        << "#det"   << "\t"
        << "side"   << "\t"
        << "strip"  << "\t"
        << "nPoints"<< "\t"
        << "itcptL" << "\t"
        << "slopeL" << "\t"
        << "itcptR" << "\t"
        << "slopeR" << endl;
}

void FillC1Parameters(int det, int side, int strip, int nPoints, double itcptL, double slopeL, double itcptR, double slopeR)
{
    fOStreamCalParameters
        << det     << "\t"
        << side    << "\t"
        << strip   << "\t"
        << nPoints << "\t"
        << itcptL  << "\t"
        << slopeL  << "\t"
        << itcptR  << "\t"
        << slopeR  << endl;
    fC1Parameters[det][side][strip][0][0] = itcptL;
    fC1Parameters[det][side][strip][0][1] = slopeL;
    fC1Parameters[det][side][strip][1][0] = itcptR;
    fC1Parameters[det][side][strip][1][1] = slopeR;
}

void StartWriteCPParameters()
{
    fGlobalParDataName = fCPParDataName;
    fOStreamCalParameters.open(fCPParDataName);
    fOStreamCalParameters
        << "#det"   << "\t"
        << "side"   << "\t"
        << "strip"  << "\t"
        << "x00" << "\t"
        << "x01" << "\t"
        << "x02" << "\t"
        << "x03" << "\t"
        << "x04" << "\t"
        << "x10" << "\t"
        << "x11" << "\t"
        << "x12" << "\t"
        << "x13" << "\t"
        << "x14" << endl;
}

void FillCPParameters(int det, int side, int strip, double x[2][5])
{
    fOStreamCalParameters
        << det     << "\t"
        << side    << "\t"
        << strip   << "\t"
        << x[0][0] << "\t"
        << x[0][1] << "\t"
        << x[0][2] << "\t"
        << x[0][3] << "\t"
        << x[0][4] << "\t"
        << x[1][0] << "\t"
        << x[1][1] << "\t"
        << x[1][2] << "\t"
        << x[1][3] << "\t"
        << x[1][4] << endl;
    fCPParameters[det][side][strip][0][0] = x[0][0];
    fCPParameters[det][side][strip][0][1] = x[0][1];
    fCPParameters[det][side][strip][0][2] = x[0][2];
    fCPParameters[det][side][strip][0][3] = x[0][3];
    fCPParameters[det][side][strip][0][4] = x[0][4];
    fCPParameters[det][side][strip][1][0] = x[1][0];
    fCPParameters[det][side][strip][1][1] = x[1][1];
    fCPParameters[det][side][strip][1][2] = x[1][2];
    fCPParameters[det][side][strip][1][3] = x[1][3];
    fCPParameters[det][side][strip][1][4] = x[1][4];
}

void StartWriteCEParameters()
{
    fGlobalParDataName = fCEParDataName;
    fOStreamCalParameters.open(fCEParDataName);
    fOStreamCalParameters
        << "#det"   << "\t"
        << "side"   << "\t"
        << "strip"  << "\t"
        << "nPoints"<< "\t"
        << "itcptL" << "\t"
        << "slopeL" << "\t"
        << "itcptR" << "\t"
        << "slopeR" << endl;
}

void FillCEParameters(int det, int side, int strip, int nPoints, double itcptL, double slopeL, double itcptR, double slopeR)
{
    fOStreamCalParameters
        << det     << "\t"
        << side    << "\t"
        << strip   << "\t"
        << nPoints << "\t"
        << itcptL  << "\t"
        << slopeL  << "\t"
        << itcptR  << "\t"
        << slopeR  << endl;
    fCEParameters[det][side][strip][0][0] = itcptL;
    fCEParameters[det][side][strip][0][1] = slopeL;
    fCEParameters[det][side][strip][1][0] = itcptR;
    fCEParameters[det][side][strip][1][1] = slopeR;
}

void StartWriteC2Parameters()
{
    fGlobalParDataName = fC2ParDataName;
    fOStreamCalParameters.open(fC2ParDataName);
    fOStreamCalParameters
        << "#det"    << "\t" 
        << "side"   << "\t" 
        << "strip"  << "\t" 
        << "entries"<< "\t" 
        << "b0"     << "\t" 
        << "b1"     << "\t" 
        << "b2"     << endl;
}

void FillC2Parameters(int det, int side, int strip, double entries, double b0, double b1, double b2)
{
    fOStreamCalParameters
        << det     << "\t" 
        << side    << "\t" 
        << strip   << "\t" 
        << entries << "\t" 
        << b0      << "\t" 
        << b1      << "\t" 
        << b2      << endl;
    fC2Parameters[det][side][strip][0] = b0;
    fC2Parameters[det][side][strip][1] = b1;
    fC2Parameters[det][side][strip][2] = b2;
}

void StartWriteC3Parameters(TString fileName)
{
    fGlobalParDataName = fileName;
    if (fGlobalParDataName.IsNull())
        fGlobalParDataName = fC3ParDataName;
    fOStreamCalParameters.open(fGlobalParDataName);
    fOStreamCalParameters
        << "#det"    << "\t" 
        << "side"    << "\t" 
        << "strip"   << "\t" 
        << "entries" << "\t" 
        << "itcpt"   << "\t" 
        << "slope"   << endl;
}

void FillC3Parameters(int det, int side, int strip, double entries, double itcpt, double slope)
{
    fOStreamCalParameters
        << det     << "\t" 
        << side    << "\t" 
        << strip   << "\t" 
        << entries << "\t" 
        << itcpt   << "\t" 
        << slope   << endl;
    fC3Parameters[det][side][strip][0] = itcpt;
    fC3Parameters[det][side][strip][1] = slope;
}

TString EndWriteParameters()
{
    fOStreamCalParameters.close();
    cout << fGlobalParDataName << endl;
    return fGlobalParDataName;
}

void GetC0Parameters(int calRun)
{
    TString calDataName = fC0ParDataName;
    if (calRun>=0)
        calDataName = Form("%s/stark_%04d.c0.dat" ,fDataParDir.Data(),calRun);
    int det, side, strip;
    double entries, itcpt, slope;
    e_info << "Using " << calDataName << endl;
    fIStreamCalParameters.open(calDataName);
    std::string dummyLine;
    std::getline(fIStreamCalParameters, dummyLine);
    while (fIStreamCalParameters >> det >> side >> strip >> entries >> itcpt >> slope)
    {
        fC0Parameters[det][side][strip][0] = itcpt;
        fC0Parameters[det][side][strip][1] = slope;
    }
    fIStreamCalParameters.close();
}

void GetC1Parameters(int calRun)
{
    TString calDataName = fC1ParDataName;
    if (calRun>=0)
        calDataName = Form("%s/stark_%04d.c1.dat" ,fDataParDir.Data(),calRun);
    e_info << "Using " << calDataName << endl;
    int det, side, strip;
    int nPoints;
    double slopeL, itcptL, slopeR, itcptR;
    fIStreamCalParameters.open(calDataName);
    std::string dummyLine;
    std::getline(fIStreamCalParameters, dummyLine);
    while (fIStreamCalParameters >> det >> side >> strip >> nPoints >> itcptL >> slopeL >> itcptR >> slopeR)
    {
        fC1Parameters[det][side][strip][0][0] = itcptL;
        fC1Parameters[det][side][strip][0][1] = slopeL;
        fC1Parameters[det][side][strip][1][0] = itcptR;
        fC1Parameters[det][side][strip][1][1] = slopeR;
    }
    fIStreamCalParameters.close();
}

void GetCPParameters(int calRun)
{
    TString calDataName = fCPParDataName;
    if (calRun>=0) {
        calDataName = Form("%s/stark_%04d.cp.dat" ,fDataParDir.Data(),calRun);
    }
    int det, side, strip;
    double x00, x01, x02, x03, x04, x10, x11, x12, x13, x14;
    e_info << "Using " << calDataName << endl;
    fIStreamCalParameters.open(calDataName);
    std::string dummyLine;
    std::getline(fIStreamCalParameters, dummyLine);
    while (fIStreamCalParameters >> det >> side >> strip >> x00 >> x01 >> x02 >> x03 >> x04 >> x10 >> x11 >> x12 >> x13 >> x14)
    {
        fCPParameters[det][side][strip][0][0] = x00;
        fCPParameters[det][side][strip][0][1] = x01;
        fCPParameters[det][side][strip][0][2] = x02;
        fCPParameters[det][side][strip][0][3] = x03;
        fCPParameters[det][side][strip][0][4] = x04;
        fCPParameters[det][side][strip][1][0] = x10;
        fCPParameters[det][side][strip][1][1] = x11;
        fCPParameters[det][side][strip][1][2] = x12;
        fCPParameters[det][side][strip][1][3] = x13;
        fCPParameters[det][side][strip][1][4] = x14;
    }
    fIStreamCalParameters.close();
}

void GetCEParameters(int calRun)
{
    TString calDataName = fCEParDataName;
    if (calRun>=0) {
        calDataName = Form("%s/stark_%04d.ce.dat" ,fDataParDir.Data(),calRun);
    }
    int det, side, strip;
    int nPoints;
    double slopeL, itcptL, slopeR, itcptR;
    e_info << "Using " << calDataName << endl;
    fIStreamCalParameters.open(calDataName);
    std::string dummyLine;
    std::getline(fIStreamCalParameters, dummyLine);
    while (fIStreamCalParameters >> det >> side >> strip >> nPoints >> itcptL >> slopeL >> itcptR >> slopeR)
    {
        fCEParameters[det][side][strip][0][0] = itcptL;
        fCEParameters[det][side][strip][0][1] = slopeL;
        fCEParameters[det][side][strip][1][0] = itcptR;
        fCEParameters[det][side][strip][1][1] = slopeR;
    }
    fIStreamCalParameters.close();
}

void GetC2Parameters(int calRun)
{
    TString calDataName = fC2ParDataName;
    if (calRun>=0)
        calDataName = Form("%s/stark_%04d.c2.dat" ,fDataParDir.Data(),calRun);
    int det, side, strip;
    double entries, b0, b1, b2;
    e_info << "Using " << calDataName << endl;
    fIStreamCalParameters.open(calDataName);
    std::string dummyLine;
    std::getline(fIStreamCalParameters, dummyLine);
    while (fIStreamCalParameters >> det >> side >> strip >> entries >> b0 >> b1 >> b2)
    {
        fC2Parameters[det][side][strip][0] = b0;
        fC2Parameters[det][side][strip][1] = b1;
        fC2Parameters[det][side][strip][2] = b2;
    }
    fIStreamCalParameters.close();
}

void GetC3Parameters(int calRun)
{
    TString calDataName = fC3ParDataName;
    if (calRun>=0) {
        calDataName = Form("%s/stark_%04d.c3.dat" ,fDataParDir.Data(),calRun);
    }
    GetC3Parameters(calDataName);
}

void GetC3Parameters(TString calDataName)
{
    if (calDataName.IsNull()) {
        e_error << "GetC3Parameters empty file name!!!" << endl;
        return;
    }
    int det, side, strip;
    double entries, itcpt, slope;
    {
        e_info << "Using " << calDataName << endl;
        fIStreamCalParameters.open(calDataName);
        std::string dummyLine;
        std::getline(fIStreamCalParameters, dummyLine);
        while (fIStreamCalParameters >> det >> side >> strip >> entries >> itcpt >> slope)
        {
            fC3Parameters[det][side][strip][0] = itcpt;
            fC3Parameters[det][side][strip][1] = slope;
        }
        fIStreamCalParameters.close();
    }
}

double CalibrateC0(int det, int side, int strip, double &energy)
{
    double itcpt = fC0Parameters[det][side][strip][0];
    double slope = fC0Parameters[det][side][strip][1];
    energy = slope * energy + itcpt;
    return energy;
}

double CalibrateC1(int det, int side, int strip, int lr, double &energy)
{
    double itcpt = fC1Parameters[det][side][strip][lr][0];
    double slope = fC1Parameters[det][side][strip][lr][1];
    energy = slope * energy + itcpt;
    return energy;
}

double CalibrateCE(int det, int side, int strip, int lr, double &energy)
{
    //double itcpt = fCEParameters[det][side][strip][lr][0];
    //double slope = fCEParameters[det][side][strip][lr][1];
    //if (slope>0.8&&slope<1.2)
    //    energy = slope * energy + itcpt;
    return energy;
}

double CalibrateC2(int det, int side, int strip, double &sum, double &pos)
{
    double b0 = fC2Parameters[det][side][strip][0];
    double b1 = fC2Parameters[det][side][strip][1];
    double b2 = fC2Parameters[det][side][strip][2];
    //auto sum = energyR + energyL;
    //auto pos = (energyR - energyL) / sum;
    sum = sum / (b0 + b1*pos + b2*pos*pos) * f241AmAlphaEnergy1;
    return sum;
}

double CalibrateC3(int det, int side, int strip, double &energy)
{
    double itcpt = fC3Parameters[det][side][strip][0];
    double slope = fC3Parameters[det][side][strip][1];
    energy = slope * energy + itcpt;
    return energy;
}
