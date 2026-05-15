#include "TObject.h"

class stark_energy_handler : public TObject
{
    public:
        stark_energy_handler(TString energyCalibrationName, TString positionCalibrationName="");
        virtual ~stark_energy_handler() {}

        double RestoreEnergy(int det, int side, int strip, double &energy);
        void RestoreEnergyRS(int det, int side, int strip, double energyL, double energyR, double &position, double &energySum); ///< for Resistive Strip
        double RestorePosition(int det, int side, int strip, double &position, double energySum);
        LKDrawing *RestorePositionAndDraw(int det, int side, int strip, double position, double energySum);

    private:
        double fC0Parameters[40][2][8][2];
        double fC1Parameters[40][2][8][2][2];
        double fCPParameters[40][2][8][2][5];
        double fC2Parameters[40][2][8][3];
        double fC3Parameters[40][2][8][2];
        const double fGate1Energy = 5.486;
        const double fGate0Energy = 3.1822;
};

stark_energy_handler::stark_energy_handler(TString energyCalibrationName, TString positionCalibrationName)
{
    int det, side, strip;
    string dummy;

    ifstream energyCalibrationFile(energyCalibrationName);
    if (energyCalibrationFile.is_open()==false) {
        e_error << "Cannot open " << energyCalibrationName << endl;
        return;
    }
    else {
        e_info << "Reading energy calibration file " << energyCalibrationName << endl;
    }

    double c0_itcpt, c0_slope, c1_itcptL, c1_slopeL, c1_itcptR, c1_slopeR, c2_par0, c2_par1, c2_par2, c3_itcpt, c3_slope;
    getline(energyCalibrationFile,dummy);
    while (energyCalibrationFile>>det>>side>>strip>>c0_itcpt>>c0_slope>>c1_itcptL>>c1_slopeL>>c1_itcptR>>c1_slopeR>>c2_par0>>c2_par1>>c2_par2>>c3_itcpt>>c3_slope)
    {
        fC0Parameters[det][side][strip][0] = c0_itcpt;
        fC0Parameters[det][side][strip][1] = c0_slope;
        fC1Parameters[det][side][strip][0][0] = c1_itcptL;
        fC1Parameters[det][side][strip][0][1] = c1_slopeL;
        fC1Parameters[det][side][strip][1][0] = c1_itcptR;
        fC1Parameters[det][side][strip][1][1] = c1_slopeR;
        fC2Parameters[det][side][strip][0] = c2_par0;
        fC2Parameters[det][side][strip][1] = c2_par1;
        fC2Parameters[det][side][strip][2] = c2_par2;
        fC3Parameters[det][side][strip][0] = c3_itcpt;
        fC3Parameters[det][side][strip][1] = c3_slope;
    }

    if (positionCalibrationName.IsNull())
        return;

    /////////////////////////////////////////////////////////////////////

    ifstream positionCalibrationFile(positionCalibrationName);
    if (positionCalibrationFile.is_open()==false) {
        e_error << "Cannot open " << positionCalibrationName << endl;
        return;
    }
    else {
        e_info << "Reading position calibration file " << positionCalibrationName << endl;
    }

    double x00 ,x01 ,x02 ,x03 ,x04 ,x10 ,x11 ,x12 ,x13 ,x14;
    getline(positionCalibrationFile,dummy);
    while (positionCalibrationFile>>det>>side>>strip >> x00 >> x01 >> x02 >> x03 >> x04 >> x10 >> x11 >> x12 >> x13 >> x14)
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
}

double stark_energy_handler::RestoreEnergy(int det, int side, int strip, double &energy)
{
    energy = fC0Parameters[det][side][strip][1] * energy + fC0Parameters[det][side][strip][0];
    return energy;
}

void stark_energy_handler::RestoreEnergyRS(int det, int side, int strip, double energyL, double energyR, double &position, double &energySum)
{
    energyL = fC1Parameters[det][side][strip][0][1] * energyL + fC1Parameters[det][side][strip][0][0];
    energyR = fC1Parameters[det][side][strip][1][1] * energyR + fC1Parameters[det][side][strip][1][0];
    //lk_debug << "c1lr " << energyL << " " << energyR << endl;
    energySum = energyL + energyR;
    position = (energyR - energyL) / energySum;
    //lk_debug << "c1pe " << position << " " << energySum << endl;
    energySum = energySum / (fC2Parameters[det][side][strip][0] + fC2Parameters[det][side][strip][1]*position + fC2Parameters[det][side][strip][2]*position*position) * fGate1Energy;
    //lk_debug << "c2pe " << position << " " << energySum << endl;
    energySum = fC3Parameters[det][side][strip][1] * energySum + fC3Parameters[det][side][strip][0];
    //lk_debug << "c3pe " << position << " " << energySum << endl;
    //position = RestorePosition(det, side, strip, position, energySum);
    //lk_debug << "pppp " << position << " " << energySum << endl;
}

LKDrawing *stark_energy_handler::RestorePositionAndDraw(int det, int side, int strip, double position, double energySum)
{
    LKDrawing* drawing = new LKDrawing();
    drawing -> Add(new TMarker(position,energySum,47));
    double x01 = fCPParameters[det][side][strip][0][1];
    double x02 = fCPParameters[det][side][strip][0][2];
    double x03 = fCPParameters[det][side][strip][0][3];
    double x11 = fCPParameters[det][side][strip][1][1];
    double x12 = fCPParameters[det][side][strip][1][2];
    double x13 = fCPParameters[det][side][strip][1][3];
    drawing -> Add(new TMarker(x01,fGate0Energy,21));
    drawing -> Add(new TMarker(x02,fGate0Energy,21));
    drawing -> Add(new TMarker(x03,fGate0Energy,21));
    drawing -> Add(new TMarker(x11,fGate1Energy,21));
    drawing -> Add(new TMarker(x12,fGate1Energy,21));
    drawing -> Add(new TMarker(x13,fGate1Energy,21));
    auto f1 = new TF1("f1x1",Form("%f * (x-%f) + %f",(fGate1Energy-fGate0Energy)/(x11-x01),x01,fGate0Energy),-1,1); f1 -> SetLineWidth(1);
    auto f2 = new TF1("f1x2",Form("%f * (x-%f) + %f",(fGate1Energy-fGate0Energy)/(x12-x02),x02,fGate0Energy),-1,1); f2 -> SetLineWidth(1);
    auto f3 = new TF1("f1x3",Form("%f * (x-%f) + %f",(fGate1Energy-fGate0Energy)/(x13-x03),x03,fGate0Energy),-1,1); f3 -> SetLineWidth(1);
    drawing -> Add(f1);
    drawing -> Add(f2);
    drawing -> Add(f3);
    double x1 = (x11-x01)/(fGate1Energy-fGate0Energy) * (energySum-fGate0Energy) + x01;
    double x2 = (x12-x02)/(fGate1Energy-fGate0Energy) * (energySum-fGate0Energy) + x02;
    double x3 = (x13-x03)/(fGate1Energy-fGate0Energy) * (energySum-fGate0Energy) + x03;
    drawing -> Add(new TMarker(x1,energySum,20));
    drawing -> Add(new TMarker(x2,energySum,20));
    drawing -> Add(new TMarker(x3,energySum,20));
    if (position>x2) position =  0.5/(x3-x2)*(position-x2);
    else             position = -0.5/(x1-x2)*(position-x2);
    drawing -> Add(new TParameter<double>("position",position));
    return drawing;
}

double stark_energy_handler::RestorePosition(int det, int side, int strip, double &position, double energySum)
{
    double x01 = fCPParameters[det][side][strip][0][1];
    double x02 = fCPParameters[det][side][strip][0][2];
    double x03 = fCPParameters[det][side][strip][0][3];
    double x11 = fCPParameters[det][side][strip][1][1];
    double x12 = fCPParameters[det][side][strip][1][2];
    double x13 = fCPParameters[det][side][strip][1][3];
    double x1 = (x11-x01)/(fGate1Energy-fGate0Energy) * (energySum-fGate0Energy) + x01;
    double x2 = (x12-x02)/(fGate1Energy-fGate0Energy) * (energySum-fGate0Energy) + x02;
    double x3 = (x13-x03)/(fGate1Energy-fGate0Energy) * (energySum-fGate0Energy) + x03;
    if (position>x2) position =  0.5/(x3-x2)*(position-x2);
    else             position = -0.5/(x1-x2)*(position-x2);
    return position;
}
