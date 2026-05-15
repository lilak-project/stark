#include "evaluate_rutherford_crosssection.C"

bool initialized_AdjustYValue = false;
double fExpEnergies[] = {4.41,  5.838, 8.13, 7.77, 9.36};
TGraphErrors *fGraphData[5];
TGraphErrors *fGraphRatio[5];

void GetArgonGraph()
{
	double ep;
	double angle;
	double diff_ov_ruth;
	double diff;
	double diff_ruth;

    fGraphData[3]  = new TGraphErrors();
    fGraphData[4]  = new TGraphErrors();
    fGraphRatio[3]  = new TGraphErrors();
    fGraphRatio[4]  = new TGraphErrors();

	std::ifstream file("data_x/argon_scattering_data.txt");
	std::string line;
	std::getline(file, line);
	while (file >> ep >> angle >> diff_ov_ruth >> diff >> diff_ruth)
	{
		int i = -1;
		for (auto j : {3,4}) { if (abs(fExpEnergies[j]-ep)<0.01) { i = j; break; } }
		if (i<0) continue;

		fGraphData[i]  -> SetPoint(fGraphData[i]->GetN(), angle, diff);
        double ratio = diff/EvalRutherford_DT(fExpEnergies[i], angle);
		fGraphRatio[i]  -> SetPoint(fGraphData[i]->GetN(), angle, ratio);
        //lk_debug << i << " " << fGraphData[i]->GetN()-1 << " " << angle << " " << ratio << endl;
	}
	file.close();
    initialized_AdjustYValue = true;
}

void AdjustYValue(int energyIndex, TGraphErrors* graphData, TGraphErrors* graphRatio, double &scaleFactor, double &xTarget, bool is99=false)
{
    if (initialized_AdjustYValue==false)
        GetArgonGraph();

    double energy = 4.41;
    if (energyIndex==4) energy = 4.41;
    else if (energyIndex==5) energy = 5.84;
    else if (energyIndex==8) energy = 8.13;

	// Target X value
	//double xTarget = 45.5;
	xTarget = 48;
    if (energyIndex==8) xTarget = 47;
    //if (energyIndex==8) xTarget = 55;

	// Get interpolated Y values
	double y3 = fGraphData[3]->Eval(xTarget); // Interpolated Y value from fGraphData[3]
	double y4 = fGraphData[4]->Eval(xTarget); // Interpolated Y value from fGraphData[4]
    //lk_debug << endl;
    //lk_debug << "y3=" << y3 << ",  y4=" << y4 << endl;

    bool usePrev = false;

    double ry4 = EvalRutherford_DT(9.36, xTarget);
    double ry3 = EvalRutherford_DT(7.77, xTarget);

	double ratioy4 = y4/ry4;
	double ratioy3 = y3/ry3;

	double a = (ratioy4-ratioy3)/(9.36-7.77);
	double b = ratioy3-a*7.77;
    //lk_debug << "a=" << a << ",  b=" << b << endl;
    //lk_debug << "ry3=" << ry3 << ",  ry4=" << ry4 << " " << EvalRutherford_DT(9.36, xTarget) << endl;
    //lk_debug << "ratio3=" << ratioy3 << ",  ratio4=" << ratioy4 << endl;

	double yRatioCentralOld, yCentralOld, yCentralNew, yRatioCentralNew;
    if (1)
        for (int i = 0; i < graphRatio->GetN(); i++)
        {
            double rnew = EvalRutherford_DT(energy, xTarget);
            yRatioCentralNew = (energy * a + b);
            double theta, theta2, y, yratio;
            graphRatio->GetPoint(i, theta2, yratio);
            //if (is99) lk_debug << theta2 << " " << xTarget << " y=" << y << " ratio=" << yratio << endl;
            if (fabs(theta2 - xTarget) < 1e-6) { // Check if X is close to the target
                yRatioCentralOld = yratio;
                break;
            }
        }
    else
        for (int i = 0; i < graphData->GetN(); i++)
        {
            double rnew = EvalRutherford_DT(energy, xTarget);
            yCentralNew = rnew*(energy * a + b);
            double theta, theta2, y, yratio;
            graphData->GetPoint(i, theta, y);
            //if (is99) lk_debug << theta2 << " " << xTarget << " y=" << y << " ratio=" << yratio << endl;
            if (fabs(theta2 - xTarget) < 1e-6) { // Check if X is close to the target
                yCentralOld = y;
                break;
            }
        }

    //scaleFactor = yCentralNew / yCentralOld;
    //lk_debug << xTarget << " " << scaleFactor << " new=" << yCentralNew << " old=" <<  yCentralOld << endl;

    //scaleFactor = yRatioCentralNew / yRatioCentralOld;
    scaleFactor = yRatioCentralNew / yRatioCentralOld;
    if (scaleFactor>10000||isinf(scaleFactor))
    {
        if (energyIndex==4) scaleFactor = 1.39;
        if (energyIndex==5) scaleFactor = 1.45;
        if (energyIndex==8) scaleFactor = 1.20;
    }
    if (is99) lk_debug << xTarget << " " << scaleFactor << " new=" << yRatioCentralNew << " old=" <<  yRatioCentralOld << endl;

    for (int i = 0; i < graphData->GetN(); i++) {
        double theta, y;
        graphData->GetPoint(i, theta, y);
        //if (fabs(theta - xTarget) < 1e-6) continue; // Skip the already updated point
        double yNew = y * scaleFactor;
        graphData->SetPoint(i, theta, yNew);
        //std::cout << "Scaled Y value at X = " << x << " from " << y << " to " << yNew << std::endl;
    }

    for (int i = 0; i < graphRatio->GetN(); i++) {
        double theta, y, ex, ey;
        graphRatio->GetPoint(i, theta, y);
        ex = graphRatio->GetErrorX(i);
        ey = graphRatio->GetErrorY(i);
        if (fabs(theta - xTarget) < 1e-6) {} // Skip the already updated point
        double yNew = y * scaleFactor;
        //if (is99) lk_debug << theta << " " << y << " * " << scaleFactor << " = " << yNew << endl;
        double eyNew = sqrt(ey*ey+yNew*0.075*yNew*0.075);
        graphRatio->SetPoint(i, theta, yNew);
        //if (fabs(theta - xTarget) < 1e-6) { lk_debug << "ratio " << yNew << " " << " " << endl; }
        if (fabs(theta - xTarget) < 1e-6) {
            double xx, yy;
            graphData->GetPoint(i, xx, yy);
            //lk_debug << "======== " << y << " " << yCentralOld/ry0 << endl;
        }
        graphRatio->SetPointError(i, ex, eyNew);
        //std::cout << "Ratio Scaled Y value at X = " << x << " from " << y << " to " << yNew  << "and " << ey << " to " << eyNew << std::endl;
    }

    //graphRatio->Print();
}
