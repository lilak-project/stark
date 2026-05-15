#include "set_draw_and_fit_cross_section.h"
#include "find_scaling_factor.C"
#include "functions.h"
#include "FitParInfo.h"

void draw_and_fit_cross_section(int energyIndex=8, int optionIndex=6, bool drawProjectionFits=true)
{
    if (Initialize(energyIndex, optionIndex)==false)
        return;
    return;

    auto groupX = fTop -> CreateGroup("crosssection",!fPickup);
    auto subXFinal = groupX -> CreateGroup("x_final");
    auto subXEachD = groupX -> CreateGroup("x_each");
    auto subREachD = groupX -> CreateGroup("r_each");
    auto subETEachD = groupX -> CreateGroup("et_each");

    if (fSelectBeamEnergy>0)
    {
        for (auto idArray : {fIDArray12, fIDArray16})
            for (auto pairID : idArray) {
                AnalyzeWithBeamEnergyAndPairID(fSelectBeamEnergy,pairID,drawProjectionFits);
                //break;
            }
    }
    else {
        double dBeamEnergy = 0.1;
        AnalyzeWithBeamEnergyAndPairID(fBeamEnergy);
        AnalyzeWithBeamEnergyAndPairID(fBeamEnergy+dBeamEnergy); fGraphQE -> Sort();
        AnalyzeWithBeamEnergyAndPairID(fGraphQE->Eval(0)); fGraphQE -> Sort();
        AnalyzeWithBeamEnergyAndPairID(fGraphQE->Eval(0)); fGraphQE -> Sort();
    }

    DrawXPair(99);
    for (auto idArray : {fIDArray12, fIDArray16}) {
        for (auto pairID : idArray) {
            DrawXPair(pairID);
            auto drawX = fTop -> FindDrawing(Form("drawX_%d",pairID));
            auto drawR = fTop -> FindDrawing(Form("drawR_%d",pairID));
            auto drawET = fTop -> FindDrawing(Form("et_%d",pairID));
            if (drawX!=nullptr) {
                auto drawX2 = subXEachD -> CreateDrawing();
                drawX -> CopyTo(drawX2);
            }
            if (drawR!=nullptr) {
                auto drawR2 = subREachD -> CreateDrawing();
                drawR -> CopyTo(drawR2);
            }
            if (drawET!=nullptr) {
                auto drawET2 = subETEachD -> CreateDrawing();
                drawET -> CopyTo(drawET2);
            }
        }
    }

    if (0)
        DrawWriteDrawings(true);
    else {
        DrawWriteDrawings(false);
        WriteParameters();
    }
}

LKDrawingGroup* AnalyzeWithBeamEnergyAndPairID(double beamEnergy, int pairID, bool addToTop)
{
    InitFitParameters2(fEnergyIndex, pairID);
    bool setPair = (pairID>=0);
    TString cutString = "fIsInGate";
    if (setPair)
        cutString = Form("fIsInGate&&fPairID==%d",pairID);

    TString mainTitle = Form("(%.3f MeV)%s",beamEnergy,(setPair?Form(" Pair-%d",pairID):""));
    TString nameHET  = MakeName(pairID,"hist","ET");
    TString nameHQT  = MakeName(pairID,"hist","QT");
    TString message  = TString("== ") + pairID + " ";

    LKDrawingGroup* groupAna = nullptr;
    TString groupName = Form("pair%d",pairID);
    groupAna = fTop -> CreateGroup(groupName,addToTop);
    auto sub_summ = groupAna;
    auto sub_fits = groupAna;
    if (fShowAll) {
        sub_summ = groupAna -> CreateGroup(groupName+"_summary",!fPickup);
        sub_fits = groupAna -> CreateGroup(groupName+"_fits");
        int dx, dy;
        LKPainter::GetPainter() -> GetSizeResize(dx, dy, 1200, 800, 0.8);
        sub_fits -> SetCanvasSize(dx,dy);
    }

    TString selectQT = Form("GetQValue(%f):180-fTheta*2>>%s",beamEnergy,nameHQT.Data());
    auto histQT = (TH2D*) fSave -> FindHist(nameHQT);
    if (histQT==nullptr||histQT==0) {
        e_info << "Filling " << nameHQT << endl;
        histQT = (bnnC*bnnQ).NewH2(nameHQT,mainTitle+";#theta_{CoM.} (deg);Q-value (MeV)");
        fTree -> Draw(selectQT,cutString,"goff");
    }
    auto drawQT = sub_summ -> CreateDrawing("qt");
    drawQT -> Add(histQT);
    if (histQT -> GetEntries()<fHistEntriesCut) {
        e_warning << message << endl;
        e_warning << "Skipping due to number of entries cut " << fHistEntriesCut << endl;
        return groupAna;
    }

    TString selectET = Form("fKeyEnergy:fTheta>>%s",nameHET.Data());
    auto histET = (TH2D*) fSave -> FindHist(nameHET);
    if (histET==nullptr||histET==0) {
        e_info << "Filling " << nameHET << endl;
        histET = (bnnL*bnnE).NewH2(nameHET,mainTitle+";#theta_{Lab.} (deg);E_{tot}");
        fTree -> Draw(selectET,cutString,"goff");
    }
    auto drawET = sub_summ -> CreateDrawing(Form("et_%d",pairID));
    drawET -> Add(histET);

    TString message2 = FitET(beamEnergy,pairID);

    e_info << message << endl;
    e_info << message2 << endl;

    return groupAna;
}


void DrawWriteDrawings(bool justDraw)
{
    bool write_top = true;
    bool write_short = true;
    bool write_each = true;
    bool write_x_canvas = true;
    bool write_graph_data = true;

    auto drawGGAR = fTop -> FindDrawing("ggar");
    //auto fit = new TF1("fit","pol1",0,180);
    //auto fit = new TF1("fit","(x<[0])*[1]*(x-[0])**2",0,20);
    //fit -> SetParameter(0,5);
    //fit -> SetParameter(1,1);
    auto fit = new TF1("fit","(x>=45)*((x<[2])*[3]*(x-[2])**2)+(x<45)*pol1(0)",0,180);
    fit -> SetParameter(2,80);
    fit -> SetParameter(3,1);
    fGraphParameterGGAR -> Fit(fit,"RQN0");
    drawGGAR -> Add(fit,"l");
    if (fFitGGAR==nullptr) {
        auto fileFitGGAR = new TFile(Form("%s/ggar_e%d.root",fDataXPath.Data(),fEnergyIndex),"recreate");
        fit -> Write();
    }

    auto subXFinal = fTop -> FindGroup("x_final");
    auto subXEachD = fTop -> FindGroup("x_each");
    auto subREachD = fTop -> FindGroup("r_each");
    auto subETEachD = fTop -> FindGroup("et_each");
    auto draw1 = fTop -> FindDrawing("draw1");
    //fTop -> Draw("v:wx=1200:wy=750");
    //fTop -> Draw("v:wx=1200:wy=500");
    if (fPickup) fTop -> Draw();
    else fTop -> Draw("v");
    if (justDraw) return;
    fTop -> WriteFitParameterFile(fTagBnn);

    if (subXFinal==nullptr) return;

    ///////////////////////////////////////////////////////////////////////////////////
    if (write_x_canvas && draw1!=nullptr) {
        auto cvs = draw1 -> GetCanvas();
        if (cvs!=nullptr) {
            cvs -> SaveAs(Form(fDataXPath+"/x_e%d.png",fEnergyIndex));
            cvs -> GetCanvas() -> SaveAs(Form(fDataXPath+"/x_e%d.eps",fEnergyIndex));
        }
    }

    ///////////////////////////////////////////////////////////////////////////////////
    auto graphData = subXFinal -> FindGraph("expdata_99");
    auto graphDataScaled = subXFinal -> FindGraph("expdata_scaled_99");
    auto graphDOR = subXFinal -> FindGraph("expdata_ov_rutherford_99");
    auto graphSOR = subXFinal -> FindGraph("expdata_scaled_ov_rutherford_99");
    auto graphFOR = subXFinal -> FindGraph("fresco_ov_rutherford_99");
    if (write_graph_data) {
        graphData -> SaveAs(Form(fDataXPath+"/cross_section_e%d.txt",fEnergyIndex));
        graphDataScaled -> SaveAs(Form(fDataXPath+"/cross_section_scaled_e%d.txt",fEnergyIndex));
        graphDOR -> SaveAs(Form(fDataXPath+"/cross_section_over_rutherford_e%d.txt",fEnergyIndex));
        graphSOR -> SaveAs(Form(fDataXPath+"/cross_section_scaled_over_rutherford_e%d.txt",fEnergyIndex));
        graphFOR -> SaveAs(Form(fDataXPath+"/fresco_kd_over_rutherford_e%d.txt",fEnergyIndex));
    }

    ///////////////////////////////////////////////////////////////////////////////////
    if (write_top) {
        //e_info << "Writting " << endl;
        //e_cout << "    " << fNameFileX << endl;
        //auto file = new TFile(fNameFileX,"recreate");
        //auto file = new TFile(fNameFileX,"recreate");
        //fTop -> Write(".flat");
        fTop -> WriteFile();

        //if (0)
        for (auto sub : {subXFinal, subXEachD, subREachD, subETEachD}) {
            auto cvs = sub -> GetCanvas();
            if (cvs!=nullptr) {
                sub -> Draw();
                cvs -> Write(sub -> GetName());
                graphData -> Write();
                fGraphRF -> Write();
                fGraphFC[0] -> Write();
                graphDOR -> Write();
                graphFOR -> Write();
                cvs -> SaveAs(Form(fDataXPath+"/figure_x%d.%s.eps",fEnergyIndex,TString(sub->GetName()).Data()));
                cvs -> SaveAs(Form(fDataXPath+"/figure_x%d.%s.png",fEnergyIndex,TString(sub->GetName()).Data()));
            }
        }

        //auto cvs = subXFinal -> GetCanvas();
        //if (cvs!=nullptr) {
        //    cvs -> Write("cvs");
        //    graphData -> Write();
        //    fGraphRF -> Write();
        //    fGraphFC[0] -> Write();
        //    graphDOR -> Write();
        //    graphFOR -> Write();
        //    cvs -> SaveAs(Form(fDataXPath+"/figure_x%d.eps",fEnergyIndex));
        //    cvs -> SaveAs(Form(fDataXPath+"/figure_x%d.png",fEnergyIndex));
        //}
    }

    ///////////////////////////////////////////////////////////////////////////////////
    if (write_short) {
        TString nameFileXShort = Form(fDataXPath+"/x%d.short.root",fEnergyIndex);
        cout << nameFileXShort << endl;
        auto fileShort = new TFile(nameFileXShort,"recreate");
        for (TString name : {"expdata_99","expdata_scaled_99","fresco_0","rutherford","expdata_ov_rutherford_99","expdata_scaled_ov_rutherford_99","fresco_ov_rutherford_99"}) {
            auto obj = subXFinal -> FindGraph(name);
            name.ReplaceAll("_99","");
            name = name + "_e" + fEnergyIndex;
            if (obj!=nullptr)
                obj -> Write(name);
        }
    }

    ///////////////////////////////////////////////////////////////////////////////////
    if (write_each) {
        auto subXEachD = fTop -> FindGroup("x_each");
        if (subXEachD!=nullptr)
        {
            subXEachD -> Draw();
            auto cvsXAllDet = subXEachD -> GetCanvas();
            if (cvsXAllDet!=nullptr)
                cvsXAllDet -> SaveAs(fDataXPath+Form("/x_each_detectors_e%d.png",fEnergyIndex));
        }
    }
}

void WriteParameters()
{
    for (auto ring : {16,12})
    {
        auto backgroundType = fBackgroundType16;
        auto idArray = fIDArray16;
        if (ring==12) {
            backgroundType = fBackgroundType12;
            idArray = fIDArray12;
        }

        TString name_fit_paramter = Form(fDataXPath+"/fit_parameters_%s.ring%d.%s.dat",fAnaName.Data(),ring,((backgroundType==0)?"linear_bg":"exp_bg"));
        e_info << "Writting " << endl;
        e_cout << "    " << name_fit_paramter << endl;
        fFileFitParameterOut.open(name_fit_paramter);

        if (backgroundType==0)
        {
            for (auto idArray : {idArray})
            {
                for (auto pairID : idArray)
                {
                    fFileFitParameterOut << endl;
                    for (auto pjbin : fProjectionBins)
                    {
                        int    mxrbin = fHistParameters[pairID][pjbin][0];
                        double cntylm = fHistParameters[pairID][pjbin][1]*0.01;
                        int    pickup = fHistParameters[pairID][pjbin][2];
                        int    remove = fHistParameters[pairID][pjbin][3];
                        if (fFitParameters[pairID][pjbin][0][0]>=0) {
                            fFileFitParameterOut << "========== " << pairID << " " << pjbin << " " << cntylm << " " << mxrbin << endl;
                            fFileFitParameterOut << pickup << " " << remove << " " << fFitRange[pairID][pjbin][0] << " " << fFitRange[pairID][pjbin][1] << endl;
                            fFileFitParameterOut << "amplitude  " << fFitParameters[pairID][pjbin][0][0] << " " << fFitParameters[pairID][pjbin][0][1] << " " << fFitParameters[pairID][pjbin][0][2] << endl;
                            fFileFitParameterOut << "mean       " << fFitParameters[pairID][pjbin][1][0] << " " << fFitParameters[pairID][pjbin][1][1] << " " << fFitParameters[pairID][pjbin][1][2] << endl;
                            fFileFitParameterOut << "sigma      " << fFitParameters[pairID][pjbin][2][0] << " " << fFitParameters[pairID][pjbin][2][1] << " " << fFitParameters[pairID][pjbin][2][2] << endl;
                            fFileFitParameterOut << "y-itercept " << fFitParameters[pairID][pjbin][3][0] << " " << fFitParameters[pairID][pjbin][3][1] << " " << fFitParameters[pairID][pjbin][3][2] << endl;
                            fFileFitParameterOut << "slope      " << fFitParameters[pairID][pjbin][4][0] << " " << fFitParameters[pairID][pjbin][4][1] << " " << fFitParameters[pairID][pjbin][4][2] << endl;
                            fFileFitParameterOut << "q-Boundary " << fFitParameters[pairID][pjbin][5][0] << " " << fFitParameters[pairID][pjbin][5][1] << " " << fFitParameters[pairID][pjbin][5][2] << endl;
                        }
                    }
                }
            }
        }
        else if (backgroundType==1)
        {
            for (auto idArray : {idArray})
            {
                for (auto pairID : idArray)
                {
                    fFileFitParameterOut << endl;
                    for (auto pjbin : fProjectionBins)
                    {
                        int    mxrbin = fHistParameters[pairID][pjbin][0];
                        double cntylm = fHistParameters[pairID][pjbin][1]*0.01;
                        int    pickup = fHistParameters[pairID][pjbin][2];
                        int    remove = fHistParameters[pairID][pjbin][3];
                        if (fFitParameters[pairID][pjbin][1][0]>=0) {
                            fFileFitParameterOut << "========== " << pairID << " " << pjbin << " " << cntylm << " " << mxrbin << endl;
                            fFileFitParameterOut << pickup << " " << remove << " " << fFitRange[pairID][pjbin][0] << " " << fFitRange[pairID][pjbin][1] << endl;
                            fFileFitParameterOut << "q-Boundary " << fFitParameters[pairID][pjbin][0][0] << " " << fFitParameters[pairID][pjbin][0][1] << " " << fFitParameters[pairID][pjbin][0][2] << endl;
                            fFileFitParameterOut << "amplitude  " << fFitParameters[pairID][pjbin][1][0] << " " << fFitParameters[pairID][pjbin][1][1] << " " << fFitParameters[pairID][pjbin][1][2] << endl;
                            fFileFitParameterOut << "mean       " << fFitParameters[pairID][pjbin][2][0] << " " << fFitParameters[pairID][pjbin][2][1] << " " << fFitParameters[pairID][pjbin][2][2] << endl;
                            fFileFitParameterOut << "sigma      " << fFitParameters[pairID][pjbin][3][0] << " " << fFitParameters[pairID][pjbin][3][1] << " " << fFitParameters[pairID][pjbin][3][2] << endl;
                            fFileFitParameterOut << "expb_amp   " << fFitParameters[pairID][pjbin][4][0] << " " << fFitParameters[pairID][pjbin][4][1] << " " << fFitParameters[pairID][pjbin][4][2] << endl;
                            fFileFitParameterOut << "expb_amp   " << fFitParameters[pairID][pjbin][5][0] << " " << fFitParameters[pairID][pjbin][5][1] << " " << fFitParameters[pairID][pjbin][5][2] << endl;
                            fFileFitParameterOut << "const_bg   " << fFitParameters[pairID][pjbin][6][0] << " " << fFitParameters[pairID][pjbin][6][1] << " " << fFitParameters[pairID][pjbin][6][2] << endl;
                        }
                    }
                }
            }
        }
        fFileFitParameterOut.close();
    }

    fParFile -> cd();
    fParTree -> Write();
    e_info << "Writting " << endl;
    e_cout << "    " << fParFile->GetName() << endl;
}


LKDrawingGroup* DrawCrossSection(int selPair)
{
    bnnT = ((selPair<12)?bnnSplitX12Com:bnnSplitX16Com);
    auto nn = bnnSplitX12Com.n()+bnnSplitX16Com.n();
    int dx = 0.1*(bnnSplitX12Com.x2()-bnnSplitX16Com.x1());
    LKBinning bnnTAll(nn,bnnSplitX16Com.x1()-dx,bnnSplitX12Com.x2()+dx);
    LKBinning bnnT100(100,bnnSplitX16Com.x1()-dx,bnnSplitX12Com.x2()+dx);
    if (selPair==99)
        bnnT = bnnTAll;
    auto arrayValue = bnnTAll.MakeArrayD();
    auto arrayError = bnnTAll.MakeArrayD();

    LKDrawingGroup* groupDCS = new LKDrawingGroup();
    LKDrawing *drawDCS = groupDCS -> CreateDrawing("dcs");
    auto graphFinalX0ForScaling = NewGraph(Form("expdata_fscaling_%d",selPair));
    auto graphData = NewGraph(Form("expdata_%d",selPair));
    auto graphFinalX1 = NewGraph(Form("expdata_woSolidAngle_%d",selPair)); 
    auto graphSASum = NewGraph(Form("solid_angle_sum_%d",selPair));
    auto histCountA = bnnTAll.NewH1(Form("histCountA_%d",selPair),";#theta_{com};raw count");
    auto histCount = bnnT.NewH1(Form("histCount_%d",selPair),";#theta_{com};raw count");
    auto histCount12 = bnnSplitX12Com.NewH1(Form("histCount12_%d",selPair),";#theta_{com};raw count");
    auto histCount16 = bnnSplitX16Com.NewH1(Form("histCount16_%d",selPair),";#theta_{com};raw count");
    drawDCS -> Add(graphFinalX0ForScaling);
    drawDCS -> Add(graphData);
    drawDCS -> Add(histCountA);
    drawDCS -> Add(histCount);
    drawDCS -> Add(histCount12);
    drawDCS -> Add(histCount16);
    drawDCS -> Add(graphSASum);

    vector<TString> titlesER = {"total", "count_stat", "count_stat1", "count_stat2", "edge_syst", "range_syst", "fit_syst", "solid"};
    TGraphErrors* graphER[10] = {0};
    for (auto i=0; i<titlesER.size(); ++i) {
        graphER[i] = NewGraph(Form("graphER%d",i));
        graphER[i] -> SetMarkerStyle(20+i);
        graphER[i] -> SetMarkerColor(i+1);
        graphER[i] -> SetLineColor(i+1);
    }

    vector<int> totalIDArray;
    for (auto ring : {16,12})
    {
        vector<int> idArray;
        if (selPair>=0 && selPair<99) {
            if (ring==12) continue;
            idArray.push_back(selPair);
        }
        else {
            if (ring==16) idArray = fIDArray16;
            if (ring==12) idArray = fIDArray12;
        }
        for (auto id : idArray)
            totalIDArray.push_back(id);

        bnnT.Reset();
        while (bnnT.Next())
        {
            auto ttaIndex = bnnT.GetItIndex();
            //auto theta_com = bnnT.GetItCenter();
            if (selPair!=99 && selPair>=12) ttaIndex += bnnSplitX12.n();
            double total_signal = 0;
            double total_background = 0;
            double total_count_edge_error = 0;
            double total_count_fit_error = 0;
            double total_count_rng_error = 0;
            double total_rwcnt = 0;
            double total_solid = 0;
            double total_solid1 = 0;
            double total_solid2 = 0;
            //double total_solid3 = 0;
            //double total_solid4 = 0;
            double theta_com  = 0;
            double theta1_com = 0;
            for (auto pairID : idArray)
            {
                if (theta_com ==0) theta_com  = fDataRange[pairID][ttaIndex][0];
                if (theta1_com==0) theta1_com = fDataRange[pairID][ttaIndex][1];
            }
            //if (((theta_com>45 && theta_com<80) || (theta_com>100 && theta_com<120))==false) continue;
            double targetDensity = fNumAtomsPerUnitArea;
            double beamCount     = fBeamCount;
            //double beamCountError = 0.1*fBeamCount;
            for (auto pairID : idArray)
            {
                double count            = fDataValues[pairID][ttaIndex][0][0];
                double count_bg         = fDataValues[pairID][ttaIndex][0][1];
                double count_fit_error  = fDataValues[pairID][ttaIndex][0][2];
                double count_rng_error  = fDataValues[pairID][ttaIndex][0][3];
                double raw              = fDataValues[pairID][ttaIndex][0][4];
                //
                double solid_angle0 = fSolidAngleRatio*fDataValues[pairID][ttaIndex][3][0];
                double solid_angle1 = fSolidAngleRatio*fDataValues[pairID][ttaIndex][3][1];
                double solid_angle2 = fSolidAngleRatio*fDataValues[pairID][ttaIndex][3][2];
                total_signal += count;
                total_background += count_bg;
                total_rwcnt += raw;
                total_solid += solid_angle0;
                total_solid1 += solid_angle1;
                total_solid2 += solid_angle2;
                //total_solid3 += solid_angle3;
                //total_solid4 += solid_angle4;
                total_count_fit_error += count_fit_error;
                total_count_rng_error += count_rng_error;
            }
            double solidAngleError1 = abs(total_solid1 - total_solid);
            double solidAngleError2 = abs(total_solid2 - total_solid);
            //double solidAngleError3 = abs(total_solid3 - total_solid);
            //double solidAngleError4 = abs(total_solid4 - total_solid);

            //if (fEnergyIndex==5 && abs(theta_com-102)<0.1) total_count_edge_error = total_signal*0.2;
            //if (fEnergyIndex==4 && (abs(theta_com-99)<0.1 || abs(theta_com-99)<117)) total_count_edge_error = total_signal*0.2;

            auto total_countX = total_signal/total_solid/targetDensity/beamCount/fmbarn;
            if (ShouldIncludeThetaBin(theta_com)==false) continue;
            /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            if (total_countX>0) {
                double ex = abs(theta_com-theta1_com);

                double sqrt_count = sqrt(total_signal + total_background);
                double sqrt_count_error = sqrt(TMath::Power(sqrt(total_count_rng_error),2) + TMath::Power(sqrt(total_count_edge_error),2));

                double sqrt_count1 = sqrt(total_signal);
                double sqrt_count2 = sqrt(total_background);
                double sqrt_count_error1 = sqrt(total_count_edge_error*total_count_edge_error);
                double sqrt_count_error2 = sqrt(total_count_rng_error*total_count_rng_error);
                double sqrt_count_error3 = sqrt(total_count_fit_error*total_count_fit_error);
                double sqrt_solidAngle = sqrt(
                        +solidAngleError1*solidAngleError1
                        +solidAngleError2*solidAngleError2
                        //+solidAngleError3*solidAngleError3
                        //+solidAngleError4*solidAngleError4
                        );
                //  Partial derivatives
                double dValue_dCount = 1.0 / (total_solid * targetDensity * beamCount * fmbarn);
                double dValue_dBeamCount = -total_signal / (total_solid * targetDensity * TMath::Power(beamCount, 2) * fmbarn);
                double dValue_dSolidAngle = -total_signal / (TMath::Power(total_solid, 2) * targetDensity * beamCount * fmbarn);
                // Error propagation
                double errorSquared0 = 0;
                errorSquared0 += TMath::Power(dValue_dCount      * sqrt_count,      2);
                errorSquared0 += TMath::Power(dValue_dCount      * sqrt_count_error,2);
                errorSquared0 += TMath::Power(dValue_dSolidAngle * sqrt_solidAngle, 2);
                double ey0 = sqrt(errorSquared0);

                double error_percentage[10] = {0};
                error_percentage[0] = 100*ey0/total_countX;
                error_percentage[1] = 100*sqrt(TMath::Power(dValue_dCount      * sqrt_count,      2))/total_countX;
                error_percentage[2] = 100*sqrt(TMath::Power(dValue_dCount      * sqrt_count1,      2))/total_countX;
                error_percentage[3] = 100*sqrt(TMath::Power(dValue_dCount      * sqrt_count2,      2))/total_countX;
                error_percentage[4] = 100*sqrt(TMath::Power(dValue_dCount      * sqrt_count_error1,2))/total_countX;
                error_percentage[5] = 100*sqrt(TMath::Power(dValue_dCount      * sqrt_count_error2,2))/total_countX;
                error_percentage[6] = 100*sqrt(TMath::Power(dValue_dCount      * sqrt_count_error3,2))/total_countX;
                error_percentage[7] = 100*sqrt(TMath::Power(dValue_dSolidAngle * sqrt_solidAngle,  2))/total_countX;
                if (selPair==99) lk_debug << theta_com << endl;
                for (auto i=0; i<titlesER.size(); ++i) {
                    if (selPair==99) lk_debug << titlesER[i] << ": " << " ratio=" << error_percentage[i] << endl;
                    graphER[i] -> SetPoint(graphER[i]->GetN(),theta_com,error_percentage[i]);
                }
                ////////////////////////////////////////////////////////////////////////////
                double errorSquared1 = 0;
                errorSquared1 += TMath::Power(dValue_dCount     * sqrt_count,      2);
                double ey1 = sqrt(errorSquared1);
                ////////////////////////////////////////////////////////////////////////////
                graphData -> SetPoint(graphData->GetN(),theta_com,total_countX);
                graphData -> SetPointError(graphData->GetN()-1,ex,ey0);
                if (theta_com>=fThetaCoMScaleRange1&&theta_com<fThetaCoMScaleRange2) {
                    graphFinalX0ForScaling -> SetPoint     (graphFinalX0ForScaling->GetN(),theta_com,total_countX);
                    graphFinalX0ForScaling -> SetPointError(graphFinalX0ForScaling->GetN()-1,ex,ey0);
                }
                graphFinalX1 -> SetPoint(graphFinalX1->GetN(),theta_com,total_countX);
                graphFinalX1 -> SetPointError(graphFinalX1->GetN()-1,ex,ey1);
                histCount -> Fill(theta_com,total_rwcnt);
                histCount12 -> Fill(theta_com,total_rwcnt);
                histCount16 -> Fill(theta_com,total_rwcnt);
                if (arrayValue!=nullptr) arrayValue -> SetAt(total_countX,ttaIndex);
                if (arrayError!=nullptr) arrayError -> SetAt(ey0,ttaIndex);
                graphSASum -> SetPoint(graphSASum->GetN(),theta_com,total_solid);
                graphSASum -> SetPointError(graphSASum->GetN()-1,ex,sqrt_solidAngle*sqrt_solidAngle);
            }
            /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        }
    }

    TF1* fitFC0 = new TF1(Form("fitFC0_%d",selPair), ScaleFunction1, 0, 180, 1);
    double scaleFC = 1;
    double ymax = 400;
    if (fEnergyIndex==4) ymax = 1200;
    if (fEnergyIndex==5) ymax = 700;
    if (selPair<12) ymax = 100;
    if (fUseSecondPeak) {
        ymax = ymax/10;
        if (fEnergyIndex==5) ymax = 8;
        if (fEnergyIndex==8) ymax = 20;
    }
    FitAmplitude1To2(1,graphFinalX0ForScaling,fGraphFC[0],fitFC0,scaleFC);
    double scaleDataToFC = 1./scaleFC;
    if (fGlobalBeamNormalizationFactor>0)
        scaleDataToFC = fGlobalBeamNormalizationFactor;
    if (isinf(scaleDataToFC))
    {
        if (fEnergyIndex==4) scaleDataToFC = 1.39;
        if (fEnergyIndex==5) scaleDataToFC = 1.45;
        if (fEnergyIndex==8) scaleDataToFC = 1.20;
    }

    auto graphDataScaled = NewGraph(Form("expdata_scaled_%d",selPair));

    auto numPointsX0 = graphData -> GetN();
    double x,y,ex,ey;
    for (auto iPoint=0; iPoint<numPointsX0; ++iPoint)
    {
        graphData -> GetPoint(iPoint,x,y);
        if (y<=0) continue;
        ex = graphData -> GetErrorX(iPoint);
        ey = graphData -> GetErrorY(iPoint);
        graphDataScaled -> SetPoint(iPoint,x,y*scaleDataToFC);
        graphDataScaled -> SetPointError(iPoint,ex,ey*scaleDataToFC);
    }
    //auto frameX = (bnnT*LKBinning(100,0,ymax)).NewH2(Form("frameX_%d",selPair),Form("%d;",selPair)+fCrossSectionTitle);
    auto frameX = (bnnTAll*LKBinning(100,0,ymax)).NewH2(Form("frameX_%d",selPair),fMainTitle+(selPair==99?"":Form(" (pair %d)",selPair))+";"+fCrossSectionTitle);

    graphFinalX0ForScaling -> SetMarkerStyle(24);
    graphFinalX0ForScaling -> SetMarkerSize(1.3);
    graphData -> SetFillColor(kYellow);
    graphData -> SetFillStyle(1001);
    graphDataScaled -> SetFillColor(kCyan-10);
    graphDataScaled -> SetFillStyle(3002);

    drawDCS -> Add(frameX);
    if (graphDataScaled->GetN()>0)
        drawDCS -> Add(graphDataScaled);
    drawDCS -> Add(new TParameter<double>(Form("scaleDataToFC_%d",selPair),scaleDataToFC));
    drawDCS -> Add(new TParameter<double>(Form("scaleFC_%d",selPair),scaleFC));

    auto graphDOR = NewGraph(Form("expdata_ov_rutherford_%d",selPair),20,0.8,kOrange-3);
    auto graphSOR = NewGraph(Form("expdata_scaled_ov_rutherford_%d",selPair),25,1,kBlue);
    auto graphFOR = NewGraph(Form("fresco_ov_rutherford_%d",selPair),25,1,kRed);
    graphFOR -> SetName(Form("fresco_ov_rutherford_%d",selPair));
    graphFOR -> SetLineColor(kMagenta);
    auto frameDOR = (bnnT*LKBinning(100,0,2)).NewH2(Form("frameDOR_%d",selPair),";#theta_{CoM.} (deg);Data/Rutherford  (Fresco/Rutherford)");
    frameDOR -> SetStats(0);

    //for (auto ring : {12:16})
    {
        bnnT.Reset();
        while (bnnT.Next())
        {
            auto ttaIndex = bnnT.GetItIndex();
            if (selPair!=99 && selPair>=12) ttaIndex += bnnSplitX12.n();
            double theta_com = 0;
            LKBinning bnnT2;
            for (auto pairID : totalIDArray) {
                if (theta_com==0) {
                    theta_com = fDataRange[pairID][ttaIndex][0];
                    if (theta_com!=0) {
                        bnnT2 = ((selPair<12)?bnnSplitX12Com:bnnSplitX16Com);
                        break;
                    }
                }
            }
            if (arrayValue->At(ttaIndex)>0) {
                double vDOR = arrayValue->At(ttaIndex)/fGraphRF->Eval(theta_com);
                double fSOR = graphDataScaled->Eval(theta_com)/fGraphRF->Eval(theta_com);
                double eDOR = arrayError->At(ttaIndex)/fGraphRF->Eval(theta_com);
                double eSOR = arrayError->At(ttaIndex)*scaleDataToFC/fGraphRF->Eval(theta_com);
                graphDOR -> SetPoint(graphDOR->GetN(),theta_com,vDOR);
                graphDOR -> SetPointError(graphDOR->GetN()-1,0.5*bnnT2.GetWX(),eDOR);
                graphSOR -> SetPoint(graphSOR->GetN(),theta_com,fSOR);
                graphSOR -> SetPointError(graphSOR->GetN()-1,0.5*bnnT2.GetWX(),eSOR);
            }
            //double vFOR = fGraphFC[0]->Eval(theta_com)/fGraphRF->Eval(theta_com);
            //graphFOR -> SetPoint(graphFOR->GetN(),theta_com,vFOR);
        }
    }

    bnnT100.Reset();
    while (bnnT100.Next())
    {
        double theta_com = bnnT100.GetItCenter();
        //double vFOR = fGraphFC[0]->Eval(theta_com)/fGraphRF->Eval(theta_com);
        double energy = 4.41;
        if (fEnergyIndex==4) energy = 4.41;
        else if (fEnergyIndex==5) energy = 5.84;
        else if (fEnergyIndex==8) energy = 8.13;
        double vFOR = fGraphFC[0]->Eval(theta_com)/EvalRutherford_DT(energy, theta_com);
        graphFOR -> SetPoint(graphFOR->GetN(),theta_com,vFOR);
    }

    drawDCS -> Add(frameDOR);
    drawDCS -> Add(graphDOR);
    drawDCS -> Add(graphSOR);
    drawDCS -> Add(graphFOR);

    sort(totalIDArray.begin(),totalIDArray.end());
    TString pairs = "det: ";
    auto id1 = totalIDArray[0];
    auto id2 = totalIDArray[0];
    auto n = totalIDArray.size();
    for (auto i=1; i<n; ++i)
    {
        auto id = totalIDArray[i];
        if (id>id2+1) {
            if (id2==id1) pairs = pairs + Form("%d,",id1);
            else pairs = pairs + Form("%d-%d,",id1,id2);
            id1 = id;
        }
        id2 = id;
    }
    if (id2==id1) pairs = pairs + Form("%d,",id1);
    else pairs = pairs + Form("%d-%d,",id1,id2);
    pairs = pairs(0,pairs.Sizeof()-2);

    TString bgtype = "";
    bgtype = bgtype + (fBackgroundType12==0?"12:lin-bg, ":"12:exp-bg, ");
    bgtype = bgtype + (fBackgroundType16==0?"16:lin-bg"  :"16:expbg");

    TString scalingTitle = Form("Data*%.2f (%.0f<#theta<%.0f)",scaleDataToFC,fThetaCoMScaleRange1,fThetaCoMScaleRange2);
    if (fUseScalingWithPrevData) {
        double thetaPick = 0;;
        AdjustYValue(fEnergyIndex, graphDataScaled, graphSOR, scaleDataToFC, thetaPick, (selPair==99));
        scalingTitle = Form("Data*%.2f (#theta=%.0f)",scaleDataToFC,thetaPick);
    }

    auto drawX = groupDCS -> CreateDrawing(Form("x_%d",selPair));
    if (fUseSecondPeak) {
        drawX -> Add(frameX,"stat0",".");
        drawX -> Add(graphDataScaled,"p",scalingTitle);
        drawX -> Add(graphData,"drawx","Data");
        drawX -> SetCreateLegend(((selPair<12)?1:1),0.4,0.06);
        drawX -> SetGridx();
        drawX -> SetGridy();
    }
    else {
        drawX -> Add(frameX,"stat0",".");
        drawX -> Add(graphData,"same5->pz","Data");
        drawX -> Add(graphDataScaled,"same5->pz",scalingTitle);
        //if (selPair==99) graphDataScaled->Print();
        drawX -> AddLegendLine(Form("#sigma_{Target} = %.2e",fNumAtomsPerUnitArea));
        drawX -> AddLegendLine(Form("E_{Beam} = %.3f",fBeamEnergy));
        drawX -> AddLegendLine(Form("#Beam = %.3e",fBeamCount));
        drawX -> AddLegendLine(Form("#atoms/A=%.2e",fNumAtomsPerUnitArea));
        drawX -> AddLegendLine(Form("#Omega ratio= %.2f",fSolidAngleRatio));
        drawX -> AddLegendLine(pairs);
        drawX -> AddLegendLine(bgtype);
        drawX -> Add(fGraphFC[0],"samel","Fresco KD");
        drawX -> Add(fGraphRF,"samel","Rutherford");
        drawX -> Add(graphFinalX0ForScaling,"p",".");
        drawX -> SetCreateLegend((selPair<12)?1:0,0.4,0.06);
        drawX -> SetGridx();
        drawX -> SetGridy();
    }

    auto drawR = groupDCS -> CreateDrawing(Form("xr_%d",selPair));
    drawR -> Add(frameDOR,"",".");
    drawR -> Add(graphDOR,"p","Data/Rutherford");
    drawR -> Add(graphSOR,"p","Data_Scaled/Rutherford");
    drawR -> Add(graphFOR,"l","Fresco_KD/Rutherford");
    //drawR -> SetCreateLegend(1,0.6,0.06);
    //drawR -> SetCreateLegend(2,0.4,0.05);
    drawR -> SetCreateLegend((selPair<12)?1:2,0.4,0.05);
    drawR -> SetGridx();
    drawR -> SetGridy();

    auto drawER = groupDCS -> CreateDrawing(Form("er_%d",selPair));
    drawER -> SetCreateFrame(Form("frame_er_%d",selPair),"","y0");
    drawER -> SetCreateLegend();
    for (auto i=0; i<titlesER.size(); ++i) {
        if (graphER[i]->GetN()>0) {
            graphER[i] -> Sort();
            if (drawER -> GetEntries()==0)
                drawER -> Add(graphER[i],"pl",titlesER[i]);
            else
                drawER -> Add(graphER[i],"pl",titlesER[i]);
        }
    }

    return groupDCS;
}

bool ShouldIncludeThetaBin(double theta_com)
{
    if ((theta_com>fFinalThetaRange[0][0]&&theta_com<fFinalThetaRange[0][1]) || (theta_com>fFinalThetaRange[1][0]&&theta_com<fFinalThetaRange[1][1]))
        return true;
    return false;
}

TString FitET(double beamEnergy, int pairID)
{
    bool setPair = (pairID>=0);
    TString cutString = "fIsInGate";
    if (setPair)
        cutString = Form("fIsInGate&&fPairID==%d",pairID);

    bool backgroundType = (pairID<12?fBackgroundType12:fBackgroundType16);

    TString mainTitle = Form("(%.3f MeV)%s",beamEnergy,(setPair?Form(" Pair-%d",pairID):""));
    TString nameHET  = MakeName(pairID,"hist","ET");
    TString nameHQT  = MakeName(pairID,"hist","QT");
    TString nameGET  = MakeName(pairID,"graph","ET");
    TString nameGET2  = MakeName(pairID,"graph","ET2");
    TString nameGSR  = MakeName(pairID,"graph","SR");
    TString nameGBD  = MakeName(pairID,"graph","BD");
    TString nameGPV  = MakeName(pairID,"graph","PV");
    TString nameHPV  = MakeName(pairID,"frame","PV");
    TString nameFBD  = MakeName(pairID,"fitB","BD");
    TString namedEE1 = MakeName(pairID,"hist","dEE1");
    TString namedEE2 = MakeName(pairID,"hist","dEE2");
    TString nameHET2 = MakeName(pairID,"hist","ET2");
    TString groupName = Form("pair%d",pairID);
    TString message;

    vector<int> projectionBinsInUse;
    if (fPickup)
    {
        bool pickupExist = false;
        TString pickup_numbers;
        for (auto pjbin : fProjectionBins) {
            int pickup = fHistParameters[pairID][pjbin][2];
            if (pickup>=1) {
                pickup_numbers = pickup_numbers + pjbin + " ";
                pickupExist = true;
                projectionBinsInUse.push_back(pjbin);
            }
            if (pickup==2)
                fHistParameters[pairID][pjbin][2] = 0;
        }
        if (pickupExist==false) {
            e_info << message << "pickup not found in this pair. Skipping..." << endl;
            return "";
        }
        e_info << "Found pickup " << pickup_numbers << endl;
    }
    else {
        for (auto pjbin : fProjectionBins)
            projectionBinsInUse.push_back(pjbin);
    }

    auto groupAna = fTop -> FindGroup(groupName);
    auto sub_summ = groupAna -> FindGroup(groupName+"_summary");
    auto sub_fits = groupAna -> FindGroup(groupName+"_fits");
    //auto sub_parv = groupAna -> CreateGroup(groupName+"_parv",false);
    auto sub_parv = sub_fits;
    auto drawET = sub_summ -> FindDrawing(Form("et_%d",pairID));

    auto histQT = sub_summ -> FindHist(nameHQT);
    auto histET = sub_summ -> FindHist(nameHET);
    auto entriesQT = histQT -> GetEntries();

    auto graphET = NewGraph(nameGET,20,0.6,kBlack,1,1,kRed);
    auto graphET2 = NewGraph(nameGET2,20,0.6,kBlack,1,1,kRed);
    auto graphQTBoundary = NewGraph(nameGBD,0,0,1,1,2,kGreen);
    auto fitQTBoundary = new TF1(nameFBD,"pol1",bnnC.x1(),80);
    fitQTBoundary -> SetLineColor(kGreen);
    fitQTBoundary -> SetLineWidth(1);

    auto histdEE1 = (bnnE*bnnE).NewH2(namedEE1,mainTitle+Form(" X6-X6  %s;E_{tot};dE",(setPair?Form("(Pair-%d)",pairID):"")));
    auto histdEE2 = (bnnE*bnnE).NewH2(namedEE2,mainTitle+Form(" X6-CSD %s;E_{tot};dE",(setPair?Form("(Pair-%d)",pairID):"")));
    TString selectdEE1 = Form("fdEOhmic:fKeyEnergy>>%s",namedEE1.Data());
    TString selectdEE2 = Form("fdEJunction:fKeyEnergy>>%s",namedEE2.Data());
    TString cutStringdEE1 = "fIsInGate && fPairID<4";
    TString cutStringdEE2 = "fIsInGate && (fPairID>=4 && fPairID<12)";
    if (fSave!=nullptr) {
        auto histdEE1 = (TH2D*) fSave -> FindHist(namedEE1);
        auto histdEE2 = (TH2D*) fSave -> FindHist(namedEE2);
    }
    else if (setPair==false) {
        fTree -> Draw(selectdEE1,cutStringdEE1,"goff");
        fTree -> Draw(selectdEE2,cutStringdEE2,"goff");
    }
    else if (pairID>=0&&pairID<4) {
        selectdEE1 = Form("fdEOhmic:fKeyEnergy>>%s",namedEE1.Data());
        cutStringdEE1 = Form("fIsInGate&&fPairID==%d",pairID);
        fTree -> Draw(selectdEE1,cutStringdEE1,"goff");
    }
    else if (pairID>=4&&pairID<12) {
        selectdEE2 = Form("fdEJunction:fKeyEnergy>>%s",namedEE2.Data());
        cutStringdEE2 = Form("fIsInGate&&fPairID==%d",pairID);
        fTree -> Draw(selectdEE2,cutStringdEE2,"goff");
    }
    auto lg1 = new TLegend(); lg1 -> AddEntry(histdEE1,cutStringdEE1);
    auto lg2 = new TLegend(); lg2 -> AddEntry(histdEE2,cutStringdEE2);
    auto drawdEE1 = sub_summ -> CreateDrawing("dee1",false);
    auto drawdEE2 = sub_summ -> CreateDrawing("dee2",false);
    drawdEE1 -> Add(histdEE1);
    drawdEE2 -> Add(histdEE2);
    drawdEE1 -> Add(lg1);
    drawdEE2 -> Add(lg2);
    drawdEE1 -> AddOption("legend_dx",0.7);
    drawdEE2 -> AddOption("legend_dx",0.7);
    drawdEE1 -> SetLegendCorner(1);
    drawdEE2 -> SetLegendCorner(1);
    drawdEE1 -> SetLogz();
    drawdEE2 -> SetLogz();

    auto graphSignalRatio = NewGraph(nameGSR);

    double meanAll, sigmaAll;
    int bin0 = bnnC.FindBin(40 );
    int bin1 = bnnC.FindBin(84 );
    int bin2 = bnnC.FindBin(96 );
    int bin3 = bnnC.FindBin(120);
    double continuityThresholdRatio = 0.03;

    LKBinning bnnET(histET);

    auto histET2 = (TH2D*) histET -> Clone(nameHET2);

    drawET -> Add(graphET,"p5");
    drawET -> Add(graphET2,"p5");
    drawET -> SetOptStat(10);
    drawET -> SetStatBottomLeftCorner();

    bnnT = ((pairID<12)?bnnSplitX12:bnnSplitX16);
    bnnET.SetProjectionBinningValues(bnnT);
    //auto drawGrid = bnnET.CreateProjectionYGrid();
    //drawET -> Add(drawGrid);
    if (fGraphK[0]!=nullptr) drawET -> Add(fGraphK[0]);
    if (fGraphK[1]!=nullptr) drawET -> Add(fGraphK[1]);

    auto arrayX = bnnT.MakeArrayD();
    auto arrayXOff0 = bnnT.MakeArrayD();
    auto arrayXOff1 = bnnT.MakeArrayD();
    auto arrayXOff2 = bnnT.MakeArrayD();
    auto arrayXOff3 = bnnT.MakeArrayD();
    auto arrayX11 = bnnT.MakeArrayD();
    auto arrayXRaw = bnnT.MakeArrayD();
    auto arrayX1 = bnnT.MakeArrayD();
    auto arrayX2 = bnnT.MakeArrayD();

    double fitRange1 = bnnET.y1();
    double fitRange2 = 15;
    if (fEnergyIndex==4 && pairID<12) { fitRange1 = 4; fitRange2 = 15; }
    if (fEnergyIndex==5 && pairID<12) { fitRange1 = 4; fitRange2 = 18; }
    if (fEnergyIndex==8 && pairID<12) { fitRange1 = 15; fitRange2 = 25; }

    //TString formula = "(x>[0])*(exp(-(x-[1]))*(x>[0]?TMath::Power(x-[0],[2]):0)+[3]) + gaus(4) + [7]*exp(-0.5*((x-[8])/[6])**2)";
    //TString formula1 = "gaus(0)";
    //TString formula2 = "gaus(0)";
    //TString formula3 = "(x>[0])*(exp(-(x-[1]))*(x>[0]?TMath::Power(x-[0],[2]):0)+[3])";
    //TString formula_short = "(x>[0])*((x-[0])^{[2]}*e^{-(x-[1])}+[3]) + gaus(4)";

    //FitParInfo init_parameter[10];
    //init_parameter[0].SetNameValue("v0_x_low"  ,    2,  1.5,   2.2,  false);
    //init_parameter[1].SetNameValue("v1_exp_pos",    5,    0,    10,  false);
    //init_parameter[2].SetNameValue("v2_x_pow"  ,    1,  0.5,     4,  false);
    //init_parameter[3].SetNameValue("v3_const"  ,    4,    0,    20,  false);
    //init_parameter[4].SetNameValue("v4_gaus0_a", 1000,    0, 1.E+4,  false);
    //init_parameter[5].SetNameValue("v5_gaus0_m",   10,    0,    20,  false);
    //init_parameter[6].SetNameValue("v6_gaus0_s",  0.5,    0,     1,  false);
    //init_parameter[7].SetNameValue("v7_gaus1_a", 1000,    0, 1.E+4,  false);
    //init_parameter[8].SetNameValue("v8_gaus1_m",   10,    0,    20,  false);

    FitParInfo dummy_parameter;

    TGraphAsymmErrors* graphParameter[2][fNumParameters+1];
    for (auto i=0; i<fNumParameters+1; ++i) {
        graphParameter[0][i] = new TGraphAsymmErrors();
        graphParameter[0][i] -> SetName(nameGPV+i+0);
        graphParameter[0][i] -> SetMarkerStyle(25);
        graphParameter[0][i] -> SetMarkerColor(kRed);
        graphParameter[0][i] -> SetMarkerSize(1.2);
        graphParameter[1][i] = new TGraphAsymmErrors();
        graphParameter[1][i] -> SetName(nameGPV+i+1);
        graphParameter[1][i] -> SetMarkerStyle(20);
    }

    auto drawGGAR = sub_parv -> CreateDrawing("ggar");
    drawGGAR -> Add(fGraphParameterGGAR,"p");
    drawGGAR -> SetCreateFrame(Form("frame_fitet_%d",pairID),"ggar","y0");

    //for (auto i : {1,2,3,6})
    //for (auto i : {0,1,2,3,6,4,7})
    //for (auto i : {0,1,2,3,4,5,6,7,8})
    for (auto i=0; i<fNumParameters+1; ++i)
    {
        auto drawGP = sub_parv -> CreateDrawing();
        drawGP -> Add(graphParameter[0][i],"lp");
        drawGP -> Add(graphParameter[1][i],"lp");
        TString title = Form("%s",fit_par_energy[i].fName.Data());
        if (i==fNumParameters) title = "GGAR";
        drawGP -> SetCreateFrame(Form("frame_gp_%d_%d",pairID,i),title,"y0");
    }

    for (auto trial : {0,1}) // trial=0 for pre-fit and collecting background parameter, trial=1 for main fitting
    {
        double parameterAverage[10] = {0};
        int countProjection = 0;
        bnnET.ResetNextProjection();
        while (auto histE = bnnET.NextProjectionY(histET2))
        {
            histE -> GetXaxis() -> SetRangeUser(0,20);
            histE -> SetFillColor(18);
            histE -> SetLineColor(kGray+1);

            bool projAll = false;
            LKDrawing* drawQ = nullptr;
            int pjbin = bnnET.GetCurrentProjectionIt();
            double pjCenter = bnnET.GetCurrentProjectionCenter();
            double pjBinWidth = bnnET.GetCurrentProjectionBinWidth();

            //int countRebin = 0;
            //double continuity = LKSAM::GetSAM()->ContinuityIndex(histE,continuityThresholdRatio);
            //double cntylm = 0.6;
            int mxrbin = 3;
            //if (cntylm==0) cntylm = 0.6;
            //while (continuity<cntylm && countRebin<mxrbin) {
            //    countRebin++;
            //    histE -> RebinX(2);
            //    continuity = LKSAM::GetSAM()->ContinuityIndex(histE,continuityThresholdRatio);
            //}
            histE -> RebinX(2);
            histE -> RebinX(2);
            int    countRebin = 2;
            double continuity = LKSAM::GetSAM()->ContinuityIndex(histE,continuityThresholdRatio);

            FitParInfo fit_par_pair[10];
            for (auto i=0; i<fNumParameters; ++i) { fit_par_pair[i].CopyFrom(fit_par_energy[i]); }

            auto binMax = histE -> GetMaximumBin();
            double mean0K = 5;
            double mean1K = 5;
            if (fGraphK[0]!=nullptr) mean0K = fGraphK[0] -> Eval(pjCenter);
            if (fGraphK[1]!=nullptr) mean1K = fGraphK[1] -> Eval(pjCenter);

            double amp = histE -> GetBinContent(binMax);
            fit_par_pair[4].SetValue(amp,amp*0.5,amp*2);
            fit_par_pair[5].SetValue(mean0K,mean0K-0.5,mean0K+1);
            fit_par_pair[7].SetValue(0.1*amp,0,0.8*amp);
            if (fEnergyIndex==4 && pairID>=12) {
                fit_par_pair[7].FixValue(0);
                fit_par_pair[8].FixValue(0);
            }
            else if (mean1K<2.5) {
                fit_par_pair[7].FixValue(0);
                fit_par_pair[8].FixValue(0);
            }
            else {
                double amp1 = GetMaxValueAround(histE,mean1K,2);
                fit_par_pair[7].SetValue(amp1,0,amp1);
                fit_par_pair[8].SetValue(mean1K,mean1K-0.5,mean1K+0.5);
            }

            TString nameFit = Form("fitE_%s",histE->GetName());
            auto fit = new TF1(nameFit,fFormula,fitRange1,fitRange2);
            fit -> SetLineWidth(1);
            fit -> SetNpx(1000);
            for (auto i=0; i<fNumParameters; ++i) fit_par_pair[i].SetPar(fit,i);

            TF1* fitSave = nullptr;
            if (fSaveFit!=nullptr) fitSave = (TF1*) fSaveFit -> Get(nameFit);
            bool usingSaveFit = false;
            if (fitSave!=nullptr) usingSaveFit = true;
            if (usingSaveFit)
            {
                if (trial==0) break;
                //fitSave -> SetName(nameFit+2);
                double fit_range1, fit_range2;
                fitSave -> GetRange(fit_range1, fit_range2);
                fit -> SetRange(fit_range1, fit_range2);
                if (TString(fitSave->GetExpFormula())==TString(fit->GetExpFormula()))
                {
                    for (auto i=0; i<fNumParameters; ++i) {
                        dummy_parameter.PullPar(fitSave,i);
                        dummy_parameter.SetPar(fit,i);
                    }
                }
            }

            histE -> Fit(fit,"Q0BR","");
            if (fFitGGAR!=nullptr) {
                auto amp1 = fit -> GetParameter(4) * fFitGGAR -> Eval(pjCenter);
                if (amp1<0.001) amp1 = 0;
                if (amp1<0.001) {
                    fit -> FixParameter(7, amp1);
                }
                else {
                    fit -> SetParameter(7, amp1);
                    fit -> SetParLimits(7, 0.2*amp1, 2*amp1);
                    //fit -> SetParLimits(7,0,20);
                }
                histE -> Fit(fit,"Q0BR","");
            }

            for (auto i=0; i<fNumParameters; ++i) fit_par_pair[i].PullPar(fit,i);

            double ggar = fit_par_pair[7].fValue/fit_par_pair[4].fValue;
            //if (trial==1 && pairID>=12) fGraphParameterGGAR -> SetPoint(fGraphParameterGGAR->GetN(),pjbin,ggar);
            if (trial==1) fGraphParameterGGAR -> SetPoint(fGraphParameterGGAR->GetN(),pjCenter,ggar);
            for (auto i=0; i<fNumParameters+1; ++i)
            {
                auto iPoint = graphParameter[trial][i] -> GetN();
                if (i==fNumParameters) {
                    graphParameter[trial][i] -> SetPoint(iPoint,pjbin,ggar);
                }
                else {
                    graphParameter[trial][i] -> SetPoint(iPoint,pjbin,fit_par_pair[i].fValue);
                    //if (i==1 || i==2) lk_debug << trial << " " << mean0K << " " << pjbin << " " << fit_par_pair[i].fValue << endl;
                    graphParameter[trial][i] -> SetPointEYlow (iPoint, abs(fit_par_pair[i].fValue-fit_par_pair[i].fLimit1));
                    graphParameter[trial][i] -> SetPointEYhigh(iPoint, abs(fit_par_pair[i].fValue-fit_par_pair[i].fLimit2));
                }
            }

            if (trial==0) {
                for (auto i=0; i<fNumParameters; ++i)
                    parameterAverage[i] += fit -> GetParameter(i);
                //parameterAverage[1] += fit -> GetParameter(1);
                //parameterAverage[2] += fit -> GetParameter(2);
                //parameterAverage[3] += fit -> GetParameter(3);
                //parameterAverage[6] += fit -> GetParameter(6);
                countProjection++;
                if (mean0K<6)
                //if (mean0K<5)
                    break;
                continue;
            }

            auto fit1 = new TF1(nameFit+"_guas",fFormula1,fitRange1,fitRange2);
            for (auto i : {0,1,2}) fit_par_pair[i+4].SetPar(fit1,i);
            fit1 -> SetLineWidth(1);
            fit1 -> SetLineColor(kBlue);
            fit1 -> SetNpx(1000);
            auto fit2 = new TF1(nameFit+"_gaus2",fFormula2,fitRange1,fitRange2);
            for (auto i : {0,1,2}) {
                if (i==2) fit_par_pair[6].SetPar(fit2,i);
                else fit_par_pair[i+7].SetPar(fit2,i);
            }
            fit2 -> SetLineWidth(1);
            fit2 -> SetLineColor(kMagenta);
            fit2 -> SetNpx(1000);
            auto fit3 = new TF1(nameFit+"_bg",fFormula3,fitRange1,fitRange2);
            for (auto i : {0,1,2,3}) fit_par_pair[i].SetPar(fit3,i);
            fit3 -> SetLineWidth(1);
            fit3 -> SetLineColor(kGreen);
            fit3 -> SetNpx(1000);

            TF1 *fitSignal = fit1;
            if (fUseSecondPeak) fitSignal = fit2;

            TF1 *fitError[2][2] = {0};
            for (auto i : {0,1})
            {
                TF1 *fitClone = fit;
                if (i==1) fitClone = fitSignal;
                for (int j : {0,1})
                {
                    fitError[i][j] = (TF1*) fitClone -> Clone();
                    auto numParameters = fitClone -> GetNpar();
                    for (auto iPar=0; iPar<numParameters; ++iPar) {
                        auto error = fitClone -> GetParError(iPar);
                        if (j==0) fitError[i][j] -> SetParameter(iPar,fitClone->GetParameter(iPar)-error);
                        if (j==1) fitError[i][j] -> SetParameter(iPar,fitClone->GetParameter(iPar)+error);
                    }
                }
            }

            double meanFit0 = fit_par_pair[5].fValue;
            double meanFit1 = fit_par_pair[8].fValue;
            double sigmaPercent = fit_par_pair[6].fValue;
            double sigmaFit0 = meanFit0*0.01*sigmaPercent;
            double sigmaFit1 = meanFit1*0.01*sigmaPercent;

            double meanSignal = fitSignal -> GetParameter(1);
            double sigmaSignal = meanSignal*0.01*fitSignal -> GetParameter(2);

            graphET -> SetPoint(graphET->GetN(),pjCenter,meanFit0);
            graphET -> SetPointError(graphET->GetN()-1,0.5*pjBinWidth,fQProjTestSigmaFactor*sigmaFit0);

            graphET2 -> SetPoint(graphET2->GetN(),pjCenter,meanFit1);
            graphET2 -> SetPointError(graphET2->GetN()-1,0.5*pjBinWidth,fQProjTestSigmaFactor*sigmaFit1);

            double qSigRange1 = meanSignal-fQProjSigmaFactor*sigmaSignal;
            double qSigRange2 = meanSignal+fQProjSigmaFactor*sigmaSignal;
            double qSigRgBin1 = histE -> GetXaxis() -> FindBin(qSigRange1);
            double qSigRgBin2 = histE -> GetXaxis() -> FindBin(qSigRange2);
            double total = fit -> Integral(qSigRange1,qSigRange2);
            double signal = fitSignal -> Integral(qSigRange1,qSigRange2);

            double totalOff0 = fitError[0][0] -> Integral(qSigRange1,qSigRange2);
            double signalOff0 = fitError[1][0] -> Integral(qSigRange1,qSigRange2);

            double totalOff1 = fitError[0][1] -> Integral(qSigRange1,qSigRange2);
            double signalOff1 = fitError[1][1] -> Integral(qSigRange1,qSigRange2);

            double qSigRange1Off2 = fit_par_pair[5].fValue-fQProjSigmaFactor2*sigmaSignal;
            double qSigRange2Off2 = fit_par_pair[5].fValue+fQProjSigmaFactor2*sigmaSignal;
            double totalOff2 = fit -> Integral(qSigRange1Off2,qSigRange2Off2);
            double signalOff2 = fitSignal -> Integral(qSigRange1Off2,qSigRange2Off2);

            double qSigRange1Off3 = fit_par_pair[5].fValue-fQProjSigmaFactor3*sigmaSignal;
            double qSigRange2Off3 = fit_par_pair[5].fValue+fQProjSigmaFactor3*sigmaSignal;
            double totalOff3 = fit -> Integral(qSigRange1Off3,qSigRange2Off3);
            double signalOff3 = fitSignal -> Integral(qSigRange1Off3,qSigRange2Off3);

            auto hmax = histE -> GetMaximum();
            auto lineQ1 = new TLine(fitRange1,0,fitRange1,hmax*1.05); lineQ1 -> SetLineStyle(1);
            auto lineQ2 = new TLine(fitRange2,0,fitRange2,hmax*1.05); lineQ2 -> SetLineStyle(1);
            auto lineQ3 = new TLine(qSigRange1,0,qSigRange1,amp*0.5); lineQ3 -> SetLineStyle(2); lineQ3 -> SetLineColor(kBlue);
            auto lineQ4 = new TLine(qSigRange2,0,qSigRange2,amp*0.5); lineQ4 -> SetLineStyle(2); lineQ4 -> SetLineColor(kBlue);
            //auto lineQC = new TLine(mean,0,mean,amp); lineQC -> SetLineColor(kRed); lineQC -> SetLineStyle(2);
            auto mk0 = new TMarker(mean0K,0,23);
            auto mk1 = new TMarker(mean1K,0,23);

            double valueRaw = 0;
            double signalRatio = 0;
            double signalRatio1 = 0;
            double signalRatio2 = 0;
            if (projAll) {}
            else {
                if (total==0) e_warning << "entries, signal, total : " << entriesQT << ", " << signal << ", " << total << endl;
                signalRatio = signal/total;
                double signalRatioOff0 = signalOff0/totalOff0;
                double signalRatioOff1 = signalOff1/totalOff1;
                double signalRatioOff2 = signalOff2/totalOff2;
                double signalRatioOff3 = signalOff3/totalOff3;
                graphSignalRatio -> SetPoint(graphSignalRatio->GetN(),pjCenter,signalRatio);
                if (total<=0||signalRatio<fSignalRatioCut) {
                    e_warning << pairID << ", " << pjCenter << ": signal ratio is " << signalRatio << " = " << signal << " / " << total << endl;
                    //fit1 -> Print("value");
                    signalRatio = 0;
                }
                bool goodData = (signalRatio>fSignalRatioCut);
                if (fEnergyIndex==4&&(pairID==16||pairID==17||pairID==18||pairID==19)) {
                    double qedge = fQTEdge -> Eval(pjCenter);
                    if (qedge>qSigRange1)
                        goodData = false;
                }
                if (goodData) {
                    valueRaw = histE -> Integral(qSigRgBin1,qSigRgBin2);
                    auto valueRaw11 = histE -> Integral(-1,1);
                    auto value = signalRatio*valueRaw;
                    auto valueOff0 = signalRatioOff0*valueRaw;
                    auto valueOff1 = signalRatioOff1*valueRaw;
                    auto valueOff2 = signalRatioOff2*valueRaw;
                    auto valueOff3 = signalRatioOff3*valueRaw;
                    auto error = (1-signalRatio)*histE -> Integral(qSigRgBin1,qSigRgBin2);
                    int idx = pjbin;
                    arrayX -> SetAt(value,idx);
                    arrayXOff0 -> SetAt(valueOff0,idx);
                    arrayXOff1 -> SetAt(valueOff1,idx);
                    arrayXOff2 -> SetAt(valueOff2,idx);
                    arrayXOff3 -> SetAt(valueOff3,idx);
                    arrayX11 -> SetAt(valueRaw11,idx);
                    arrayXRaw -> SetAt(valueRaw,idx);
                    arrayX1 -> SetAt(qSigRgBin1,idx);
                    arrayX2 -> SetAt(qSigRgBin2,idx);
                }
            }

            if (projAll) {;}
            else if (fShowAll) drawQ = sub_fits -> CreateDrawing(Form("q_%d_%d_%d",pairID,pjbin,trial));
            if (drawQ!=nullptr) {
                drawQ -> Add(histE);
                drawQ -> Add(fit3);
                drawQ -> Add(fit2);
                drawQ -> Add(fit1);
                drawQ -> Add(fit);
                drawQ -> SetFitObjects(histE,fit);
                //drawQ -> Add(lineQ1);
                //drawQ -> Add(lineQ2);
                drawQ -> Add(lineQ3);
                drawQ -> Add(lineQ4);
                drawQ -> Add(mk0);
                drawQ -> Add(mk1);
                //drawQ -> Add(lineQC);
                drawQ -> AddOption("pave_dx",0.350);
                //drawQ -> AddOption("pave_line_dy",0.070);
                drawQ -> AddOption("pave_line_dy",0.040);
                drawQ -> SetOptFit(1);
                auto lg = new TLegend();
                lg -> AddEntry(histE,Form("#Rebin(<%d) = %d",mxrbin,countRebin));
                lg -> AddEntry(histE,Form("Continuity = %.2f",continuity));
                lg -> AddEntry(histE,Form("Sig. Ratio = %.2f",signalRatio));
                lg -> AddEntry(histE,Form("signal = %.2f",signal/total*valueRaw));
                lg -> AddEntry(histE,Form("inrange = %.2f",valueRaw));
                //lg -> AddEntry(histE,Form("Signal ratio = %.2f",signalRatio));
                drawQ -> Add(lg);
                if (pairID<12&&fEnergyIndex==8)
                    drawQ -> SetStatCorner(1);
                drawQ -> SetLegendBelowStats();
                //drawQ -> SetLegendCorner(0);
            }
        } // end of projection

        if (trial==0 && countProjection>0) {
            for (auto i=0; i<fNumParameters; ++i)
                parameterAverage[i] = parameterAverage[i]/countProjection;
            //fit_par_energy[0].SetValue(parameterAverage[0],0,fit_par_energy[0].fLimit2);
            //fit_par_energy[1].SetValue(parameterAverage[1],0,fit_par_energy[1].fLimit2);
            //fit_par_energy[2].FixValue(parameterAverage[2]);
            //fit_par_energy[2].SetValue(parameterAverage[2],0,parameterAverage[2]*2);
            //fit_par_energy[3].SetValue(parameterAverage[3],0,parameterAverage[3]*5);
            //fit_par_energy[6].FixValue(parameterAverage[6]);
            //fit_par_energy[6].SetValue(parameterAverage[6],0.5*parameterAverage[6],parameterAverage[6]);
            //if (fEnergyIndex==4 && pairID<12){
            //    fit_par_energy[1].FixValue(parameterAverage[1]);
            //    fit_par_energy[2].FixValue(parameterAverage[2]);
            //}
            //else {
                fit_par_energy[0].FixValue(parameterAverage[0]);
                //fit_par_energy[1].FixValue(parameterAverage[1]);
                //fit_par_energy[1].SetValue(parameterAverage[1],0.5*parameterAverage[1],2*parameterAverage[1]);
                double p11 = parameterAverage[1]-1; if (p11<0) p11 = 0;
                double p12 = parameterAverage[1]+1;
                fit_par_energy[1].SetValue(parameterAverage[1],p11,p12);
                fit_par_energy[2].FixValue(parameterAverage[2]);
                //fit_par_energy[3].FixValue(parameterAverage[3]);
                //fit_par_energy[2].SetValue(parameterAverage[2],0.5*parameterAverage[2],2*parameterAverage[2]);
                fit_par_energy[3].SetValue(parameterAverage[3],0.5*parameterAverage[3],2*parameterAverage[3]);
            //}
            //fit_par_energy[1].SetValue(parameterAverage[1],0.8*parameterAverage[1],1.2*parameterAverage[1]);
            //fit_par_energy[2].SetValue(parameterAverage[2],0.8*parameterAverage[2],1.2*parameterAverage[2]);
            //fit_par_energy[3].SetValue(parameterAverage[3],0.8*parameterAverage[3],1.2*parameterAverage[3]);
            continue;
        }
    }

    auto num12Rings = 1;//setPair?1:fIDArray12.size();
    auto num16Rings = 1;//setPair?1:fIDArray16.size();

    bnnET.ResetNextProjection();
    while (bnnET.NextProjection())
    {
        int pjbin = bnnET.GetCurrentProjectionIt();
        double pjCenter = ThetaLabToCom(bnnET.GetCurrentProjectionCenter());
        double count = arrayX -> At(pjbin);
        double countOff0 = arrayXOff0 -> At(pjbin);
        double countOff1 = arrayXOff1 -> At(pjbin);
        double countOff2 = arrayXOff2 -> At(pjbin);
        double countOff3 = arrayXOff3 -> At(pjbin);
        double count_fit_error = 0;
        if (countOff0>0&&countOff1>0) count_fit_error = abs(countOff1 - countOff0);
        double count_rng_error = 0;
        if (countOff2>0&&countOff3>0) count_rng_error = abs(countOff2 - countOff3);
        double count11 = arrayX11 -> At(pjbin);
        double raw = arrayXRaw -> At(pjbin);
        double xsum1 = arrayX1 -> At(pjbin);
        double xsum2 = arrayX2 -> At(pjbin);
        if (count==0) continue;
        double x1 = ThetaLabToCom(bnnET.GetCurrentProjectionLowEdge());
        double x2 = ThetaLabToCom(bnnET.GetCurrentProjectionUpEdge());
        double ex = 0.5*(x2-x1);
        double theta_com = pjCenter;
        double theta_lab = (180. - theta_com)/2.;
        double theta1_com = x1;
        double theta2_com = x2;
        double theta1_lab = (180. - theta1_com)/2.;
        double theta2_lab = (180. - theta2_com)/2.;
        double solid_angle0 = 1;
        double solid_angle1 = 1;
        double solid_angle2 = 1;
        double solid_angle3 = 1;
        double solid_angle4 = 1;
        double theta1_com_tg1 = fGraphThetaError[0] -> Eval(theta1_com);
        double theta2_com_tg1 = fGraphThetaError[0] -> Eval(theta2_com);
        double theta1_com_tg2 = fGraphThetaError[1] -> Eval(theta1_com);
        double theta2_com_tg2 = fGraphThetaError[1] -> Eval(theta2_com);

        double theta1_lab_rad = TMath::DegToRad()*theta1_lab;
        double theta2_lab_rad = TMath::DegToRad()*theta2_lab;
        double cos1 = TMath::Cos(theta1_com*TMath::DegToRad());
        double cos2 = TMath::Cos(theta2_com*TMath::DegToRad());

        if (theta_lab<42.) {
            double dphi = 0.52359878;
            double solid_anglex = abs( dphi * (TMath::Cos(theta1_com*TMath::DegToRad())-TMath::Cos(theta2_com*TMath::DegToRad())) );

            solid_angle0 = CalculatePairSolidAngle(TMath::DegToRad()*theta1_com,     TMath::DegToRad()*theta2_com, num12Rings);
            solid_angle1 = CalculatePairSolidAngle(TMath::DegToRad()*theta1_com_tg1, TMath::DegToRad()*theta2_com_tg1, num12Rings);
            solid_angle2 = CalculatePairSolidAngle(TMath::DegToRad()*theta1_com_tg2, TMath::DegToRad()*theta2_com_tg2, num12Rings);
            //solid_angle3 = CalculatePairSolidAngle(TMath::DegToRad()*theta1_com_si1, TMath::DegToRad()*theta2_com_si1, num12Rings);
            //solid_angle4 = CalculatePairSolidAngle(TMath::DegToRad()*theta1_com_si2, TMath::DegToRad()*theta2_com_si2, num12Rings);
        }
        else {
            double dphi = 0.390075;
            double solid_anglex = abs( dphi * (TMath::Cos(theta1_com*TMath::DegToRad())-TMath::Cos(theta2_com*TMath::DegToRad())) );

            solid_angle0 = Calculate16RingSolidAngle(TMath::DegToRad()*theta1_com,     TMath::DegToRad()*theta2_com, num16Rings);
            solid_angle1 = Calculate16RingSolidAngle(TMath::DegToRad()*theta1_com_tg1, TMath::DegToRad()*theta2_com_tg1, num16Rings);
            solid_angle2 = Calculate16RingSolidAngle(TMath::DegToRad()*theta1_com_tg2, TMath::DegToRad()*theta2_com_tg2, num16Rings);
            //solid_angle3 = Calculate16RingSolidAngle(TMath::DegToRad()*theta1_com_si1, TMath::DegToRad()*theta2_com_si1, num16Rings);
            //solid_angle4 = Calculate16RingSolidAngle(TMath::DegToRad()*theta1_com_si2, TMath::DegToRad()*theta2_com_si2, num16Rings);
        }

        double valueCorrected = 0;
        int ttaIndex = pjbin;
        if (pairID>=12) ttaIndex += bnnSplitX12.n();
        if (pairID>=0) {
            valueCorrected = count/solid_angle0/fNumAtomsPerUnitArea/fBeamCount/fmbarn;
            fDataRange [pairID][ttaIndex][0] = theta_com;
            fDataRange [pairID][ttaIndex][1] = theta1_com;
            fDataRange [pairID][ttaIndex][2] = theta2_com;
            //
            fDataValues[pairID][ttaIndex][0][0] = count;
            fDataValues[pairID][ttaIndex][0][1] = abs(raw-count);
            fDataValues[pairID][ttaIndex][0][2] = count_fit_error;
            fDataValues[pairID][ttaIndex][0][3] = count_rng_error;
            fDataValues[pairID][ttaIndex][0][4] = raw;
            fDataValues[pairID][ttaIndex][3][0] = solid_angle0;
            fDataValues[pairID][ttaIndex][3][1] = solid_angle1;
            fDataValues[pairID][ttaIndex][3][2] = solid_angle2;
            //
            //fDataValues[pairID][ttaIndex][1][0] = count;
            //fDataValues[pairID][ttaIndex][1][1] = abs(raw-count);
            //fDataValues[pairID][ttaIndex][1][2] = count_fit_error;
            //fDataValues[pairID][ttaIndex][1][3] = count_rng_error;
            //
        }
    }

    return message;
}

void DrawXPair(int pairID)
{
    LKDrawingGroup* sub;
    if (pairID==99) sub = fTop -> FindGroup("x_final");
    else sub = fTop -> FindGroup(Form("pair%d_summary",pairID));
    if (pairID!=99) {
        if (sub==nullptr) {
            e_warning << "sub do not exist" << endl;
            return;
        }
        else {
            if (sub->FindDrawing(Form("et_%d",pairID))==nullptr) {
                e_warning << "et do not exist in " << pairID << endl;
                return;
            }
        }
    }

    auto groupDCS = DrawCrossSection(pairID);
    auto drawX0 = groupDCS -> FindDrawing(Form("x_%d",pairID));
    auto drawR0 = groupDCS -> FindDrawing(Form("xr_%d",pairID));
    auto drawER0 = groupDCS -> FindDrawing(Form("er_%d",pairID));
    auto graphSASum = (TGraphErrors*) groupDCS -> FindObject(Form("solid_angle_sum_%d",pairID));
    auto histCountA = (TH1D*) groupDCS -> FindObject(Form("histCountA_%d",pairID));
    auto histCount = (TH1D*) groupDCS -> FindObject(Form("histCount_%d",pairID));
    auto histCount12 = (TH1D*) groupDCS -> FindObject(Form("histCount12_%d",pairID));
    auto histCount16 = (TH1D*) groupDCS -> FindObject(Form("histCount16_%d",pairID));
    if (fGlobalBeamNormalizationFactor==0)
        fGlobalBeamNormalizationFactor = ((TParameter<double>*) groupDCS -> FindObject(Form("scaleDataToFC_%d",99))) -> GetVal();

    TString mainTitle = Form("(%.3f MeV)%s",fSelectBeamEnergy,(pairID!=99?Form(" Pair-%d",pairID):""));

    auto draw1 = sub -> CreateDrawing(Form("drawX_%d",pairID));
    drawX0 -> CopyTo(draw1);
    draw1 -> SetMainTitle(mainTitle);
    fGraphData[3] -> SetMarkerStyle(21);
    fGraphData[4] -> SetMarkerStyle(21);
    fGraphData[3] -> SetMarkerSize(1.5);
    fGraphData[4] -> SetMarkerSize(1.5);
    fGraphData[3] -> SetMarkerColor(kMagenta);
    fGraphData[4] -> SetMarkerColor(kGreen+2);
    fGraphData[3] -> SetLineWidth(2);
    fGraphData[4] -> SetLineWidth(2);
    fGraphData[3] -> SetLineColor(kMagenta);
    fGraphData[4] -> SetLineColor(kGreen+2);
    if (!fUseSecondPeak) {
        if (pairID==99) draw1 -> Add(fGraphData[3],"samep","7.77");
        if (pairID==99) draw1 -> Add(fGraphData[4],"samep","9.36");
    }

    auto drawC = sub -> CreateDrawing(Form("drawC_%d",pairID),fAllDraw1);
    double max00 = 1.15*histCount->GetMaximum();
    double max12 = 1.15*histCount12->GetMaximum();
    double max16 = 1.15*histCount16->GetMaximum();
    double max = max00;
    if (max12>max) max = max12;
    if (max16>max) max = max16;
    histCountA -> SetMaximum(max);
    drawC -> Add(histCountA);
    //drawC -> Add(histCount);
    drawC -> Add(histCount12,"hist");
    drawC -> Add(histCount16,"hist");
    drawC -> SetMainTitle(mainTitle);

    auto drawR = sub -> CreateDrawing(Form("drawR_%d",pairID));
    drawR0 -> CopyTo(drawR);
    drawR -> SetMainTitle(mainTitle);
    fGraphRatio[3] -> SetMarkerStyle(21);
    fGraphRatio[4] -> SetMarkerStyle(21);
    fGraphRatio[3] -> SetMarkerSize(1.5);
    fGraphRatio[4] -> SetMarkerSize(1.5);
    fGraphRatio[3] -> SetMarkerColor(kMagenta);
    fGraphRatio[4] -> SetMarkerColor(kGreen+2);
    fGraphRatio[3] -> SetLineWidth(2);
    fGraphRatio[4] -> SetLineWidth(2);
    fGraphRatio[3] -> SetLineColor(kMagenta);
    fGraphRatio[4] -> SetLineColor(kGreen+2);
    drawR -> Add(fGraphRatio[3],"samep","7.77");
    drawR -> Add(fGraphRatio[4],"samep","9.36");

    auto drawS = sub -> CreateDrawing(Form("drawS_%d",pairID),fAllDraw1);
    drawS -> Add(drawS->MakeGraphFrame(graphSASum,";theta_com;solid_angle"));
    drawS -> Add(graphSASum,"p");
    drawS -> SetMainTitle(mainTitle);
    
    auto drawER = sub -> CreateDrawing(Form("drawE_%d",pairID));
    drawER0 -> CopyTo(drawER);
}

bool Initialize(int energyIndex, int optionIndex)
{
    fEnergyIndex = energyIndex;

    fUseScalingWithPrevData = true;

    fELARK = new ELARK();
    fELARK -> AddPar("config_stark.mac");
    fELARK -> Init();

    fPar = new LKParameterContainer("config_fit_etheta_com.mac");
    fQProjTestSigmaFactor = fPar -> GetParDouble("q_proj_test_sigma_factor");
    fQProjSigmaFactor     = fPar -> GetParDouble("q_proj_sigma_factor");
    fQProjSigmaFactor2    = 2.;
    fQProjSigmaFactor3    = 3.;
    fHistQEntryCut        = fPar -> GetParDouble("histQ_entry_cut");
    fNumAtomsPerUnitArea  = fPar -> GetParDouble(Form("e%d/num_atoms_per_area",fEnergyIndex));
    fBeamEnergy           = fPar -> GetParDouble(Form("e%d/beam_energy",fEnergyIndex));
    fBeamCount            = fPar -> GetParDouble(Form("e%d/beam_count",fEnergyIndex));
    fThetaCoMScaleRange1  = fPar -> GetParDouble("theta_com_scale_range",0);
    fThetaCoMScaleRange2  = fPar -> GetParDouble("theta_com_scale_range",1);
    fSolidAngleRatio      = fPar -> GetParDouble("solid_angle_ratio");
    fAllDraw1             = fPar -> GetParBool("draw_all_draw1");
    if (fPar -> CheckPar(Form("e%d/global_beam_normalization_factor",fEnergyIndex)))
        fGlobalBeamNormalizationFactor = fPar -> GetParDouble(Form("e%d/global_beam_normalization_factor",fEnergyIndex));
    //fHistEntriesCut = 1000;

    if (fBinWidthForFit==0)  fBinWidthForFit       = fPar -> GetParInt("bin_width_for_fit");
    if (fBackgroundType12<0) fBackgroundType12 = fPar -> GetParInt("background_type_12");
    if (fBackgroundType16<0) fBackgroundType16 = fPar -> GetParInt("background_type_16");
    bnnE = fPar -> GetBinning("binning_e");
    bnnQ = fPar -> GetBinning("binning_q");
    bnnL = fPar -> GetBinning("binning_ttaLab");
    bnnC = fPar -> GetBinning("binning_ttaCoM");
    bnnT0 = bnnC;
    bnnT0.SetW(fBinWidthForFit);
    bnnE.SetName("bnnE");
    bnnQ.SetName("bnnQ");
    bnnL.SetName("bnnL");
    bnnC.SetName("bnnC");
    bnnT0.SetName("bnnT");
    bnnE.Print();
    bnnQ.Print();
    bnnL.Print();
    bnnC.Print();
    bnnT0.Print();
    bnnSplitX12.SetName("bnnSplitX12");
    bnnSplitX16.SetName("bnnSplitX16");
    bnnSplitX12 = fPar -> GetBinning(Form("e%d/binng_Split_x_12",fEnergyIndex));//LKBinning(4,32,40);
    bnnSplitX16 = fPar -> GetBinning(Form("e%d/binng_Split_x_16",fEnergyIndex));//LKBinning(8,51.5,71.5);
    bnnSplitX12Com.SetName("bnnSplitX12Com");
    bnnSplitX16Com.SetName("bnnSplitX16Com");
    bnnSplitX12Com = LKBinning(bnnSplitX12.n(),ThetaLabToCom(bnnSplitX12.x2()),ThetaLabToCom(bnnSplitX12.x1()));
    bnnSplitX16Com = LKBinning(bnnSplitX16.n(),ThetaLabToCom(bnnSplitX16.x2()),ThetaLabToCom(bnnSplitX16.x1()));
    bnnSplitX12Com.Print();
    bnnSplitX16Com.Print();

    fTagBnn = "_";
    fTagBnn += LKMisc::RemoveTrailing0(Form("%d",bnnSplitX16Com.n ()),true) + "_";
    fTagBnn += LKMisc::RemoveTrailing0(Form("%f",bnnSplitX16Com.x1()),true).ReplaceAll(".","p") + "_";
    fTagBnn += LKMisc::RemoveTrailing0(Form("%f",bnnSplitX16Com.x2()),true).ReplaceAll(".","p") + "__";
    fTagBnn += LKMisc::RemoveTrailing0(Form("%d",bnnSplitX12Com.n ()),true) + "_";
    fTagBnn += LKMisc::RemoveTrailing0(Form("%f",bnnSplitX12Com.x1()),true).ReplaceAll(".","p") + "_";
    fTagBnn += LKMisc::RemoveTrailing0(Form("%f",bnnSplitX12Com.x2()),true).ReplaceAll(".","p");

    fMaxPoints = bnnT0.n()+2;
    if (fMaxPoints<bnnSplitX12.n()+bnnSplitX16.n()+2) fMaxPoints=bnnSplitX12.n()+bnnSplitX16.n()+2;

    fSelectBeamEnergy = fBeamEnergy;
    fBeamEnergyString = Form("%f",fSelectBeamEnergy);
    fBeamEnergyString.ReplaceAll(".","");
    while (fBeamEnergyString[fBeamEnergyString.Sizeof()-2]=='0')
        fBeamEnergyString = fBeamEnergyString(0,fBeamEnergyString.Sizeof()-2);

    fFinalThetaRange[0][0] = fPar -> GetParDouble(Form("e%d/final_theta_range",fEnergyIndex),0);
    fFinalThetaRange[0][1] = fPar -> GetParDouble(Form("e%d/final_theta_range",fEnergyIndex),1);
    fFinalThetaRange[1][0] = fPar -> GetParDouble(Form("e%d/final_theta_range",fEnergyIndex),2);
    fFinalThetaRange[1][1] = fPar -> GetParDouble(Form("e%d/final_theta_range",fEnergyIndex),3);

    for (auto i : {0,1}) {
        auto fileK = new TFile(Form("%s/et_kinematic_line.e%d.%d.root",fDataXPath.Data(),fEnergyIndex,i));
        if (fileK->IsZombie()) {
            e_error << "Need " << fileK -> GetName() << " !!" << endl;
        }
        else {
            fGraphK[i] = (TGraph*) fileK -> Get("Graph");
            fGraphK[i] -> SetName("GraphK");
            fGraphK[i] -> SetLineColor(kRed);
        }
    }

    vector<int> exIDs = fPar -> GetParVInt(Form("e%d/exclude_pair_ids",fEnergyIndex));
    for (auto id : exIDs) e_info << "Excluding detector-" << id << endl;


    if (fIDArray12.size()==0)
        for (auto pairID=0; pairID<12; ++pairID)
            if (LKMisc::ValueIsInArray(pairID,exIDs)==false)
                fIDArray12.push_back(pairID);

    if (fIDArray16.size()==0)
        for (auto pairID=12; pairID<fNumPairs; ++pairID)
            if (LKMisc::ValueIsInArray(pairID,exIDs)==false)
                fIDArray16.push_back(pairID);

    fProjectionBins.push_back(0);
    for (auto bin=1; bin<=bnnT0.n(); bin++)
        fProjectionBins.push_back(bin);

    const int maxNumPar = 8;
    fHistParameters = new int**[fNumPairs];
    fFitParameters = new double***[fNumPairs];
    fFitRange = new double**[fNumPairs];
    for (int iPair=0; iPair<fNumPairs; ++iPair) {
        fHistParameters[iPair] = new int*[bnnT0.n()+1];
        fFitParameters[iPair] = new double**[bnnT0.n()+1];
        fFitRange[iPair] = new double*[bnnT0.n()+1];
        for (int pjbin=0; pjbin<bnnT0.n()+1; ++pjbin)
        {
            fHistParameters[iPair][pjbin] = new int[10];
            fHistParameters[iPair][pjbin][0] = -1;
            fHistParameters[iPair][pjbin][1] = -1;
            fHistParameters[iPair][pjbin][2] = 0;
            fHistParameters[iPair][pjbin][3] = 0;
            //fHistParameters[iPair][pjbin][4] = -1;
            //fHistParameters[iPair][pjbin][5] = -1;
            //fHistParameters[iPair][pjbin][6] = -1;
            //fHistParameters[iPair][pjbin][7] = -1;
            //fHistParameters[iPair][pjbin][8] = -1;
            //fHistParameters[iPair][pjbin][9] = -1;
            fFitRange[iPair][pjbin] = new double[2];
            fFitRange[iPair][pjbin][0] = bnnQ.x1();
            fFitRange[iPair][pjbin][1] = bnnQ.x2();
            fFitParameters[iPair][pjbin] = new double*[maxNumPar];
            for (int iPar=0; iPar<maxNumPar; ++iPar) {
                fFitParameters[iPair][pjbin][iPar] = new double[3];
                for (int iVE=0; iVE<3; ++iVE) {
                    fFitParameters[iPair][pjbin][iPar][iVE] = -1;
                }
            }
            bool backgroundType = (iPair<12?fBackgroundType12:fBackgroundType16);
            if (backgroundType==0) {
                fFitParameters[iPair][pjbin][0][0] = 100;
                fFitParameters[iPair][pjbin][0][1] = 0;
                fFitParameters[iPair][pjbin][0][2] = 1000;
                fFitParameters[iPair][pjbin][1][0] = 0.1;
                fFitParameters[iPair][pjbin][1][1] = -0.8;
                fFitParameters[iPair][pjbin][1][2] = 0.8;
                fFitParameters[iPair][pjbin][2][0] = 0.2;
                fFitParameters[iPair][pjbin][2][1] = 0.1;
                fFitParameters[iPair][pjbin][2][2] = 0.3;
                fFitParameters[iPair][pjbin][3][0] = 10;
                fFitParameters[iPair][pjbin][3][1] = 0;
                fFitParameters[iPair][pjbin][3][2] = 1000;
                fFitParameters[iPair][pjbin][4][0] = -1;
                fFitParameters[iPair][pjbin][4][1] = -200;
                fFitParameters[iPair][pjbin][4][2] = 0;
                fFitParameters[iPair][pjbin][5][0] = -1;
                fFitParameters[iPair][pjbin][5][1] = -1;
                fFitParameters[iPair][pjbin][5][2] = 0.5;
            }
            if (backgroundType==1) {
                fFitParameters[iPair][pjbin][0][0] = -1;
                fFitParameters[iPair][pjbin][0][1] = -1;
                fFitParameters[iPair][pjbin][0][2] = 0.5;
                fFitParameters[iPair][pjbin][1][0] = 100;
                fFitParameters[iPair][pjbin][1][1] = 0;
                fFitParameters[iPair][pjbin][1][2] = 1000;
                fFitParameters[iPair][pjbin][2][0] = 0.1;
                fFitParameters[iPair][pjbin][2][1] = -0.8;
                fFitParameters[iPair][pjbin][2][2] = 0.8;
                fFitParameters[iPair][pjbin][3][0] = 0.2;
                fFitParameters[iPair][pjbin][3][1] = 0.1;
                fFitParameters[iPair][pjbin][3][2] = 0.3;
                fFitParameters[iPair][pjbin][4][0] = 0.1;
                fFitParameters[iPair][pjbin][4][1] = 0;
                fFitParameters[iPair][pjbin][4][2] = 0.4;
                fFitParameters[iPair][pjbin][5][0] = 5;
                fFitParameters[iPair][pjbin][5][1] = 0;
                fFitParameters[iPair][pjbin][5][2] = 10;
                fFitParameters[iPair][pjbin][6][0] = 1;
                fFitParameters[iPair][pjbin][6][1] = 0;
                fFitParameters[iPair][pjbin][6][2] = 10;
            }
        }
    }

    fDataValues = new double***[fNumPairs];
    fDataRange = new double**[fNumPairs];
    for (int iPair=0; iPair<fNumPairs; ++iPair) {
        fDataValues[iPair] = new double**[fMaxPoints];
        fDataRange[iPair] = new double*[fMaxPoints];
        for (int iPoint=0; iPoint<fMaxPoints; ++iPoint) {
            fDataValues[iPair][iPoint] = new double*[fNumVariables];
            fDataRange[iPair][iPoint] = new double[3];
            fDataRange[iPair][iPoint][0] = 0;
            fDataRange[iPair][iPoint][1] = 0;
            fDataRange[iPair][iPoint][2] = 0;
            for (int iVar=0; iVar<fNumVariables; ++iVar) {
                fDataValues[iPair][iPoint][iVar] = new double[fNumVE];
                for (int iVE=0; iVE<fNumVE; ++iVE) {
                    fDataValues[iPair][iPoint][iVar][iVE] = 0;
                }
            }
        }
    }

    fGraphRF = new TGraph();
    fGraphRF -> SetName("rutherford");
    fGraphRF -> SetLineColor(kBlack);
    TString rutherfordFileName  = fPar -> GetParString(Form("e%d/rutherford_file", fEnergyIndex));
    std::ifstream fileRF(rutherfordFileName);
    std::string line;
    std::getline(fileRF, line);
    double value[10];
    while (fileRF >> value[0] >> value[1] >> value[2] >> value[3] >> value[4] >> value[5] >> value[6] >> value[7] >> value[8] >> value[9]) {
        if (value[1]>2000)
            continue;
        fGraphRF -> SetPoint(fGraphRF->GetN(),value[0],value[1]);
    }

    fNumFresco = fPar -> GetParN(Form("e%d/fresco_file" ,fEnergyIndex));
    for (auto ifc=0; ifc<fNumFresco; ++ifc) {
        TString frescoName = fPar -> GetParString(Form("e%d/fresco_file" ,fEnergyIndex),ifc);
        e_info << "Fresco file-" << ifc << ": " << frescoName << endl;
        fGraphFC[ifc] = new TGraph();
        fGraphFC[ifc] -> SetName(Form("fresco_%d",ifc));
        fGraphFC[ifc] -> SetTitle(frescoName);
        //fGraphFC[ifc] -> SetLineColor(ifc+2);
        if (ifc==0) fGraphFC[ifc] -> SetLineColor(kRed);
        if (ifc==1) fGraphFC[ifc] -> SetLineColor(kBlue);
        ifstream file2(frescoName);
        for (int i=0; i<13; ++i) { std::getline(file2, line); }
        double xfc, yfc;
        while (file2 >> xfc >> yfc) {
            if (yfc>2000)
                continue;
            fGraphFC[ifc] -> SetPoint(fGraphFC[ifc]->GetN(),xfc,yfc);
        }
    }

    fQTEdge = (TGraph*) (new TFile(fDataXPath+"/qt_edge.6789.root","read"))->Get("Graph");
    fQTEdge -> SetLineColor(kRed);
    fQTEdge -> SetLineWidth(2);

    //TString path = "data_reco";
    TString path = "data_ana";
    TString cut_pid_name = "all";
    if (optionIndex==0) cut_pid_name = "cutX_ProtonAndStop";
    else if (optionIndex==1) cut_pid_name = "cutX_ProtonAndERange";
    else if (optionIndex==2) cut_pid_name = "cutX_Stopped";
    else if (optionIndex==3) cut_pid_name = "cutX_ProtonInERange";
    else if (optionIndex==4) cut_pid_name = "cutX_ProtonInERange1";
    else if (optionIndex==5) cut_pid_name = "cutX_ProtonInERange2";
    else if (optionIndex==6) cut_pid_name = "cutX_ProtonCut";
    else if (optionIndex==7) cut_pid_name = "cutX_ProtonCutWideE";

    TString tag = "ana_xmulth";

    fNameFileX = Form(fDataXPath+"/x%d_%s.%s.root",fEnergyIndex,cut_pid_name.Data(),tag.Data());

    fAnaName = Form("KO2421_e%d_%s",fEnergyIndex,cut_pid_name.Data());
    e_info << "analysis name : " << fAnaName << endl;

    TString topName = fAnaName;
    topName.ReplaceAll("KO2421","anaX");
    fTop = new LKDrawingGroup(topName);
    fSave = new LKDrawingGroup(Form("data_lilak/anaX_e%d_cutX_ProtonCut.root",fEnergyIndex),"xprint");
    if (fSave->GetEntries()==0) {
        e_warning << "Failed!" << endl;
        e_warning << "Trying path " << fDataXPath << endl;
        fSave = new LKDrawingGroup(Form("%s/anaX_e%d_cutX_ProtonCut.root",fDataXPath.Data(),fEnergyIndex),"xprint");
        if (fSave->GetEntries()!=0)
            e_info << "Good!" << endl;
        else {
            e_warning << "No save!" << endl;
            fSave = nullptr;
        }
    }

    //fSaveFit = new LKDrawingGroup(fTop->GetSaveFileName("","","fit_parameters"));
    fSaveFit = nullptr;

    if (0)
    {
        TString recoFileName = Form("%s/%s_NT360.%s.root",path.Data(),fAnaName.Data(),tag.Data());
        //TString recoFileName = Form("%s/%s_NT180.ana_eve.root",path.Data(),fAnaName.Data());
        //recoFileName = "/home/ejungwoo/lilak/ko2421/macros/data_ana/KO2421_e4_cutX_pidAll_NT360.ana_xmulth.root";
        e_info << recoFileName << endl;
        auto fileReco = new TFile(recoFileName.Data(),"read");
        if (fileReco->IsOpen()==true) {
            e_info << "Openning " << recoFileName << endl;
            auto runHeader = (LKParameterContainer*) fileReco -> Get("RunHeader");
            fMainName = runHeader -> GetParString("MainName");
            fMainName = fMainName(0,9);
            fMainName.ReplaceAll("_e","-E");
        }
        if (fileReco->IsOpen()==false)
            e_error << "Cannot open " << recoFileName << endl;

        fMainTitle = Form("[%s] (%.3f MeV)",fAnaName.Data(),fSelectBeamEnergy);
        if (fileReco!=nullptr) fTree = (TTree*) fileReco -> Get("event");
    }

    if (fSave==nullptr && fTree==nullptr ) {
        e_error << "No save no tree" << endl;
        return false;
    }

    fGraphQE = new TGraph();
    fGraphQE -> SetMarkerStyle(20);

    auto file_tta_error = new TFile(fDataXPath+"/theta_conversion.root","read");
    fGraphThetaError[0] = (TGraph*) file_tta_error -> Get("tta_to_tta12_a0_tg0");
    fGraphThetaError[1] = (TGraph*) file_tta_error -> Get("tta_to_tta12_a0_tg1");
    fGraphThetaError[2] = (TGraph*) file_tta_error -> Get("tta_to_tta12_a1_tg0");
    fGraphThetaError[3] = (TGraph*) file_tta_error -> Get("tta_to_tta12_a1_tg1");

    //auto file_tta_error2 = new TFile(fDataXPath+"/theta_conversion.2.root","read");
    //file_tta_error2 -> Print();
    //fGraphThetaError[4] = (TGraph*) file_tta_error2 -> Get("tta_to_tta12_tg1");
    //fGraphThetaError[5] = (TGraph*) file_tta_error2 -> Get("tta_to_tta12_tg2");
    //fGraphThetaError[6] = (TGraph*) file_tta_error2 -> Get("tta_to_tta12_si1");
    //fGraphThetaError[7] = (TGraph*) file_tta_error2 -> Get("tta_to_tta12_si2");

    fUseSavedLinearBG = fPar -> GetParBool("use_saved_linear_bg_fit_parameter");
    for (auto ring : {12,16})
    {
        bool backgroundType = (ring==12?fBackgroundType12:fBackgroundType16);
        TString name_fit_paramter = Form(fDataXPath+"/fit_parameters_%s.%s.dat",fAnaName.Data(),"exp_bg");
        //if (backgroundType==0)
            //name_fit_paramter = Form(fDataXPath+"/fit_parameters_%s.dat",fAnaName.Data());
            //name_fit_paramter = Form(fDataXPath+"/fit_parameters_%s.%s.dat",fAnaName.Data(),"linear_bg");
        name_fit_paramter = Form(fDataXPath+"/fit_parameters_%s.ring%d.%s.dat",fAnaName.Data(),ring,((backgroundType==0)?"linear_bg":"exp_bg"));

        if (fUseSavedLinearBG)
            name_fit_paramter = Form(fDataXPath+"/fit_parameters_%s.dat",fAnaName.Data());
        e_info << "Reading " << name_fit_paramter << " for " << ring << endl;

        TString dummy;
        int pairID, pjbin, mxrbin, pickup, remove;
        double cntylm, range1, range2;
        fFileFitParameterIn.open(name_fit_paramter);
        while(fFileFitParameterIn >> dummy >> pairID >> pjbin >> cntylm >> mxrbin >> pickup >> remove >> range1 >> range2)
        {
            //lk_debug << dummy << " "<<  pairID << " "<<  pjbin << " "<<  cntylm << " "<<  mxrbin << " "<<  pickup << " "<<  remove << " "<<  range1 << " "<<  range2 << endl;
            if (ring==12 && pairID>=12) continue;
            if (ring==16 && pairID< 12) continue;
            double** fitpv = fFitParameters[pairID][pjbin];
            if (pickup==1) {
                fPickup = true;
                e_note << Form("Found pickup at pair=%d, pjbin=%d",pairID,pjbin) << endl;
                e_cout << "    continuity: " << cntylm << endl;
                e_cout << "    max rebin : " << mxrbin << endl;
                e_cout << "    pickup    : " << pickup << endl;
                e_cout << "    remove    : " << remove << endl;
            }
            if (pairID>=0) {
                fHistParameters[pairID][pjbin][0] = mxrbin;
                fHistParameters[pairID][pjbin][1] = cntylm;
                fHistParameters[pairID][pjbin][2] = pickup;
                fHistParameters[pairID][pjbin][3] = remove;
                if (range1<-1) range1 = -1;
                //if (ring==16 && range1>-0.8) range1 = -0.8;
                fFitRange[pairID][pjbin][0] = range1;
                fFitRange[pairID][pjbin][1] = range2;
            }

            int numPar = 6;
            if (!fUseSavedLinearBG&&backgroundType==0) numPar = 6;
            else if (!fUseSavedLinearBG&&backgroundType==1) numPar = 7;
            for (auto i=0; i<numPar; ++i)
            {
                TString parName;
                double value;
                double value1;
                double value2;
                fFileFitParameterIn >> parName >> value >> value1 >> value2;
                if (pairID>=0) {
                    fitpv[i][0] = value; 
                    fitpv[i][1] = value1; 
                    fitpv[i][2] = value2;
                }
                if (pickup==1) e_cout << "    " << std::left << std::setw(10) << parName << ": " << std::left << std::setw(10) << value << " (" << value1 << ", " << value2 << ")" << endl;
            }
            if (backgroundType==1)
            {
                if (fUseSavedLinearBG)
                {
                    double v50 = fitpv[5][0];
                    double v51 = fitpv[5][1];
                    double v52 = fitpv[5][2];
                    double v00 = fitpv[0][0];
                    double v01 = fitpv[0][1];
                    double v02 = fitpv[0][2];
                    double v10 = fitpv[1][0];
                    double v11 = fitpv[1][1];
                    double v12 = fitpv[1][2];
                    double v20 = fitpv[2][0];
                    double v21 = fitpv[2][1];
                    double v22 = fitpv[2][2];
                    fitpv[0][0] = v50;
                    fitpv[0][1] = v51;
                    fitpv[0][2] = v52;
                    fitpv[1][0] = v00;
                    fitpv[1][1] = v01;
                    fitpv[1][2] = v02;
                    fitpv[2][0] = v10;
                    fitpv[2][1] = v11;
                    fitpv[2][2] = v12;
                    fitpv[3][0] = v20;
                    fitpv[3][1] = v21;
                    fitpv[3][2] = v22;
                    fitpv[4][0] = 0.2;
                    fitpv[4][1] = 0;
                    fitpv[4][2] = 2;
                    fitpv[5][0] = 5;
                    fitpv[5][1] = 0;
                    fitpv[5][2] = 20;
                    fitpv[6][0] = 0;
                    fitpv[6][1] = -5;
                    fitpv[6][2] = 5;
                    if (fitpv[0][1]<-1) fitpv[0][1] = -1; // qrange
                    if (fitpv[2][0]>1) fitpv[2][0] = 1; // mean
                    if (fitpv[3][0]>0.35)  {
                        fitpv[3][0] = 0.35; // sigma
                    }
                    fitpv[3][1] = 0.1; // sigma
                    fitpv[3][2] = 0.5; // sigma
                }
            }
        }
        fFileFitParameterIn.close();
    }

    fParFile = new TFile(Form("%s/parameters_e%d.root",fDataXPath.Data(),fEnergyIndex),"recreate");
    fParTree = new TTree("par","");
    //fParTree -> Branch("pairID",&pairID);
    //fParTree -> Branch("pjbin" ,&pjbin);
    //if (fBackgroundType16==0) {
    //    fParTree -> Branch("amp"   ,&bValue[0]);
    //    fParTree -> Branch("mean"  ,&bValue[1]);
    //    fParTree -> Branch("sigma" ,&bValue[2]);
    //    fParTree -> Branch("itcpt" ,&bValue[3]);
    //    fParTree -> Branch("slope" ,&bValue[4]);
    //    fParTree -> Branch("qbndr" ,&bValue[5]);
    //}
    //if (fBackgroundType16==1) {
    //    fParTree -> Branch("pair"  ,&bIndex[0]);
    //    fParTree -> Branch("pj"    ,&bIndex[1]);
    //    fParTree -> Branch("theta" ,&bIndex[2]);
    //    fParTree -> Branch("qb"    ,&bValue[0]);
    //    fParTree -> Branch("amp"   ,&bValue[1]);
    //    fParTree -> Branch("mean"  ,&bValue[2]);
    //    fParTree -> Branch("sigma" ,&bValue[3]);
    //    fParTree -> Branch("e_amp" ,&bValue[4]);
    //    fParTree -> Branch("e_off" ,&bValue[5]);
    //    fParTree -> Branch("const" ,&bValue[6]);
    //}

    fGraphParameterGGAR = new TGraph();
    fGraphParameterGGAR -> SetMarkerStyle(20);
    fGraphParameterGGAR -> SetMarkerSize(0.8);

    if (fEnergyIndex==5||fEnergyIndex==8)
    {
        TString nameGGAR = Form("%s/ggar_e%d.root",fDataXPath.Data(),fEnergyIndex);
        auto fileFitGGAR = new TFile(nameGGAR);
        if (fileFitGGAR->IsOpen()) fFitGGAR = (TF1*) fileFitGGAR -> Get("fit");
        else e_error << "Cannot find " << nameGGAR << endl;
        //fFitGGAR = nullptr; // for writting new fit
    }

    return true;
}

double GetFirstBoundary(TH1D* histQ, double thresholdFactor)
{
    auto nbins = histQ -> GetXaxis() -> GetNbins();
    auto threshold = histQ->GetMaximum()*thresholdFactor;
    for (auto bin=1; bin<=nbins; ++bin) {
        double value = histQ -> GetBinContent(bin);
        if (value>threshold) {
            return (histQ -> GetXaxis() -> GetBinCenter(bin));
        }
    }
    return -999.999;
}

TGraphErrors* FitAmplitude1To2(int idx, TGraphErrors* graph1, TGraph* graph2, TF1* fit, double &scale)
{
    fGraphFix[idx] = graph2;
    graph1 -> Fit(fit,"RQN0");
    scale = fit -> GetParameter(0);

    auto numPoints = graph1 -> GetN();
    auto graph3 = NewGraph(Form("scaled_%s",graph1->GetName()));
    double x,y,ex,ey;
    for (auto iPoint=0; iPoint<numPoints; ++iPoint)
    {
        graph1 -> GetPoint(iPoint,x,y);
        ex = graph1 -> GetErrorX(iPoint);
        ey = graph1 -> GetErrorY(iPoint);
        graph3 -> SetPoint(iPoint,x,y/scale);
        graph3 -> SetPointError(iPoint,ex,ey/scale);
    }
    return graph3;
}

TString MakeName(int pairID, TString header, TString tag)
{
    return ( header+tag+"_"+fBeamEnergyString+((pairID>=0)?Form("_pair%d",pairID):"") );
}

void InitFitParameters2(int energyIndex, int pairID)
{
    fit_par_energy[0].SetNameValue("v0_x_low", 1, 0, 2, false);
    fit_par_energy[1].SetNameValue("v1_exp_pos", 6, 4, 9, false);
    fit_par_energy[2].SetNameValue("v2_x_pow", 1, 0, 3, false);
    fit_par_energy[3].SetNameValue("v3_const", 4, 0, 20, false);
    if (pairID>=4)
        fit_par_energy[3].SetNameValue("v3_const", 4, 0, 40, false);
    fit_par_energy[4].SetNameValue("v4_gaus0_a", 1000, 0, 1.E+3, false);
    fit_par_energy[5].SetNameValue("v5_gaus0_m", 10, 0, 20, false);
    fit_par_energy[6].SetNameValue("v6_gaus0_s", 7, 3, 15, false);
    fit_par_energy[7].SetNameValue("v7_gaus1_a", 1000, 0, 1.E+2, false);
    fit_par_energy[8].SetNameValue("v8_gaus1_m", 10, 0, 20, false);
    if (pairID>=4&&pairID<12) {
        fit_par_energy[7].SetNameValue("v7_gaus1_a", 0,0,0,true);
        fit_par_energy[8].SetNameValue("v8_gaus1_m", 0,0,0,true);
    }
    if (energyIndex==5 && (pairID<12)) {
        //fit_par_energy[0].SetNameValue("v0_x_low", 1, 1, 2, false);
        fit_par_energy[1].SetNameValue("v1_exp_pos", 6, 0, 9, false);
        fit_par_energy[2].SetNameValue("v2_x_pow", 6, 4, 8, false);
        fit_par_energy[3].SetNameValue("v3_const", 4, 0, 20, false);
    }
    if (energyIndex==4 && (pairID<12)) {
        fit_par_energy[1].SetNameValue("v1_exp_pos", 6, 0, 9, false);
        fit_par_energy[2].SetNameValue("v2_x_pow", 2, 2, 5, false);
    }
}

double CalculatePairSolidAngle(double theta1, double theta2, int numDetectors)
{
    double dphi = (2*TMath::Pi()/12) * numDetectors;
    double solid_angle = abs( dphi * (TMath::Cos(theta1)-TMath::Cos(theta2)) );
    return solid_angle;
}

double Calculate16RingSolidAngle(double theta1, double theta2, int numDetectors)
{
    double dphi = (2*TMath::Pi()/16) * numDetectors;
    double solid_angle = abs( dphi * (TMath::Cos(theta1)-TMath::Cos(theta2)) );
    return solid_angle;
}
