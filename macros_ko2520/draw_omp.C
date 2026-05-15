#include "yield.C"

double find_minimum(TGraph* graph, double x1, double x2, int itNumber=3, bool EvalWithSpline=true);
double theta_com_to_lab(double theta_com);

void draw_omp()
{
    bool setLogy = true;
    bool showRing = false;
    bool saveFigures = true;
    bool findPeak = false;
    bool showMin = false;
    bool showBeta = true;
    bool showStatErrorCut = false;
    bool showExpectedDataPoints = false;
    vector<int> NaArray = {21,23,24,25};
    //vector<int> NaArray = {21,25};
    int n_mse = 3;
    int mse[]    = {10,7,5};    // minimum statistical error
    double cs_mse[] = {0.5,1,2}; // value of cross section where it satistify the minimum statistical error
    //int mse[]    = {7,5,3};    // minimum statistical error
    //double cs_mse[] = {1,2,5.5}; // value of cross section where it satistify the minimum statistical error

    if (0) {
        setLogy = true;
        showRing = false;
        saveFigures = true;
        findPeak = false;
        showMin = false;
        showBeta = false;
        showStatErrorCut = false;
        showExpectedDataPoints = false;
        NaArray.clear();
        NaArray.push_back(21);
        NaArray.push_back(23);
        NaArray.push_back(24);
        NaArray.push_back(25);
    }
    else if (1) {
        setLogy = true;
        showRing = true;
        saveFigures = true;
        findPeak = false;
        showMin = false;
        showBeta = false;
        showStatErrorCut = true;
        showExpectedDataPoints = false;
        NaArray.clear();
        //NaArray.push_back(21);
        NaArray.push_back(25);
    }

    double theta_lab_range[3][2] = {{50,70}, {30,40}, {9,27}};
    double theta_com_range[3][2] = {{40,80}, {100,120}, {126,162}};

    auto top = new LKDrawingGroup("top");
    auto gd_lab = top -> CreateGroup("lab",false);
    auto gd_com = top -> CreateGroup("com");

    auto top1 = new LKDrawingGroup("top1");
    auto draw_a_com = top1 -> CreateDrawing("draw_a_com");
    auto draw_a_lab = top1 -> CreateDrawing("draw_a_lab");

    double legend_width = 0.45;
    draw_a_com -> SetCreateLegend(0,legend_width,0.07);
    draw_a_lab -> SetCreateLegend(1,legend_width,0.07);
    if (setLogy) {
        draw_a_com -> SetLogy();
        draw_a_lab -> SetLogy();
    }

    LKBinning bnnx_lab(100,0,70);
    LKBinning bnnx_com(100,25,180);
    LKBinning bnnx_com_x[3];
    bnnx_com_x[0] = LKBinning(100,40,80);
    bnnx_com_x[1] = LKBinning(100,100,120);
    bnnx_com_x[2] = LKBinning(100,126,162);
    LKBinning bnny(100,0,30);
    if (setLogy)
        bnny = LKBinning(100,0.1,2000);

    auto hist_com = (bnnx_com*bnny).NewH2("hist_com",";#theta_{c.m.};d#rho/d#Omega (mb/sr)");
    hist_com -> SetStats(0);
    draw_a_com -> Add(hist_com,"",".");

    auto hist_lab = (bnnx_lab*bnny).NewH2("hist_lab",";#theta_{lab};Cross section");
    hist_lab -> SetStats(0);
    draw_a_lab -> Add(hist_lab,"",".");

    for (auto hist : {hist_com,hist_lab}) {
        hist -> GetXaxis() -> SetLabelFont(132);
        hist -> GetYaxis() -> SetLabelFont(132);
        hist -> GetXaxis() -> SetTitleFont(132);
        hist -> GetYaxis() -> SetTitleFont(132);
        hist -> GetXaxis() -> SetLabelSize(0.045);
        hist -> GetYaxis() -> SetLabelSize(0.045);
        hist -> GetXaxis() -> SetTitleSize(0.050);
        hist -> GetYaxis() -> SetTitleSize(0.050);
        hist -> GetXaxis() -> SetTitleOffset(1.12);
        hist -> GetYaxis() -> SetTitleOffset(1.2);
    }

    map<int, double> df;
    df[21] = 0.505;
    df[23] = 0.487;
    df[24] = 0.377;
    df[25] = 0.310;

    LKParameter colors("colors_4.mac");
    if (NaArray.size()==2) colors = LKParameter("colors_2.mac");
    LKParameter colors_ring("colors_ring.mac");
    lk_debug << colors_ring.GetColor() << endl;

    int count = 0;
    for (auto ANa : NaArray)
    {
        auto draw_com = gd_com -> CreateDrawing("com",true);
        auto draw_lab = gd_lab -> CreateDrawing("lab",true);
        draw_lab -> SetCanvasMargin(0.12,0.05,0.13,0.05);
        draw_com -> SetCanvasMargin(0.12,0.05,0.13,0.05);
        draw_lab -> Add(hist_lab,".");
        draw_com -> Add(hist_com,".");
        if (showBeta) draw_com -> AddLegendLine(Form("#beta_{2} = %.3f",df[ANa]));
        if (showBeta) draw_lab -> AddLegendLine(Form("#beta_{2} = %.3f",df[ANa]));
        draw_com -> SetCreateLegend(0,0.32,0.07);
        draw_lab -> SetCreateLegend(1,0.32,0.07);
        //draw_com -> SetLegendTransparent();
        //draw_lab -> SetLegendTransparent();
        if (setLogy) {
            draw_com -> SetLogy();
            draw_lab -> SetLogy();
        }

        if (showRing)
        {
            for (auto ring : {0,1,2})
            {
                double x1 = theta_com_range[ring][0];
                double x2 = theta_com_range[ring][1];
                TGraph *graph = new TGraph();
                graph -> SetPoint(0,x1,bnny.x1());
                graph -> SetPoint(1,x1,bnny.x2());
                graph -> SetPoint(2,x2,bnny.x2());
                graph -> SetPoint(3,x2,bnny.x1());
                graph -> SetPoint(4,x1,bnny.x1());
                graph -> SetFillColor(colors_ring.GetColor());
                draw_com -> Add(graph,"samefl",".");
                auto tt = new TText((x1+x2)/2.,bnny.Lerp(0.00005),Form("R%d",ring+1));
                tt -> SetTextFont(132);
                tt -> SetTextAlign(22);
                tt -> SetTextSize(0.05);
                draw_com -> Add(tt);
            }
        }

        for (auto iType : {0})
        {
            TString fileNameOMP = Form("data_omp/cx_avg_drhbc_%dNa.txt",ANa);
            if (iType==1) fileNameOMP = Form("data_omp/cx_rchb_%dNa.txt",ANa);
            auto graph_com = new TGraph();
            auto graph_lab = new TGraph();
            graph_com -> SetName(Form("graph_%dNa_%d",ANa,iType));
            graph_lab -> SetName(Form("graph_%dNa_%d",ANa,iType));
            cout << fileNameOMP << endl;
            ifstream file(fileNameOMP);
            std::string line;
            std::getline(file, line);
            std::getline(file, line);
            double theta_com, y;

            while (file >> theta_com >> y)
            {
                if ((bnnx_com*bnny).IsInside(theta_com,y))
                {
                    auto theta_lab = theta_com_to_lab(theta_com);
                    graph_com -> SetPoint(graph_com->GetN(),theta_com,y);
                    graph_lab -> SetPoint(graph_lab->GetN(),theta_lab,y);
                }
            }
            file.close();

            auto graph_cmd = new TGraphErrors();
            graph_cmd -> SetName(Form("graph_%dNa_%d_edata",ANa,iType));
            if (iType==0 && showExpectedDataPoints)
            {
                for (auto ring : {0,1,2})
                {
                    double x1 = theta_com_range[ring][0];
                    double x2 = theta_com_range[ring][1];
                    LKBinning bnnx_com(x1,x2);
                    bnnx_com.SetWX(4);
                    //bnnx_com.Print();
                    bnnx_com.Reset();
                    while (bnnx_com.Next()) {
                        auto data_x = bnnx_com.GetItCenter();
                        auto data_y = graph_com -> Eval(data_x,0,"S");
                        double yield_4 = calculate_yield(4, data_x-2,data_x+2,data_y);
                        double error_y_percentage = sqrt(yield_4)/yield_4*100;
                        cout << "x=" << data_x << " y=" << data_y << ", ep=" << error_y_percentage << " %, e=" << data_y*0.01*error_y_percentage << endl;
                        graph_cmd -> SetPoint(graph_cmd->GetN(),data_x,data_y);
                        graph_cmd -> SetPointError(graph_cmd->GetN()-1,2,data_y*0.01*error_y_percentage);
                    }
                    //double xmin = find_minimum(graph_com,x1,x2);
                    //double yatm = graph_com -> Eval(xmin,0,"S");
                    //draw_com -> AddLegendLine(Form("y_{%d}=%f",ring+1,yatm));
                }
                //for (auto bnnx : {bnnx_com1,bnnx_com2,bnnx_com3})
                //{
                //    if (bnnx.IsInside(theta_com))
                //    {
                //        double ee = 0.1;
                //        double yh = y + y * (TMath::Power(10,ee) - 1);
                //        double yl = y - y * (1 - TMath::Power(10,-ee));
                //        graph_com_s -> SetPoint(graph_com_s->GetN(),theta_com,y);
                //        graph_com_s -> SetPointError(graph_com_s->GetN()-1, 0, 0, yl, yh);
                //        cout << ANa << " " << theta_com << " " << yl << " " << yh << endl;
                //    }
                //    draw_com -> Add(graph_com_s,"3",".");
                //}
            }

            auto title_reaction = Form("^{%d}Na+p",ANa);
            //auto title = Form("Deformed ^{%d}Na+p",ANa);
            //if (iType==1) title = Form("Spherical ^{%d}Na+p",ANa);
            auto title = "Deformed";
            if (iType==1) title = "Spherical";
            graph_com -> SetLineColor(colors.GetColor(count));
            graph_lab -> SetLineColor(colors.GetColor(count));
            graph_com -> SetLineWidth(2);
            graph_lab -> SetLineWidth(2);
            graph_cmd -> SetMarkerStyle(24);
            //graph_cmd -> SetLineWidth(2);
            if (iType==1) {
                graph_com -> SetLineStyle(2);
                graph_lab -> SetLineStyle(2);
            }

            double xp_com = find_minimum(graph_com, 60, 80);
            double xp_lab = find_minimum(graph_lab, 60, 80);

            draw_a_com -> Add(graph_com,"l",title);
            draw_a_lab -> Add(graph_lab,"l",title);
            draw_com -> Add(graph_com,"l",title);
            draw_lab -> Add(graph_lab,"l",title);
            if (showExpectedDataPoints)
                draw_com -> Add(graph_cmd,"p","Expected");

            if (showMin)
                for (auto ring : {0,1,2})
                {
                    double x1 = theta_com_range[ring][0];
                    double x2 = theta_com_range[ring][1];
                    double xmin = find_minimum(graph_com,x1,x2);
                    double yatm = graph_com -> Eval(xmin,0,"S");
                    draw_com -> AddLegendLine(Form("y_{%d}=%f",ring+1,yatm));
                }

            if (findPeak) {
                draw_com -> AddLegendLine(Form("%f",xp_com));
                draw_lab -> AddLegendLine(Form("%f",xp_lab));
            }

            //auto tt = new TLatex(bnnx_com.Lerp(1.-legend_width),bnny.Lerp(0.95),title_reaction);
            //tt -> SetTextAlign(33);
            //tt -> SetTextFont(132);
            //tt -> SetTextSize(0.08);
            //draw_com -> Add(tt);

            //auto pv = new TPaveText(bnnx_com.x1()+0.03*bnnx_com.d(),bnny.x2()-0.18*bnny.d(),bnnx_com.x1()+0.3*bnnx_com.d(),bnny.x2());
            //auto pv = new TPaveText(bnnx_com.x1()+0.03*bnnx_com.d(),bnny.x2()-0.18*bnny.d(),bnnx_com.x1()+0.3*bnnx_com.d(),bnny.x2());
            //auto pv = new TPaveText(0.03,0.82,0.3,1,"NDC");
            auto pv = new TPaveText(0.15,0.80,0.4,0.95,"NDC");
            pv -> AddText(title_reaction);
            pv -> SetTextAlign(22);
            pv -> SetTextFont(132);
            pv -> SetTextSize(0.08);
            pv -> SetBorderSize(1);
            draw_com -> Add(pv,"same auto");
        }

        if (showStatErrorCut)
        {
            for (auto i=0; i<n_mse; ++i)
            {
                int error = mse[i];
                double value = cs_mse[i];
                auto line = new TLine(bnnx_com.x1(),value,bnnx_com.Lerp(0.91),value);
                if (i==0) {}//line -> SetLineStyle(1);
                if (i==1) line -> SetLineStyle(9);
                if (i==2) line -> SetLineStyle(2);
                if (i==3) line -> SetLineStyle(4);
                draw_com -> Add(line,"samel");
                auto tt = new TText(bnnx_com.Lerp(0.91),value,Form(" %.d %s",error,"%"));
                tt -> SetTextAlign(11);
                tt -> SetTextFont(132);
                tt -> SetTextSize(0.035);
                draw_com -> Add(tt);
                if (i==n_mse-1) {
                    ///auto t2 = new TLatex(bnnx_com.Lerp(0),value,"error_{stat}");
                    auto t2 = new TLatex(bnnx_com.Lerp(0.025),2.2,"statistical error");
                    t2 -> SetTextAlign(11);
                    t2 -> SetTextFont(132);
                    t2 -> SetTextSize(0.035);
                    //draw_com -> Add(t2);

                    //auto t3 = new TLatex(bnnx_com.Lerp(0.03),0.55,"stat. error (4 days)");
                    auto t3 = new TLatex(bnnx_com.Lerp(0.03),0.55,"stat. error (32 hours)");
                    t3 -> SetTextAlign(11);
                    t3 -> SetTextFont(132);
                    t3 -> SetTextSize(0.035);
                    draw_com -> Add(t3);
                }
            }
            //auto line 
            //calculate_yield(4,theta_com_range[0][0],theta_com_range[0][1]);
            //calculate_yield(4,theta_com_range[1][0],theta_com_range[1][1]);
            //calculate_yield(4,theta_com_range[2][0],theta_com_range[2][1]);
        }

        if (0)
        if (ANa==23)
        {
            TString fileNamePre = Form("data_pre/x_%dNapp_8MeV.csv",ANa);
            auto graph = new TGraph(fileNamePre,"%lg, %lg");
            graph -> SetMarkerStyle(20);
            graph -> Print();
            draw_com -> Add(graph,"p","J. Hellstrom (1969)");
        }

        count++;
    }

    //top1 -> SetCanvasSize(1200,500);
    //top1 -> Draw();

    for (auto group : {gd_com})
    { 
        //group -> SetCanvasDivision(4,1);
        group -> SetCanvasDivision(2,2);
        if (NaArray.size()==2) {
            //group -> SetCanvasDivision(1,2);
            group -> SetCanvasSize(500,250,0);
            group -> SetCanvasDivision(2,1);
        }
        group -> Draw();
        if (saveFigures)
            //group -> GetCanvas() -> SaveAs(Form("figures/fig_%s_%d%d%s%s.png",group->GetName(),group->GetDivX(),group->GetDivY(),(setLogy?".logy":""),(showRing?".ring":"")));
            group -> GetCanvas() -> SaveAs(Form("figures/fig_%s_%d%d%s%s.pdf",group->GetName(),group->GetDivX(),group->GetDivY(),(setLogy?".logy":""),(showRing?".ring":"")));
    }

    //top -> Draw();
    //top -> Save();
}

double theta_com_to_lab(double theta_com) { return (180. - theta_com)/2.; }

double find_minimum(TGraph* graph, double x1, double x2, int itNumber, bool EvalWithSpline)
{
    LKBinning bnn(20,x1,x2);
    double yMin = DBL_MAX;
    double xAtYMin = 0;

    TString splineOption = "";
    if (EvalWithSpline) splineOption = "S";

    int count = 0;
    while (count++<itNumber)
    {
        bnn.Reset();
        while (bnn.Next())
        {
            auto x = bnn.GetItValue();
            auto y = graph -> Eval(x,0,"S");
            if (y<yMin) {
                xAtYMin = x;
                yMin = y;
            }
        }
        auto ww = bnn.w();
        bnn = LKBinning(10,xAtYMin-ww,xAtYMin+ww);
    }

    return xAtYMin;
}
