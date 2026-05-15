void geometry()
{
    auto top = new LKDrawingGroup("SiPosition");
    auto group = top -> CreateGroup();
    auto drawSetup = group -> CreateDrawing("drawSetup");
    group -> CreateDrawing();

    double radius[3] = {95, 102, 126.5};
    double distance[3] = {138.5, 138.5, 77.5};
    double detectorWidth = 75.0;
    double detectorThickness = 4.;
    double targetRadius = 11.;
    double targetFB = 2;
    double siPositionError = 1; //mm

    double z1 = -20;
    double z2 = 220;
    int lstyle = 2;

    //int detectors[] = {0,1,2};
    int detectors[] = {1,2};

    auto hist = new TH2D("hist",";z;y",100,z1,z2,100,-20,155);
    hist -> SetStats(0);
    drawSetup -> Add(hist);

    auto line = new TLine(z1,0,z2,0);
    line -> SetLineStyle(lstyle);
    drawSetup -> Add(line);

    TVector3 posDet[3][4];
    for (auto iDet : detectors)
    {
        posDet[iDet][0] = TVector3(0,radius[iDet]+0.5*detectorThickness,distance[iDet]-0.5*detectorWidth);
        posDet[iDet][1] = TVector3(0,radius[iDet]+0.5*detectorThickness,distance[iDet]+0.5*detectorWidth);
        posDet[iDet][2] = TVector3(0,radius[iDet]-0.5*detectorThickness,distance[iDet]-0.5*detectorWidth);
        posDet[iDet][3] = TVector3(0,radius[iDet]-0.5*detectorThickness,distance[iDet]+0.5*detectorWidth);
    }

    // drawing
    //for (auto iVertex : {0,1,2})
    for (auto iVertex : {0})
    {
        double pm = iVertex;
        if (iVertex==2) pm = -1;
        TVector3 posVertex(0,pm*targetRadius,0);
        //drawSetup -> Add(new TMarker(posVertex.Z(),posVertex.Y(),20));
        auto graphVertex = new TGraphErrors();
        graphVertex -> SetPoint(0,posVertex.Z(),posVertex.Y());
        graphVertex -> SetPointError(0,0,targetRadius);
        graphVertex -> SetMarkerStyle(20);
        drawSetup -> Add(graphVertex,"samepl");

        for (auto iDet : detectors)
        {
            if (iVertex==0)
            {
                auto graphDet = new TGraph();
                graphDet -> SetPoint(0,posDet[iDet][0].Z(),posDet[iDet][0].Y());
                graphDet -> SetPoint(1,posDet[iDet][1].Z(),posDet[iDet][1].Y());
                graphDet -> SetPoint(2,posDet[iDet][3].Z(),posDet[iDet][3].Y());
                graphDet -> SetPoint(3,posDet[iDet][2].Z(),posDet[iDet][2].Y());
                graphDet -> SetPoint(4,posDet[iDet][0].Z(),posDet[iDet][0].Y());
                drawSetup -> Add(graphDet,"l");

                double yc = 0.5*(posDet[iDet][0] + posDet[iDet][2]).Y();
                auto liney = new TLine(posDet[iDet][1].Z(),yc,z2,yc);
                liney -> SetLineStyle(2);
                drawSetup -> Add(liney);

                auto txy = new TLatex(z2,(iDet==0?yc-7:yc+7),Form("y=%.1f",yc));
                txy -> SetTextFont(42);
                txy -> SetTextSize(0.04);
                txy -> SetTextAlign(32);
                drawSetup -> Add(txy);

                if (iDet!=0)
                {
                    double offy = 0;
                    if (iDet==1) offy = -20;
                    auto line0 = new TLine(posDet[iDet][0].Z(),0+offy,posDet[iDet][0].Z(),posDet[iDet][0].Y());
                    auto line1 = new TLine(posDet[iDet][1].Z(),0+offy,posDet[iDet][1].Z(),posDet[iDet][1].Y());
                    line0 -> SetLineStyle(lstyle);
                    line1 -> SetLineStyle(lstyle);
                    drawSetup -> Add(line0);
                    drawSetup -> Add(line1);

                    int pmy = 1;
                    if (iDet==1) pmy = -1;
                    auto tx0 = new TLatex(posDet[iDet][0].Z(),8*pmy,Form(" z=%.1f",posDet[iDet][0].Z()));
                    auto tx1 = new TLatex(posDet[iDet][1].Z(),8*pmy,Form(" z=%.1f",posDet[iDet][1].Z()));
                    for (auto tx : {tx0,tx1}) {
                        tx -> SetTextFont(42);
                        tx -> SetTextSize(0.04);
                        tx -> SetTextAlign(12);
                        //if (iDet==1) tx -> SetTextAlign(32);
                        drawSetup -> Add(tx);
                    }
                }
            }

            if (0)
            //if (iDet==0)
            {
                auto graph1 = new TGraph();
                drawSetup -> Add(graph1);
                graph1 -> SetPoint(0,posVertex.Z(),posVertex.Y());
                graph1 -> SetPoint(1,posDet[iDet][0].Z(),posDet[iDet][0].Y());
                auto graph2 = new TGraph();
                drawSetup -> Add(graph2,"samep");
                graph2 -> SetPoint(0,posVertex.Z(),posVertex.Y());
                graph2 -> SetPoint(1,posDet[iDet][1].Z(),posDet[iDet][1].Y());
            }

            if (iVertex==0 && (iDet==1 || iDet==2))
            {
                for (auto j : {0,1}) {
                    TVector3 posDetIn = 0.5*(posDet[iDet][0] + posDet[iDet][2]);
                    if (j==1) posDetIn = 0.5*(posDet[iDet][1] + posDet[iDet][3]);
                    auto graph1 = new TGraph();
                    drawSetup -> Add(graph1);
                    graph1 -> SetPoint(0,posVertex.Z(),posVertex.Y());
                    graph1 -> SetPoint(1,posDetIn.Z(),posDetIn.Y());
                    TVector3 v3(posDetIn.Z()-posVertex.Z(),posDetIn.Y()-posVertex.Y(),0);
                    TLatex* tx = new TLatex(posDetIn.Z(),posDetIn.Y()+10,Form("#theta=%.1f^{#circ}",v3.Phi()*TMath::RadToDeg()));
                    tx -> SetTextFont(42);
                    tx -> SetTextSize(0.04);
                    tx -> SetTextAlign(12);
                    if (j==1) tx -> SetTextAlign(32);
                    drawSetup -> Add(tx);
                }
            }
        }
    }

    TGraph* graph_tta_from_tgpos_error[2][2] = {0};
    for (auto iAxis : {0,1}) {
        for (auto pm : {0,1}) {
            graph_tta_from_tgpos_error[iAxis][pm] = new TGraph();
            graph_tta_from_tgpos_error[iAxis][pm] -> SetMarkerStyle(20);
            graph_tta_from_tgpos_error[iAxis][pm] -> SetName(Form("tta_to_tta12_a%d_tg%d",iAxis,pm));
        }
    }

    for (int iAxis : {0,1})
    {
        for (auto iDet : detectors)
        {
            int color = iDet + 2;
            if (iDet==0) color = kRed;
            if (iDet==1) color = kBlue;
            if (iDet==2) color = kGreen+1;

            int nn = 20;
            //int nn = 10;
            double z0 = posDet[iDet][0].Z();
            double z1 = posDet[iDet][1].Z();
            double dz = (z1-z0)/nn;
            for (int iz=0; iz<=nn; ++iz)
            {
                double z = z0 + iz*dz;
                TVector3 posDetIn = 0.5*(posDet[iDet][0] + posDet[iDet][2]);
                posDetIn.SetZ(z);

                auto mk = new TMarker(posDetIn.Z(),posDetIn.Y(),20);
                mk -> SetMarkerColor(color);
                mk -> SetMarkerSize(0.8);
                drawSetup -> Add(mk);

                double ttalab[3] = {0};
                double ttacom[3] = {0};
                for (auto iVertex : {0,1,2})
                {
                    double pm = iVertex;
                    if (iVertex==2) pm = -1;
                    TVector3 posVertex;
                    if      (iAxis==0) posVertex = TVector3(0,pm*targetRadius,0);
                    else if (iAxis==1) posVertex = TVector3(0,0,pm*targetFB);
                    ttalab[iVertex] = (posDetIn-posVertex).Theta()*TMath::RadToDeg();
                    ttacom[iVertex] = 180. - ttalab[iVertex] * 2.;
                }

                cout << iz << setw(10) << z << setw(10) << ttalab[0] << setw(10) << ttalab[1] << setw(10) << ttalab[2] << endl;

                graph_tta_from_tgpos_error[iAxis][0] -> SetPoint(graph_tta_from_tgpos_error[iAxis][0]->GetN(),ttacom[0],ttacom[1]);
                graph_tta_from_tgpos_error[iAxis][1] -> SetPoint(graph_tta_from_tgpos_error[iAxis][1]->GetN(),ttacom[0],ttacom[2]);

                for (auto iSiPos : {0,1,2})
                {
                    double pm = iSiPos;
                    if (iSiPos==2) pm = -1;
                    TVector3 posVertex(0,0,pm*siPositionError/*mm*/); // shifting vertex instead of si position for convinience
                    ttalab[iSiPos] = (posDetIn-posVertex).Theta()*TMath::RadToDeg();
                    ttacom[iSiPos] = 180. - ttalab[iSiPos] * 2.;
                }
            }
        }
        graph_tta_from_tgpos_error[iAxis][0] -> Sort();
        graph_tta_from_tgpos_error[iAxis][1] -> Sort();

        auto drawtg = group -> CreateDrawing(Form("draw_pos%d",iAxis));
        if (iAxis==0) drawtg -> SetCreateFrame("drawtg",Form(";#theta_{com};#theta_{com,1,2} (pos. error y #pm %.1f mm)",targetRadius));
        if (iAxis==1) drawtg -> SetCreateFrame("drawtg",Form(";#theta_{com};#theta_{com,1,2} (pos. error z #pm %.1f mm)",targetFB));
        drawtg -> Add(graph_tta_from_tgpos_error[iAxis][0],"pl");
        drawtg -> Add(graph_tta_from_tgpos_error[iAxis][1],"pl");
        drawtg -> SetGridx();
        drawtg -> SetGridy();
    }

    group -> SetCanvasSize(1800,1200);
    group -> Draw("v");
    group -> GetCanvas() -> SaveAs("PositionErrorAndTheta.png");
    group -> GetCanvas() -> SaveAs("PositionErrorAndTheta.eps");
    //drawSetup -> Draw();

    TString fileName = "data/theta_conversion.root";
    auto file_tta = new TFile(fileName,"recreate");
    graph_tta_from_tgpos_error[0][0] -> Write();
    graph_tta_from_tgpos_error[0][1] -> Write();
    graph_tta_from_tgpos_error[1][0] -> Write();
    graph_tta_from_tgpos_error[1][1] -> Write();
    cout << fileName << endl;
    file_tta -> ls();

    //top -> Save();
}
