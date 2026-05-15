void draw_detectors()
{
    auto elark = new ELARK();
    elark -> AddPar("config_stark.mac");
    elark -> Init();
    auto drawing = elark -> NewHistPlaneDrawing("elark");
    auto cvs = LKPainter::GetPainter() -> CanvasFull("",0.8);
    cvs -> SetMargin(0.02,0.02,0.11,0.08);
    drawing -> SetCanvas(cvs);
    drawing -> Draw();
    cvs -> SaveAs("data_x/detectors_phi_ring.eps");

    auto numDetectors = elark -> GetNumSiDetectors();
    double solid_angle_sum = 0.;
    for (auto det=0; det<numDetectors; ++det) {
        auto detector = elark -> GetSiDetector(det);

        auto pairID = elark -> FindEPairDetectorID(det,0);
        detector = elark -> GetSiDetector(pairID);
        double phi1 = TMath::DegToRad()*detector -> GetPhi1();
        double phi2 = TMath::DegToRad()*detector -> GetPhi2();
        double tta1 = TMath::DegToRad()*detector -> GetTheta1();
        double tta2 = TMath::DegToRad()*detector -> GetTheta2();
        double solid_angle = abs((phi2-phi1)*(TMath::Cos(tta1)-TMath::Cos(tta2)));

        //cout << pairID << " ->  phi=(" << phi1 << ", " << phi2 << "), theta=" << " (" << tta1 << ", " << tta2 << ") -> " << solid_angle << endl;
        //cout << pairID << " ->  phi=(" << phi2-phi1 << "), theta=" << " (" << tta1 << ", " << tta2 << ") -> " << solid_angle << endl;
        solid_angle_sum += solid_angle;
    }
    cout << solid_angle_sum << endl;

    for (auto det : {0,15})
    {
        auto detector = elark -> GetSiDetector(det);
        auto pairID = elark -> FindEPairDetectorID(det,0);
        detector = elark -> GetSiDetector(pairID);
        double tta1 = detector -> GetTheta1();
        double tta2 = detector -> GetTheta2();
        cout << det << " "  << tta1 << " " << tta2 << endl;
    }
}
