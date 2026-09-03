void print_detectors()
{
    LKSiliconMapping mapping;
    if (!mapping.Load("mapping_ko2520")) {
        cerr << "Failed to load mapping_ko2520" << endl;
        return;
    }

    cout << "Number of detectors: " << mapping.GetNumDetectors() << endl;
    cout << "idx  det#  type  ring   dE/E phi#"
         << "       radius       z-center          phi          phi1          phi2" << endl;
    for (auto iDet=0; iDet<mapping.GetNumDetectors(); ++iDet) {
        auto det = mapping.GetDetectorByVectorIndex(iDet);
        if (det == nullptr) {
            cerr << "Missing detector at vector index " << iDet << endl;
            continue;
        }
        printf("%3d %5d %5s %5s %5s %4d %12.4f %14.4f %12.4f %12.4f %12.4f\n",
               det->detIndex, det->detNumber, det->detType.Data(), det->ringType.Data(),
               det->dEE.Data(), det->phiNumber, det->ringRadius, det->ringZ,
               det->phi, det->phi1, det->phi2);
    }

    cout << endl << "12-ring dE-E pairs (matched by phi#)" << endl;
    cout << "phi#   dE idx/det#/phi       E idx/det#/phi" << endl;
    for (auto phiNumber=0; phiNumber<12; ++phiNumber) {
        const LKSiliconMapping::DetectorInfo *detDE = nullptr;
        const LKSiliconMapping::DetectorInfo *detE = nullptr;
        for (auto iDet=0; iDet<mapping.GetNumDetectors(); ++iDet) {
            auto det = mapping.GetDetectorByVectorIndex(iDet);
            if (det == nullptr || det->ringNumber != 12 || det->phiNumber != phiNumber)
                continue;
            if (det->dEE == "dE") detDE = det;
            if (det->dEE == "E")  detE = det;
        }
        if (detDE != nullptr && detE != nullptr) {
            printf("%4d   %3d/%4d/%7.1f       %3d/%4d/%7.1f\n",
                   phiNumber, detDE->detIndex, detDE->detNumber, detDE->phi,
                   detE->detIndex, detE->detNumber, detE->phi);
        }
        else {
            printf("%4d   MISSING dE or E detector\n", phiNumber);
        }
    }
}
