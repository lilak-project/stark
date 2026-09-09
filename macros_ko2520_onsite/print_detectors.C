/// Dump the full silicon mapping.
///   option : "all" (summary + detector table + channel table),
///            "summary", "detector", "channel", "group" (channels listed under each detector)
void print_detectors(TString mappingPath = "mapping_ko2520_0904", TString option = "all")
{
    LKSiliconMapping mapping;
    if (!mapping.Load(mappingPath)) {
        cerr << "Failed to load " << mappingPath << endl;
        return;
    }

    mapping.Print(option);

    cout << endl;
    mapping.PrintDetectorChannels();

    cout << endl << "== 12-ring dE-E pairs (matched by phi#) ======================================" << endl;
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
