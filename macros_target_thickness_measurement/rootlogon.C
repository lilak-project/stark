{
    TString message = "libs: ";

    // Resolve LILAK build path
    TString lilakBuild = TString(gSystem->Getenv("LILAK_PATH")) + "/build";
    if (lilakBuild == "/build" || !gSystem->AccessPathName(
            (TString(gSystem->Getenv("HOME")) + "/Research/lilak/build/libLILAK.dylib")))
        lilakBuild = TString(gSystem->Getenv("HOME")) + "/Research/lilak/build";

    TString libName = lilakBuild + "/libLILAK";
    int loadv = gSystem->Load(libName);
    if (loadv == 0 || loadv == 1) {
        message = message + "LILAK ";
        gROOT->ProcessLine("#include \"LKCompiled.h\"");

        // Add LILAK include paths so headers are found in macros
        TString src = lilakBuild + "/../source";
        gSystem->AddIncludePath(" -I" + src + "/tool/drawing");
        gSystem->AddIncludePath(" -I" + src + "/tool");
        gSystem->AddIncludePath(" -I" + src + "/base");
        gSystem->AddIncludePath(" -I" + lilakBuild);
    }
    else {
        cout << "Error while loading " << libName << endl;
    }

    cout << message << endl;

    for (TString fileName : {"README.lilak","readme.lilak"})
    {
        string line;
        ifstream readme(fileName);
        if (readme.is_open()) {
            while (getline(readme, line))
                cout << line << endl;
            break;
        }
    }
}
