TGraphErrors *NewGraph(TString name, int mst=-1, double msz=-1, int mcl=-1, int lst=-1, int lsz=-1, int lcl=-1)
{
    auto graph = new TGraphErrors();
    graph -> SetName(name);
    if (mst<=0) mst = 20;
    if (msz<=0) msz = 1;
    if (mcl<0) mcl = kBlack;
    graph -> SetMarkerStyle(mst);
    graph -> SetMarkerSize(msz);
    graph -> SetMarkerColor(mcl);
    if (lst<0) lst = 1;
    if (lsz<0) lsz = 1;
    if (lcl<0) lcl = mcl;
    graph -> SetLineStyle(lst);
    graph -> SetLineWidth(lsz);
    graph -> SetLineColor(lcl);
    graph -> SetFillStyle(0);
    return graph;
}

double GetMaxValueAround(TH1D *histE, double mean1K, int binRange)
{
    auto bin0 = histE -> GetXaxis() -> FindBin(mean1K);
    auto binMax = histE -> GetXaxis() -> GetNbins();
    double valueMax = 0;
    for (auto binOff=-binRange; binOff<=binRange; ++binOff)
    {
        int bin = bin0 + binOff;
        if (bin<1) continue;
        if (bin>binMax) continue;
        double value = histE -> GetBinContent(bin);
        if (value>valueMax) valueMax = value;
    }
    return valueMax;
}

