#include "TFile.h"
#include "TF1.h"
#include "TH1.h"
#include "TObjArray.h"

#include <algorithm>
#include <fstream>
#include <iomanip>
#include <map>
#include <tuple>
#include <vector>

namespace {
using Key=std::tuple<int,int,int>;
struct Row { int det,side,strip; double entries,intercept,slope; };

void collect(TObject *object,std::map<Key,TH1*> &result)
{
    if (auto h=dynamic_cast<TH1*>(object)) {
        int d,s,t;
        if (h->GetDimension()==1 &&
            sscanf(h->GetName(),"stage1_energy_det%d_side%d_strip%d",&d,&s,&t)==3)
            result[{d,s,t}]=h;
        return;
    }
    if (auto a=dynamic_cast<TObjArray*>(object))
        for (auto item:*a) collect(item,result);
}
}

void correct_run77_c0()
{
    const char *inputName="data_calibration/run_0077/ko2520_0077.c0.par";
    const char *outputName="data_calibration/run_0077/ko2520_0077.c0.corrected.par";
    TFile histFile("ko2520_0077.si_calibration_stage_1.root","READ");
    std::map<Key,TH1*> histograms;
    collect(histFile.Get("top"),histograms);

    std::ifstream input(inputName);
    std::string header;
    std::getline(input,header);
    std::vector<Row> rows;
    Row row;
    while (input>>row.det>>row.side>>row.strip>>row.entries>>row.intercept>>row.slope)
        rows.push_back(row);

    std::vector<double> goodSlopes;
    for (const auto &r:rows)
        if (r.slope>0.0055 && r.slope<0.0075 && TMath::Abs(r.intercept)<0.7)
            goodSlopes.push_back(r.slope);
    std::sort(goodSlopes.begin(),goodSlopes.end());
    auto median=goodSlopes[goodSlopes.size()/2];

    std::ofstream output(outputName);
    output<<"# det side strip entries intercept slope\n"<<std::setprecision(12);
    std::ofstream report("data_calibration/run_0077/ko2520_0077.c0.correction_report.txt");
    report<<"# det side strip old_intercept old_slope gd_adc new_intercept new_slope\n";
    int corrected=0;
    for (auto &r:rows) {
        // Also recover active QQQ5 junction strips for which the automatic
        // calibration left a zero slope.  Without this, those strips vanish
        // completely from the calibrated spectrum.
        bool missingActiveQQQ = r.det>=1 && r.det<=4 && r.side==0
            && r.entries>0 && r.slope<=0;
        bool suspicious=missingActiveQQQ || (r.slope>0 &&
            (r.slope<0.0045 || r.slope>0.0085 || TMath::Abs(r.intercept)>0.7));
        auto found=histograms.find({r.det,r.side,r.strip});
        if (suspicious && found!=histograms.end()) {
            auto h=found->second;
            auto low=3.1822/(median*1.30);
            auto high=3.1822/(median*0.72);
            int first=h->FindBin(low),last=h->FindBin(high),maximum=first;
            for (int bin=first+1;bin<=last;++bin)
                if (h->GetBinContent(bin)>h->GetBinContent(maximum)) maximum=bin;
            auto peak=h->GetBinCenter(maximum);
            if (h->GetBinContent(maximum)>=2 && peak>0) {
                TF1 fit("gd_refit","gaus",peak-35,peak+35);
                fit.SetParameters(h->GetBinContent(maximum),peak,10);
                h->Fit(&fit,"Q0NR");
                auto mean=fit.GetParameter(1);
                if (mean>low && mean<high) peak=mean;
                auto oldIntercept=r.intercept,oldSlope=r.slope;
                r.intercept=0;
                r.slope=3.1822/peak;
                report<<r.det<<" "<<r.side<<" "<<r.strip<<" "
                    <<oldIntercept<<" "<<oldSlope<<" "<<peak<<" 0 "<<r.slope<<"\n";
                ++corrected;
            }
        }
        output<<r.det<<" "<<r.side<<" "<<r.strip<<" "<<r.entries<<" "
              <<r.intercept<<" "<<r.slope<<"\n";
    }
    Info("correct_run77_c0","median slope %.8f; corrected %d channels; wrote %s",
        median,corrected,outputName);
}
