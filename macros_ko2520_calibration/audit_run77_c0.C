#include "TFile.h"
#include "TH1.h"
#include "TObjArray.h"

#include <vector>

namespace {
void collect(TObject *object, std::vector<TH1*> &histograms)
{
    if (auto hist=dynamic_cast<TH1*>(object)) {
        if (hist->GetDimension()==1 && TString(hist->GetName()).BeginsWith("stage6_energy_"))
            histograms.push_back(hist);
        return;
    }
    if (auto array=dynamic_cast<TObjArray*>(object))
        for (auto item:*array) collect(item,histograms);
}
double integral(TH1 *h,double lo,double hi)
{
    return h->Integral(h->FindBin(lo),h->FindBin(hi));
}
}

void audit_run77_c0()
{
    TFile file("ko2520_0077.si_calibration_stage_6.root");
    std::vector<TH1*> histograms;
    collect(file.Get("top"),histograms);
    printf("# channel entries Gd_counts Am_counts Am_over_Gd\n");
    for (auto h:histograms) {
        auto gd=integral(h,3.1822-0.16,3.1822+0.16);
        auto am=integral(h,5.486-0.16,5.486+0.16);
        if (h->GetEntries()>=100 && (gd==0 || am/(gd+1.)>2.5))
            printf("%s %.0f %.0f %.0f %.3f\n",h->GetName(),h->GetEntries(),gd,am,am/(gd+1.));
    }
}
