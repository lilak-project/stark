#ifndef SKCHANNELFILTERTASK_HH
#define SKCHANNELFILTERTASK_HH

#include "LKTask.h"
#include "TClonesArray.h"
#include <vector>

class LKRun;

/// Evaluate every rule against the original SiChannel array, then compact once.
class SKChannelFilterTask : public LKTask
{
    public:
        SKChannelFilterTask();
        bool Init() override;
        void Exec(Option_t*) override;

    private:
        TClonesArray *fSiChannelArray = nullptr; //!
        TClonesArray *fFitDataArray = nullptr; //!
        LKRun *fEventRun = nullptr; //!
        double fNoiseLevel = -1;
        double fSignalToNoiseRatio = -1;
        double fSignalThreshold = 0;
        bool fBipolar = false;
        bool fFilterTime = false;
        double fTimeMin = 0; ///< Inclusive lower hit-time bound, in waveform TB.
        double fTimeMax = 0; ///< Inclusive upper hit-time bound, in waveform TB.
        bool fFilterEnergy = false;
        double fEnergyMin = 0; ///< Inclusive lower bound on stored channel energy.
        double fEnergyMax = 0; ///< Inclusive upper bound; ADC amplitude after LKSetSiChannelTask.
        struct EnergyRangeRule {
            int address[4]; // cobo, asad, aget, channel; -1 matches any value
            double minimum;
            double maximum;
        };
        // More specified address fields win; ties use the last matching rule.
        std::vector<EnergyRangeRule> fEnergyRanges; //! Parsed runtime rules
        bool fFilterFiredChannelCount = false;
        int fFiredChannelCountMin = 0; ///< Inclusive count before any channel removal.
        int fFiredChannelCountMax = 0;
        // Event veto in matching cobo/asad/aget; -1 is wildcard.
        // Use PFData flags when available, otherwise scan SiChannel waveforms at the signal rail.
        std::vector<std::vector<int>> fSaturationEventRules;
        // cobo, asad (-1 wildcard), first/last AGET (inclusive), junction/ohmic counts.
        // Count original fired channels over the union of the selected AGETs.
        std::vector<std::vector<int>> fSideMultiplicityRules;
        int fMinChannels = 1;
        // Each rule: cobo, asad, aget, comparison code, count.
        std::vector<std::vector<int>> fCrosstalkRules;

    ClassDef(SKChannelFilterTask, 7)
};
#endif
