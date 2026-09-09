#include "SKChannelFilterTask.h"
#include "LKRun.h"
#include "LKSiChannel.h"
#include "LKPulseFitData.h"
#include <algorithm>
#include <cmath>
#include <map>
#include <regex>
#include <string>
#include <sstream>
#include <stdexcept>
#include <tuple>

ClassImp(SKChannelFilterTask)

SKChannelFilterTask::SKChannelFilterTask() : LKTask("SKChannelFilterTask", "Filter silicon channels before reconstruction") {}

bool SKChannelFilterTask::Init()
{
    fEventRun = dynamic_cast<LKRun*>(fRun);
    if (!fEventRun || !fRun->GetBranchA("SiChannel", "LKSiChannel")) return false;
    fSiChannelArray = fRun->KeepBranchA("SiChannel", "LKSiChannel");
    fFitDataArray = fRun->GetBranchA("PFData", false);
    fPar->UpdatePar(fNoiseLevel, "SKChannelFilterTask/Filter/NoiseLevel");
    fPar->UpdatePar(fSignalToNoiseRatio, "SKChannelFilterTask/Filter/SignalToNoiseRatio");
    fPar->UpdatePar(fSignalThreshold, "SKChannelFilterTask/Filter/SignalThreshold");
    fPar->UpdatePar(fBipolar, "SKChannelFilterTask/Filter/Bipolar");
    fPar->UpdatePar(fMinChannels, "SKChannelFilterTask/Filter/MinChannels");
    if (!std::isfinite(fNoiseLevel) || !std::isfinite(fSignalToNoiseRatio) ||
        !std::isfinite(fSignalThreshold) || fSignalThreshold < 0 || fMinChannels < 0 ||
        (fNoiseLevel < 0 && fNoiseLevel != -1) ||
        (fSignalToNoiseRatio < 0 && fSignalToNoiseRatio != -1)) {
        lk_error << "Invalid Filter thresholds (use -1 to disable noise/SNR cuts)" << endl;
        return false;
    }
    fSideMultiplicityRules.clear();
    const TString sideParameter = "SKChannelFilterTask/Filter/SideMultiplicity";
    for (int i=0; i<fPar->GetEntriesFast(); ++i) {
        auto par = fPar->GetParameter(i);
        if (!par || par->IsCommentOut() || TString(par->GetName()) != sideParameter) continue;
        std::string raw = par->GetRaw().Data();
        std::replace(raw.begin(), raw.end(), ',', ' ');
        std::istringstream input(raw);
        std::vector<int> rule;
        try {
            auto integer = [&](std::string token) {
                TString expanded(token);
                fPar->ReplaceVariables(expanded);
                token = expanded.Data();
                size_t used = 0;
                int value = std::stoi(token, &used);
                if (used != token.size()) throw std::invalid_argument("expected integers");
                return value;
            };
            if (raw.find('=') == std::string::npos) {
                // Legacy positional form: cobo asad agetMin agetMax junction ohmic.
                for (std::string token; input >> token;) rule.push_back(integer(token));
            } else {
                std::map<std::string, std::string> fields;
                for (std::string token; input >> token;) {
                    auto equals = token.find('=');
                    if (equals == std::string::npos || equals == 0 || equals+1 == token.size())
                        throw std::invalid_argument("use key=value without spaces around =");
                    auto key = token.substr(0, equals);
                    if (key != "cobo" && key != "asad" && key != "aget" && key != "junction" && key != "ohmic")
                        throw std::invalid_argument("unknown field: " + key);
                    if (!fields.emplace(key, token.substr(equals+1)).second)
                        throw std::invalid_argument("duplicate field: " + key);
                }
                if (fields.size() != 5)
                    throw std::invalid_argument("require cobo, asad, aget, junction, ohmic");
                const auto &aget = fields.at("aget");
                auto colon = aget.find(':');
                int first = integer(aget.substr(0, colon));
                int last = colon == std::string::npos ? first : integer(aget.substr(colon+1));
                rule = {integer(fields.at("cobo")), integer(fields.at("asad")), first, last,
                        integer(fields.at("junction")), integer(fields.at("ohmic"))};
            }
            if (rule.size() != 6 || rule[0] < -1 || rule[1] < -1 ||
                rule[2] < 0 || rule[3] < rule[2] || rule[4] < 0 || rule[5] < 0)
                throw std::invalid_argument("expected cobo asad agetMin agetMax junctionCount ohmicCount; valid inclusive AGET range and nonnegative counts");
            fSideMultiplicityRules.push_back(rule);
        } catch (const std::exception &error) {
            lk_error << "Invalid " << sideParameter << ": " << raw << " (" << error.what() << ")" << endl;
            return false;
        }
    }
    const TString countParameter = "SKChannelFilterTask/Filter/FiredChannelCountRange";
    fFilterFiredChannelCount = fPar->CheckPar(countParameter);
    if (fFilterFiredChannelCount) {
        std::string value = fPar->GetParRaw(countParameter).Data();
        std::replace(value.begin(), value.end(), ',', ' ');
        std::istringstream input(value);
        std::string low, high, extra;
        try {
            if (!(input >> low >> high) || (input >> extra))
                throw std::invalid_argument("expected two nonnegative integers: min max");
            auto count = [](const std::string &token) {
                size_t used = 0;
                const int result = std::stoi(token, &used);
                if (used != token.size() || result < 0)
                    throw std::invalid_argument("counts must be nonnegative integers");
                return result;
            };
            fFiredChannelCountMin = count(low);
            fFiredChannelCountMax = count(high);
            if (fFiredChannelCountMin > fFiredChannelCountMax)
                throw std::invalid_argument("require min <= max");
        } catch (const std::exception &error) {
            lk_error << "Invalid " << countParameter << ": " << value << " (" << error.what() << ")" << endl;
            return false;
        }
    }
    const TString timeParameter = "SKChannelFilterTask/Filter/TimeRange";
    fFilterTime = fPar->CheckPar(timeParameter);
    if (fFilterTime) {
        if (fPar->GetParN(timeParameter) != 2 ||
            !fPar->GetParString(timeParameter, 0).IsFloat() ||
            !fPar->GetParString(timeParameter, 1).IsFloat()) {
            lk_error << timeParameter << " requires two numeric values: min max (inclusive TB)" << endl;
            return false;
        }
        fTimeMin = fPar->GetParDouble(timeParameter, 0);
        fTimeMax = fPar->GetParDouble(timeParameter, 1);
        if (!std::isfinite(fTimeMin) || !std::isfinite(fTimeMax) || fTimeMin > fTimeMax) {
            lk_error << "Invalid " << timeParameter << ": require finite min <= max" << endl;
            return false;
        }
    }
    const TString energyParameter = "SKChannelFilterTask/Filter/EnergyRange";
    fFilterEnergy = false;
    fEnergyRanges.clear();
    // Read every occurrence, including legacy entries retained for repeated names.
    // Two values set the default; six values set an address-specific override.
    for (int i=0; i<fPar->GetEntriesFast(); ++i) {
        auto par = fPar->GetParameter(i);
        if (!par || par->IsCommentOut() || TString(par->GetName()) != energyParameter) continue;
        // Preserve address separators: whole-line arithmetic evaluation can collapse -1 lists.
        std::string value = par->GetRaw().Data();
        std::replace(value.begin(), value.end(), ',', ' ');
        std::istringstream input(value);
        std::vector<std::string> tokens;
        for (std::string token; input >> token;) {
            TString expanded(token);
            fPar->ReplaceVariables(expanded);
            tokens.emplace_back(expanded.Data());
        }
        try {
            if (tokens.size() != 2 && tokens.size() != 6)
                throw std::invalid_argument("expected min max OR cobo asad aget channel min max");
            const size_t offset = tokens.size()-2;
            auto number = [](const std::string &token) {
                size_t used = 0;
                double result = std::stod(token, &used);
                if (used != token.size() || !std::isfinite(result))
                    throw std::invalid_argument("bounds must be finite numbers");
                return result;
            };
            const double minimum = number(tokens[offset]);
            const double maximum = number(tokens[offset+1]);
            if (minimum > maximum) throw std::invalid_argument("require min <= max");
            if (offset == 0) {
                fFilterEnergy = true;
                fEnergyMin = minimum;
                fEnergyMax = maximum;
            } else {
                EnergyRangeRule rule;
                for (int j=0; j<4; ++j) {
                    size_t used = 0;
                    rule.address[j] = std::stoi(tokens[j], &used);
                    if (used != tokens[j].size() || rule.address[j] < -1)
                        throw std::invalid_argument("addresses must be integers >= -1");
                }
                rule.minimum = minimum;
                rule.maximum = maximum;
                fEnergyRanges.push_back(rule);
            }
        } catch (const std::exception &error) {
            lk_error << "Invalid " << energyParameter << ": " << value << " (" << error.what() << ")" << endl;
            return false;
        }
    }
    fSaturationEventRules.clear();
    for (int i=0; i<fPar->GetEntriesFast(); ++i) {
        auto par = fPar->GetParameter(i);
        if (!par || par->IsCommentOut() || TString(par->GetName()) != "SKChannelFilterTask/Filter/SaturationEvent") continue;
        std::string value = par->GetRaw().Data();
        std::replace(value.begin(), value.end(), ',', ' ');
        std::istringstream input(value);
        std::vector<int> address;
        try {
            for (std::string token; input >> token;) {
                TString expanded(token);
                fPar->ReplaceVariables(expanded);
                token = expanded.Data();
                size_t used = 0;
                int id = std::stoi(token, &used);
                if (used != token.size() || id < -1)
                    throw std::invalid_argument("addresses must be integers >= -1");
                address.push_back(id);
            }
            if (address.size() != 3)
                throw std::invalid_argument("expected cobo asad aget");
            fSaturationEventRules.push_back(address);
        } catch (const std::exception &error) {
            lk_error << "Invalid SKChannelFilterTask/Filter/SaturationEvent: " << value << " (" << error.what() << ")" << endl;
            return false;
        }
    }
    if (!fSaturationEventRules.empty() && !fFitDataArray) {
        lk_info << "SaturationEvent: PFData absent; checking SiChannel waveforms (negative: ADC 0, positive: ADC 4095)" << endl;
    }
    fCrosstalkRules.clear();
    const std::regex syntax(R"(^\s*(-?\d+)\s+(-?\d+)\s+(-?\d+)\s*(>=|<=|==|!=|>|<|=)?\s*(\d+)\s*$)");
    // Iterate actual entries: repeated names are intentionally cumulative, including legacy entries.
    for (int i=0; i<fPar->GetEntriesFast(); ++i) {
        auto par = fPar->GetParameter(i);
        if (!par || par->IsCommentOut()) continue;
        TString name = par->GetName();
        if (name != "SKChannelFilterTask/Filter/Crosstalk" && name != "SKChannelFilterTask/Filter/Crosstask") continue;
        std::string value = par->GetValue().Data();
        std::replace(value.begin(), value.end(), ',', ' ');
        std::smatch match;
        if (!std::regex_match(value, match, syntax)) {
            lk_error << "Invalid " << name << ": " << value << endl;
            return false;
        }
        try {
            int cobo = std::stoi(match[1]), asad = std::stoi(match[2]), aget = std::stoi(match[3]);
            int count = std::stoi(match[5]);
            if (cobo < -1 || asad < -1 || aget < -1) return false;
            std::string op = match[4];
            // A bare count is the maximum allowed run length, like >count.
            int comparison = (op.empty() || op == ">") ? 0 : op == ">=" ? 1 :
                             op == "<" ? 2 : op == "<=" ? 3 : op == "!=" ? 5 : 4;
            fCrosstalkRules.push_back({cobo, asad, aget, comparison, count});
        } catch (const std::exception&) {
            lk_error << "Crosstalk integer out of range: " << value << endl;
            return false;
        }
    }
    return true;
}

void SKChannelFilterTask::Exec(Option_t*)
{
    const int n = fSiChannelArray->GetEntriesFast();
    // Inspect original channels before any channel or PFData compaction.
    if (!fSaturationEventRules.empty()) {
        std::vector<bool> saturated(n, false);
        if (fFitDataArray) {
            for (int i=0; i<fFitDataArray->GetEntriesFast(); ++i) {
                auto fit = static_cast<LKPulseFitData*>(fFitDataArray->At(i));
                if (fit && fit->fIsSaturated && fit->fHitIndex >= 0 && fit->fHitIndex < n)
                    saturated[fit->fHitIndex] = true;
            }
        }
        for (int i=0; i<n; ++i) {
            auto channel = static_cast<LKSiChannel*>(fSiChannelArray->At(i));
            if (!channel) continue;
            bool selected = false;
            for (const auto &rule : fSaturationEventRules) {
                if ((rule[0] != -1 && rule[0] != channel->GetCobo()) ||
                    (rule[1] != -1 && rule[1] != channel->GetAsad()) ||
                    (rule[2] != -1 && rule[2] != channel->GetAget())) continue;
                selected = true;
                break;
            }
            if (!selected) continue;
            if (!fFitDataArray) {
                // Older reconstructed files may not have set fInverted on junction channels.
                // LKSetSiChannelTask normal polarity: side 0 is negative, side 1 is positive.
                const bool inverted = channel->GetInverted() || channel->GetSide() == 0;
                const int rail = inverted ? 0 : 4095;
                auto waveform = channel->GetWaveformY();
                for (int tb=0; tb<512; ++tb) {
                    if (waveform[tb] == rail) { saturated[i] = true; break; }
                }
            }
            if (!saturated[i]) continue;
            fEventRun->DoNotFillCurrentEvent();
            fEventRun->StopCurrentEvent();
            lk_info << "Event rejected: saturated channel " << channel->GetCobo() << "/"
                    << channel->GetAsad() << "/" << channel->GetAget() << "/" << channel->GetChan() << endl;
            return;
        }
    }
    if (fFilterFiredChannelCount) {
        int fired = 0;
        for (int i=0; i<n; ++i) {
            auto channel = static_cast<LKSiChannel*>(fSiChannelArray->At(i));
            if (channel && std::isfinite(channel->GetEnergy()) && channel->GetEnergy() > fSignalThreshold)
                ++fired;
        }
        // This is an event cut on the original multiplicity, independent of channel cuts,
        // board grouping and adjacency. Do not reduce the count using other rejection masks.
        if (fired < fFiredChannelCountMin || fired > fFiredChannelCountMax) {
            fEventRun->DoNotFillCurrentEvent();
            fEventRun->StopCurrentEvent();
            lk_info << "Event rejected: " << fired << " fired channels outside ["
                    << fFiredChannelCountMin << ", " << fFiredChannelCountMax << "]" << endl;
            return;
        }
    }
    for (const auto &rule : fSideMultiplicityRules) {
        int counts[2] = {0, 0};
        for (int i=0; i<n; ++i) {
            auto channel = static_cast<LKSiChannel*>(fSiChannelArray->At(i));
            if (!channel || !std::isfinite(channel->GetEnergy()) || channel->GetEnergy() <= fSignalThreshold)
                continue;
            if ((rule[0] != -1 && channel->GetCobo() != rule[0]) ||
                (rule[1] != -1 && channel->GetAsad() != rule[1]) ||
                channel->GetAget() < rule[2] || channel->GetAget() > rule[3]) continue;
            ++counts[channel->GetSide() ? 1 : 0];
        }
        if (counts[0] != rule[4] || counts[1] != rule[5]) {
            fEventRun->DoNotFillCurrentEvent();
            fEventRun->StopCurrentEvent();
            lk_info << "Event rejected: selected junction/ohmic multiplicity "
                    << counts[0] << "/" << counts[1] << ", required " << rule[4] << "/" << rule[5] << endl;
            return;
        }
    }
    std::vector<bool> reject(n, false);
    using Group = std::tuple<int,int,int>;
    std::map<Group, std::map<int, std::vector<int>>> groups;
    for (int i=0; i<n; ++i) {
        auto channel = static_cast<LKSiChannel*>(fSiChannelArray->At(i));
        if (!channel) { reject[i] = true; continue; }
        const double signedNoise = channel->GetNoiseScale();
        const double noise = std::abs(signedNoise);
        const double energy = channel->GetEnergy();
        reject[i] = (fNoiseLevel >= 0 && (!std::isfinite(noise) || noise >= fNoiseLevel)) ||
                    (fBipolar && std::signbit(signedNoise));
        // A positive signal with zero baseline variance has infinite SNR.
        if (fSignalToNoiseRatio >= 0 &&
            (!std::isfinite(noise) || !std::isfinite(energy) || energy <= 0 ||
             (noise > 0 && energy/noise < fSignalToNoiseRatio))) reject[i] = true;
        bool filterEnergy = fFilterEnergy;
        double energyMin = fEnergyMin, energyMax = fEnergyMax;
        int bestSpecificity = -1;
        const int address[] = {channel->GetCobo(), channel->GetAsad(), channel->GetAget(), channel->GetChan()};
        for (const auto &rule : fEnergyRanges) {
            bool matches = true;
            int specificity = 0;
            for (int j=0; j<4; ++j) {
                if (rule.address[j] == -1) continue;
                ++specificity;
                if (rule.address[j] != address[j]) { matches = false; break; }
            }
            if (matches && specificity >= bestSpecificity) {
                bestSpecificity = specificity;
                filterEnergy = true;
                energyMin = rule.minimum;
                energyMax = rule.maximum;
            }
        }
        // Resolve overrides before applying the cut: an override can widen the default.
        if (filterEnergy && (!std::isfinite(energy) || energy < energyMin || energy > energyMax))
            reject[i] = true;
        const double time = channel->GetTime();
        if (fFilterTime && (!std::isfinite(time) || time < fTimeMin || time > fTimeMax))
            reject[i] = true;
        // Build from ALL original signals, even those rejected by another rule.
        if (std::isfinite(energy) && energy > fSignalThreshold)
            groups[{channel->GetCobo(),channel->GetAsad(),channel->GetAget()}][channel->GetChan()].push_back(i);
    }
    for (const auto &rule : fCrosstalkRules) {
        for (const auto &group : groups) {
            if ((rule[0] != -1 && rule[0] != std::get<0>(group.first)) ||
                (rule[1] != -1 && rule[1] != std::get<1>(group.first)) ||
                (rule[2] != -1 && rule[2] != std::get<2>(group.first))) continue;
            auto first = group.second.begin();
            while (first != group.second.end()) {
                auto end = std::next(first);
                int previous = first->first, count = 1;
                while (end != group.second.end() && end->first == previous+1) {
                    previous = end->first; ++count; ++end;
                }
                const int limit = rule[4];
                const bool bad = rule[3] == 0 ? count > limit : rule[3] == 1 ? count >= limit :
                                 rule[3] == 2 ? count < limit : rule[3] == 3 ? count <= limit :
                                 rule[3] == 4 ? count == limit : count != limit;
                if (bad)
                    for (auto it=first; it!=end; ++it)
                        for (int index : it->second) reject[index] = true;
                first = end;
            }
        }
    }
    std::vector<int> newIndex(n, -1);
    int kept = 0;
    for (int i=0; i<n; ++i) {
        if (reject[i]) fSiChannelArray->RemoveAt(i);
        else newIndex[i] = kept++;
    }
    fSiChannelArray->Compress();
    // Clear dangling pairs and rebuild using surviving channels and their new indices.
    for (int i=0; i<kept; ++i) {
        auto channel = static_cast<LKSiChannel*>(fSiChannelArray->At(i));
        channel->SetPairArrayIndex(-1);
        channel->SetPairChannel(nullptr);
        channel->SetEnergy2(-99.9);
    }
    for (int i=0; i<kept; ++i) {
        auto channel = static_cast<LKSiChannel*>(fSiChannelArray->At(i));
        if (!channel->IsPairedChannel() || channel->GetDirection() != 0) continue;
        for (int j=0; j<kept; ++j) {
            auto pair = static_cast<LKSiChannel*>(fSiChannelArray->At(j));
            if (i != j && channel->IsPair(pair)) {
                channel->SetPairArrayIndex(j);
                channel->SetPairChannel(pair);
                channel->SetEnergy2(pair->GetEnergy());
                break;
            }
        }
    }
    if (fFitDataArray) {
        for (int i=0; i<fFitDataArray->GetEntriesFast(); ++i) {
            auto fit = static_cast<LKPulseFitData*>(fFitDataArray->At(i));
            if (!fit) continue;
            int old = fit->fHitIndex;
            if (old < 0 || old >= n || newIndex[old] < 0) fFitDataArray->RemoveAt(i);
            else fit->fHitIndex = newIndex[old];
        }
        fFitDataArray->Compress();
    }
    if (kept < fMinChannels) {
        fEventRun->DoNotFillCurrentEvent();
        fEventRun->StopCurrentEvent();
    }
    lk_info << kept << "/" << n << " si-channels kept" << endl;
}
