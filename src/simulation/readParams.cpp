#include "SimulationParams.h"
#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>
#include <cctype>
#include <cmath>
#include <stdexcept>

// トリム関数（先頭・末尾の空白除去）
static inline void trim(std::string &s) {
    auto notspace = [](int ch){ return !std::isspace(ch); };
    s.erase(s.begin(), std::find_if(s.begin(), s.end(), notspace));
    s.erase(std::find_if(s.rbegin(), s.rend(), notspace).base(), s.end());
}

static inline std::string toLower(std::string s) {
    std::transform(s.begin(), s.end(), s.begin(), [](unsigned char c){ return std::tolower(c); });
    return s;
}

static double parseDoubleStrict(const std::string& value, const std::string& key) {
    std::size_t used = 0;
    double parsed = 0.0;
    try {
        parsed = std::stod(value, &used);
    } catch (const std::exception&) {
        throw std::runtime_error("invalid finite numeric value for " + key + ": " + value);
    }
    if (used != value.size() || !std::isfinite(parsed)) {
        throw std::runtime_error("invalid finite numeric value for " + key + ": " + value);
    }
    return parsed;
}

static int parseIntStrict(const std::string& value, const std::string& key) {
    std::size_t used = 0;
    int parsed = 0;
    try {
        parsed = std::stoi(value, &used);
    } catch (const std::exception&) {
        throw std::runtime_error("invalid integer value for " + key + ": " + value);
    }
    if (used != value.size()) {
        throw std::runtime_error("invalid integer value for " + key + ": " + value);
    }
    return parsed;
}

void SimulationParams::setUpstreamAcousticProfile(const std::string& profile) {
    if (profile != "modern" && profile != "legacy_main")
        throw std::runtime_error("unknown upstreamAcousticProfile: " + profile);
    upstreamAcousticProfile = profile;
    const SimulationParams defaults;
    r_inlet = profile == "legacy_main" ? 0.15 : defaults.r_inlet;
    L_inlet = profile == "legacy_main" ? 0.50 : defaults.L_inlet;
    L_sub = profile == "legacy_main" ? 0.25 : defaults.L_sub;
    N_sub = profile == "legacy_main" ? 3 : defaults.N_sub;
}

bool SimulationParams::loadFromFile(const fs::path& filename, std::string& err) {
    std::ifstream ifs(filename);
    if (!ifs) { err = "Cannot open parameter file"; return false; }

    std::string line;
    std::vector<std::string> tokens;
    auto nextLine = [&]() -> std::string {
        while (std::getline(ifs, line)) {
            trim(line);
            if (line.empty() || line[0] == '#') continue;
            return line;
        }
        return "";
    };

    try {
        nmode  = std::stoi(nextLine());
        nsurfz = std::stoi(nextLine());
        nstep  = std::stoi(nextLine());
        nwrite = std::stoi(nextLine());
        dt     = std::stod(nextLine());
        zetaL   = std::stod(nextLine());
        zetaR   = std::stod(nextLine());

        std::istringstream iss(nextLine());
        iss >> kc1 >> kc2 >> kc3;

        ncont = std::stoi(nextLine());

        freqFile  = nextLine();
        modeFile  = nextLine();
        surfFile  = nextLine();
        inputDir  = nextLine();
        resultDir = nextLine();

        // The historical format contains this first iforce slot before the
        // acoustic constants and a second, effective slot below.  Keep
        // consuming it so existing files remain readable, but do not let an
        // accidental duplicate silently obscure which value drives the run.
        const int legacyIforce = std::stoi(nextLine());
        ps     = std::stod(nextLine());
        rho    = std::stod(nextLine());
        mu     = std::stod(nextLine());
        mass   = std::stod(nextLine());
        c_sound= std::stod(nextLine());

        iforce  = std::stoi(nextLine());
        forcef  = std::stod(nextLine());
        famp    = std::stod(nextLine());

        modeSelection.clear();

        std::set<std::string> namedParameters;
        std::string glottalProfile = "modern";
        std::string upstreamProfile = "modern";
        std::map<std::string, bool> legacyOverrides;

        // Optional trailing values, in order:
        // forceDirection (0/1), contactReferenceFrequencyHz,
        // flowBlendLengthMm, flowSeparationAreaRatio.
        // A frequency can be supplied without forceDirection because it is
        // unambiguously not 0 or 1 in ordinary use.
        // Named values may be placed anywhere in this trailing section without
        // shifting the historical positional values.
        std::vector<std::string> optional;
        for (std::string value = nextLine(); !value.empty(); value = nextLine()) {
            const auto equal = value.find('=');
            if (equal == std::string::npos) {
                optional.push_back(value);
                continue;
            }

            std::string key = value.substr(0, equal);
            std::string rhs = value.substr(equal + 1);
            trim(key);
            trim(rhs);
            key = toLower(key);
            key.erase(std::remove_if(key.begin(), key.end(), [](char c) {
                return c == '_' || c == '-';
            }), key.end());
            if (!namedParameters.insert(key).second) {
                throw std::runtime_error(key + " was specified more than once");
            }
            if (key == "modeselection") {
                for (char& c : rhs) {
                    if (c == ',' || c == ';' || c == '[' || c == ']') c = ' ';
                }
                std::istringstream selected(rhs);
                int modeNumber = 0;
                while (selected >> modeNumber) modeSelection.push_back(modeNumber);
                if (modeSelection.empty() || (!selected.eof() && selected.fail())) {
                    throw std::runtime_error("invalid modeSelection");
                }
            } else {
                if (rhs.empty()) throw std::runtime_error(key + " must not be empty");
                if (key == "glottalmethodprofile") glottalProfile = toLower(rhs);
                else if (key == "upstreamacousticprofile") upstreamProfile = toLower(rhs);
                else if (key == "legacyarea" || key == "legacyflowcoefficients"
                    || key == "legacywallviscosity" || key == "legacyseparation"
                    || key == "legacypostseparation" || key == "legacysurfacepressureload"
                    || key == "legacycontact" || key == "legacyintegrator") {
                    const int enabled = parseIntStrict(rhs, key);
                    if (enabled != 0 && enabled != 1) throw std::runtime_error(key + " must be 0 or 1");
                    legacyOverrides[key] = enabled != 0;
                }
                else if (key == "writedeformationvtu") {
                    const int enabled = parseIntStrict(rhs, key);
                    if (enabled != 0 && enabled != 1) throw std::runtime_error(key + " must be 0 or 1");
                    writeDeformationVtu = enabled != 0;
                }
                else if (key == "vocaltractlengthm") L_vt = parseDoubleStrict(rhs,key);
                else if (key == "leftfrequencyfile") leftFrequencyFile = rhs;
                else if (key == "rightfrequencyfile") rightFrequencyFile = rhs;
                else if (key == "leftmodefile") leftModeFile = rhs;
                else if (key == "rightmodefile") rightModeFile = rhs;
                else if (key == "leftsurfacenasfile" || key == "leftsurfacefile") leftSurfaceNasFile = rhs;
                else if (key == "rightsurfacenasfile" || key == "rightsurfacefile") rightSurfaceNasFile = rhs;
                else if (key == "initialgapmm") initialGapMm = parseDoubleStrict(rhs, key);
                else if (key == "pressureramptimesec") {
                    pressureRampTimeSec = parseDoubleStrict(rhs, key);
                    pressureRampTimeSecExplicit = true;
                }
                else if (key == "pressureprofile") pressureProfile = toLower(rhs);
                else if (key == "pressureinitialpa") {
                    pressureInitialPa = parseDoubleStrict(rhs, key);
                    pressureInitialPaExplicit = true;
                }
                else if (key == "pressuretargetpa") {
                    pressureTargetPa = parseDoubleStrict(rhs, key);
                    pressureTargetPaExplicit = true;
                }
                else if (key == "pressurefinalpa") {
                    pressureFinalPa = parseDoubleStrict(rhs, key);
                    pressureFinalPaExplicit = true;
                }
                else if (key == "pressurestartsec") pressureStartSec = parseDoubleStrict(rhs, key);
                else if (key == "pressurerisetimesec") {
                    pressureRiseTimeSec = parseDoubleStrict(rhs, key);
                    pressureRiseTimeSecExplicit = true;
                }
                else if (key == "pressureholdtimesec") {
                    pressureHoldTimeSec = parseDoubleStrict(rhs, key);
                    pressureHoldTimeSecExplicit = true;
                }
                else if (key == "pressurefalltimesec") {
                    pressureFallTimeSec = parseDoubleStrict(rhs, key);
                    pressureFallTimeSecExplicit = true;
                }
                else if (key == "apcontactdiagnosticenabled") {
                    const int enabled = parseIntStrict(rhs, key);
                    if (enabled != 0 && enabled != 1)
                        throw std::runtime_error("apContactDiagnosticEnabled must be 0 or 1");
                    apContactDiagnosticEnabled = enabled != 0;
                }
                else if (key == "apcontactexcludedrowsperend")
                    apContactExcludedRowsPerEnd = parseIntStrict(rhs, key);
                else if (key == "apcontactexcludeddistancemm")
                    apContactExcludedDistanceMm = parseDoubleStrict(rhs, key);
                else throw std::runtime_error("unknown named parameter: " + key);
            }
        }
        glottalMethods.configure(glottalProfile, legacyOverrides);
        setUpstreamAcousticProfile(upstreamProfile);
        std::size_t optionalIndex = 0;
        if (!optional.empty()) {
            const double first = std::stod(optional[0]);
            if (std::abs(first) < 1.0e-12 || std::abs(first - 1.0) < 1.0e-12) {
                forceDirection = static_cast<int>(first);
                optionalIndex = 1;
            }
        }
        if (optionalIndex < optional.size()) contactReferenceFrequencyHz = std::stod(optional[optionalIndex++]);
        if (optionalIndex < optional.size()) flowBlendLengthMm = std::stod(optional[optionalIndex++]);
        if (optionalIndex < optional.size()) flowSeparationAreaRatio = std::stod(optional[optionalIndex++]);
        if (optionalIndex != optional.size()) throw std::runtime_error("too many optional parameter values");
        if (legacyIforce != iforce) {
            std::cerr << "[Parameters] legacy iforce=" << legacyIforce
                      << " differs from effective iforce=" << iforce
                      << "; using the latter.\n";
        }
    } catch (const std::exception& exception) {
        err = std::string("Parse error: ") + exception.what();
        return false;
    }

    return true;
}

bool SimulationParams::validate(std::string& err) const {
    if (upstreamAcousticProfile != "modern" && upstreamAcousticProfile != "legacy_main") {
        err = "unknown upstreamAcousticProfile"; return false;
    }
    if (upstreamAcousticProfile == "legacy_main"
        && (r_inlet != 0.15 || L_inlet != 0.50 || L_sub != 0.25 || N_sub != 3)) {
        err = "legacy_main requires its complete upstream circuit preset"; return false;
    }
    if (glottalMethods.any() && iforce != 0) {
        err = "legacy fluid methods require iforce=0"; return false;
    }
    if (glottalMethods.any() && zetaL != zetaR) {
        err = "legacy comparisons require equal left/right damping"; return false;
    }
    if (nmode <= 0) { err = "nmode must be > 0"; return false; }
    std::set<int> selectedModes;
    for (int modeNumber : modeSelection) {
        if (modeNumber < 1 || modeNumber > nmode) {
            err = "modeSelection entries must be between 1 and nmode";
            return false;
        }
        if (!selectedModes.insert(modeNumber).second) {
            err = "modeSelection must not contain duplicate mode numbers";
            return false;
        }
    }
    if (nstep <= 0) { err = "nstep must be > 0"; return false; }
    if (dt <= 0.0)  { err = "dt must be > 0"; return false; }
    if (nwrite <= 0){ err = "nwrite must be > 0"; return false; }
    if (ncont < 0)  { err = "ncont must be >= 0"; return false; }
    if (N_sub < 1)  { err = "N_sub must be >= 1"; return false; }
    if (N_vt < 1)   { err = "N_vt must be >= 1"; return false; }
    if (!std::isfinite(L_vt) || L_vt < 0.) { err = "vocalTractLengthM must be finite and >=0"; return false; }
    if (iforce < 0 || iforce > 1) {
        err = "iforce must be 0 (flow) or 1 (prescribed force)";
        return false;
    }
    if (forceDirection < 0 || forceDirection > 1) {
        err = "forceDirection must be 0 (x) or 1 (opposing y)";
        return false;
    }
    if (contactReferenceFrequencyHz <= 0.0) {
        err = "contactReferenceFrequencyHz must be > 0";
        return false;
    }
    if (flowBlendLengthMm <= 0.0) {
        err = "flowBlendLengthMm must be > 0";
        return false;
    }
    if (!std::isfinite(flowSeparationAreaRatio)
        || flowSeparationAreaRatio < 1.0) {
        err = "flowSeparationAreaRatio must be finite and >= 1.0";
        return false;
    }
    if (!std::isfinite(initialGapMm)) {
        err = "initialGapMm must be finite";
        return false;
    }
    if (!std::isfinite(pressureRampTimeSec) || pressureRampTimeSec <= 0.0) {
        err = "pressureRampTimeSec must be finite and > 0";
        return false;
    }
    if (pressureProfile != "legacy_cosine" && pressureProfile != "linear_sweep"
        && pressureProfile != "linear_cycle") {
        err = "unknown pressureProfile"; return false;
    }
    const bool hasSweepFields = pressureInitialPaExplicit || pressureTargetPaExplicit
        || pressureFinalPaExplicit || pressureStartSec != 0.0
        || pressureRiseTimeSecExplicit || pressureHoldTimeSecExplicit
        || pressureFallTimeSecExplicit;
    if (pressureProfile == "legacy_cosine") {
        if (hasSweepFields) {
            err = "pressure sweep parameters require linear_sweep or linear_cycle"; return false;
        }
    } else {
        if (iforce != 0) {
            err = "time-dependent pressure requires iforce=0"; return false;
        }
        if (pressureRampTimeSecExplicit) {
            err = "pressureRampTimeSec is exclusive to legacy_cosine"; return false;
        }
        if (!pressureInitialPaExplicit || !pressureTargetPaExplicit
            || !pressureRiseTimeSecExplicit) {
            err = "pressureInitialPa, pressureTargetPa and pressureRiseTimeSec are required"; return false;
        }
        if (!std::isfinite(pressureInitialPa) || !std::isfinite(pressureTargetPa)
            || !std::isfinite(pressureFinalPa) || !std::isfinite(pressureStartSec)
            || !std::isfinite(pressureRiseTimeSec) || !std::isfinite(pressureHoldTimeSec)
            || !std::isfinite(pressureFallTimeSec)
            || pressureInitialPa < 0.0 || pressureTargetPa < 0.0
            || pressureStartSec < 0.0 || pressureRiseTimeSec <= 0.0) {
            err = "invalid pressure sweep values or durations"; return false;
        }
        if (pressureProfile == "linear_sweep") {
            if (pressureFinalPaExplicit || pressureHoldTimeSecExplicit
                || pressureFallTimeSecExplicit) {
                err = "linear_sweep does not accept final/hold/fall options"; return false;
            }
        } else {
            if (!pressureFinalPaExplicit || !pressureFallTimeSecExplicit
                || pressureFinalPa < 0.0 || pressureHoldTimeSec < 0.0
                || pressureFallTimeSec <= 0.0
                || pressureTargetPa < std::max(pressureInitialPa, pressureFinalPa)) {
                err = "linear_cycle requires valid peak/final/hold/fall values"; return false;
            }
        }
        const double endSec = pressureStartSec + pressureRiseTimeSec
            + (pressureProfile == "linear_cycle"
                ? pressureHoldTimeSec + pressureFallTimeSec : 0.0);
        if (!std::isfinite(endSec)) {
            err = "pressure sweep schedule has nonfinite end time"; return false;
        }
        const double largestPressure = std::max(pressureInitialPa,
            std::max(pressureTargetPa, pressureFinalPa));
        // The empirical upstream loss alpha1 must not become negative.
        if (0.185 - 2.5e-5 * largestPressure < 0.0) {
            err = "pressure sweep makes alpha1 negative (empirical model limit)";
            return false;
        }
    }
    if (apContactExcludedRowsPerEnd < 0) {
        err = "AP contact row counts must be >= 0";
        return false;
    }
    if (!std::isfinite(apContactExcludedDistanceMm)
        || apContactExcludedDistanceMm < 0.0) {
        err = "AP contact distances must be finite and >= 0";
        return false;
    }
    if (apContactExcludedRowsPerEnd > 0 && apContactExcludedDistanceMm > 0.0) {
        err = "set only one of apContactExcludedRowsPerEnd and apContactExcludedDistanceMm";
        return false;
    }
    if (leftFrequencyFile.empty() || rightFrequencyFile.empty()
        || leftModeFile.empty() || rightModeFile.empty()
        || leftSurfaceNasFile.empty() || rightSurfaceNasFile.empty()) {
        err = "left/right frequency, mode, and surface NAS files must not be empty";
        return false;
    }
    // 追加チェック（例: ファイル/ディレクトリ存在確認を入れるならここ）
    return true;
}

double SimulationParams::pressureCommandAt(double time) const {
    if (!std::isfinite(time) || time < 0.0)
        throw std::runtime_error("pressureCommandAt requires nonnegative finite time");
    if (pressureProfile == "legacy_cosine") {
        const double rampTime = pressureRampTimeSecExplicit ? pressureRampTimeSec : 0.05;
        const double ratio = time < rampTime
            ? 0.5 * (1.0 - std::cos(std::acos(-1.0) * time / rampTime)) : 1.0;
        return ps * ratio;
    }
    const auto interpolate = [time](double from, double to, double start, double duration) {
        if (time <= start) return from;
        if (time >= start + duration) return to;
        const double fraction = (time - start) / duration;
        return from + (to - from) * fraction;
    };
    if (pressureProfile == "linear_sweep")
        return interpolate(pressureInitialPa, pressureTargetPa,
                           pressureStartSec, pressureRiseTimeSec);
    if (pressureProfile == "linear_cycle") {
        const double peakTime = pressureStartSec + pressureRiseTimeSec;
        const double fallStart = peakTime + pressureHoldTimeSec;
        if (time < peakTime)
            return interpolate(pressureInitialPa, pressureTargetPa,
                               pressureStartSec, pressureRiseTimeSec);
        if (time < fallStart) return pressureTargetPa;
        return interpolate(pressureTargetPa, pressureFinalPa,
                           fallStart, pressureFallTimeSec);
    }
    throw std::runtime_error("unknown pressure profile");
}

void SimulationParams::print(std::ostream& os) const {
    os << "SimulationParams:\n";
    os << "  glottalMethodProfile = " << glottalMethods.profile << "\n";
    os << "  nmode   = " << nmode << "\n";
    os << "  modeSelection = ";
    if (modeSelection.empty()) {
        os << "all modes 1..nmode";
    } else {
        for (std::size_t i = 0; i < modeSelection.size(); ++i) {
            if (i) os << ",";
            os << modeSelection[i];
        }
    }
    os << "\n";
    os << "  nsurfz  = " << nsurfz << "\n";
    os << "  nstep   = " << nstep << "\n";
    os << "  nwrite  = " << nwrite << "\n";
    os << "  dt      = " << dt << " [s]\n";
    os << "  zetaL    = " << zetaL << "\n";
    os << "  zetaR    = " << zetaR << "\n";
    os << "  kc1     = " << kc1 << "\n";
    os << "  kc2     = " << kc2 << "\n";
    os << "  mass    = " << mass << "\n";
    os << "  iforce  = " << iforce << "\n";
    os << "  forcef  = " << forcef << "\n";
    os << "  famp    = " << famp << "\n";
    os << "  forceDirection = " << forceDirection << "\n";
    os << "  contactReferenceFrequencyHz = " << contactReferenceFrequencyHz << "\n";
    os << "  flowBlendLengthMm = " << flowBlendLengthMm << "\n";
    os << "  flowSeparationAreaRatio = " << flowSeparationAreaRatio << "\n";
    os << "  initialGapMm = " << initialGapMm << " [mm]\n";
    os << "  pressureProfile = " << pressureProfile << "\n";
    if (pressureProfile != "legacy_cosine") {
        os << "  pressureInitialPa = " << pressureInitialPa << "\n"
           << "  pressureTargetPa = " << pressureTargetPa << "\n"
           << "  pressureFinalPa = " << pressureFinalPa << "\n"
           << "  pressureStartSec = " << pressureStartSec << "\n"
           << "  pressureRiseTimeSec = " << pressureRiseTimeSec << "\n"
           << "  pressureHoldTimeSec = " << pressureHoldTimeSec << "\n"
           << "  pressureFallTimeSec = " << pressureFallTimeSec << "\n";
    }
    os << "  pressureRampTimeSec = " << pressureRampTimeSec
       << (pressureRampTimeSecExplicit ? " [s] (explicit)\n" : " [s] (legacy omitted option)\n");
    os << "  apContactDiagnosticEnabled = " << apContactDiagnosticEnabled << "\n";
    os << "  apContactExcludedRowsPerEnd = " << apContactExcludedRowsPerEnd << "\n";
    os << "  apContactExcludedDistanceMm = " << apContactExcludedDistanceMm << " [mm]\n";
    os << "  ps      = " << ps << " [Pa]\n";
    os << "  rho     = " << rho << " [kg/m^3]\n";
    os << "  mu      = " << mu << " [Pa·s]\n";
    os << "  inputDir  = " << inputDir.string() << "\n";
    os << "  resultDir = " << resultDir.string() << "\n";
    os << "  freqFile  = " << freqFile.string() << "\n";
    os << "  modeFile  = " << modeFile.string() << "\n";
    os << "  surfFile  = " << surfFile.string() << "\n";
    os << "  leftFrequencyFile  = " << leftFrequencyFile.string() << "\n";
    os << "  rightFrequencyFile = " << rightFrequencyFile.string() << "\n";
    os << "  leftModeFile       = " << leftModeFile.string() << "\n";
    os << "  rightModeFile      = " << rightModeFile.string() << "\n";
    os << "  leftSurfaceNasFile  = " << leftSurfaceNasFile.string() << "\n";
    os << "  rightSurfaceNasFile = " << rightSurfaceNasFile.string() << "\n";
}
