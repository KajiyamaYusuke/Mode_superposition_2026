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
                if (key == "leftfrequencyfile") leftFrequencyFile = rhs;
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
                else if (key == "diagnosticoutputintervalsteps") {
                    diagnosticOutputIntervalSteps = parseIntStrict(rhs, key);
                }
                else if (key == "perturbationenabled") {
                    const int enabled = parseIntStrict(rhs, key);
                    if (enabled != 0 && enabled != 1) {
                        throw std::runtime_error("perturbationEnabled must be 0 or 1");
                    }
                    perturbationEnabled = enabled != 0;
                }
                else if (key == "perturbationtimesec") perturbationTimeSec = parseDoubleStrict(rhs, key);
                else if (key == "perturbationmodeindex") perturbationModeIndex = parseIntStrict(rhs, key);
                else if (key == "perturbationamplitudeatprobemm") {
                    perturbationAmplitudeAtProbeMm = parseDoubleStrict(rhs, key);
                }
                else if (key == "perturbationpattern") perturbationPattern = toLower(rhs);
                else if (key == "apcontactdiagnosticenabled") {
                    const int enabled = parseIntStrict(rhs, key);
                    if (enabled != 0 && enabled != 1)
                        throw std::runtime_error("apContactDiagnosticEnabled must be 0 or 1");
                    apContactDiagnosticEnabled = enabled != 0;
                }
                else if (key == "apcontactanalysisrowsperend")
                    apContactAnalysisRowsPerEnd = parseIntStrict(rhs, key);
                else if (key == "apcontactanalysisdistancemm")
                    apContactAnalysisDistanceMm = parseDoubleStrict(rhs, key);
                else if (key == "apcontactexcludedrowsperend")
                    apContactExcludedRowsPerEnd = parseIntStrict(rhs, key);
                else if (key == "apcontactexcludeddistancemm")
                    apContactExcludedDistanceMm = parseDoubleStrict(rhs, key);
                else if (key == "contactdetailstartsec")
                    contactDetailStartSec = parseDoubleStrict(rhs, key);
                else if (key == "contactdetailendsec")
                    contactDetailEndSec = parseDoubleStrict(rhs, key);
                else if (key == "diagnosticloadmode") diagnosticLoadMode = toLower(rhs);
                else if (key == "diagnosticreferencewindowstartsec")
                    diagnosticReferenceWindowStartSec = parseDoubleStrict(rhs, key);
                else if (key == "diagnosticreferencewindowendsec")
                    diagnosticReferenceWindowEndSec = parseDoubleStrict(rhs, key);
                else if (key == "fixednodeidsfile") fixedNodeIdsFile = rhs;
                else throw std::runtime_error("unknown named parameter: " + key);
            }
        }
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
    if (diagnosticOutputIntervalSteps <= 0) {
        err = "diagnosticOutputIntervalSteps must be > 0";
        return false;
    }
    if (!std::isfinite(perturbationTimeSec) || perturbationTimeSec < 0.0) {
        err = "perturbationTimeSec must be finite and >= 0";
        return false;
    }
    if (perturbationModeIndex < 1 || perturbationModeIndex > nmode) {
        err = "perturbationModeIndex must be between 1 and nmode";
        return false;
    }
    if (!std::isfinite(perturbationAmplitudeAtProbeMm)
        || perturbationAmplitudeAtProbeMm < 0.0) {
        err = "perturbationAmplitudeAtProbeMm must be finite and >= 0";
        return false;
    }
    if (perturbationPattern != "symmetric_opening") {
        err = "perturbationPattern must be symmetric_opening";
        return false;
    }
    if (apContactAnalysisRowsPerEnd < 0 || apContactExcludedRowsPerEnd < 0) {
        err = "AP contact row counts must be >= 0";
        return false;
    }
    if (!std::isfinite(apContactAnalysisDistanceMm)
        || !std::isfinite(apContactExcludedDistanceMm)
        || apContactAnalysisDistanceMm < 0.0
        || apContactExcludedDistanceMm < 0.0) {
        err = "AP contact distances must be finite and >= 0";
        return false;
    }
    if (apContactAnalysisRowsPerEnd > 0 && apContactAnalysisDistanceMm > 0.0) {
        err = "set only one of apContactAnalysisRowsPerEnd and apContactAnalysisDistanceMm";
        return false;
    }
    if (apContactExcludedRowsPerEnd > 0 && apContactExcludedDistanceMm > 0.0) {
        err = "set only one of apContactExcludedRowsPerEnd and apContactExcludedDistanceMm";
        return false;
    }
    if ((contactDetailStartSec < 0.0) != (contactDetailEndSec < 0.0)
        || (contactDetailStartSec >= 0.0 && contactDetailEndSec < contactDetailStartSec)) {
        err = "contact detail interval must be disabled with two negative values or have end >= start";
        return false;
    }
    if (diagnosticLoadMode != "live" && diagnosticLoadMode != "freeze_contact"
        && diagnosticLoadMode != "freeze_fluid" && diagnosticLoadMode != "freeze_all") {
        err = "diagnosticLoadMode must be live, freeze_contact, freeze_fluid, or freeze_all";
        return false;
    }
    if (!std::isfinite(diagnosticReferenceWindowStartSec)
        || !std::isfinite(diagnosticReferenceWindowEndSec)
        || diagnosticReferenceWindowStartSec < 0.0
        || diagnosticReferenceWindowEndSec <= diagnosticReferenceWindowStartSec) {
        err = "diagnostic reference window must be finite with 0 <= start < end";
        return false;
    }
    if (diagnosticLoadMode != "live"
        && diagnosticReferenceWindowEndSec > perturbationTimeSec) {
        err = "diagnostic reference window must end no later than perturbationTimeSec for frozen loads";
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

void SimulationParams::print(std::ostream& os) const {
    os << "SimulationParams:\n";
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
    os << "  pressureRampTimeSec = " << pressureRampTimeSec
       << (pressureRampTimeSecExplicit ? " [s] (explicit)\n" : " [s] (legacy omitted option)\n");
    os << "  diagnosticOutputIntervalSteps = " << diagnosticOutputIntervalSteps << "\n";
    os << "  perturbationEnabled = " << perturbationEnabled << "\n";
    os << "  perturbationTimeSec = " << perturbationTimeSec << " [s]\n";
    os << "  perturbationModeIndex = " << perturbationModeIndex << "\n";
    os << "  perturbationAmplitudeAtProbeMm = "
       << perturbationAmplitudeAtProbeMm << " [mm]\n";
    os << "  perturbationPattern = " << perturbationPattern << "\n";
    os << "  apContactDiagnosticEnabled = " << apContactDiagnosticEnabled << "\n";
    os << "  apContactAnalysisRowsPerEnd = " << apContactAnalysisRowsPerEnd << "\n";
    os << "  apContactAnalysisDistanceMm = " << apContactAnalysisDistanceMm << " [mm]\n";
    os << "  apContactExcludedRowsPerEnd = " << apContactExcludedRowsPerEnd << "\n";
    os << "  apContactExcludedDistanceMm = " << apContactExcludedDistanceMm << " [mm]\n";
    os << "  contactDetailIntervalSec = [" << contactDetailStartSec << ", "
       << contactDetailEndSec << "]\n";
    os << "  diagnosticLoadMode = " << diagnosticLoadMode << "\n";
    os << "  diagnosticReferenceWindowSec = ["
       << diagnosticReferenceWindowStartSec << ", "
       << diagnosticReferenceWindowEndSec << "]\n";
    os << "  fixedNodeIdsFile = " << fixedNodeIdsFile.string() << "\n";
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
