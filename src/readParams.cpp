#include "SimulationParams.h"
#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>
#include <cctype>

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

        // Optional trailing values, in order:
        // forceDirection (0/1), contactReferenceFrequencyHz, flowBlendLengthMm.
        // A frequency can be supplied without forceDirection because it is
        // unambiguously not 0 or 1 in ordinary use.
        std::vector<std::string> optional;
        for (std::string value = nextLine(); !value.empty(); value = nextLine()) optional.push_back(value);
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
        if (optionalIndex != optional.size()) throw std::runtime_error("too many optional parameter values");
        if (legacyIforce != iforce) {
            std::cerr << "[Parameters] legacy iforce=" << legacyIforce
                      << " differs from effective iforce=" << iforce
                      << "; using the latter.\n";
        }
    } catch (...) {
        err = "Parse error (check file format)";
        return false;
    }

    return true;
}

bool SimulationParams::validate(std::string& err) const {
    if (nmode <= 0) { err = "nmode must be > 0"; return false; }
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
    // 追加チェック（例: ファイル/ディレクトリ存在確認を入れるならここ）
    return true;
}

void SimulationParams::print(std::ostream& os) const {
    os << "SimulationParams:\n";
    os << "  nmode   = " << nmode << "\n";
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
    os << "  ps      = " << ps << " [Pa]\n";
    os << "  rho     = " << rho << " [kg/m^3]\n";
    os << "  mu      = " << mu << " [Pa·s]\n";
    os << "  inputDir  = " << inputDir.string() << "\n";
    os << "  resultDir = " << resultDir.string() << "\n";
    os << "  freqFile  = " << freqFile.string() << "\n";
    os << "  modeFile  = " << modeFile.string() << "\n";
    os << "  surfFile  = " << surfFile.string() << "\n";
}
