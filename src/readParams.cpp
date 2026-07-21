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

        iforce = std::stoi(nextLine());
        ps     = std::stod(nextLine());
        rho    = std::stod(nextLine());
        mu     = std::stod(nextLine());
        mass   = std::stod(nextLine());
        c_sound= std::stod(nextLine());

        iforce  = std::stoi(nextLine());
        forcef  = std::stod(nextLine());
        famp    = std::stod(nextLine());
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
    os << "  ps      = " << ps << " [Pa]\n";
    os << "  rho     = " << rho << " [kg/m^3]\n";
    os << "  mu      = " << mu << " [Pa·s]\n";
    os << "  inputDir  = " << inputDir.string() << "\n";
    os << "  resultDir = " << resultDir.string() << "\n";
    os << "  freqFile  = " << freqFile.string() << "\n";
    os << "  modeFile  = " << modeFile.string() << "\n";
    os << "  surfFile  = " << surfFile.string() << "\n";
}
