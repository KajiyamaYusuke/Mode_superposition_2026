#include "Geometry.h"
#include <fstream>
#include <sstream>
#include <iostream>
#include <algorithm>
#include <stdexcept>
#include <vector>
#include <cmath>
#include <set>
#include <iomanip>
#include <unordered_map>
#include <unordered_set>
#include <queue>
#include <limits>
#include <numeric>

inline bool isSamePointRounded(double x1, double y1, double z1,
                               double x2, double y2, double z2,
                               double eps = 1e-3) {
    // 2点間のユークリッド距離の2乗が eps の2乗より小さければ「同じ点」とみなす
    double dx = x1 - x2;
    double dy = y1 - y2;
    double dz = z1 - z2;
    return (dx * dx + dy * dy + dz * dz) < (eps * eps);
}

static std::string trim(const std::string& s) {
    size_t b = s.find_first_not_of(" \t\r\n");
    size_t e = s.find_last_not_of(" \t\r\n");
    if (b == std::string::npos) return "";
    return s.substr(b, e - b + 1);
}

static std::vector<std::string> splitComma(const std::string& line) {
    std::vector<std::string> tokens;
    std::stringstream ss(line);
    std::string t;

    while (std::getline(ss, t, ',')) {
        tokens.push_back(trim(t));
    }

    return tokens;
}

static std::vector<double> makeLevels(std::vector<double> values, double tol) {
    std::sort(values.begin(), values.end());

    std::vector<double> levels;
    std::vector<double> cluster;

    for (double v : values) {
        if (cluster.empty()) {
            cluster.push_back(v);
        } else if (std::abs(v - cluster.back()) <= tol) {
            cluster.push_back(v);
        } else {
            double sum = 0.0;
            for (double c : cluster) sum += c;
            levels.push_back(sum / cluster.size());

            cluster.clear();
            cluster.push_back(v);
        }
    }

    if (!cluster.empty()) {
        double sum = 0.0;
        for (double c : cluster) sum += c;
        levels.push_back(sum / cluster.size());
    }

    return levels;
}

static int nearestLevel(const std::vector<double>& levels, double v, double tol) {
    int best = -1;
    double bestDist = std::numeric_limits<double>::max();

    for (int i = 0; i < static_cast<int>(levels.size()); ++i) {
        double d = std::abs(levels[i] - v);
        if (d < bestDist) {
            bestDist = d;
            best = i;
        }
    }

    if (bestDist > tol) return -1;
    return best;
}

static std::vector<int> traceAllComponentsAsOneRow(
    const std::unordered_map<int, std::set<int>>& adj,
    const std::vector<int>& allNodes,
    const std::vector<Point>& points)
{
    std::vector<int> result;

    if (allNodes.empty()) return result;

    std::unordered_set<int> allSet(allNodes.begin(), allNodes.end());
    std::unordered_set<int> visitedAll;

    std::vector<std::vector<int>> components;

    // ------------------------
    // 1. 連結成分を抽出
    // ------------------------
    for (int start : allNodes) {
        if (visitedAll.count(start)) continue;

        std::vector<int> comp;
        std::vector<int> stack;

        stack.push_back(start);
        visitedAll.insert(start);

        while (!stack.empty()) {
            int cur = stack.back();
            stack.pop_back();

            comp.push_back(cur);

            auto it = adj.find(cur);
            if (it == adj.end()) continue;

            for (int nb : it->second) {
                if (!allSet.count(nb)) continue;

                if (!visitedAll.count(nb)) {
                    visitedAll.insert(nb);
                    stack.push_back(nb);
                }
            }
        }

        components.push_back(comp);
    }

    if (components.empty()) return result;

    if (components.size() > 1) {
        std::cout << "[trace WARN] components=" << components.size()
                  << " sizes=";
        for (const auto& c : components) {
            std::cout << c.size() << " ";
        }
        std::cout << "\n";
    }

    auto dist2 = [&](int a, int b) {
        double dx = points[a].x - points[b].x;
        double dy = points[a].y - points[b].y;
        double dz = points[a].z - points[b].z;
        return dx * dx + dy * dy + dz * dz;
    };

    // ------------------------
    // 2. 各成分を1本のポリラインとして並べる
    // ------------------------
    std::vector<std::vector<int>> orderedComps;

    for (const auto& comp : components) {
        std::unordered_set<int> compSet(comp.begin(), comp.end());

        std::vector<int> endpoints;

        for (int p : comp) {
            int deg = 0;

            auto it = adj.find(p);
            if (it != adj.end()) {
                for (int nb : it->second) {
                    if (compSet.count(nb)) ++deg;
                }
            }

            if (deg <= 1) {
                endpoints.push_back(p);
            }
        }

        int start = comp[0];

        if (!endpoints.empty()) {
            start = endpoints[0];

            // とりあえず x が小さい端点から開始
            for (int p : endpoints) {
                if (points[p].x < points[start].x) {
                    start = p;
                }
            }
        } else {
            // 閉ループなど。x最小点から開始
            for (int p : comp) {
                if (points[p].x < points[start].x) {
                    start = p;
                }
            }
        }

        std::vector<int> ordered;
        std::unordered_set<int> visited;

        int prev = -1;
        int cur = start;

        while (cur >= 0) {
            ordered.push_back(cur);
            visited.insert(cur);

            auto it = adj.find(cur);
            if (it == adj.end()) break;

            int next = -1;
            double bestD2 = std::numeric_limits<double>::max();

            for (int nb : it->second) {
                if (!compSet.count(nb)) continue;
                if (nb == prev) continue;
                if (visited.count(nb)) continue;

                double d = dist2(cur, nb);
                if (d < bestD2) {
                    bestD2 = d;
                    next = nb;
                }
            }

            prev = cur;
            cur = next;
        }

        // 取りこぼした孤立点・分岐点があれば末尾に追加
        if (ordered.size() != comp.size()) {
            for (int p : comp) {
                if (!visited.count(p)) {
                    ordered.push_back(p);
                }
            }
        }

        orderedComps.push_back(ordered);
    }

    // ------------------------
    // 3. 成分同士を端点距離が近い順につなぐ
    // ------------------------
    std::vector<bool> used(orderedComps.size(), false);

    // 最初の成分：先頭点の x が一番小さいもの
    int currentComp = -1;
    double minX = std::numeric_limits<double>::max();

    for (int c = 0; c < static_cast<int>(orderedComps.size()); ++c) {
        if (orderedComps[c].empty()) continue;

        int p0 = orderedComps[c].front();
        if (points[p0].x < minX) {
            minX = points[p0].x;
            currentComp = c;
        }
    }

    if (currentComp < 0) return result;

    result = orderedComps[currentComp];
    used[currentComp] = true;

    while (true) {
        int last = result.back();

        int bestComp = -1;
        bool reverseBest = false;
        double bestD2 = std::numeric_limits<double>::max();

        for (int c = 0; c < static_cast<int>(orderedComps.size()); ++c) {
            if (used[c]) continue;
            if (orderedComps[c].empty()) continue;

            int front = orderedComps[c].front();
            int back  = orderedComps[c].back();

            double dFront = dist2(last, front);
            double dBack  = dist2(last, back);

            if (dFront < bestD2) {
                bestD2 = dFront;
                bestComp = c;
                reverseBest = false;
            }

            if (dBack < bestD2) {
                bestD2 = dBack;
                bestComp = c;
                reverseBest = true;
            }
        }

        if (bestComp < 0) break;

        if (reverseBest) {
            std::reverse(orderedComps[bestComp].begin(),
                         orderedComps[bestComp].end());
        }

        result.insert(result.end(),
                      orderedComps[bestComp].begin(),
                      orderedComps[bestComp].end());

        used[bestComp] = true;
    }

    return result;
}

struct gridpoint {
    double x, y, z;
};

void Geometry::loadFromVTK(const std::string& filename) {
    std::ifstream file(filename);
    if (!file) throw std::runtime_error("Cannot open VTK file: " + filename);

    std::string line;
    std::vector<double> pointsBuffer;
    std::vector<int> connectivityBuffer;
    std::vector<int> offsetsBuffer;
    std::vector<int> typesBuffer;

    while (std::getline(file, line)) {
        // Points
        if (line.find("<Points>") != std::string::npos) {
            while (std::getline(file, line) && line.find("</Points>") == std::string::npos) {
                if (line.find("<DataArray") != std::string::npos) {
                    while (std::getline(file, line) && line.find("</DataArray>") == std::string::npos) {
                        std::istringstream ss(line);
                        double x, y, z;
                        while (ss >> y >> z >> x) {
                            pointsBuffer.push_back(x);
                            pointsBuffer.push_back(y);
                            pointsBuffer.push_back(z);
                        }
                    }
                }
            }
        }

        // Cells
        else if (line.find("<Cells>") != std::string::npos) {
            while (std::getline(file, line) && line.find("</Cells>") == std::string::npos) {
                if (line.find("connectivity") != std::string::npos) {
                    while (std::getline(file, line) && line.find("</DataArray>") == std::string::npos) {
                        std::istringstream ss(line);
                        int id;
                        while (ss >> id) connectivityBuffer.push_back(id);
                    }
                } else if (line.find("offsets") != std::string::npos) {
                    while (std::getline(file, line) && line.find("</DataArray>") == std::string::npos) {
                        std::istringstream ss(line);
                        int id;
                        while (ss >> id) offsetsBuffer.push_back(id);
                    }
                } else if (line.find("types") != std::string::npos) {
                    while (std::getline(file, line) && line.find("</DataArray>") == std::string::npos) {
                        std::istringstream ss(line);
                        int id;
                        while (ss >> id) typesBuffer.push_back(id);
                    }
                }
            }
        }
    }

    // 座標格納
    nPoints = pointsBuffer.size() / 3;
    points.resize(nPoints);
    for (int i = 0; i < nPoints; ++i) {
        points[i].x = pointsBuffer[i*3];
        points[i].y = pointsBuffer[i*3+1];
        points[i].z = pointsBuffer[i*3+2];
    }

    // 接続情報格納
    if (!offsetsBuffer.empty()) {
        // 従来の処理 (offsets がある場合)
        nCells = static_cast<int>(offsetsBuffer.size());
        connect.resize(nCells);
        int start = 0;
        for (int i = 0; i < nCells; ++i) {
            int end = offsetsBuffer[i];
            connect[i].assign(connectivityBuffer.begin() + start,
            connectivityBuffer.begin() + end);
            start = end;
        }

        offsets = offsetsBuffer;
    } 
    else if (!typesBuffer.empty()) {
        // (offsets がなく、types がある場合)
        nCells = static_cast<int>(typesBuffer.size());
        connect.resize(nCells);
        offsets.resize(nCells);
        
        int connectivity_index = 0; // connectivityBuffer内の現在位置
        for (int i = 0; i < nCells; ++i) {
            int numPointsPerCell = 0;

            // types[i] の値を見て、頂点数を決定する
            if (typesBuffer[i] == 10) { // VTK_TETRA (四面体)
                numPointsPerCell = 4;
            } 
            else if (typesBuffer[i] == 5) { // VTK_TRA Triangle (三角形)
                numPointsPerCell = 3;
            }
            else if (typesBuffer[i] == 8) { // VTK_QUAD (四角形)
                numPointsPerCell = 4;
            }
            else if (typesBuffer[i] == 12) { // VTK_CUBE (六面体)
                numPointsPerCell = 8;
            }
            else {
            // 未対応のセルタイプ
                throw std::runtime_error("Unsupported cell type found: " + std::to_string(typesBuffer[i]));
            }

            // connectivityBufferから頂点数だけコピー
            int start = connectivity_index;
            int end = start + numPointsPerCell;
            if(end > connectivityBuffer.size()) {
                throw std::runtime_error("Connectivity array size mismatch.");
            }

            connect[i].assign(connectivityBuffer.begin() + start,
            connectivityBuffer.begin() + end);
            
            connectivity_index = end;
            offsets[i] = connectivity_index; // offsets配列も再構築しておく
        }
    }
    else {
        // offsets も types もない場合
        nCells = 0;
        connect.clear();
            offsets.clear();
    }
    types = typesBuffer;

    std::cout << "Geometry: Loaded " << nPoints << " points, "
              << nCells << " cells." << std::endl;


}

void Geometry::surfExtract(const std::string &surfaceFile, int nsurfz_param) {
    std::ifstream ifs(surfaceFile);
    if (!ifs) {
        std::cerr << "Cannot open surface file: " << surfaceFile << std::endl;
        return;
    }

    nsurfz = nsurfz_param;

    std::vector<double> xvec, yvec;
    std::string line;

    int nsurfl_file = 0;
    ifs >> nsurfl_file;
    nsurfl = nsurfl_file;
    surflx.resize(nsurfl);
    surfly.resize(nsurfl);

    while (std::getline(ifs, line)) {
        if (line.empty() || line[0]=='#') continue;
        std::istringstream iss(line);
        double x, y;
        iss >> y >> x;
        xvec.push_back(x);
        yvec.push_back(y);
    }

    surflx = xvec;
    surfly = yvec;

    nsurfl = static_cast<int>(xvec.size());
    surfp.assign(nsurfl, std::vector<int>(nsurfz, -1));

    if (points.empty()) {std ::cout << "だめですよー" << std::endl; 
        return;}

    // z方向インデックス生成
    zmax = points[0].z;
        for (const auto& p : points)
            if (p.z > zmax) zmax = p.z;

    //std::cout<<"[DEBUG] zmax="<<zmax<<std::endl;

    double dz = zmax / double(nsurfz - 1);
    surflz.assign(nsurfz, 0.0);
    for (int j = 0; j < nsurfz; ++j) {
        surflz[j] = j * dz;
    }

    surfl.clear();
    for (size_t ci = 0; ci < types.size(); ++ci) {
        if (types[ci] == 9) {            // 9 = quad
            const auto &cell = connect[ci]; 
            for (int j = 0; j < 4; ++j) {
                surfl.push_back(cell[j]);
            }
        }
    }

    std::sort(surfl.begin(), surfl.end());
    surfl.erase(std::unique(surfl.begin(), surfl.end()), surfl.end());

    // nos を更新
    int nos = static_cast<int>(surfl.size());

    for (int i = 0; i < nsurfl; ++i) {
        for (int j = 0; j < nsurfz; ++j) {
            int foundIndex = -1;
            for (int k = 0; k < nos; ++k) {
                int pid = surfl[k]; // 0-based node index

                double sx = surflx[i];
                double sy = surfly[i]; 
                double sz = surflz[j];  

                if (isSamePointRounded(points[pid].x, points[pid].y, points[pid].z, sx, sy, sz)) {
                    foundIndex = pid;
                        break;
                    }
            }
            if (foundIndex == -1) {
                //std::cerr << "surfExtract: Warning - couldn't find point for ("
                //        << surflx[i] << ", " << surfly[i] << ", " << surflz[j] << ")\n";
            }
            surfp[i][j] = foundIndex;
        }
    }

    std::ofstream ofs("../output/surfp_output.csv");
    if (!ofs) {
        std::cerr << "Failed to open surfp_output.csv for writing.\n";
        return;
    }

    ofs << "# i j node_id x y z\n";
    for (int i = 0; i < nsurfl; ++i) {
        for (int j = 0; j < nsurfz; ++j) {
            int nid = surfp[i][j];
            if (nid < 0) continue;
            const auto& g = points[nid];
            ofs << i << ", " << j << ", "
                << nid << ", " << g.x << ", " << g.y << ", " << g.z << "\n";
        }
    }


    //int i = 1;  // 調べたい行

   for (int i = 0; i < nsurfl; ++i) {
        for (int j = 0; j < nsurfz; ++j) {
            int pid = surfp[i][j];
            if (pid >= 0)
                std::cout << "i=" << std::setw(2) << i << " j=" << std::setw(2) << j
                        << " pid=" << std::setw(4) << pid
                        << " (x,y,z)=( " << std::setw(8) <<points[pid].x << ", " 
                                        << std::setw(8) <<points[pid].y << ", " 
                                        << std::setw(8) <<points[pid].z << " )\n";
        }
    }  


}

// ------------------------
// surfArea
// ------------------------
void Geometry::surfArea() {

    /* nsurfl = static_cast<int>(surfp.size());
    nsurfz = static_cast<int>(surfp[0].size()); */

    xsup = points[0].x;
    for (const auto &p : points) {
        if (p.x > xsup) xsup = p.x;
    }

    ymid.assign(nsurfz, 0.0);

    double ymidconst = points[0].y;
    for (const auto &p : points) {
        if (p.y > ymidconst) ymidconst = p.y;
    }

    for (int j = 0; j < nsurfz; ++j) {
        double ymax = -1e9;
        for (int i = 0; i < nsurfl; ++i) {
            int idx = surfp[i][j];  // 2D→1Dインデックス変換
            if (points[idx].y > ymax) {
                ymax = points[idx].y;
            }
        }
        //ymid[j] = ymax ;  // j列における最大y値　
        ymid[j] = ymidconst;
    }



    double min_diff = 1e9;
    nxsup = 0;
    double prev_x = -1e9; // 前回のX座標を記憶するための変数
    
    for (int i = 0; i < nsurfl; ++i) {
        int nid = surfp[i][0];
        
        // メッシュが途切れている場合は終了
        if (nid < 0) {
            break; 
        }

        double curr_x = points[nid].x;

        // 【超重要】X座標が前回から進んでいない（dx ≒ 0）場合は、
        // 平坦な管（声帯の出口）に到達した証拠なので、迷わず強制終了する！
        if (i > 0 && std::abs(curr_x - prev_x) < 1e-4) {
            break;
        }

        double diff = std::abs(curr_x - xsup);

        // 差が明確に縮まっている間は更新
        if (diff < min_diff - 1e-5) {
            min_diff = diff;
            nxsup = i ; // 要素数として扱うため i + 1
        } 
        else {
            break; // 遠ざかり始めたら終了
        }

        prev_x = curr_x; 
    }


    if (nsurfl < 2 || nsurfz < 2) return; // 十分なサーフェスがない場合


    sarea.assign(nxsup, std::vector<double>(nsurfz-1, 0.0));

    for (int i = 1; i < nxsup - 1; ++i) {
        for (int j = 1; j < nsurfz - 1; ++j) {
            int pid_left  = surfp[i-1][j];
            int pid_right = surfp[i+1][j];
            int pid_down  = surfp[i][j-1];
            int pid_up    = surfp[i][j+1];

            if (pid_left < 0 || pid_right < 0 || pid_down < 0 || pid_up < 0) {
                continue;
            }

            const auto& pL = points[pid_left];
            const auto& pR = points[pid_right];
            const auto& pD = points[pid_down];
            const auto& pU = points[pid_up];

            double ti_x = 0.5 * (pR.x - pL.x);
            double ti_y = 0.5 * (pR.y - pL.y);
            double ti_z = 0.5 * (pR.z - pL.z);

            double tj_x = 0.5 * (pU.x - pD.x);
            double tj_y = 0.5 * (pU.y - pD.y);
            double tj_z = 0.5 * (pU.z - pD.z);

            double cx = ti_y * tj_z - ti_z * tj_y;
            double cy = ti_z * tj_x - ti_x * tj_z;
            double cz = ti_x * tj_y - ti_y * tj_x;

            sarea[i][j] = std::sqrt(cx * cx + cy * cy + cz * cz);
        }
    }

    std::ofstream fsA("../output/SurfArea.dat");
    
    for (int i = 1; i < nxsup-1; i++){
        for (int j = 1; j < nsurfz-1; ++j){
            fsA << sarea[i][j] << " ";
        }
        fsA << std::endl;
    }
}

void Geometry::surfExtractFromNAS(
    const std::string& nasFile,
    int nsurfl_param,
    int nsurfz_param)
{
    std::ifstream ifs(nasFile);
    if (!ifs) {
        throw std::runtime_error("Cannot open NAS file: " + nasFile);
    }

    nsurfl = nsurfl_param;
    nsurfz = nsurfz_param;

    struct RawQuad {
        int eid;
        int pid;
        int nasNode[4];
    };

    struct Quad {
        int eid;
        int pid;
        int n[4]; // VTU index
    };

    // NAS node ID -> NAS座標
    std::map<int, Point> grid;

    // NAS ID のまま保持する CQUAD4
    std::vector<RawQuad> rawQuads;

    // VTU index に変換後の CQUAD4
    std::vector<Quad> quads;

    std::string line;

    int gridLineCount = 0;
    int quadLineCount = 0;

    // ------------------------
    // 1. NASを読む
    // ------------------------
    while (std::getline(ifs, line)) {
        line = trim(line);

        if (line.empty()) continue;
        if (line[0] == '$') continue;

        if (line.rfind("GRID", 0) == 0) {
            ++gridLineCount;

            auto tokens = splitComma(line);

            if (tokens.size() >= 6) {
                try {
                    int id = std::stoi(tokens[1]);

                    double gx = std::stod(tokens[3]);
                    double gy = std::stod(tokens[4]);
                    double gz = std::stod(tokens[5]);

                    // 重要:
                    // ここは既存コードに合わせている。
                    // NAS raw (gx, gy, gz) -> internal Point(x,y,z) = (gz, gx, gy)
                    grid[id] = Point{gz, gx, gy};
                }
                catch (...) {
                    std::cerr << "[NAS WARN] Failed to parse GRID line: "
                              << line << "\n";
                }
            }
        }
        else if (line.rfind("CQUAD4", 0) == 0) {
            ++quadLineCount;

            std::replace(line.begin(), line.end(), ',', ' ');
            std::istringstream iss(line);

            std::string tag;
            int eid, pid;
            int n1, n2, n3, n4;

            iss >> tag >> eid >> pid >> n1 >> n2 >> n3 >> n4;

            if (!iss) {
                std::cerr << "[NAS WARN] Failed to parse CQUAD4 line: "
                          << line << "\n";
                continue;
            }

            RawQuad rq;
            rq.eid = eid;
            rq.pid = pid;
            rq.nasNode[0] = n1;
            rq.nasNode[1] = n2;
            rq.nasNode[2] = n3;
            rq.nasNode[3] = n4;

            rawQuads.push_back(rq);
        }
    }

    // std::cout << "[NAS DEBUG] gridLineCount = " << gridLineCount << "\n";
    // std::cout << "[NAS DEBUG] quadLineCount = " << quadLineCount << "\n";
    // std::cout << "[NAS DEBUG] grid.size() = " << grid.size() << "\n";
    // std::cout << "[NAS DEBUG] rawQuads.size() = " << rawQuads.size() << "\n";

    // ------------------------
    // 2. RawQuad の NAS node ID を VTU index に変換
    // ------------------------
    int quadOkCount = 0;
    int quadGridMissingCount = 0;
    int quadVtuMissingCount = 0;

    for (const auto& rq : rawQuads) {
        Quad q;
        q.eid = rq.eid;
        q.pid = rq.pid;

        bool ok = true;

        for (int a = 0; a < 4; ++a) {
            int nasId = rq.nasNode[a];

            auto it = grid.find(nasId);
            if (it == grid.end()) {
                ++quadGridMissingCount;
                ok = false;
                break;
            }

            const Point& gp = it->second;

            int vtuIndex = -1;

            for (int k = 0; k < static_cast<int>(points.size()); ++k) {
                if (isSamePointRounded(points[k].x, points[k].y, points[k].z,
                                       gp.x, gp.y, gp.z)) {
                    vtuIndex = k;
                    break;
                }
            }

            if (vtuIndex < 0) {
                ++quadVtuMissingCount;

                if (quadVtuMissingCount <= 10) {
                    std::cout << "[NAS WARN] VTU point not found for NAS node "
                              << nasId
                              << " mapped=("
                              << gp.x << ", "
                              << gp.y << ", "
                              << gp.z << ")\n";
                }

                ok = false;
                break;
            }

            q.n[a] = vtuIndex;
        }

        if (ok) {
            quads.push_back(q);
            ++quadOkCount;
        }
    }

    // std::cout << "[NAS DEBUG] quadOkCount = " << quadOkCount << "\n";
    // std::cout << "[NAS DEBUG] quadGridMissingCount = "
    //           << quadGridMissingCount << "\n";
    // std::cout << "[NAS DEBUG] quadVtuMissingCount = "
    //           << quadVtuMissingCount << "\n";
    // std::cout << "[NAS DEBUG] quads.size() = " << quads.size() << "\n";

    if (quads.empty()) {
        std::cerr << "[NAS ERROR] No valid CQUAD4 elements were converted. "
                  << "Check coordinate mapping grid[id] = {gz,gx,gy} "
                  << "and isSamePointRounded tolerance.\n";
        return;
    }

    // ------------------------
    // 3. 節点を一意に集める
    // ------------------------
    std::set<int> uniqueNodes;

    for (const auto& q : quads) {
        for (int a = 0; a < 4; ++a) {
            uniqueNodes.insert(q.n[a]);
        }
    }

    // std::cout << "[NAS DEBUG] uniqueNodes.size() = "
    //           << uniqueNodes.size() << "\n";

    if (uniqueNodes.empty()) {
        std::cerr << "[NAS ERROR] uniqueNodes is empty.\n";
        return;
    }

    // ------------------------
    // 4. zレベルを作る
    // ------------------------
    std::vector<double> zValues;
    zValues.reserve(uniqueNodes.size());

    for (int pid : uniqueNodes) {
        zValues.push_back(points[pid].z);
    }

    
const double tolZCluster = 1.0e-3;  // zレベルを作るための許容値
const double tolZAssign  = 1.0e-2;  // 既存節点をzレベルに割り当てる許容値

    std::vector<double> zLevels = makeLevels(zValues, tolZCluster);

    // std::cout << "[NAS DEBUG] zValues.size() = "
    //           << zValues.size() << "\n";
    // std::cout << "[NAS DEBUG] zLevels.size() = "
    //           << zLevels.size() << "\n";

    if (!zLevels.empty()) {
        std::cout << "[NAS DEBUG] zLevels front/back = "
                  << zLevels.front()
                  << " / "
                  << zLevels.back()
                  << "\n";
    }

    if (static_cast<int>(zLevels.size()) != nsurfz_param) {
        std::cerr << "[WARN] detected z levels = "
                  << zLevels.size()
                  << ", expected nsurfz = "
                  << nsurfz_param
                  << "\n";
    }

    nsurfz = static_cast<int>(zLevels.size());

    if (nsurfz == 0) {
        std::cerr << "[NAS ERROR] nsurfz became 0.\n";
        return;
    }

    auto getZLevel = [&](int pid) {
        return nearestLevel(zLevels, points[pid].z, tolZAssign);
    };

    // ------------------------
    // 5. zレベルごとに i方向接続を作る
    //
    // CQUAD4の節点順に依存せず、
    // 各quad内で同じzレベルに属する2点を接続する。
    // ------------------------
    std::vector<std::unordered_map<int, std::set<int>>> adjByZ(nsurfz);

    int sameZPairCount = 0;
    int badQuadZGroupCount = 0;
    int invalidZNodeCount = 0;

    for (const auto& q : quads) {
        int ids[4] = {q.n[0], q.n[1], q.n[2], q.n[3]};

        std::map<int, std::vector<int>> groups;

        for (int a = 0; a < 4; ++a) {
            int pid = ids[a];
            int jz = getZLevel(pid);

            if (jz < 0) {
                ++invalidZNodeCount;
                continue;
            }

            groups[jz].push_back(pid);
        }

        for (const auto& kv : groups) {
            int jz = kv.first;
            const auto& v = kv.second;

            if (v.size() == 2) {
                int a = v[0];
                int b = v[1];

                adjByZ[jz][a].insert(b);
                adjByZ[jz][b].insert(a);

                ++sameZPairCount;
            }
            else {
                ++badQuadZGroupCount;

                if (badQuadZGroupCount <= 20) {
                    std::cout << "[NAS WARN] quad eid=" << q.eid
                              << " pid=" << q.pid
                              << " zLevel=" << jz
                              << " has " << v.size()
                              << " nodes\n";
                }
            }
        }
    }

    // std::cout << "[NAS DEBUG] sameZPairCount = "
    //           << sameZPairCount << "\n";
    // std::cout << "[NAS DEBUG] badQuadZGroupCount = "
    //           << badQuadZGroupCount << "\n";
    // std::cout << "[NAS DEBUG] invalidZNodeCount = "
    //           << invalidZNodeCount << "\n";

    if (sameZPairCount == 0) {
        std::cerr << "[NAS ERROR] No same-z pairs were found. "
                  << "Cannot build i-direction polylines.\n";
        return;
    }


    std::vector<std::vector<int>> nodesByZ(nsurfz);

    int assignedZNodeCount = 0;
    int unassignedZNodeCount = 0;

    for (int pid : uniqueNodes) {
        int jz = getZLevel(pid);

        if (jz >= 0) {
            nodesByZ[jz].push_back(pid);
            ++assignedZNodeCount;
        } else {
            ++unassignedZNodeCount;

            if (unassignedZNodeCount <= 20) {
                std::cout << "[Z ASSIGN WARN] pid=" << pid
                        << " z=" << points[pid].z
                        << "\n";
            }
        }
    }

// std::cout << "[NAS DEBUG] assignedZNodeCount = "
//           << assignedZNodeCount << "\n";

// std::cout << "[NAS DEBUG] unassignedZNodeCount = "
//           << unassignedZNodeCount << "\n";



    // ------------------------
    // 6. 各zレベルでポリラインをたどる
    // ------------------------
    std::vector<std::vector<int>> rowsByZ(nsurfz);

    for (int j = 0; j < nsurfz; ++j) {
        std::cout << "[NAS DEBUG] adjByZ[" << j << "].size() = "
                  << adjByZ[j].size() << "\n";

        rowsByZ[j] = traceAllComponentsAsOneRow(adjByZ[j], nodesByZ[j], points);

        std::cout << "z-level j=" << j
                  << " nodes=" << rowsByZ[j].size()
                  << "\n";
    }

    // ------------------------
    // 7. zレベル間で向きをそろえる
    // ------------------------
    auto dist2 = [&](int p, int q) {
        double dx = points[p].x - points[q].x;
        double dy = points[p].y - points[q].y;
        double dz = points[p].z - points[q].z;
        return dx * dx + dy * dy + dz * dz;
    };

    for (int j = 1; j < nsurfz; ++j) {
        if (rowsByZ[j].empty() || rowsByZ[j - 1].empty()) continue;

        int a0 = rowsByZ[j][0];
        int a1 = rowsByZ[j].back();

        int b0 = rowsByZ[j - 1][0];
        int b1 = rowsByZ[j - 1].back();

        double same = dist2(a0, b0) + dist2(a1, b1);
        double flip = dist2(a0, b1) + dist2(a1, b0);

        if (flip < same) {
            std::reverse(rowsByZ[j].begin(), rowsByZ[j].end());
        }
    }

    // ------------------------
    // 8. surfp[i][j] を構築
    // ------------------------
    int detectedNsurfl = 0;

    for (const auto& row : rowsByZ) {
        detectedNsurfl =
            std::max(detectedNsurfl, static_cast<int>(row.size()));
    }

    if (detectedNsurfl != nsurfl_param) {
        std::cerr << "[WARN] detected nsurfl = "
                  << detectedNsurfl
                  << ", expected nsurfl = "
                  << nsurfl_param
                  << "\n";
    }

    nsurfl = detectedNsurfl;

    if (nsurfl == 0) {
        std::cerr << "[NAS ERROR] nsurfl became 0.\n";
        return;
    }

    surfp.assign(nsurfl, std::vector<int>(nsurfz, -1));

    for (int j = 0; j < nsurfz; ++j) {
        for (int i = 0; i < static_cast<int>(rowsByZ[j].size()); ++i) {
            surfp[i][j] = rowsByZ[j][i];
        }
    }

    // ------------------------
    // 9. surfp検査
    // ------------------------
    int missing = 0;
    int badZOrder = 0;
    int largeDz = 0;
    int largeDi = 0;

    const double dzLimit = 1.0;
    const double diLimit = 5.0;

    auto dist3 = [&](int p, int q) {
        double dx = points[p].x - points[q].x;
        double dy = points[p].y - points[q].y;
        double dz = points[p].z - points[q].z;
        return std::sqrt(dx * dx + dy * dy + dz * dz);
    };

    for (int i = 0; i < nsurfl; ++i) {
        for (int j = 0; j < nsurfz; ++j) {
            if (surfp[i][j] < 0) {
                ++missing;
            }
        }
    }

    // j方向 = z方向
    for (int i = 0; i < nsurfl; ++i) {
        for (int j = 0; j < nsurfz - 1; ++j) {
            int p0 = surfp[i][j];
            int p1 = surfp[i][j + 1];

            if (p0 < 0 || p1 < 0) continue;

            double z0 = points[p0].z;
            double z1 = points[p1].z;
            double dz = std::abs(z1 - z0);

            if (z1 < z0) {
                ++badZOrder;

                std::cout << "[BAD Z ORDER] i=" << i
                          << " j=" << j
                          << " z0=" << z0
                          << " z1=" << z1
                          << "\n";
            }

            if (dz > dzLimit) {
                ++largeDz;

                std::cout << "[LARGE DZ] i=" << i
                          << " j=" << j
                          << " dz=" << dz
                          << " z0=" << z0
                          << " z1=" << z1
                          << "\n";
            }
        }
    }

    // i方向 = 表面に沿う方向
    for (int i = 0; i < nsurfl - 1; ++i) {
        for (int j = 0; j < nsurfz; ++j) {
            int p0 = surfp[i][j];
            int p1 = surfp[i + 1][j];

            if (p0 < 0 || p1 < 0) continue;

            double d = dist3(p0, p1);

            if (d > diLimit) {
                ++largeDi;

                std::cout << "[LARGE I STEP] i=" << i
                          << " j=" << j
                          << " dist=" << d
                          << "\n";
            }
        }
    }

    std::cout << "[surfp check] missing=" << missing
              << " badZOrder=" << badZOrder
              << " largeDz=" << largeDz
              << " largeDi=" << largeDi
              << "\n";

    // ------------------------
    // 10. zmax更新
    // ------------------------
    zmax = -1.0e100;

    for (int i = 0; i < nsurfl; ++i) {
        for (int j = 0; j < nsurfz; ++j) {
            int pid = surfp[i][j];
            if (pid < 0) continue;

            zmax = std::max(zmax, points[pid].z);
        }
    }

    if (zmax < -1.0e90) {
        std::cerr << "[WARN] zmax could not be computed from surfp.\n";
    }

    // ------------------------
    // 11. CSV出力
    // ------------------------
    std::ofstream ofs("../output/surfp_output.csv");
    if (!ofs) {
        std::cerr << "Failed to open surfp_output.csv for writing.\n";
        return;
    }

    ofs << "# i j node_id x y z\n";

    for (int i = 0; i < nsurfl; ++i) {
        for (int j = 0; j < nsurfz; ++j) {
            int nid = surfp[i][j];
            if (nid < 0) continue;

            const auto& g = points[nid];

            ofs << i << ", " << j << ", "
                << nid << ", "
                << g.x << ", "
                << g.y << ", "
                << g.z << "\n";
        }
    }

    std::cout << "Surface extracted by NAS connectivity: "
              << "(" << nsurfl << " x " << nsurfz << ") grid.\n";
}
//     }   */


// ------------------------
// デバッグ用
// ------------------------
void Geometry::print() const {
   std::cout<< "nxsup=" <<nxsup<<std::endl;
    std::cout<< "xsup=" <<xsup<<std::endl;
    std::cout<< "ymid=" <<ymid[0]<<std::endl;
    std::cout<< "zmax=" <<zmax<<std::endl;
    std::cout<< "nsurfl=" <<nsurfl<<std::endl;
    std::cout<< "nsurfz=" <<nsurfz<<std::endl;

/*     for (int i = 1; i < nxsup-1; ++i) {            // Fortran: 2..nxsup-1
        for (int j = 1; j < nsurfz-1; ++j) {  // Fortran: 2..nsurfz-1
            std::cout<<sarea[i][j]<<" ";
        }  
        std::cout<<std::endl; 
    }  */  
}
