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

inline bool isSamePointRounded(double x1, double y1, double z1,
                               double x2, double y2, double z2,
                               double scale = 1e2) {
    return (std::llround(x1*scale) == std::llround(x2*scale)) &&
           (std::llround(y1*scale) == std::llround(y2*scale)) &&
           (std::llround(z1*scale) == std::llround(z2*scale));
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
        // ★ 新しい処理 (offsets がなく、types がある場合) ★
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
                std::cerr << "surfExtract: Warning - couldn't find point for ("
                        << surflx[i] << ", " << surfly[i] << ", " << surflz[j] << ")\n";
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



    double min_diff = 1e3;
    nxsup = 0;
    for (int i = 0; i < nsurfl; ++i) {
        double diff = std::abs(points[surfp[i][0]].x - xsup);
        if (diff < min_diff) {
            min_diff = diff;
            nxsup = i+1;
        }
    }



    if (nsurfl < 2 || nsurfz < 2) return; // 十分なサーフェスがない場合


    sarea.assign(nxsup, std::vector<double>(nsurfz-1, 0.0));

    for (int i = 1; i < nxsup-1; ++i) {            // Fortran: 2..nxsup-1
        for (int j = 1; j < nsurfz-1; ++j) {  // Fortran: 2..nsurfz-1
            int pid_left  = surfp[i-1][j];
            int pid_right = surfp[i+1][j];
            int pid_down  = surfp[i][j-1];
            int pid_up    = surfp[i][j+1];

            if (pid_left < 0 || pid_right < 0 || pid_down < 0 || pid_up < 0) continue;

            double dx = 0.5 * (points[pid_right].x - points[pid_left].x);
            double dy = 0.5 * (points[pid_right].y - points[pid_left].y);
            double ds = std::sqrt(dx*dx + dy*dy);
            double dz = 0.5 * (points[pid_up].z - points[pid_down].z);
            sarea[i][j] = ds * dz;
        }   
    }     
}

void Geometry::surfExtractFromNAS(const std::string& nasFile, int nsurfl_param, int nsurfz_param) {
    std::ifstream ifs(nasFile);
    if (!ifs) throw std::runtime_error("Cannot open NAS file: " + nasFile);

    nsurfl = nsurfl_param;
    nsurfz = nsurfz_param;

    std::map<int, Point> grid;
    std::vector<std::vector<int>> quads;

    std::string line;
    while (std::getline(ifs, line)) {
        if (line.rfind("GRID", 0) == 0) {
            std::vector<std::string> tokens;
            std::stringstream ss(line);
            std::string token;

            while (std::getline(ss, token, ',')) {
                tokens.push_back(token);
            }

            // "GRID,1,,0.00000,0.00000,0.00000"
            // → tokens = ["GRID", "1", "", "0.00000", "0.00000", "0.00000"]
            if (tokens.size() >= 6) {
                try {
                    int id = std::stoi(tokens[1]);
                    double x = std::stod(tokens[3]);
                    double y = std::stod(tokens[4]);
                    double z = std::stod(tokens[5]);
                    grid[id] = {z, x, y};
                } catch (const std::exception& e) {
                    std::cerr << "Failed to parse GRID line: " << line << "\n";
                }
            }
        }

        else if (line.rfind("CQUAD4", 0) == 0) {
            std::replace(line.begin(), line.end(), ',', ' ');
            std::istringstream iss(line);
            std::string tag;
            int eid, pid;
            int n1, n2, n3, n4;
            iss >> tag >> eid >> pid >> n1 >> n2 >> n3 >> n4;
            quads.push_back({n1, n2, n3, n4});
        }
    }

    // すべての節点IDを一意に集める
    std::set<int> uniqueNodes;
    for (auto& q : quads)
        uniqueNodes.insert(q.begin(), q.end());
    
    std::cout<<"uniquenodes="<<uniqueNodes.size()<<"\n";

    // 座標付きで並べ替え用にベクトル化
    struct NodeRef { double x, y, z; };
    std::vector<NodeRef> nodes;
    for (int nid : uniqueNodes) {
        const auto& p = grid[nid];
        nodes.push_back({p.x, p.y, p.z});
    }

    // x → z の順でソート（i: x方向, j: z方向）
    const double eps = 1e-3;

    std::sort(nodes.begin(), nodes.end(), [&](const NodeRef &a, const NodeRef &b){

        if (std::abs(a.x - b.x) > eps) 
            return a.x < b.x; // 主キー: x (i方向)
        else
            return a.z < b.z;                               // 副キー: z (j方向)

    });

    // 2D格子化
    if (nodes.size() != nsurfl * nsurfz) {
        std::cerr << "Warning: node count (" << nodes.size() 
                  << ") != nsurfl*nsurfz (" << nsurfl * nsurfz << ")\n";
    }

    surfp.assign(nsurfl, std::vector<int>(nsurfz, -1));

    int idx = 0;
    for (int i = 0; i < nsurfl; ++i) {
        for (int j = 0; j < nsurfz; ++j) {
            if (idx >= (int)nodes.size()) break;

            const auto& n = nodes[idx];  // ソート済みの座標
            int vtuIndex = -1;

            // VTU上のインデックスを探索
            for (int k = 0; k < (int)points.size(); ++k) {
                if (isSamePointRounded(points[k].x, points[k].y, points[k].z,
                                    n.x, n.y, n.z)) {
                    vtuIndex = k;
                    break;
                }
            }

            if (vtuIndex == -1) {
                std::cerr << "Warning: couldn't find VTU point for (" 
                        << n.x << ", " << n.y << ", " << n.z << ") \n";
            }

            surfp[i][j] = vtuIndex;
            ++idx;
        }
    }

    zmax = points[0].z;
    for (const auto& p : points)
        if (p.z > zmax) zmax = p.z;



    std::cout << "Surface extracted: " << nodes.size() 
              << " nodes → (" << nsurfl << " x " << nsurfz << ") grid.\n";

    // surfpの内容を出力して確認
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

/*    for (int i = 0; i < nsurfl; ++i) {
        for (int j = 0; j < nsurfz; ++j) {
            int pid = surfp[i][j];
            if (pid >= 0)
                std::cout << "i=" << std::setw(2) << i << " j=" << std::setw(2) << j
                        << " pid=" << std::setw(4) << pid
                        << " (x,y,z)=( " << std::setw(8) <<points[pid].x << ", " 
                                        << std::setw(8) <<points[pid].y << ", " 
                                        << std::setw(8) <<points[pid].z << " )\n";
        }
    }   */
}

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
