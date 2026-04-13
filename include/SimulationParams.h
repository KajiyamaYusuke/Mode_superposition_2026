
#pragma once

#if __has_include(<filesystem>)
  #include <filesystem>
  namespace fs = std::filesystem;
#elif __has_include(<experimental/filesystem>)
  #include <experimental/filesystem>
  namespace fs = std::experimental::filesystem;
#else
  #error "No filesystem support"
#endif

#include <string>
#include <iostream>

/// Simulation parameters (読み取り専用の設定群)
struct SimulationParams {
    // デフォルトコンストラクタ（メンバ初期化子で既に初期値を定義）
    SimulationParams() = default;

    // --- 数値パラメータ ---
    int    nmode   = 20;       // モード数
    int    nsurfz  = 4;        // spanwise 分割数
    int    ncont   = 0;
    int    nstep   = 10000;    // 総ステップ数
    int    nwrite  = 100;      // 出力間隔（ステップ）
    double dt      = 1e-5;     // 時間刻み [s]
    double zeta    = 0.0;      // 減衰比

    // contact stiffness (例)。意味はプロジェクトに合わせて調整してください
    double kc1     = 1e6;
    double kc2     = 1e6;
    double kc3     = 1e6;

    double mass    = 1.0;      // 単位はコード内で統一（例: kg）

    // --- 物理パラメータ（単位をコメント） ---
    int    iforce  = 1;        // 力の種類フラグ（Fortran互換）
    double forcef  = 0.0;
    double famp    = 0.0;
    double ps      = 101325.0; // 静圧 [Pa]
    double rho     = 1.225;    // 密度 [kg/m^3] (空気の初期値)
    double mu      = 1.81e-5;  // 動粘性係数 [Pa·s]（参考値）
    double c_sound = 340.0;

    // --- 音響管（Vocal Tract / Subglottal）パラメータ ---
    double L_inlet = 0.50;  // 吸気管の長さ [m]
    double r_inlet = 0.15;  // 吸気管の半径 [m]

    double L_sub   = 0.25;  // 声門下管の長さ [m]
    double r_sub   = 0.0125;// 声門下管の半径 [m] (2.5cm / 2)
    int    N_sub   = 3;     // 声門下管のセクション数 (Nsecgに対応)

    double L_vt    = 0.25;   // 声道の長さ [m]
    double A_vt    = 0.0125;  // 声道の断面積 [m^2]
    int    N_vt    = 10;    // 声道のセクション数 (Nsecpに対応)

    // --- ファイル/ディレクトリパス ---
    fs::path inputDir  = ".";
    fs::path resultDir = "./results";
    fs::path freqFile  = "no_mem_freq.txt";
    fs::path modeFile  = "no_mem_mode.vtk";
    fs::path surfFile  = "surface.txt";

    // --- IO / 検証 ---
    // filename を読み込み、エラー文字列は err に格納して false を返す
    bool loadFromFile(const fs::path& filename, std::string& err);

    // パラメータの整合性チェック。問題があれば err に入れて false を返す
    bool validate(std::string& err) const;

    // デバッグ出力（標準出力か、指定した ostream に出す）
    void print(std::ostream& os = std::cout) const;
};


