
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
#include <bits/stdc++.h>
#include <iostream>

/// Simulation parameters (読み取り専用の設定群)
struct SimulationParams {
    // デフォルトコンストラクタ（メンバ初期化子で既に初期値を定義）
    SimulationParams() = default;

    // --- 数値パラメータ ---
    int    nmode   = 20;       // モード数
    // Empty: use modes 1..nmode. Otherwise use only these one-based mode numbers.
    std::vector<int> modeSelection;
    int    nsurfz  = 4;        // spanwise 分割数
    int    ncont   = 0;
    int    nstep   = 10000;    // 総ステップ数
    int    nwrite  = 100;      // 出力間隔（ステップ）
    double dt      = 1e-5;     // 時間刻み [s]
    double zetaL    = 0.0;
    double zetaR    = 0.0;      // 減衰比


    // contact stiffness (例)。意味はプロジェクトに合わせて調整してください
    double kc1     = 1e6;
    double kc2     = 1e6;
    double kc3     = 1e6;

    double mass    = 1.0;      // 単位はコード内で統一（例: kg）

    // --- 物理パラメータ（単位をコメント） ---
    int    iforce  = 0;        // 力の種類フラグ（Fortran互換）
    double forcef  = 0.0;
    double famp    = 0.0;
    // Prescribed-force direction for iforce=1:
    // 0 = common x direction, 1 = opposing y directions (gap direction).
    int    forceDirection = 0;
    // Temporary fixed reference for the legacy reduced contact coefficients.
    // It deliberately does not follow the participating structural modes.
    double contactReferenceFrequencyHz = 40.0;
    // Physical distance over which wall pressure blends after separation.
    double flowBlendLengthMm = 0.5;
    // Flow separates at the first downstream section satisfying A/Amin.
    double flowSeparationAreaRatio = 1.0;
    // Diagnostic controls.  pressureRampTimeSec is the documented value for
    // new runs; legacy files that omit it retain the historical 0.05 s ramp.
    double initialGapMm = 0.0;
    double pressureRampTimeSec = 0.10;
    bool pressureRampTimeSecExplicit = false;
    int diagnosticOutputIntervalSteps = 5;
    bool perturbationEnabled = false;
    double perturbationTimeSec = 0.20;
    int perturbationModeIndex = 5;  // user-facing, one-based source mode
    double perturbationAmplitudeAtProbeMm = 0.001;
    std::string perturbationPattern = "symmetric_opening";
    // AP-end/contact diagnostics.  AP rows are surface-grid j indices; a
    // positive distance selects rows whose reference z lies within that many
    // millimetres of either end.  Row and distance selectors are exclusive.
    bool apContactDiagnosticEnabled = false;
    int apContactAnalysisRowsPerEnd = 3;
    double apContactAnalysisDistanceMm = 0.0;
    int apContactExcludedRowsPerEnd = 0;
    double apContactExcludedDistanceMm = 0.0;
    double contactDetailStartSec = -1.0;
    double contactDetailEndSec = -1.0;
    std::string diagnosticLoadMode = "live";
    double diagnosticReferenceWindowStartSec = 0.25;
    double diagnosticReferenceWindowEndSec = 0.30;
    fs::path fixedNodeIdsFile;
    double ps      = 1325.0; // 静圧 [Pa]
    double rho     = 1.225;    // 密度 [kg/m^3] (空気の初期値)
    double mu      = 1.81e-5;  // 動粘性係数 [Pa·s]（参考値）
    double c_sound = 340.0;

    // --- 音響管（Vocal Tract / Subglottal）パラメータ ---
    double r_inlet = 8.0 * 1e-2;    // 8.0 cm
    double L_inlet = 17.5 * 1e-2;   // 17.5 cm

    double L_sub   = 0.15;  // 声門下管の長さ [m]
    double r_sub   = 0.0125;// 声門下管の半径 [m] (2.5cm / 2)
    int    N_sub   = 10;     // 声門下管のセクション数 (Nsecgに対応)

    double L_vt    = 0 * 1e-2;   // 声道の長さ [m]
    double r_vt    = 1.25 * 1e-2;   // 声道の断面積 [m^2]
    int    N_vt    = 10;    // 声道のセクション数 (Nsecpに対応)

    // --- ファイル/ディレクトリパス ---
    fs::path inputDir  = ".";
    fs::path resultDir = "./results";
    fs::path freqFile  = "no_mem_freq.txt";
    fs::path modeFile  = "no_mem_mode.vtk";
    fs::path surfFile  = "surface.txt";

    // Per-fold model inputs. Relative paths are resolved from the directory
    // containing the parameter file.
    fs::path leftFrequencyFile  = "M5_test/M5_freq_T3_d2_b12c3.txt";
    fs::path rightFrequencyFile = "M5_test/M5_freq_T3_d2_b2c2.txt";
    fs::path leftModeFile       = "M5_test/M5_mode_T3_b12c3.vtu";
    fs::path rightModeFile      = "M5_test/M5_mode_T3_b2c2.vtu";
    fs::path leftSurfaceNasFile  = "M5_test/M5_surface_T3_d2.nas";
    fs::path rightSurfaceNasFile = "M5_test/M5_surface_T3_d2.nas";

    // --- IO / 検証 ---
    // filename を読み込み、エラー文字列は err に格納して false を返す
    bool loadFromFile(const fs::path& filename, std::string& err);

    // パラメータの整合性チェック。問題があれば err に入れて false を返す
    bool validate(std::string& err) const;

    // デバッグ出力（標準出力か、指定した ostream に出す）
    void print(std::ostream& os = std::cout) const;
};
