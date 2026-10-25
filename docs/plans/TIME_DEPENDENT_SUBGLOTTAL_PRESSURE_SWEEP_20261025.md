# Codex実装仕様：声門下駆動圧の線形増圧・減圧・往復掃引と圧力依存係数の毎ステップ更新

作成：2026-10-25  
対象：`KajiyamaYusuke/Mode_superposition_2026` の `main`（確認時 SHA: `3fa95c6bcdaf8587885ae3710c8c8143d4c351da`）。**実装着手時にはorigin/mainを再取得し、現時点の実コードを優先する。**  
目的：M5等の3D声帯モデルで、連成計算を1回走らせながら入口の声門下「指定駆動圧」を時間に対し増圧／減圧／往復させる。経験式の係数も同じ指定圧力を使って各ステップ更新する。  
範囲：**実装と最小限のソフトウェア動作確認のみ**。物理妥当性の再検証、発振閾圧の自動判定、パラメータ探索、フル掃引計算の自動実行は今回は不要。

## 0. 結論：計算順序

ユーザーの「圧力更新 → 経験係数計算 → 抵抗等の計算」は**基本的に正しい**。ただし抵抗は2種類あるので区別する。

1. **時刻 `t_n=n*dt` の変形後声帯形状（既存の `q_n`）**から、50断面の位置、`A(x,t_n)`、`G3(x,t_n)`、`Amin`、剥離点を更新（既存の `Simulation::run` → `updateChannelSections()`）。これは圧力スケジュールとは独立で、既存の先行順序を変えなくてよい。
2. **指定駆動圧 `Pcmd(t_n)` を確定**（線形増圧／減圧／往復。旧プロファイルなら従来のcosineランプ）。`lungPressure_` はこの値である。**声門直下の音響変動圧 `Pg(t)` とは異なる**。
3. 新規掃引プロファイルでは、**同じ `Pcmd(t_n)` から** `alpha1,alpha2,beta` を更新し、`R2, Rr` を再計算。**`legacy_cosine` では従来どおり係数は `params.ps` で初期化した固定値を維持**し、既存実行を変えない。
4. **形状由来**の `Lg(t_n),Rv(t_n)` をA/G3と剥離位置から計算（modern）または既存のlegacy係数関数で計算。**`beta(t_n)` を使うため legacy計算より前にステップ3を行うこと。**
5. `Rk(t_n)=beta(t_n)*rho/Amin(t_n)^2` を算出し、声門下・声道の音響状態とともに声門流量 `Ug(t_n)` を既存NR式で更新。`R2(t_n)` は同じステップの上流管更新に、`Rr(t_n)` は同じステップの放射境界に使用する。
6. 更新した `Ug, Pg, Pd` から壁面圧力・節点流体力を計算、一般化力へ射影。既存の接触反復とNewmarkで `q_{n+1}` を確定。**構造・接触・Newmark・メッシュ幾何・剥離則・圧力分布式は変更しない。**
7. `t_n, Pcmd, alpha1,alpha2,beta,R2,Rr` の小容量CSVを指定間隔で保存。

因果関係：**形状 → A/G3 → 幾何係数、指定圧力 → 経験係数 → 流体係数、両者＋音響履歴 → Ug/Pg/Pd → 壁面力 → 次の形状**。すべて同じ `t_n` のパラメータを使用する。圧力スケジュールだけを `calcFlowStep()` で更新する設計は禁止（legacy `FlowModel::advance()` の `LegacyGlottal::coefficients(... beta_ ...)` が先に実行されるため）。

## 1. 現行コードの確認済み事実と変更理由

- `src/fluid/flowModel.cpp::initialize` にて、`alpha1_=-2.5e-5*params.ps+0.185`、`alpha2_=1.6e-3*params.ps+0.6`、`beta_=1.125e-4*params.ps+0.1375`、`R2_`、声道有効時の `Rr_` を**一度だけ**計算している。
- `calcFlowStep` は、指定がなければ0.05 s（明示`pressureRampTimeSec` があればその長さ）のcosine圧力ランプで `lungPressure_=sp_->ps*rampRatio_` を設定している。
- `FlowModel::advance` のlegacy分岐では `calcFlowStep` に入る**前**に `LegacyGlottal::coefficients(rho_,mu_,beta_,...)` が呼ばれる。更新順序を誤ると旧流体係数だけ1ステップ遅れ／固定になる。
- `calcFlowStep` は `R2_` を上流気道、`beta_` を `Rk1`、`Rr_` を放射端に使用。 `Lu,Cu,Lui,Cui,La,Ca,Lr` は音響管の幾何等によるもので、圧力掃引で更新しない。
- `src/forces/legacySurfacePressureLoad.cpp` の旧圧力経路は `flowModel.lungPressure()` を参照する。これも当該時刻の `Pcmd` であること。
- `src/simulation/simulation.cpp` の `SIMULATION_RUN_DIR` と `writeDeformationVtu=0` は既存機能。既存スカラー出力やカレントの `params_used.txt` を壊さない。

## 2. 圧力プロファイルの外部仕様（入力paramの末尾に名前付き指定）

### 2.1 対応プロファイル

`pressureProfile = legacy_cosine | linear_sweep | linear_cycle`

- **legacy_cosine（デフォルト）**：既存と完全同一の `Ps=params.ps`・従来cosine立上げ。経験係数は初期化時の `params.ps` に固定。既存`pressureRampTimeSec` は従来どおり有効。未指定パラメータとの互換性が最優先。
- **linear_sweep**：`P0` を `tStart` まで維持、`T_rise` 秒で `Ptarget` へ線形遷移し、その後は `Ptarget` で保持。`Ptarget<P0` を指定すれば線形減圧。変化時間は名称が `Rise` でも増減共通。
- **linear_cycle**：`P0` を開始時刻まで保持、`Ppeak` へ `T_rise` 秒で増圧、`T_hold` 秒保持、`Pfinal` へ `T_fall` 秒で減圧し、その後は `Pfinal` 保持。 `Ppeak>=max(P0,Pfinal)` を必須とする。

名前付きキー（大文字小文字と既存key正規化規則に従う）：

| キー | 意味 | 新プロファイルでの扱い |
| --- | --- | --- |
| `pressureProfile` | 上記モード | 省略時`legacy_cosine` |
| `pressureInitialPa` | `P0` [Pa] | 必須 |
| `pressureTargetPa` | `Ptarget/Ppeak` [Pa] | 必須 |
| `pressureFinalPa` | `Pfinal` [Pa] | `linear_cycle`のみ必須 |
| `pressureStartSec` | 線形遷移開始 `tStart` [s] | 省略時0 |
| `pressureRiseTimeSec` | `T_rise` [s] | 必須、正 |
| `pressureHoldTimeSec` | `T_hold` [s] | cycleのみ。省略時0 |
| `pressureFallTimeSec` | `T_fall` [s] | cycleのみ必須、正 |

**既存の位置引数 `ps` は削除・変更しない。** 新プロファイルでは上記の明示的な `pressureInitialPa/TargetPa` が指定圧力を決め、`ps` は旧形式との互換入力として残す。入力不備はエラーにし、`ps` を目標圧へ**暗黙流用しない**。 `pressureRampTimeSec` は`legacy_cosine`専用。線形プロファイルと同時指定した場合は警告して無視するのではなく明確に拒否する（誤設定防止）。

### 2.2 厳密な時間関数（端点・範囲）

`r(t;a,b)=clamp((t-a)/(b-a),0,1)`（有効な `b>a` を検証）。

`linear_sweep`：
- `t<tStart`：`Pcmd=P0`
- `tStart<=t<tStart+T_rise`：`Pcmd=P0+(Ptarget-P0)*r(t;tStart,tStart+T_rise)`
- それ以降：`Pcmd=Ptarget`

`linear_cycle`：
- `t<tStart`：`Pcmd=P0`
- `tStart<=t<tStart+T_rise`：`P0→Ppeak` の線形補間
- `tStart+T_rise<=t<tStart+T_rise+T_hold`：`Ppeak`
- `tStart+T_rise+T_hold<=t<tStart+T_rise+T_hold+T_fall`：`Ppeak→Pfinal` の線形補間
- 以降：`Pfinal`

時刻 `t=n*dt` を使い、`n+1`の先取り・余分な圧力ランプの二重適用はしない。`0`秒付近および各転換点で有限かつ連続の圧力にする。`P0=Ptarget` は一定圧力として許容してよいが、時間長は正とする。

入力例（**線形増圧**）：

~~~ini
pressureProfile = linear_sweep
pressureInitialPa = 0
pressureTargetPa = 2000
pressureStartSec = 0
pressureRiseTimeSec = 1.0
writeDeformationVtu = 0
~~~

入力例（**線形減圧**）：

~~~ini
pressureProfile = linear_sweep
pressureInitialPa = 2000
pressureTargetPa = 0
pressureStartSec = 0
pressureRiseTimeSec = 1.0
writeDeformationVtu = 0
~~~

入力例（**増圧→保持→減圧**）：

~~~ini
pressureProfile = linear_cycle
pressureInitialPa = 0
pressureTargetPa = 2000
pressureFinalPa = 0
pressureStartSec = 0
pressureRiseTimeSec = 1.0
pressureHoldTimeSec = 0.2
pressureFallTimeSec = 1.0
writeDeformationVtu = 0
~~~

上記cycleの全時間を観察するには `dt=1e-5 s` なら **`nstep>=220000`** が必要（最終端点 `t=2.2 s` 自体の保存が必要ならさらに1ステップ）。計算時間や容量負荷のため、Codexは本番フル計算を自動実行しない。

## 3. 毎ステップ更新する経験式と係数

`P=Pcmd(t_n)` [Pa]（音響変動 `Pg` ではない）から

~~~text
alpha1 = 0.185   - 2.5e-5   * P
alpha2 = 0.600   + 1.6e-3   * P
beta   = 0.1375  + 1.125e-4 * P
Asub   = pi * r_sub^2
Avt    = pi * r_vt^2
R2     = alpha1/(Asub^2) * sqrt(rho*mu*c_sound)
Rr     = alpha2*rho*c_sound/(9*pi^2*Avt)       [hasVocalTract の場合のみ]
Rk     = beta*rho/(Amin^2)                     [声門が開放している場合]
~~~

- `alpha1, alpha2, beta` は各ステップ `Pcmd` から**全部再計算**。 `R2, Rr` は更新済みalphaを使用し、`Rr=0` は声道なしの場合維持する。 `Rk` は`calcFlowStep`で最新betaと最新Aminから計算される既存式を維持。
- modernの `Lg,Rv` は最新 `A,G3,xSep` から従来通り再積分。legacyの `Lg,Rv` は既存の面積ベース則のままで変更しない。legacy `LegacyGlottal::coefficients()` の`Rk`に渡すbetaは**必ず今回の更新値**。
- `rho,mu,c_sound, r_sub,r_vt` 等の定数と `Lu,Cu,Lui,Cui,La,Ca,Lr` は初期化時のまま。気道の音響状態 `Pu/Uu/Pd/Ud` は既存更新則を継続し、掃引転換点でリセットしない。
- 現行経験式の出典は「3D流体計算との比較から決めた経験的な圧力補正式」と説明し、今回は各時刻の圧力へ適用する**準定常近似**と記録する。元の有効圧力範囲は未確認であるため「全圧力で妥当」と記さない。係数を勝手に一定値へクランプしたり、圧力で`Pg`を代用したりしない。
- 非有限・負の指定圧力は入力エラー。掃引端点で `alpha1 < 0`、`alpha2 <= 0`、`beta <= 0` となる条件は不適切な係数として入力時に拒否（現行式なら`alpha1`が高圧で負になり得る）。**経験式の適用範囲の保証ではなく符号の最小限の検査**と注記する。`alpha1=0`はモデル上許容可能だがゼロ抵抗となることをログに明示してよい。

## 4. 実装設計：変更位置・依存順序の保証

### 4.1 `SimulationParams.h`／`readParams.cpp`

- `PressureProfile`列挙もしくは既存命名規則に合う文字列、関連数値、明示入力フラグを追加。既存の位置パラメータと読み込み順を変えない。
- 名前付きパラメータに新キーを追加。重複、余計なフィールド、無効なプロファイル、期間`<=0`、負/NaN/Inf、順序の矛盾、`linear_cycle`の非ピーク設定を拒否。
- 未指定時の`legacy_cosine`の**従来の`pressureRampTimeSec`解釈（未指定0.05秒、指定時のみその値）を厳密に維持**。
- `iforce=1`と時間依存圧力プロファイルの同時指定は拒否または明示的に無効扱い（推奨：入力エラー）。独立加振で掃引値を使用したと誤認させない。

### 4.2 `FlowModel.h`／`flowModel.cpp`

- 小さな`pressureCommandAt(time)`と`preparePressureStep(time)`（仮名）を用意し、`preparePressureStep`は**`FlowModel::advance()`の先頭で分岐前に1回だけ**呼ぶ。
- legacy_cosineは従来と同じ`rampRatio_`と`lungPressure_`の計算を移設する。係数の再計算をしない。旧固定圧力結果を保持。
- linear_sweep/linear_cycleは各時刻の`Pcmd`を計算して`lungPressure_`に代入し、その値でalpha群、R2、Rrを更新。
- `calcFlowStep()`内部の現行cosine再計算を除き、**先頭で確定した`lungPressure_`を上流気道更新で使用**する。legacy/new/closedのどの分岐でも圧力と係数を一貫させる。
- `R2,Rr`算出用の`Asub,Avt,rho*mu*c_sound`は初期化由来の値を再利用する（計算式を再定義しない）。声道なしでは`Rr=0`。
- `rampRatio_`は旧プロファイルの互換性のため維持。新プロファイルでは、曖昧な全体圧力比を発明せず、各活動線形区間の進行度（0→1）または0/1の定義を文書化。**圧力の正本は`lungPressure_`**とする。
- `FlowModel::acousticCoefficients()`は最新値を返す。初期化直後の値も意味が分かるようにする（旧は`params.ps`値、新は`Pcmd(0)`値）。

### 4.3 `Simulation::run()`／出力

- 各ステップの `FlowModel::advance`実行後、`pressure_schedule.csv`に追跡可能な最小限の列を記録：
  `step,time_s,pressure_command_pa,alpha1,alpha2,beta,R2,Rr`。単位をヘッダやREADMEに明記。**nwriteまたは新しい軽量な間引き設定**を使い、200kステップでもファイルサイズを抑える。時刻0・終了近傍の値は失わないこと。
- `runDir/txt/manifest.txt`、`params_used.txt`に圧力プロファイル、全パラメータ、係数を動的更新する方式、実際の設定開始/終了時刻、圧力端点を記録。
- 既存の`acoustic_coefficients.csv`は初期値のスナップショットであり、時変計算全体の一定値と誤読させないよう `initial` 等と明記。旧CSVリーダとの互換性が必要ならファイル構造は維持し、manifestで説明を補う。
- 既存の`acoustic_trace.csv`と`lungPressure()`が新方式でも正しい値を返すこと。追加の全節点診断は不要。`writeDeformationVtu=0`を実行例で案内。
- 本番長時間計算、M5の素材設定や声門下管の分割数、音響モデル、接触設定、Newmark、VTU書式、既存データパスを勝手に変更しない。

## 5. 最低限の実装チェック（物理検証ではない）

複数条件の本計算や発振閾圧の解析は**しない**。最小限：
1. `icpx` Releaseビルド＋既存CTestを実行。コンパイル不具合がなければ十分。
2. `legacy_cosine`を未指定／明示指定した短い計算で、変更前と`Pcmd,R2,Rr,beta,Ug,Pg`が一致するか確認（可能なら既存出力のbitwise一致）。**既存の全流体係数固定**が守られること。
3. プロファイル関数の数値単体試験：増圧 `0→2000 Pa/1s` では`t=0,0.5,1s`で`0,1000,2000 Pa`、減圧 `2000→0`で`2000,1000,0 Pa`、cycle `0→2000`（1 s）→保持（0.2 s）→`0`（1 s）で`t=0,0.5,1,1.1,1.2,1.7,2.2`の境界・中点を確認。**長時間の連成再計算は不要**。
4. `Pcmd=0,1000,2000 Pa`で係数式が参照値と一致すること、`Pcmd`変更と同じステップでlegacy `Rk`/modern `Rk`、`R2,Rr`へ反映されること。声道`L_vt=0`の場合の`Rr=0`と、声道ありの更新も確認。
5. 新しい入力の不正例（負圧、期間0、負保持時間、欠損、未知名、線形プロファイルと旧`pressureRampTimeSec`同時指定）を拒否すること。CSVの`step,time`と実際の駆動圧一致を確認。
6. 最大でも短時間スモークテスト。大量VTUや新規GB級バックアップ、フル掃引の自動実行は行わず、Cドライブの容量を圧迫しない。

## 6. 最終報告

- 修正したファイル、実装したスケジュール、デフォルト互換性、更新順序、各係数が変更される箇所を説明。
- 最小テストの結果（ビルド/CTest、圧力時間関数、modern/legacyの係数確認）。未実施の物理検証は「未実施」と明記。
- 実行例の`param.txt`断片と、`./simulation ../input/param.txt`の従来コマンド。
- 出力する圧力は**入口に指定した駆動圧**であり、`Pg`の実測値や音響圧力履歴ではないことを明記。
- 経験係数の圧力依存を時刻ごとに適用するのは**準定常近似**であり、3D比較で同定した経験式の適用圧力範囲は別途確認事項。今回は物理検証なし。
- 既存の未コミット変更と本番入力・過去出力を保護。変更用branch／PRで実装し、合格確認までmainへ直接force pushしない。
