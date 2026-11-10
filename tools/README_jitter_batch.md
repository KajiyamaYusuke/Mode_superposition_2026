# 主ピーク対応・流量下降Jitterの一括解析（既定3%）

`batch_analyze_airflow_jitter.py` は既存の `analyze_airflow_metrics.py` の
`acf_guided_extract` → `peak_associated_extract` → `all_event_metrics` を呼びます。
既定の主指標は **`peak-associated-closure/main_3pct`**。`--flow-threshold-pct`で下降閾値を変更できます。
ACF・主ピーク選択・10%明瞭さ判定・MAD床・交差の選択規則・Jitter/RAP/PPQ5の計算式は共用し、
3%の既存イベント時刻と数値を維持します。
CPPやソルバは実行しません。

## 実行

リポジトリルートで、まず入力を検証します。dry-runは検出器を実行せず、ファイルを生成しません。

```bash
python3 tools/batch_analyze_airflow_jitter.py \
  --manifest configs/jitter_batch_cases.example.csv --dry-run

python3 tools/batch_analyze_airflow_jitter.py \
  --manifest configs/jitter_batch_cases.example.csv \
  --out-dir analysis_runs/jitter_batch_20261108/batch_001
```

サンプルmanifestは既存の0.25/0.30・900 Paの回帰用2ケースです。新条件のパスや数値は含みません。
自分のCSVに行を追加して使ってください。出力先は毎回**存在しないディレクトリ**を指定します。
同名の結果は上書きしません。

パスだけでも実行できます。次の `/path/to/...` は置き換えるダミーです。

```bash
python3 tools/batch_analyze_airflow_jitter.py \
  --paths '/path/to/run A/snapshot/dat' '/path/to/run_B/snapshot/dat/airflow_vt.dat' \
  --dt 1e-5 --start 0.50 --end 1.00 \
  --out-dir analysis_runs/jitter_batch_20261108/paths_001
```

直接指定では入力順に安全なcase_idを生成します。名前から剛性比や圧力は推定せずCSVに`NA`と記録します。
`--manifest`と`--paths`は同時に指定できません。
入力ディレクトリは配下の`airflow_vt.dat`が**ちょうど1件**の場合だけ解決します。
複数runを含む親ディレクトリの自動一括列挙はしません。

別のCWDから起動する場合はスクリプトを絶対パスで指定してください。
**CSV内相対パス、CLI入力・出力の相対パスはすべてリポジトリルート基準**です。
manifestが置かれたフォルダやCWD基準ではありません。絶対パス、空白・非ASCIIを含むパスにも対応します。

## 下降閾値を変更する

このオプションはバッチCLI用です。省略すれば従来の3%、`--flow-threshold-pct 3`の明示も同じ結果です。
例えば主ピーク流量の5%への下降交差を解析するには次を使います。

```bash
python3 tools/batch_analyze_airflow_jitter.py \
  --manifest configs/jitter_batch_allRatio.csv \
  --flow-threshold-pct 5 \
  --out-dir analysis_runs/jitter_batch_20261108/allRatio_5pct_001
```

`--flow-threshold-pct 2.5`のような小数も指定できます。指定範囲は有限な **0より大きく100より小さい値**。
0%は厳密ゼロ到達の別方式（最初の到達を採用）なので、このオプションでは指定できません。
百分率はそのバッチの全ケース共通で、CSVを書き換える必要はありません。
**ピークの明瞭さの10%やMAD床は変更しません。**

非既定値では出力名・列名も実際の閾値に合わせます。

| 指定 | イベントCSV | summaryの未到達数列 | 面積との時差の基準時刻列 |
|---|---|---|---|
| 省略/3 | `events_3pct.csv` | `n_without_3pct_event` | `t_Ug_3pct_s` |
| 5 | `events_5pct.csv` | `n_without_5pct_event` | `t_Ug_5pct_s` |
| 2.5 | `events_2p5pct.csv` | `n_without_2p5pct_event` | `t_Ug_2p5pct_s` |

periods/metric_terms、`has_*_event`、`Amin_at_*_mm2`も同じ接尾辞になります。
3%の従来ファイル名・列名を維持し、追加の`flow_threshold_pct`列と`event_method`・設定IDで実際の閾値を追跡できます。
レポートと比較図にも指定値を表示し、異なる閾値の結果を同じ比較群に混ぜません。
以降の3%という記述・出力例は省略時の設定です。

## 主ピーク選択の感度検証（A/B/C）

主解析の既定は引き続きA（`--peak-selector prominence`）です。
次の指定で、同じACF領域・下降枝処理・Jitter計算のまま選択条件を切り替えられます。

| オプション値 | 方式 | 候補と同値規則 |
|---|---|---|
| `prominence` | A・既定 | 10%判定通過候補のprominence最大。従来の早い時刻・ID規則。 |
| `clear-flow` | B・主要感度 | 同じ判定通過候補の実標本流量最大。同値は早い補間時刻、次にID。 |
| `raw-flow` | C・補助感度のみ | 明瞭さ判定で落ちた実raw極大も含む流量最大。同値はBと同じ。 |

単一のバッチでB方式を試す場合は、通常のコマンドへ`--peak-selector clear-flow`を追加します。
主ピークが変わると次主ピークまでの下降枝も再探索します。閾値の数値だけを差し替える処理ではありません。
選択法はsummary/inputsの`peak_selector`と設定IDに保存し、比較群も分離します。

保存済みAとの回帰確認を含む一括比較は次です。

```bash
python3 tools/batch_analyze_airflow_jitter.py \
  --manifest configs/jitter_batch_allRatio.csv \
  --compare-peak-selectors --comparison-thresholds 3 5 \
  --baseline-root analysis_runs/jitter_batch_20261108 \
  --out-dir analysis_runs/jitter_batch_20261108/peak_selection_comparison_next
```

10%も比較する場合は`--comparison-thresholds 3 5 10`とします。
参照先は`<baseline-root>/allRatio_001`（3%）または`allRatio_<値>pct_001`です。
同じcase_id・流量SHA256・dt・窓の保存済みAを要求し、全イベントID・親ID・時刻・全間隔・
Jitter/RAP/PPQ5が完全一致しなければ停止します。
一括比較は`--flow-threshold-pct`ではなく`--comparison-thresholds`を使います。
`--peak-selector`は既定のままで、A/B/Cすべてを評価します。

比較出力は`comparison.csv`、`main_peak_changes.csv`、`events_diff.csv`、`region_effects.csv`、
`single_region_impacts.csv`、全イベント/周期/差項、raw候補と拒否理由、境界診断、レポートと最小限の図です。
`single_region_impacts.csv`は1領域だけを観測された変更先ピークへ置換し、全列を再計算した診断で、
同時変更の効果と別です。領域ごとの効果を加算して全体差とみなしたり、物理的因果の寄与と呼びません。

イベント差分は先に全列で時刻差2 ps以内を対応付け、親の移動も区別します。
残りを同じ探索領域で比較し、時刻変更・追加・欠落を記録します。
この対応処理は比較表だけに用い、計算対象のイベント時刻や個数は変更しません。
ACF境界診断も採否ゲートには使いません。Cを新しい主解析として採用する機能ではありません。

## CSV manifest

```csv
case_id,flow_path,run_path,area_path,stiffness_ratio,pressure_pa,label,dt_s,start_s,end_s,notes
r0p25,analysis_runs/manual_run_archive/r0p25_P900_20261104T064537_9fba0102/snapshot/dat/airflow_vt.dat,,,0.25,900,r=0.25,1e-5,0.50,1.00,reference
r0p30,,analysis_runs/manual_run_archive/r0p30_P900_20261104T071153_3b5765f4/snapshot,,0.30,900,r=0.30,1e-5,0.50,1.00,reference
```

| 列 | 指定方法 |
|---|---|
| `case_id` | 必須・重複禁止。ASCII英数字で開始し、英数字、`_`、`-`、単一の`.`で構成。空白・`..`・パス区切り不可。 |
| `flow_path` / `run_path` | 各行でどちらか一方だけ必須。前者は既存の`airflow_vt.dat`、後者はrun/snapshot/dat等のディレクトリ。 |
| `area_path` | 任意。同じdatディレクトリの`area.dat`だけを許可。省略時もその隣接ファイルだけ自動検出。 |
| `stiffness_ratio`, `pressure_pa` | 任意の数値メタデータ。空欄/`NA`は不明。材料条件への変換・推測なし。 |
| `label` | 表示名。省略時case_id。 |
| `dt_s` | 1 solver stepの秒数。空欄/`NA`は`--dt`、既定`1e-5`。保存間隔は実stepのstrideから検証。 |
| `start_s`, `end_s` | 半開解析窓。空欄/`NA`は`--start`/`--end`、既定`[0.50,1.00)`秒。 |
| `notes` | 任意メモ。検出・計算には不使用。 |

余分な列、壊れたCSV、重複/危険なIDなどの構造エラーは開始前に拒否します。
入力不存在・不正値・非有限流量・不等間隔step・負流量・時間被覆不足などの行エラーは、そのケースの失敗として残します。
各runの全入力を事前検証して一覧を表示してから、manifestの順に解析します。

## 共通の交差選択規則（以下は既定3%）

100 ms/20 msの局所ACFは探索スケールにだけ使います。
10%prominence・10%上下振幅・既存MAD床を満たす極大から、領域内最大prominenceを主ピークにします。
副ピークは保存し、主3%列に混入させません。

主ピークのrawサンプルから次の主ピークのrawサンプル直前までの実波形を調べます。
閾値は`0.03 × Upeak`（0流量基準）。最初の厳密な0到達までの最後の実下降交差を選びます。
0未到達なら、その枝内の最初の実測全体最小点までの最後の実下降交差を選びます。
時刻は実2サンプル間の線形補間。未到達時に最小点や仮イベントを代入しません。

全採用イベントの隣接間隔を使用します。欠落をまたぐ長間隔、短間隔を削除・分割せず、
全`N−1`のJitter差項・`N−2`のRAP項・`N−4`のPPQ5項を既存式で集計します。
最終枝は窓末尾で右打ち切りです。
詳細の正本は前回の`analysis_runs/jitter_all_events_20261107/peak_associated_closure_20261107/selection_rules.md`です。

## 出力

```text
<out-dir>/
  report.md
  summary.csv
  inputs_resolved.csv
  failures.csv
  events_3pct.csv
  periods_3pct.csv
  metric_terms_3pct.csv
  area_alignment.csv
  peaks.csv
  cases/<case_id>/           # 完了ケースごとにatomic renameで確定する数値CSV
  figures/jitter_by_stiffness_ratio.png  # 明示ratioと定義済みJitterがある場合のみ
```

`summary.csv`は失敗も含め1ケース1行、入力順です。Jitter/RAP/PPQ5 [%]、イベント/全間隔/差項数、
平均・標準偏差・最小・最大・中央間隔 [ms]、周期CV [%]、主ピーク数・3%未到達数、複数交差枝数、
参考の0/10%到達数、面積対応の集計を記録します。
定義不能の指標は**空欄**で、`INSUFFICIENT_EVENTS`を付けます。
Jitterが定義できてもPPQ5の周期数が足りなければ総合statusは同じですが、定義できた値は保持し、各指標のstatusを併記します。

`event_span_fraction=(last_event_s-first_event_s)/(end_s-start_s)`は観測イベント幅だけを表します。
探索率・物理的閉鎖率・Jitter計算時間率とは呼びません。
圧力・窓・固定設定が異なる条件は比較線を分けます。不明圧力はケースごとに独立点、
不明ratioは表のみです。定義不能点を0%として描かず、失敗点で線を途切れさせます。

`events_3pct.csv`は旧単一runと同じイベントID・主ピークID・時刻と実bracketを保持します。
`periods_3pct.csv`は全`T_i`、`next_T_s=T_(i+1)`、正規化Jitter寄与を含み、末尾のnext/寄与は空欄です。
`metric_terms_3pct.csv`にはRAP/PPQ5も含む全分子・項番号・分母平均を保存します。
`peaks.csv`には主/副ピークを区別して残します。生波形コピー、大量PNG、JSONは生成しません。

`inputs_resolved.csv`には元のパス・解決済みパス、実dt/stride/fs、SHA256、commit、
既存コアとバッチコードのhash、固定設定ID/選択規則を保存します。
`failures.csv`は失敗が0件でもヘッダーを持ちます。
途中結果は`report.md`の`IN_PROGRESS`で識別でき、未完了ケースを成功とは記録しません。

## 面積は独立した補助情報

同じdatディレクトリの`area.dat`だけを使い、既存読込関数の有効48断面最小値を再利用します。
認識する形式はstep+50断面、DATの0-based列2〜49、認識可能なmm²/m²ヘッダーです。
単位不明、異なる形式、時刻不足等は`AREA_FAILED`、欠落は`AREA_NOT_AVAILABLE`とし、流量解析を継続します。
NaN断面は0へ置換しません。全内側断面が非有限の行や負値/無限大は既存読込規則で拒否します。

主ピークraw時刻〜次主ピークraw時刻未満の**面積のネイティブ標本**から最小値を求めます。
同率最小は最早時刻、断面も同率なら最小列番号を代表にし、tie・平坦部・断面切替回数を併記します。
`min_section_index`は保存された50断面の0-based番号（DAT列番号−1）、
`min_section_column_zero_based`は元DATの列番号です。

`lag_ms=1000*(t_Amin_min-t_Ug_3pct)`なので、正なら面積最小が後です。
`Amin_at_3pct_mm2`は同時刻の実標本値、または隣接するネイティブAmin標本の線形補間です。
補間方法を行ごとに明記し、端を越える外挿・波形の再サンプリング・位相シフトはしません。
3%未到達の主ピークにも面積行を残し、そのflow時刻/lagは空欄です。
右打ち切り・先頭領域・面積被覆不足はフラグで示します。

集計のlag中央値/IQRと面積最小値中央値は、3%イベントがあり面積対応に成功した枝が対象です。
近接や小面積だけで完全閉鎖とは判定しません。面積情報で3%イベントを変更しません。

## 終了コードとテスト

- `0`：全ケースを解析済み。定義不能や任意面積のみの失敗は行statusに残します。
- `1`：全ケースの流量解析/事前検証が失敗。
- `2`：部分失敗・未実行ケースあり、またはCLI/manifest構造/出力先のエラー。
- `130`：割り込み。完了ケースを残し、進行中/残りを明示。

既定はエラー後も継続（`--continue-on-error`）。`--fail-fast`は最初の失敗後を`NOT_RUN_FAIL_FAST`にします。

```bash
# バッチ用合成・入出力テスト
python3 -B -m unittest discover -s tools/tests -p test_batch_analyze_airflow_jitter.py -v

# 既存0.25/0.30との全精度回帰も実行（既存結果は非改変）
JITTER_BATCH_REAL_REGRESSION=1 python3 -B -m unittest discover \
  -s tools/tests -p test_batch_analyze_airflow_jitter.py -v
```

回帰値は0.25が61イベント/60間隔・21.925136%、0.30が58イベント/57間隔・16.986170%。
丸めた値だけでなく、旧イベント/主ピークID・全時刻・全間隔・Jitter/RAP/PPQ5の全精度一致を検証します。
人手の全イベント正解照合は未実施であり、バッチ化によって検出の物理的正確さを新たに証明するものではありません。
