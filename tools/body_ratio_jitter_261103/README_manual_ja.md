# 手動計算後の保存と解析

各条件の計算が終わった直後、次の計算で `output/` が上書きされる前に、リポジトリのルートから実行します。

```bash
python3 tools/body_ratio_jitter_261103/save_manual_run.py --label r0p50_P900
```

`--label` は各条件で変更します。英数字、`.`、`-`、`_` が使用できます。計算そのものは実行しません。標準ではリポジトリの `output/` を読み、`analysis_runs/manual_run_archive/` 内に毎回別の記録を作ります。必要なパッケージは NumPy、SciPy、Matplotlib です。

`SIMULATION_RUN_DIR` を使った場合や既存ケースを解析する場合は、出力先または `output/` を含むケースディレクトリを指定します。

```bash
python3 tools/body_ratio_jitter_261103/save_manual_run.py \
  --run-dir analysis_runs/jitter_body_ratio_261103_900_1400_1s/cases/r0p50_P900 \
  --label r0p50_P900
```

標準の解析窓は **[0.50, 1.00) 秒**です。別の計算長を使った場合は、計算範囲内の窓を明示できます。窓の端は流量出力間隔 `5 × dt` に合わせてください。

```bash
python3 tools/body_ratio_jitter_261103/save_manual_run.py \
  --run-dir output --window 0.30 0.50 --label earlier_window
```

solverが保存した `output/txt/params_used.txt` を優先します。これが無い場合だけ `--param input/param.txt` を指定できますが、計算時の値と一致することは確認できず、記録は `user_supplied_not_independently_verified` と表示されます。両方が存在する場合、`--param` の中身が `params_used.txt` と異なれば保存を中止します。保存済み `params_used.txt` がある場合は `--param` は不要です。元の `input/param.txt`、計算結果、既存アーカイブは変更しません。

保存先を変える場合は `--save-dir` を指定します。実行のたびに日時と識別子が付いたフォルダが作られ、同じラベルでも既存の記録を上書きしません。

```text
analysis_runs/manual_run_archive/
  summary.csv                 # 全アーカイブを追記。窓・条件・状態も記録
  <label>_<UTC日時>_<ID>/
    metrics.json              # パラメータ、出典、ハッシュ、判定、指標
    summary_row.csv           # その実行の一覧行
    waveforms.csv             # 秒、Ug[m³/s]、左右変位[mm]、あれば音圧・最小面積
    waveforms.png             # 全時間の流量・左右変位・音圧など
    flow_events.png           # 主窓末尾付近の流量とイベント候補
    flow_spectrum.png         # 主窓の流量スペクトル
    flow_events.csv           # 各検出方法のイベント候補
    flow_cycles.csv           # 各検出方法の周期候補
    snapshot/                 # params_used.txtと解析した元の波形ファイル
```

`summary.csv` の `amplitude_ratio_R_over_L` は、指定窓内の左右**観測点変位**の半peak-to-peak振幅 `((max−min)/2)` による右/左比です。全表面平均の振幅比ではありません。左振幅が実質ゼロなら空欄にします。`metrics.json` には半peak-to-peak、5–95パーセンタイル半幅、交流RMSも左右別に保存します。

Jitterの主信号はUgです。最大値法と自己相関・テンプレート法を比較し、`local`、`RAP`、`PPQ5`、周期CVを記録します。周期抽出に必要な条件を満たさない場合は `UNDEFINED_JITTER` と理由を保存し、0%にはしません。`VALID_*` はあくまで今回の暫定判定です。テンプレート法の閾値は最初の900 Pa波形の観察後に決めたもので、他条件の実波形に対する独立検証や閾値感度は残っています。`flow_events.csv` と `flow_cycles.csv` は候補も含むため、すべてが有効なJitterの周期列とは限りません。

流量の第1列は秒ではなくstepとして読み、`t = step × dt` とします。出力の刻みが5stepで揃い、paramの `nstep` に対応するデータが最後まであること、変位の時刻が一致し有限であることを検査します。実行途中や中断済みの結果、paramの不一致、欠落・非有限データはエラーにします。標準 `output/` に終了状態ファイルが無い場合、必要な波形が揃っていることは確認しますが、solverの終了コードまでは証明できません。

VTUはコピーしません。記録されているモード・surface入力が現在アクセスできる場合は現時点のSHA256も保存しますが、古いrunで当時使用されたファイルと現在のファイルが同一である証明にはなりません。
