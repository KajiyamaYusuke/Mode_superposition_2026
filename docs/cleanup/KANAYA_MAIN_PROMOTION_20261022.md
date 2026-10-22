# 検証済みKanaya実装のmain昇格

2026-10-22。これは既存実装の移行であり、新しい物理モデルや周期的閉鎖の達成を意味しない。

## 移植元と範囲

移植元は `4b4b6d0c854602a76a25fafe0f3493894b63e0f6`。前回報告後に未コミット実装がこのコミットへ保存され、移行開始時はcleanだった。最新版main `be2b6f2219b21f91d747d52fe0c2bc63329581d1` から独立worktreeを作り、src/include/tests/tools/docs/CMakeの52ファイルを選択的に移植した。古い実験ブランチをmainへmergeしていない。

含む機能：独立した旧流体4 cpp、6フラグ切替、R6境界tie修正、旧上流音響回路プロファイル、読取専用診断、診断OFF通常ビルド、VTU抑制、R0開口分析、面積波形比較、全体積VTU/PVDと軽量表面VTP/PVDの後処理出力。

ソース・CMake・テストの物理計算部分は検証済み実装と同じ。後処理4 Pythonファイルのパス設定だけを環境変数／CLIで指定可能にした。モデル式、NRカーネル、接触、Newmark、モード正規化、左右鏡映、入力の物理値は変更していない。取得した重複MDとWindows Zone.Identifierは実装ではないため移植から除外し、移行アーカイブには保存した。

build/output/result、実行バイナリ、生マニフェスト、巨大VTU、個人用入力、pycacheはGitに含めない。`input/param.txt` はmainの内容を維持する。元のdirtyな作業ディレクトリは更新せず保護し、mainのクリーンな本番worktreeを別に用意する。

## デフォルトと明示オプション

無指定は `glottalMethodProfile=modern` と `upstreamAcousticProfile=modern`。診断のビルドデフォルトはOFFで、旧方式を自動的に有効化しない。

```text
glottalMethodProfile = legacy_fluid_all
upstreamAcousticProfile = legacy_main
writeDeformationVtu = 0
```

上の各行は独立した明示オプション。S3は旧流体と旧上流を両方指定する。旧接触・旧積分への切替は未実装で、要求すると拒否する。不正なprofileも拒否する。既存の通常VTU出力はデフォルト有効。

## Releaseと回帰ゲート

```bash
cmake -S . -B build/release -DCMAKE_BUILD_TYPE=Release -DENABLE_KANAYA_DIAGNOSTICS=OFF
cmake --build build/release -j
ctest --test-dir build/release --output-on-failure
cmake -S . -B build/diagnostic -DCMAKE_BUILD_TYPE=Release -DENABLE_KANAYA_DIAGNOSTICS=ON
cmake --build build/diagnostic -j
ctest --test-dir build/diagnostic --output-on-failure
```

icpx 2025.1、OpenMP、O3、fast-mathなし。OFFのCTest5件、ONの6件、Python50件（移動済み参照のSHA照合1件を含む）を検査する。Pythonは既存の解析環境（NumPy2、SciencePlots、VTKを含む）を使用し、環境やNumPy差を物理コード変更で隠さない。関数形式の5テストも省略しない。

移行ゲートは80/300stepの無指定／明示modern、診断OFF、VTU抑制、旧6フラグ、旧上流プロファイル、不正値拒否、元C++の回路初期化9係数との比較を含む。最終パッケージでは20モード・1600Pa・ζ左右0.06・dt=1e-5・20,000stepのS3を1条件だけ再実行し、保存済みのAmin/Ug/圧力/一般化力/接触/音声を同じ条件で比較する。仕事の最終区間除外規約は変更しない。

本番フル検証は前回実測を踏まえ8物理コアへ固定し、OMP_PROC_BIND=close、OMP_PLACES=cores、BLAS/MKL/BLIS=1。CPUトポロジーとcpusetを再確認する。検証結果・所要時間・SHA・本番の実パス・main SHA・削除対象はGitに含めない恒久promotionアーカイブのreport.mdに保存する。

## 後処理の移動可能なパス

通常のソルバ実行に以下の環境変数は不要。過去の実験結果を解析する場合だけ明示する。

| 変数 | 意味／未指定時 |
|---|---|
| KANAYA_CASE_ROOT | 旧run一覧。未指定はrepoのoutput/kanaya_ablation |
| KANAYA_NEXT_ROOT | 新規後処理の出力。未指定はoutput/kanaya_ablation_next |
| KANAYA_ARCHIVE_ROOT | 継承済みアーカイブ。未指定はoutput/kanaya_archive |
| KANAYA_BUILD_DIR | 診断ON build。未指定はbuild/diagnostic |
| KANAYA_PARAMETER_TEMPLATE | 旧matrix実行のパラメータ雛形。未指定は保護対象先頭のinput/param.txt |
| KANAYA_PROTECTED_ROOTS | 保護対象をOSのパス区切りで指定。matrixは未指定なら現在のrepoだけ |
| KANAYA_COMPILER_ROOT | Intel compilerルート。未指定はPATHのicpxから取得 |
| KANAYA_R0_CASE | VTU exporterのR0ケースディレクトリ |
| KANAYA_VTU_SPEC_DIR | exporterの取得済みplan.md／plan_provenance.jsonの場所 |
| KANAYA_NEXT_REFERENCE_ROOT | exporterのR0開口再生検証結果を持つディレクトリ |
| KANAYA_PRISTINE_VTU | q_15001照合用の既存正規VTU |

面積比較は `compare_area_waveforms.py --r0 <R0-case> --s2 <S2-case> --s3 <S3-case> --output <new-dir>`。物理計算は実行せず、area.datだけを入力波形にする。VTU表面化は `export_r0_surface_view.py --volume-root <verified-export-dir>`。

履歴のrunにはmode/NAS/frequency、params、modal_native.binと検証用の参照出力が必要であり、これらはGitには同梱しない。移行アーカイブのsnapshot/outputと以前のprevious/next、恒久visualizationを継承する。古い履歴の生paramsに含まれる外部入力パスはアーカイブ内の原本であり、自動的に別モデルへ置換しない。再実行はSHAと同一入力を確認してから行う。

過去の`prepare`/`production`/`finish`スクリプトは特定実験の監査規約を保持するもので、全てが任意の既存ディレクトリ上で即時再開できるわけではない。新規出力先を使い、参照配置を指定する。既存ケースやアーカイブへ上書きしない。全dtの巨大VTUや新規連成計算を後処理の代用として開始しない。

## 本番運用と保護

本番のmain checkoutには、`input/param.txt` が指定するモード／周波数／NASを同じ相対位置へSHA検証付きで用意する。元のinput/param.txtを20モード・1600Paへ書き換えない。短時間の起動確認にはコピーした別paramsを使い、元入力は保持する。通常の `build/release/simulation` と入力ファイルの組み合わせを検査する。

mainへの反映、本番buildと正常終了、永続アーカイブの照合が全て完了するまで、移植元worktreeとそのローカルブランチを削除しない。削除はGitのworktree管理コマンドで指定対象だけに限定する。
