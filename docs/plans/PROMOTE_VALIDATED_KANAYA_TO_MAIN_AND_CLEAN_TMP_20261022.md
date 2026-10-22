# Codex実行指示：検証済みKanaya実装をmainへ昇格し、本番環境へ反映してtmp worktreeを整理

作成日：2026-10-22
対象：KajiyamaYusuke/Mode_superposition_2026（WSL上の実体は最初に確認）
依頼：既存の検証済み実装を正式mainへ取り込み、ユーザーが普段使う本番作業ディレクトリを更新し、検証用の/tmp worktreeとローカル一時ブランチを整理する。
**削除は全移行・検証・必要最小限の記録保存が完了した最後だけ。main未反映のソースを消さない。**

## 重要な追加方針：Cドライブ容量優先（2026-10-22）

ユーザーから、Cドライブ容量不足のため**今回の検証結果（大容量データ）を最終的に削除する可能性が高い**との連絡あり。**現時点でユーザーが即時削除を指示したわけではない。自動削除はしない。** ただし新規の1GB級アーカイブ複製は作らず、main移行前後の安全を保ったうえで削除可能な分類とコマンドを報告する。

- **必須保存：** 最終ソース（Git commit/main）、新規ヘッダ/実装/テスト/解析スクリプト、再現可能な入力設定（モデルファイル名、param、実効オプションとコミットSHA）、比較レポート、R0/S2/S3の主要な集計CSV（最大・最小面積、開口周期・周波数、エネルギー、検証成否）、変更・検証manifest。最小限の容量で再現条件を残す。
- **任意保存：** 生波形の小容量CSV（特にR0の最狭面積、流量、変位）。容量に余裕がなければ最小限の時系列だけ。
- **削除候補：** 再生成可能なVTU時系列、巨大なmodal-native履歴、重複した過去の検証output、試行の一時ログ、重複build、過去archiveの重複バックアップ。**元のモデル入力や唯一の原本・mainへの未反映コードは削除しない。**
- 前回の約1.33GBアーカイブ（previous/next）、約6.49GBの生成VTUなどは容量点検対象。既に存在するarchiveをさらに複製してはいけない。Git反映・参照検証が済んだ後なら大容量部分は保存不要と判断できるが、**削除対象一覧・推定回収容量・再生成可否をユーザーへ提示し、削除は明示対象のみ**とする。
- WSL2のデータがCドライブ上の ext4.vhdx にある場合、Linux内部で削除してもWindowsから見たVHDXサイズがすぐには縮まらないことがある。必要なら利用中プロセスを止めて安全に仮想ディスクを縮小する別工程を案内する。**作業中にVHDXを直接削除・圧縮しない。**

## 0. 最優先の前提と禁止事項

- ソースの実体は前回報告の /tmp/simulation-main-kanaya-ablation-20261021 にある。**ローカルブランチ codex/kanaya-legacy-ablation-20261021 のHEADは bc62cf8f3a2c01fb1ed0a33f6d4a6864f6c1d2d5、機能の多くは未コミット変更・未追跡ファイル。ブランチ名をcheckoutするだけでは実装は移らない。**
- 前回の保存先：/home/kajiyama/code/kanaya_validation_archive/20261022_r0_s3/previous/ と next/。S3追加コード、R0の開口解析と結果、VTUエクスポートが追加実装されていればその状態も調べる。
- 原本保護対象：/home/kajiyama/code/simulation-main と /home/kajiyama/code/Mode_superposition_2026 の両方。両者が同一Gitリポジトリのworktreeか別cloneかを git rev-parse --show-toplevel、git remote -v、git worktree list --porcelain で特定する。
- GitHub上の確認時のmainは be2b6f2219b21f91d747d52fe0c2bc63329581d1。ただし実行直前に origin/main をfetchしてSHAを再確認する。mainはその時点の最新版を基準にする。
- **禁止：** 未監査の git reset --hard / git clean -fd / git stash / git checkout -f / git push --force、原本への無条件cp -r、/tmpの直接rm -rf、検証結果・バイナリ・入力の上書き、既存の本番パラメータの勝手な変更。
- ユーザーの目的は「実装の正式反映」。旧流体方式の無条件デフォルト化、音響回路・接触・構造の新規物理変更、追加パラメータ探索は求められていない。**デフォルトは現在のmodernのまま**、旧方式や上流プロファイルは明示オプションで選べる状態を保つ。
- 参照：各worktreeのAGENTS.md、最新2件の comparison_report.md、R0/S2/S3の保存マニフェスト、provenance.json、validation.json、preservation_audit.json、前回のbitwise回帰結果、現在のworktreeのGit status。

## 1. 作業開始時の状態・保護ゲート

1. 対象3ディレクトリの存在、git remote、branch/HEAD、git status --porcelain=v1、git worktree list --porcelain、未追跡一覧、容量を記録する。実際のローカル名・pathが報告と違う場合は発見した値を使い、出典と理由を記録。
2. tmp側から git diff --binary HEAD と git diff --cached --binary を採取。未追跡ファイルは git ls-files --others --exclude-standard で一覧化し、src/include/tests/tools/docs/CMake等の実装と生成出力を分類。**未追跡ファイルはGit diffに含まれない**ため、コピー・コミット対象漏れを防ぐ。
3. tmp実装の**未コミットソースと未追跡ソース・テスト・解析ツール**、必要な小容量param・比較レポート・集計CSVのみを保存。すでにprevious/next/アーカイブに存在する同一SHAのファイルは再コピーしない。Git差分/コミットやファイルリスト・SHA256 manifestで復元可能性を確保し、**新たな大容量promotionアーカイブを作らない**。巨大なVTU・modal-native・ビルド・実行バイナリの複製を必須にしない。元のモデル入力は既存の本番フォルダに保持する。
4. 追加されたR0 VTU時系列エクスポートが実装済みか調べる。**指示書がGitHubに存在するだけなら、機能を実装済みと誤認しない。** 実装されていれば検証済みか確認し、安全に同梱できる場合だけ対象へ含める。
5. mainまたは本番環境に未コミットのユーザー変更がある場合は保護する。**既存のdirty状態を理由に無断で退避・破棄しない。** 独立したクリーンworktreeで統合作業を行う。

## 2. mainから独立統合worktreeを作成する

- 更新取得後の origin/main から新規統合ブランチを作る。例：
~~~bash
git -C /home/kajiyama/code/Mode_superposition_2026 fetch origin
git -C /home/kajiyama/code/Mode_superposition_2026 worktree add \
  -b codex/promote-kanaya-validated-20261022 \
  /home/kajiyama/code/kanaya-main-promotion-20261022 origin/main
~~~
- 上の基点ディレクトリが別Git・存在しない・dirty等の場合は適切な対象repoに読み替え、実際のコマンドとSHAを報告。既存worktreeを消さず、同名path/branchが既存なら別名を選び、無断で再利用しない。
- tmpの変更を**追跡済み差分＋未追跡のソース・ヘッダ・テスト・Pythonツールに分けて選択的に移植**する。tmpブランチの古いHEADをそのままmainにmergeするのは禁止。現行mainからの変更・追加があるため、diffを3-wayで機械適用しても各変更を目視監査し、必要なコンフリクトだけ解消。検証していない物理変更を混入しない。
- 取り込む候補：src/fluid/legacyChannelArea.cpp、legacyFlowCoefficients.cpp、legacyWallPressure.cpp、src/forces/legacySurfacePressureLoad.cpp、対応するinclude/、src/simulation/SimulationParamsおよびparser、src/fluid/flowModel.cppのupstreamAcousticProfile、必要なCMakeLists、各種テスト・tools/kanaya_ablation/解析、正式実装に必要な診断・出力制御。
- 上記は例であって、**実ソースの全変更・全未追跡を網羅して監査すること**。以前のR6境界修正、S3の読取専用観測、通常OFF診断ビルド、音響回路プロファイル、VC?等、実在ファイルを漏らさない。存在しないファイルを作ったと偽らない。
- **含めない：** build/、output/、result/、run生成物、GB級VTU、個人PC絶対パス入りの生マニフェスト、Python __pycache__、不要なバイナリ、修正していないinput/param.txt。入力の例示が必要なら新しい例示ファイルに分離して絶対パスを取り除く。
- パス依存の実行設定は安全な相対パス/引数へ整理するか再現READMEに制約を記載。入力実ファイルは本番環境に存在するか確認。
- 追加コードのデフォルトOFF/modern互換性を確認。診断・旧方式切替を含めても**mainの無指定通常実行はmodern**であることを明示。

## 3. 統合後の必須品質ゲート（CPUは実測最速を使う）

1. 差分と対象全ファイルの確認（git diff --check、git status、CMakeに追加cppが全て列挙されること、ビルド時に未定義参照がないこと、同名関数重複がないこと）。
2. icpx Release、診断通常OFFと必要なら診断ONの2構成を**別build-dir**でビルド。OpenMP有効、fast-mathなし。必要なパスを確認して手元のモデルを使用。CTest全件＋Pythonテスト全件を実行。
3. 無指定modern、明示modern、旧流体6フラグ、upstreamAcousticProfile=legacy_main、異常なprofileの拒否、VTUオプションについて短いsmoke test。最小限、前回と同じ80/300stepのbitwise回帰とS3回路初期化9係数／オラクル比較。入力・バイナリ・出力を同じ条件で比較し、異なる条件同士で数値を比較しない。
4. 本番へ反映する最終パッケージとして、**短時間の再現テストと前回の保存済み数値回帰を優先**。過去の参照CSVと一致を照合でき、実装に物理変更がなければ0.2sフル再計算を自動で追加しない。フル再実行が不可欠な場合もVTU出力OFF・ケース専用の小容量出力で**最大1条件**に限定し、終了後に必要なCSVと検証結果だけ保持する。前回の全bitwise一致がコンパイラ/ビルド環境の影響で維持できなければ誤差の原因・量を記載し、物理差を隠さない。
5. 0.2sの主計算の高速化は前回実測に従う。単独計算ならOMP=8を第一候補、実際のcpuset/NUMA/SMTを確認。OMP_PROC_BIND=close、OMP_PLACES=cores、BLAS/MKL/BLISスレッド=1を候補にし、既知の高速affinityを選ぶ。ビルド並列は利用可能CPU数とメモリに合わせる。独立比較を2本同時に走らせる場合だけ4 threads×2ジョブ等の前回ベンチを使う。**OpenMP 16を機械的に使うこと、同じCPUへの過剰バインド、巨大VTU出力による時間・ディスク浪費を避ける。**
6. 再計算は合否判定に必要な最少本数。既存R0/S2/S3の**集計・参照値を優先**し、スモーク用と検証用を別dirとし、ユーザー原本への出力上書きを禁止。重要な集計CSVとvalidation reportのみ永続保存し、容量見積もりを報告。統合前後のモード数・接触係数・質量・ζ・dt・音響回路が一致することを列挙。

**停止条件**：必要な**ソース/再現設定/小容量の検証記録**が未保存、差分移植漏れ、Release/CTest失敗、modern default変化の説明不能、物理出力差が説明不能、入力データ不足、ローカルmainの保護不能の場合はmainを変更せず報告する。失敗を隠すために物理パラメータを調整しない。

## 4. mainへの反映と「本番環境」の更新

- 品質ゲート通過後、統合worktreeでsrc/include/tests/tools/docs/CMakeだけを**選択的にgit add**し、変更を適切な数のコミットへまとめる。不要な巨大生成ファイル・作業ログ・秘密情報・機器固有パスはGitに含めない。
- 統合ブランチを通常pushし、main向けPRを作成して変更範囲と検証結果・ベースSHAを記載。レビュー不要で明確に全ゲート合格なら、ユーザーの本依頼に基づきmainへ通常mergeしてよい。保護規則やCI未通過があれば強行しない。**force pushは禁止**。コミットSHA/PR/merge後のmain SHAを確認する。
- PRツールが利用できない場合でも、mainが期待SHAから進んでいないことを確認して安全な通常merge/fast-forwardを実施できる。mainが進んでいたら再統合・必要なテストを行う。mainの既存変更を上書きすることは禁止。
- その後、本番ディレクトリを特定する。候補 /home/kajiyama/code/simulation-main と /home/kajiyama/code/Mode_superposition_2026 をそれぞれ確認し、**本番実行に使用している方**へ反映する。cleanであれば git fetch と git merge --ff-only origin/main（またはpull --ff-only）で更新。dirtyならその内容を保持し、勝手にreset/checkoutしない。必要なら別のクリーン本番worktreeを構成し、旧本番を残したままパスを明示して引き渡す。
- 本番ディレクトリで正式buildを作成・実行起動確認。既存のinput/param.txtを変更しない。通常実行 ./simulation ../input/param.txt の動作・参照モデルパスを確認し、既存出力へ上書きする本番フル実行はユーザーの既存データを保護できる場合のみ実施。
- ここまで完了したらmainに必要ファイルが実際に存在することをGitとファイルSHAで再確認。**単にGithub側に指示書MDだけが増えた状態を、本番コード反映完了と報告しない。**

## 5. tmp worktree / 一時ブランチの削除（最後に実行）

消す対象を厳密に限定：
- worktree = /tmp/simulation-main-kanaya-ablation-20261021
- local branch = codex/kanaya-legacy-ablation-20261021
- **この自動整理では消してはいけないもの：** /home/kajiyama/code/kanaya_validation_archive/**、/home/kajiyama/code/simulation-main、本番のMode_superposition_2026、main、検証用とは無関係の別worktree、ユーザーの任意の未コミット変更。古いarchive内の大容量データは**別途削除候補に分類して提示**するが、今回のtmp整理と同時に無断削除しない。

削除の必要条件：
1. mainに取り込むべき追跡・未追跡の全実装が反映済み。統合元とmainをファイルリスト・SHA／意図した差分で監査。
2. 前回および今回の**必要なソース・実効入力設定・小容量の主要比較結果・レポート**がmain/既存の永続保存先に残っている。manifestと主要SHAを確認。GB級生出力・重複VTU・modal-nativeの全量保存は削除ゲートとしない。
3. GitHub main反映と本番ローカルbuild/起動テストが完了。未確認ならtmpは残して保留。
4. tmpのgit statusとworktree情報を最終保存して、別用途の未コミットユーザーファイルがないことを確認。

削除手順（安全な状態でのみ実行すること）：
~~~bash
git -C /home/kajiyama/code/Mode_superposition_2026 worktree list --porcelain
# Gitのworktree管理コマンドで対象/tmpだけを削除する。
git -C /home/kajiyama/code/Mode_superposition_2026 worktree remove \
  /tmp/simulation-main-kanaya-ablation-20261021
# mainへ必要な全ソースが移行し、小容量の再現資料が保存されているが、
# tmpに残る未コミット生成物ゆえ通常removeが拒否された場合に限り、
# 対象パスを再検査して worktree remove --force を利用してよい。
git -C /home/kajiyama/code/Mode_superposition_2026 branch -d \
  codex/kanaya-legacy-ablation-20261021
~~~
- git branch -d が「未merge」のため失敗した場合、mainが実際に**作業内容すべてを包含すること**と既存archiveで復元可能なことを確認した後だけ、git branch -D の使用を許す。元ブランチのHEADは監査記録へ残す。
- tmpがGit管理worktreeでない場合はgit worktree removeを使わない。実体と親リポジトリ・サブディレクトリを確認し、対応する安全な処理を選ぶ。ディレクトリ丸ごとの削除を推測で行わない。
- 今回の統合作業用worktree/branchもmain反映後に、cleanかつ適切な記録を残した上で整理してよい。検証用の古いoutput/kanaya_ablationのアーカイブを消さない。

## 6. 最終成果物と報告

最終的に次を日本語で報告：
- main取り込み後のURL・commit SHA・PR番号、変更した主要ファイル、実際の本番ディレクトリとそのHEAD。
- modernがデフォルトで変わっていないことと切替オプションの動作状況。R0/S3参照と本番検証の一致状況、CTest/Pythonの件数、フル計算の実測所要時間/CPU/使用出力dir。
- 残した小容量の保存先・移行元/先の保護監査とSHA、tmpのworktree/branchを削除したか、残した場合の理由。さらにCドライブ容量圧迫の原因となる**削除候補（絶対パス/容量/再生成性/削除時の影響）**と、削除後にWSL仮想ディスクのサイズ縮小が必要かを簡潔に報告する。
- 未反映・未検証の機能（VTU exporterなど）は明示。旧成功モデル全体の再現や安定リミットサイクルの証明を、今回の本番昇格によって新たに達成したとは言わない。
- **main反映／本番環境更新／tmp削除の3工程を個別に成功・未完了・保留で報告すること。**
