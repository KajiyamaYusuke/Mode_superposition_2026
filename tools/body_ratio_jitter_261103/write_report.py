#!/usr/bin/env python3
"""Write a factual, conservative report from persisted sweep status and metrics."""
import csv
from pathlib import Path
from collections import Counter
import numpy as np

ROOT=Path(__file__).resolve().parents[2]/'analysis_runs/jitter_body_ratio_261103_900_1400_1s'

def main():
    path=ROOT/'summary/jitter_summary.csv';rows=list(csv.DictReader(path.open()));
    if len(rows)!=36:raise ValueError('expected 36 cases')
    sim=Counter(r['simulation_status'] for r in rows);qc=Counter(r['jitter_status'] for r in rows)
    lines=['# 261103 Body剛性比と流量周期ゆらぎ', '', '## 計算と判定','',f'- ケース数: {len(rows)}。ソルバ状態: {dict(sim)}。',f'- Jitter判定: {dict(qc)}。','- 左Body 12 kPa、Cover 3 kPa、右Body 1.2–11.4 kPa、900/1400 Pa、dt=1e-5秒、各100000ステップ。','- surfaceは従来のM5_surface_T3_d2.nas、VTU出力OFF。入力とバイナリのhashは各caseのstatus.jsonに記録。','- 主窓 [0.50,1.00) 秒。Ugの連続イベント間隔からlocal/RAP/PPQ5/CVを算出。未確定は0で埋めない。','- テンプレート法は候補として記録。基本周期の曖昧さやA/B不一致は未定義として扱い、定常性・イベント位置は図と周期CSVで確認する。','']
    lines.extend(['- 主指標はUg由来。単純な極大法が曖昧でもテンプレート法が品質条件を満たす場合はテンプレート周期を使用し、`primary_method` に区別する。テンプレート法の閾値は900 Pa代表波形の観察後に決めた暫定値であり、感度と独立検証が必要。',''])
    by={}
    for pressure in ('900','1400'):
        q=[r for r in rows if r['pressure_pa']==pressure and r['local']!=''];q.sort(key=lambda r:float(r['local']));by[pressure]={r['ratio']:float(r['local']) for r in q}
        lines.append(f'## {pressure} Pa\n')
        if q:
            lo,hi=q[0],q[-1];lines.append(f"算出可能 {len(q)}/18。Jitter localの最小 {float(lo['local']):.4g}%（比 {lo['ratio']}）、最大 {float(hi['local']):.4g}%（比 {hi['ratio']}）。判定状態と対応する波形を参照。")
        else:lines.append('Jitter localを信頼して算出できる条件はありません。')
        lines.append('')
    common=sorted(set(by['900'])&set(by['1400']),key=float)
    lines.append('## 同一比での圧力差\n')
    if common:
        delta=[by['1400'][r]-by['900'][r] for r in common];lines.append(f'両条件で算出可能な比は{len(common)}点。1400−900 Paの平均差は{np.mean(delta):.4g}ポイント（最大差 {max(delta):.4g}、最小差 {min(delta):.4g}）。')
    else:lines.append('両圧力で算出可能な共通比はありません。')
    lines+=['','## 解釈と制約','','自己相関、波形相似度、ピーク/テンプレート候補はJitterとは別の補助指標。高値から直ちにカオス、乱流、声質、右Body剛性だけの因果関係を結論しない。自動解析で周期境界が曖昧なケースは流量波形と周期CSVを確認する。','', '図と数値は `summary/`、各runの流量と解析は `cases/` 以下に保存。','','解析コードは合成信号で検証しているが、実波形の形状変化やサブハーモニクスによる検出誤りが残り得る。']
    out=ROOT/'summary/report.md'
    with out.open('x') as f:f.write('\n'.join(lines)+'\n')
    print(out)

if __name__=='__main__': main()
