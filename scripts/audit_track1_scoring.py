"""Execute public evaluator using synthetic truth only. Never import groundtruth/app."""
import csv, hashlib, importlib.util, itertools, json, sys, tempfile
from dataclasses import asdict
from pathlib import Path

ROOT=Path('data/submission_audit')
spec=importlib.util.spec_from_file_location('public_track1_evaluation',ROOT/'public/evaluation.py')
ev=importlib.util.module_from_spec(spec);sys.modules[spec.name]=ev;spec.loader.exec_module(ev)
def v(n): return ('chrSynthetic',n,'A','C')
A,B=v(1),v(2);X=[v(i) for i in range(3,14)]
def rows(items):
    return [ev.SubmissionRow(frozenset(s),1-i*.05,i+1) for i,s in enumerate(items)]
strategies={'A':[[A,B]],'B':[[A,B],[A],[B]],
            'C':[[A,B],[A,X[0]],[A,X[1]]],
            'D':[[A],[B]],'E':[[A,B]]+[[x] for x in X[:9]]}
truths={'AB':[A,B],'AX1':[A,X[0]],'AX2':[A,X[1]],'AX9':[A,X[8]],
        'X1X2':[X[0],X[1]],'unlisted':[X[9],X[10]],'single_A':[A]}
results={s:{t:asdict(ev.score_proband('SYNTHETIC',rows(items),frozenset(truth))) for t,truth in truths.items()} for s,items in strategies.items()}
for s in results: print(s,[(t,r['rank_points'],round(r['f_max'],6)) for t,r in results[s].items()])
assert (results['A']['AB']['rank_points'],results['A']['AB']['f_max'])==(100,1)
assert (results['D']['AB']['rank_points'],results['D']['AB']['f_max'])==(50,1)
assert results['B']['AX1']['f_max']==.5
assert abs(results['C']['AX1']['f_max']-.8)<1e-10
assert abs(results['D']['AX1']['f_max']-2/3)<1e-10
# Exact full match at rank 4 overrides the better half-credit at rank 1.
late=ev.score_proband('SYNTHETIC',rows([[A,B],[X[1]],[X[2]],[A,X[0]]]),frozenset([A,X[0]]))
assert late.rank_points==25 and results['A']['AX1']['rank_points']==50
# Strictly lower singleton tails weakly dominate A for every possible pair in this universe.
count=0
for truth in itertools.combinations([A,B]+X,2):
    sa=ev.score_proband('SYNTHETIC',rows(strategies['A']),frozenset(truth))
    se=ev.score_proband('SYNTHETIC',rows(strategies['E']),frozenset(truth))
    assert se.rank_points>=sa.rank_points and se.f_max>=sa.f_max
    count+=1
# Secondary label is informational in this implementation, not excluded.
secondary=rows([[X[0]],[A,B]]);secondary[0].finding_type='secondary'
sr=ev.score_proband('SYNTHETIC',secondary,frozenset([A,B]))
assert sr.rank_points==50 and sr.f_max==.8
# Shared EPCR forces the false variant into the same F threshold as the correct pair.
ties=rows([[A,B],[X[0]]]);ties[1].epcr=ties[0].epcr
assert ev.score_proband('SYNTHETIC',ties,frozenset([A,B])).f_max==.8
# Loader: exact chromosome spelling; sorting and stable tie breaks; row cap.
fields=['proband_id','chrom_1','pos_1','ref_1','alt_1','chrom_2','pos_2','ref_2','alt_2','epcr','finding_type','notes']
with tempfile.TemporaryDirectory() as d:
    p=Path(d)/'synthetic.csv'
    def write(n):
        with p.open('w',newline='') as f:
            w=csv.DictWriter(f,fieldnames=fields);w.writeheader()
            for i in range(n):w.writerow(dict(proband_id='SYNTHETIC',chrom_1='chrSynthetic',pos_1=i+1,ref_1='a',alt_1='c',epcr=.5 if i==0 else .8,finding_type='primary'))
    write(3);loaded=ev.load_submission(str(p))['SYNTHETIC']
    assert [next(iter(r.variants))[1] for r in loaded]==[2,3,1]
    write(11)
    try: ev.load_submission(str(p));raise AssertionError('row cap absent')
    except ValueError: pass
out=dict(strategy_definitions={k:[[z[1] for z in group] for group in vv] for k,vv in strategies.items()},results=results,
         late_full_match=asdict(late),secondary_counterexample=asdict(sr),exhaustive_pair_checks=count,
         public_sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in (ROOT/'public').iterdir() if p.is_file()})
(ROOT/'synthetic_results.json').write_text(json.dumps(out,indent=2))
print('All checks passed; singleton-tail dominance checked for',count,'synthetic causal pairs.')
