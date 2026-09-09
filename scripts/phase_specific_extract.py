"""Second retrieval pass with seeds occurring at most once in hg38 (or absent alternate seeds)."""
import concurrent.futures as cf
import json
from pathlib import Path
import subprocess
import argparse

ROOT=Path('data/phasing').resolve()
def run(lane):
    d=ROOT/f'lane{lane}'
    if not (d/'complete.json').exists():raise RuntimeError('Initial extraction not complete')
    if (d/'specific_complete.json').exists():return json.loads((d/'specific_complete.json').read_text())
    with (d/'specific_progress.log').open('w') as log:
        r=subprocess.run([str(ROOT/'phase_extract'),str(ROOT/'specific_seeds.txt'),str(d/'R1.fastq'),str(d/'R2.fastq'),str(d/'specific_R1.fastq'),str(d/'specific_R2.fastq')],stderr=log,stdout=subprocess.PIPE,text=True,check=True)
    result=json.loads(r.stdout);result['lane']=lane
    result['limitation']='Retrieval excludes seeds with >=2 exact genomic copies; not a complete error-tolerant mapping of all original reads.'
    expected=json.loads((d/'complete.json').read_text())['pairs_retained']
    if result['pairs_scanned']!=expected:raise RuntimeError('Not all initial pairs scanned')
    (d/'specific_complete.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result),flush=True)
    return result
if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('lanes',type=int,nargs='+');args=parser.parse_args()
    with cf.ThreadPoolExecutor(max_workers=3) as pool:list(pool.map(run,args.lanes))
