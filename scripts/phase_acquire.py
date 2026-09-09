"""Download all lanes in bounded batches, extract target pairs, remove only task scratch FASTQs."""
import concurrent.futures as cf
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import time
from huggingface_hub import hf_hub_download

os.environ.setdefault('HF_HUB_DISABLE_PROGRESS_BARS','1')
ROOT=Path('data/phasing').resolve()
INVENTORY=json.loads((ROOT/'dataset_inventory.json').read_text())
SIZES={r['filename']:r['size'] for r in INVENTORY['files']}

def lane(number):
    target=ROOT/f'lane{number}'
    target.mkdir(exist_ok=True)
    if (target/'complete.json').exists():return {'lane':number,'already_complete':True}
    scratch=ROOT/'scratch'/f'lane{number}'
    scratch.mkdir(parents=True,exist_ok=True)
    names=INVENTORY['lanes'][str(number)]
    def download(name):
        result=hf_hub_download('SageBio/mva-hackathon-2026-data',name,repo_type='dataset',revision=INVENTORY['revision'],local_dir=scratch)
        p=Path(result)
        if p.stat().st_size!=SIZES[name]:raise RuntimeError('size mismatch')
        print(json.dumps({'lane':number,'downloaded_mate':names.index(name)+1,'bytes':p.stat().st_size}),flush=True)
        return p
    with cf.ThreadPoolExecutor(max_workers=2) as pool:paths=list(pool.map(download,names))
    with (target/'scan_progress.log').open('w') as log:
        result=subprocess.run([str(ROOT/'phase_extract'),str(ROOT/'seeds.txt'),*(str(p) for p in paths),str(target/'R1.fastq'),str(target/'R2.fastq')],stderr=log,capture_output=False,stdout=subprocess.PIPE,text=True,check=True)
    audit=json.loads(result.stdout)
    audit.update(lane=number,files=names,revision=INVENTORY['revision'],compressed_bytes=sum(SIZES[n] for n in names),
                 method='all paired records scanned, exact 31-mer candidate retrieval in either mate; both mates retained',
                 limitation='Not an error-tolerant whole-genome mapper; seed-negative reads are not aligned.')
    (target/'complete.json').write_text(json.dumps(audit,indent=2)+'\n')
    # Delete only the newly downloaded, successful task scratch copies. Originals and retained pairs stay.
    for p in paths:p.unlink()
    print(json.dumps(audit),flush=True)
    return audit

if __name__=='__main__':
    import argparse
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--lanes',type=int,nargs='+',default=list(range(1,5)))
    parser.add_argument('--workers',type=int,default=2)
    args=parser.parse_args()
    if any(i not in range(1,5) for i in args.lanes):raise ValueError('lanes must be 1..4')
    free=shutil.disk_usage(ROOT).free
    if free<24*min(args.workers,len(args.lanes))*1024**3:raise RuntimeError('Need 24 GiB free per simultaneous lane; lower concurrency if necessary')
    with cf.ThreadPoolExecutor(max_workers=args.workers) as pool:results=list(pool.map(lane,args.lanes))
    if all((ROOT/f'lane{i}/complete.json').exists() for i in range(1,5)):
        results=[json.loads((ROOT/f'lane{i}/complete.json').read_text()) for i in range(1,5)]
        (ROOT/'acquisition_complete.json').write_text(json.dumps(results,indent=2)+'\n')
