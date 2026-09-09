"""Stage ignored replay inputs from an independently authorized local workspace.

Links reference/raw inputs and copies small historical preparation snapshots.
Never runs analysis, uploads data or copies files to a public path.
"""
import argparse,json,re,shutil
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('--authorized-workspace',type=Path,required=True);a=p.parse_args()
source=a.authorized_workspace.resolve()/'data';dest=Path('data').resolve()
if dest.exists():raise SystemExit('Destination data/ already exists; use a clean checkout for replay.')
if not source.is_dir():raise SystemExit('Authorized source data directory missing.')
dest.mkdir()
for suffix in ('','.tbi'):
    choices=list(source.glob('*.vcf.gz'+suffix))
    if len(choices)!=1:raise SystemExit('Expected one source VCF/index; choose an unambiguous authorized workspace.')
    (dest/('input.vcf.gz'+suffix)).symlink_to(choices[0])
for name in ('reference',):
    (dest/name).symlink_to(source/name,target_is_directory=True)
phase=dest/'phasing';phase.mkdir()
for name in ('seeds.txt','specific_seeds.txt','ucsc_region.json','vcf_phase_audit.json','dataset_inventory.json'):
    shutil.copy2(source/'phasing'/name,phase/name)
inventory=json.loads((phase/'dataset_inventory.json').read_text());inventory['lanes']={}
for lane in range(1,5):
    names=[r['filename'] for r in inventory['files'] if re.search(r'_L00'+str(lane)+r'_R[12]_001.fastq.gz$',r['filename'])]
    if len(names)!=2:raise SystemExit('Lane manifest incomplete')
    inventory['lanes'][str(lane)]=sorted(names)
(phase/'dataset_inventory.json').write_text(json.dumps(inventory,indent=2))
# Public-panel queries are experiment configuration, kept private here.
audit=json.loads((source/'phasing/population_panel_audit.json').read_text())
queries=[[r['pos'],r['ref'],r['alt']] for r in audit['queries']]
(phase/'panel_queries.json').write_text(json.dumps(queries,indent=2))
for name in ('chr15.fa','chr15.fa.fai'):
    (phase/name).symlink_to(source/'phasing'/name)
external=source/'phasing/bcftools-1.22'
if external.exists():(phase/'bcftools-1.22').symlink_to(external,target_is_directory=True)
print('Ignored local inputs staged. FASTQs still require authorized download. No analysis run.')
