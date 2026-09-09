"""Sensitivity analysis for overlooked phase anchors at MAPQ/BQ >=20."""
from pathlib import Path
import subprocess
import json
import phase_align as p

root=p.ROOT;exe=root/'bcftools-1.22/bcftools'
with (root/'lower_quality.log').open('w') as log:
    commands=[
        [str(exe),'mpileup','-f',str(root/'chr15.fa'),'-r',f'chr15:{p.START+1}-{p.END}',
         '-q','20','-Q','20','-a','FORMAT/AD,FORMAT/DP','-Ob','-o',str(root/'lower_quality.bcf'),str(root/'target.bam')],
        [str(exe),'call','-m','-v','-a','GQ','-Ov','-o',str(root/'lower_quality.vcf'),str(root/'lower_quality.bcf')]]
    for command in commands:subprocess.run(command,stdout=log,stderr=subprocess.STDOUT,check=True)
result=p.graph(root/'lower_quality.vcf','lower_quality',mapq=20,baseq=20)
p.whatshap(root/'lower_quality.vcf','lower_quality',mapq=20)
print(json.dumps({k:result[k] for k in ('candidate_observations','graph_phase','any_edge_candidate_connected')}))
