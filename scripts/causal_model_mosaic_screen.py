"""Exploratory selected-BAM SNV screen; no calibrated mosaic sensitivity."""
import json
import pysam

b=pysam.AlignmentFile('data/phasing/target.bam','rb')
f=pysam.FastaFile('data/phasing/chr15.fa')
v=pysam.VariantFile('data/input.vcf.gz')
called={r.pos for r in v.fetch('15',40180000,40221137)}
hits=[]; tested=0
for col in b.pileup('chr15',40180000,40221137,truncate=True,stepper='nofilter',min_base_quality=0,max_depth=10000,ignore_overlaps=False):
    p=col.reference_pos; ref=f.fetch('chr15',p,p+1); fragments={}
    for pr in col.pileups:
        r=pr.alignment; q=pr.query_position
        if q is None or r.is_duplicate or r.is_secondary or r.is_supplementary or r.mapping_quality<30 or r.query_qualities[q]<30: continue
        obs=(r.query_sequence[q],r.is_reverse)
        if r.query_name not in fragments: fragments[r.query_name]=obs
        elif fragments[r.query_name][0]!=obs[0]: fragments[r.query_name]=('N',False)
    n=sum(a!='N' for a,s in fragments.values())
    if n<20: continue
    tested+=1
    for alt in 'ACGT':
        if alt==ref: continue
        support=[s for a,s in fragments.values() if a==alt]; k=len(support)
        if k>=4 and .05<=k/n<=.30 and len(set(support))==2:
            hits.append(dict(pos=p+1,ref=ref,alt=alt,alt_fragments=k,total_fragments=n,vaf=k/n,forward=support.count(False),reverse=support.count(True),in_original_vcf=p+1 in called,context=f.fetch('chr15',p-10,p+11)))
out=dict(screen='Exploratory SNV-only low-VAF screen in selected alignments, >=4 fragments, both strands, VAF .05-.30, depth >=20, MQ/BQ >=30; no calibrated somatic error model. Not a mosaic variant caller.',tested_positions=tested,candidates=hits)
open('data/causal_models/mosaic_snv_screen.json','w').write(json.dumps(out,indent=2))
print('Positions:',tested,'candidates:',len(hits))
