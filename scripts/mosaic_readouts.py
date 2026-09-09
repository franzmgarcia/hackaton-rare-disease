"""Extract chromosome-scale downstream readouts from a single-sample variant-only VCF.

AD/(ADref+ADalt) is an unphased alternate-read fraction, NOT a phased BAF.
Depth at variant sites is not unbiased genomic coverage or absolute copy number.
"""
from array import array
from collections import Counter, defaultdict
import gzip
import json
from pathlib import Path
import numpy as np

VCF=Path('data/input.vcf.gz')
OUT=Path('data/mosaic'); OUT.mkdir(exist_ok=True)
cols={k:array('I') for k in ('chrom','pos','ref_reads','alt_reads','dp','gq')}
audit=Counter(); per_chr=defaultdict(Counter); hist={}
with gzip.open(VCF,'rt') as f:
    for line in f:
        if line.startswith('#'): continue
        x=line.rstrip('\n').split('\t'); audit['vcf_records']+=1
        chrom=x[0].removeprefix('chr')
        if chrom not in [str(i) for i in range(1,23)]: continue
        audit['autosomal_records']+=1
        if x[6]!='PASS' or len(x[3])!=1 or len(x[4])!=1 or x[4] not in 'ACGT': continue
        audit['pass_biallelic_autosomal_snvs']+=1
        if len(x)!=10: raise ValueError('Exactly one sample required')
        fields=dict(zip(x[8].split(':'),x[9].split(':')))
        try:
            dp=int(fields.get('DP','.')); gq=int(fields.get('GQ','.'))
            ad=[int(v) for v in fields.get('AD','.').split(',')]
        except ValueError: audit['missing_numeric_fields']+=1; continue
        gt=fields.get('GT','.').replace('|','/')
        c=int(chrom); pos=int(x[1]); per_chr[c]['pass_snv_records']+=1
        per_chr[c]['gt_'+gt]+=1
        if dp>0:
            key=(c,(pos-1)//1000000)
            if key not in hist: hist[key]=np.zeros(501,dtype=np.uint32)
            hist[key][min(dp,500)]+=1
        if gt not in ('0/1','1/0'): continue
        audit['pass_heterozygous_snvs']+=1
        if len(ad)!=2: audit['non_diploid_AD_shape']+=1;continue
        n=sum(ad)
        if gq<30 or not 20<=n<=100 or not 20<=dp<=100: continue
        audit['selected_heterozygous_snvs']+=1
        for key,v in zip(cols,(c,pos,ad[0],ad[1],dp,gq)): cols[key].append(v)

np.savez_compressed(OUT/'heterozygous_readouts.npz',**{k:np.asarray(v,dtype=np.uint32) for k,v in cols.items()})
bins=[]
for (c,b),counts in sorted(hist.items()):
    count=int(counts.sum()); median=int(np.searchsorted(counts.cumsum(),(count+1)//2))
    bins.append(dict(chrom=c,start=b*1000000+1,end=(b+1)*1000000,variant_sites=count,median_variant_site_dp=median))
(OUT/'variant_site_depth_bins.json').write_text(json.dumps(bins)+'\n')
(OUT/'extraction_audit.json').write_text(json.dumps(dict(audit=audit,chromosomes=per_chr,
    selection='autosomes; PASS; biallelic SNV; GT 0/1 or 1/0; GQ>=30; DP and AD sum 20..100',
    limitations=['variant-only ascertainment','diploid genotype caller ascertainment','unphased alternate fraction','no GC/mappability correction','single tissue/sample; tissue unspecified']),indent=2)+'\n')
print(json.dumps(audit))
