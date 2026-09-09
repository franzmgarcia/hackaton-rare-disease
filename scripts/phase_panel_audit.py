"""Check exact candidate alleles in a public phased panel; never upload patient data."""
import json
from pathlib import Path
import pysam

ROOT=Path('data/phasing')
URL='https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/1kGP_high_coverage_Illumina.chr15.filtered.SNV_INDEL_SV_phased_panel.vcf.gz'
tab=pysam.TabixFile(URL,index=str(ROOT/'reference_panel.vcf.gz.tbi'))
chrom='chr15' if 'chr15' in tab.contigs else '15'
header=list(tab.header)
sample_count=len(next(s for s in header if s.startswith('#CHROM')).split('\t'))-9
records=[]
for pos,ref,alt in json.loads((ROOT/'panel_queries.json').read_text()):
    rows=list(tab.fetch(chrom,pos-1,pos))
    exact=[]
    for line in rows:
        x=line.split('\t')
        if int(x[1])!=pos or x[3]!=ref or alt not in x[4].split(','):continue
        idx=x[4].split(',').index(alt)+1
        gi=x[8].split(':').index('GT')
        gts=[v.split(':')[gi].replace('|','/').split('/') for v in x[9:]]
        exact.append(dict(alternate_allele_count=sum(g.count(str(idx)) for g in gts),called_alleles=sum(a!='.' for g in gts for a in g),phased='|' in x[9].split(':')[gi]))
    records.append(dict(pos=pos,ref=ref,alt=alt,overlapping_records=len(rows),exact_matches=exact))
out=dict(source=URL,samples=sample_count,queries=records,
         limitations='Panel absence is not a phase observation. No haplotype phase is assigned to a patient-private allele by this audit.')
(ROOT/'population_panel_audit.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
