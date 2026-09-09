"""Bounded hidden-lesion audit; selected BAM is never used to call copy number.

Coordinates in JSON output are 1-based inclusive unless named start0/end0.
Genome-wide structural candidates are a screen of the supplied small-variant VCF,
not a replacement for structural-variant calling from unbiased alignments.
"""
import collections
import gzip
import hashlib
import json
from pathlib import Path
import pysam

OUT = Path('data/causal_models')
VCF = 'data/input.vcf.gz'

def record(r):
    s = next(iter(r.samples.values()))
    return dict(chrom=r.contig, pos=r.pos, end=r.stop, ref=r.ref,
                alts=r.alts, filters=list(r.filter), gt=s.get('GT'),
                ad=s.get('AD'), dp=s.get('DP'), gq=s.get('GQ'))

def main():
    OUT.mkdir(exist_ok=True)
    gene = json.loads(Path('data/reference/BUB1B.json').read_text())
    t = next(t for t in gene['Transcript'] if t['id']=='ENST00000287598')
    exons = sorted(t['Exon'], key=lambda e:e['start'])
    v = pysam.VariantFile(VCF)
    rows = []
    for r in v.fetch('15',40050000,40250000):
        d = record(r)
        d['in_gene'] = gene['start'] <= r.pos <= gene['end']
        d['canonical_exons'] = [i for i,e in enumerate(exons,1)
                                if r.pos <= e['end'] and r.stop >= e['start']]
        d['nearest_canonical_exon_boundary_bp'] = min(abs(r.pos-p) for e in exons for p in (e['start'],e['end']))
        d['any_annotated_exon_overlap'] = any(r.pos <= e['end'] and r.stop >= e['start'] for tr in gene['Transcript'] for e in tr['Exon'])
        d['promoter_heuristic_2kb'] = gene['start']-2000 <= r.pos < gene['start']
        rows.append(d)
    (OUT/'bub1b_locus_audit.json').write_text(json.dumps(rows,indent=2))
    cds = [p for e in exons for p in range(max(e['start'],t['Translation']['start']),min(e['end'],t['Translation']['end'])+1)]
    known_sites = {'regulatory_rs576524605':40117088, 'c2386_minus11':cds[2385]-11,
                   'exon9_acceptor_boundary':exons[8]['start']-1,
                   'c2535_plus2':40212650}
    sites = {name:dict(pos=p, records=[record(r) for r in v.fetch('15',p-1,p)]) for name,p in known_sites.items()}
    bam = pysam.AlignmentFile('data/phasing/target.bam','rb')
    fa = pysam.FastaFile('data/phasing/chr15.fa')
    for name,d in sites.items():
        p=d['pos']; d['reference_base']=fa.fetch('chr15',p-1,p)
        d['within_original_retrieval']=40180001 <= p <= 40245000
        d['selected_bam_bases'] = dict(collections.Counter(
            read.query_sequence[q] for read in bam.fetch('chr15',p-1,p)
            if not(read.is_duplicate or read.is_secondary or read.is_supplementary or read.is_unmapped) and read.mapping_quality>=30
            for q,rp in read.get_aligned_pairs(matches_only=True)
            if rp==p-1 and read.query_qualities[q]>=30)) if d['within_original_retrieval'] else None
    coverage=[]
    for i,e in enumerate(exons,1):
        inside=e['start']>=40180001 and e['end']<=40245000
        cov=None
        if inside:
            counts=bam.count_coverage('chr15',e['start']-1,e['end'],quality_threshold=30,
                read_callback=lambda r:not(r.is_unmapped or r.is_duplicate or r.is_secondary or r.is_supplementary) and r.mapping_quality>=30)
            depths=sorted(map(sum,zip(*counts)))
            cov=dict(median_selected_read_depth=depths[len(depths)//2],bases_ge10=sum(x>=10 for x in depths),bases=len(depths),zero_bases=depths.count(0))
        coverage.append(dict(exon=i,start=e['start'],end=e['end'],within_retrieval=inside,coverage=cov))
    clips=collections.defaultdict(set); large=[]; discord=collections.Counter(); eligible=0
    for r in bam.fetch('chr15',40180000,40221137):
        if r.is_unmapped or r.is_duplicate or r.is_secondary or r.is_supplementary or r.mapping_quality<30: continue
        eligible+=1
        cig=r.cigartuples or []
        for op,n,p in [(cig[0][0],cig[0][1],r.reference_start+1),(cig[-1][0],cig[-1][1],r.reference_end)]:
            if op==4 and n>=20: clips[(p,'left' if p==r.reference_start+1 else 'right')].add(r.query_name)
        if any(op in (1,2) and n>=20 for op,n in cig):
            large.append(dict(start=r.reference_start+1,cigar=r.cigarstring,mapq=r.mapping_quality))
        if r.is_read1 and not r.mate_is_unmapped:
            if r.next_reference_id!=r.reference_id: discord['interchromosomal_pairs']+=1
            elif abs(r.template_length)>2000: discord['same_chromosome_over_2kb_pairs']+=1
    audit=dict(sites=sites,exons=coverage,eligible_selected_alignments=eligible,
        softclip_clusters_ge3=[dict(pos=p,side=s,unique_fragments=len(q)) for (p,s),q in sorted(clips.items()) if len(q)>=3],
        indels_ge20_alignments=large,discordant_pair_counts=dict(discord),
        limits='Selected read retrieval and conservative primary alignment can omit alternate structures; no unbiased CN, SV sensitivity, or genome-wide callability inferred.')
    (OUT/'selected_bam_audit.json').write_text(json.dumps(audit,indent=2))
    # Unrestricted whole-VCF screen, including failed filters and all ALT alleles.
    bins=collections.defaultdict(list)
    annotation=OUT/'ncbiRefSeq.txt.gz'
    if annotation.exists():
        with gzip.open(annotation,'rt') as f:
            for line in f:
                a=line.rstrip().split('\t'); chrom=a[2].removeprefix('chr'); start,end=int(a[4]),int(a[5])
                for b in range(start//100000,(end-1)//100000+1): bins[(chrom,b)].append((start,end,a[12]))
    stats=collections.Counter(); events=[]
    for r in pysam.VariantFile(VCF):
        stats['records']+=1
        for alt in r.alts or ():
            symbolic=alt.startswith('<') or '[' in alt or ']' in alt
            if symbolic: stats['symbolic_or_breakend_alt']+=1
            if alt=='*': stats['spanning_deletion_star_alt']+=1; continue
            if symbolic or abs(len(alt)-len(r.ref))>=50:
                d=record(r);d['screened_alt']=alt;d['length_change']=None if symbolic else len(alt)-len(r.ref)
                d['overlapping_genes']=sorted({g for b in range((r.pos-1)//100000,(r.stop-1)//100000+1) for a,z,g in bins.get((r.contig,b),[]) if a<r.stop and z>r.pos-1})
                events.append(d); stats['structural_size_or_symbolic_alt']+=1
    (OUT/'genomewide_structural_screen.json').write_text(json.dumps(dict(stats=stats,events=events,annotation_sha256=hashlib.sha256(annotation.read_bytes()).hexdigest() if annotation.exists() else None),indent=2))
    print(json.dumps(dict(locus_records=len(rows),gene_records=sum(x['in_gene'] for x in rows),sites=sites,genomewide=stats,clip_clusters=audit['softclip_clusters_ge3'],large_indel_alignments=len(large),discordant=discord),indent=2))

if __name__=='__main__': main()
