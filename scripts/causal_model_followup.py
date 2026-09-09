"""Follow up structural-size VCF records and selected-read breakpoints, not phase."""
import collections, gzip, json, math
from pathlib import Path
import pysam
O=Path('data/causal_models')
genes=set('BUB1B BUB1 BUB3 MAD1L1 MAD2L1 MAD2L1BP TRIP13 CEP57 CENATAC SMC5 SMC6 SLF1 SLF2 NSMCE1 NSMCE2 NSMCE3 NSMCE4A CEP192 DDX11 RAD18 CDC20 TTK ZW10 ZWILCH KNTC1 CENPE CENPF KNL1 NDC80 NDE1 RABGAP1 ESCO2 SGO1 SGO2 SMC1A SMC3 RAD21 STAG2 PDS5B WAPL NCAPD2 NCAPD3 NCAPH NCAPH2 AURKA AURKB PLK1'.split())
transcripts=collections.defaultdict(list)
with gzip.open(O/'ncbiRefSeq.txt.gz','rt') as f:
    for line in f:
        a=line.rstrip().split('\t')
        transcripts[a[12]].append(dict(id=a[1],chrom=a[2].removeprefix('chr'),start=int(a[4]),end=int(a[5]),cdsstart=int(a[6]),cdsend=int(a[7]),exons=list(zip(map(int,a[9].rstrip(',').split(',')),map(int,a[10].rstrip(',').split(','))))))
events=json.loads((O/'genomewide_structural_screen.json').read_text())['events']
for d in events:
    # Trim identical prefix/suffix to exclude the unaltered VCF anchor.
    ref,alt=d['ref'],d['screened_alt']; start=d['pos']-1
    while ref and alt and ref[0]==alt[0]: ref=ref[1:];alt=alt[1:];start+=1
    while ref and alt and ref[-1]==alt[-1]: ref=ref[:-1];alt=alt[:-1]
    end=start+len(ref)
    d['affected_start0']=start; d['affected_end0']=end
    d['exon_overlaps']=[]; d['coding_overlaps']=[]; distances=[]
    for g in d['overlapping_genes']:
        for t in transcripts[g]:
            if t['chrom']!=d['chrom']: continue
            for a,b in t['exons']:
                overlap=(a<end and b>start) if end>start else a<start<b
                distances.extend((abs(start-a),abs(start-b),abs(end-a),abs(end-b)))
                if overlap:
                    d['exon_overlaps'].append(g+':'+t['id'])
                    c,z=max(a,t['cdsstart']),min(b,t['cdsend'])
                    if c<z and ((c<end and z>start) if end>start else c<start<z): d['coding_overlaps'].append(g+':'+t['id'])
    d['exon_overlaps']=sorted(set(d['exon_overlaps']));d['coding_overlaps']=sorted(set(d['coding_overlaps']))
    d['nearest_exon_boundary_bp']=min(distances) if distances else None
(O/'genomewide_structural_annotated.json').write_text(json.dumps(events,indent=2))
subset=[d for d in events if set(d['overlapping_genes'])&genes]
(O/'alternative_mechanism_candidates.json').write_text(json.dumps(dict(mechanism_genes=sorted(genes),events=subset),indent=2))
bam=pysam.AlignmentFile('data/phasing/target.bam','rb')
clusters=json.loads((O/'selected_bam_audit.json').read_text())['softclip_clusters_ge3']
for c in clusters:
    out=[]
    for r in bam.fetch('chr15',c['pos']-1,c['pos']):
        if r.is_unmapped or r.is_duplicate or r.is_secondary or r.is_supplementary or r.mapping_quality<30: continue
        cig=r.cigartuples
        if c['side']=='right' and r.reference_end==c['pos'] and cig[-1][0]==4 and cig[-1][1]>=20:
            s=r.query_sequence[-cig[-1][1]:]; q=r.query_qualities[-cig[-1][1]:]
            freq=collections.Counter(s)
            out.append(dict(reverse=r.is_reverse,clip_length=len(s),entropy_bits=-sum(n/len(s)*math.log2(n/len(s)) for n in freq.values()),median_base_quality=sorted(q)[len(q)//2],clip_sequence=s,mate_chrom=r.next_reference_name,mate_pos=r.next_reference_start+1,template_length=r.template_length,has_SA=r.has_tag('SA')))
    c['alignments']=out
(O/'softclip_followup.json').write_text(json.dumps(clusters,indent=2))
print('Genome-wide exon-overlapping events',sum(bool(d['exon_overlaps']) for d in events),'coding',sum(bool(d['coding_overlaps']) for d in events))
for d in subset: print(d['chrom'],d['pos'],d['length_change'],d['overlapping_genes'],'exons',d['exon_overlaps'],'coding',d['coding_overlaps'],'distance',d['nearest_exon_boundary_bp'])
for c in clusters: print('clip',c['pos'],[(x['reverse'],x['clip_length'],round(x['entropy_bits'],2),x['median_base_quality'],x['mate_chrom'],x['template_length']) for x in c['alignments']])
