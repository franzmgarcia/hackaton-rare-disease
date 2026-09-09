"""Review existing tail and excluded candidates; no phase or score calculation."""
import collections, gzip, json
from pathlib import Path
import pysam
from panel_cds import CODE, COMPLEMENT
O=Path('data/submission_audit'); jobs=[]; wanted=collections.defaultdict(set); panel=[]
for r in json.loads(Path('data/candidates/panel_cds.json').read_text()):
    g=json.loads(Path('data/reference',r['gene']+'.json').read_text()); d=dict(r);hits=[]; distances=[]
    for t in g['Transcript']:
        if t.get('biotype')!='protein_coding' or not t.get('Translation'):continue
        ex=sorted(t['Exon'],key=lambda e:e['start']);cs,ce=t['Translation']['start'],t['Translation']['end']
        for left,right in zip(ex,ex[1:]):
            for b in (left['end'],right['start']):distances.append(min(abs(r['pos']-b),abs(r['pos']+len(r['ref'])-1-b)))
        if r['effect']!='outside_canonical_CDS' or len(r['ref'])!=1 or len(r['alt'])!=1:continue
        coding=[p for e in ex for p in range(max(e['start'],cs),min(e['end'],ce)+1)]
        if t['strand']==-1:coding.reverse()
        if r['pos'] in coding:
            index=coding.index(r['pos']); cp=coding[index//3*3:index//3*3+3]
            h=dict(transcript=t['id'],strand=t['strand'],codon_positions=cp,coding_position=index+1,offset=index%3)
            hits.append(h);wanted['chr'+r['chrom']].update(cp);jobs.append((d,h))
    d['nearest_splice_boundary_any_coding_transcript']=min(distances) if distances else None
    d['alternate_CDS_overlaps']=hits;panel.append(d)
bases={}
with gzip.open('data/reference/hg38.fa.gz','rt') as f:
    chrom=None;offset=0
    for line in f:
        if line.startswith('>'):chrom=line[1:].split()[0];offset=0;continue
        seq=line.strip()
        if chrom in wanted:
            for p in wanted[chrom]:
                if offset<p<=offset+len(seq):bases[(chrom,p)]=seq[p-offset-1].upper()
        offset+=len(seq)
for d,h in jobs:
    seq=''.join(bases[('chr'+d['chrom'],p)] for p in h['codon_positions'])
    alt=d['alt'];ref=d['ref']
    if h['strand']==-1:seq=seq.translate(COMPLEMENT);alt=alt.translate(COMPLEMENT);ref=ref.translate(COMPLEMENT)
    assert seq[h['offset']]==ref
    new=seq[:h['offset']]+alt+seq[h['offset']+1:]
    h.update(ref_codon=seq,alt_codon=new,protein_change=CODE[seq]+str((h['coding_position']-1)//3+1)+CODE[new],effect='synonymous' if CODE[seq]==CODE[new] else 'stop_gained' if CODE[new]=='*' else 'missense')
(O/'tail_panel_review.json').write_text(json.dumps(panel,indent=2))
for d in panel:
    if d['alternate_CDS_overlaps']:print(d['gene'],d['pos'],[(h['transcript'],h['protein_change'],h['effect']) for h in d['alternate_CDS_overlaps']])
bam=pysam.AlignmentFile('data/phasing/target.bam','rb');locus=json.loads(Path('data/causal_models/bub1b_locus_audit.json').read_text())
for d in locus:
    if not d['in_gene'] or d['pos'] in (40209701,40220612):continue
    p=d['pos'];d['independent_selected_fragments']=None
    if 40180001<=p<=40245000 and len(d['ref'])==1 and all(len(a)==1 for a in d['alts']):
        fragments={}
        for r in bam.fetch('chr15',p-1,p):
            if r.is_duplicate or r.is_secondary or r.is_supplementary or r.mapping_quality<30:continue
            for q,rp in r.get_aligned_pairs(matches_only=True):
                if rp==p-1 and r.query_qualities[q]>=30:
                    b=r.query_sequence[q]
                    if r.query_name not in fragments:fragments[r.query_name]=b
                    elif fragments[r.query_name]!=b:fragments[r.query_name]='discordant'
        d['independent_selected_fragments']=dict(collections.Counter(fragments.values()))
tail=[d for d in locus if d['in_gene'] and d['pos'] not in (40209701,40220612)]
(O/'tail_bub1b_review.json').write_text(json.dumps(tail,indent=2))
print('BUB1B support',[(d['pos'],d['independent_selected_fragments']) for d in tail])
