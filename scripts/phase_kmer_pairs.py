"""Alignment-independent exact-context screen for candidate alleles on paired fragments."""
import argparse
from collections import Counter,defaultdict
import hashlib
import itertools
import json
from pathlib import Path
import ahocorasick
import mappy as mp

ROOT=Path('data/phasing');POSITIONS=(40209701,40216470,40220612,40221421)

def run(lanes):
    reference=json.loads((ROOT/'ucsc_region.json').read_text())['dna'].upper()
    rows=json.loads((ROOT/'vcf_phase_audit.json').read_text())['region']
    alternate={r['pos']:r['alt'] for r in rows if r['pos'] in POSITIONS}
    automaton=ahocorasick.Automaton()
    for pos in POSITIONS:
        i=pos-1-40180000
        for allele,base in enumerate((reference[i],alternate[pos])):
            sequence=reference[:i]+base+reference[i+1:]
            for offset in range(31):
                seed=sequence[i-offset:i-offset+31]
                for pattern,site_offset in ((seed,offset),(mp.revcomp(seed),30-offset)):
                    old=automaton.get(pattern,[])
                    automaton.add_word(pattern,old+[(pos,allele,site_offset)])
    automaton.make_automaton();edges=defaultdict(Counter);obs=defaultdict(Counter);direct=[];counts=Counter()
    for lane in lanes:
        a=mp.fastx_read(str(ROOT/f'lane{lane}/specific_R1.fastq'))
        b=mp.fastx_read(str(ROOT/f'lane{lane}/specific_R2.fastq'))
        for r1,r2 in itertools.zip_longest(a,b):
            if r1 is None or r2 is None or r1[0]!=r2[0]:raise RuntimeError('paired IDs mismatch')
            counts['pairs_screened']+=1;calls=defaultdict(set)
            for r in (r1,r2):
                for end,hits in automaton.iter(r[1].upper()):
                    start=end-30
                    for pos,allele,offset in hits:
                        if ord(r[2][start+offset])-33>=30:calls[pos].add(allele)
            counts['ambiguous_site_calls']+=sum(len(v)>1 for v in calls.values())
            unique={p:next(iter(v)) for p,v in calls.items() if len(v)==1}
            for p,allele in unique.items():obs[p][allele]+=1
            for p,q in itertools.combinations(sorted(unique),2):edges[(p,q)][f'{unique[p]}{unique[q]}']+=1
            if 40209701 in unique and 40220612 in unique:
                direct.append({'lane':lane,'fragment':hashlib.sha256((str(lane)+':'+r1[0]).encode()).hexdigest()[:24],
                               'alleles':[unique[40209701],unique[40220612]]})
    output={'lanes':lanes,'counts':counts,'site_observations':dict(obs),'pair_edges':[{'positions':k,'allele_pairs':dict(v)} for k,v in edges.items()],
            'direct_candidate_fragments':direct,'limitations':'Exact 31-base contexts, variant base Q>=30; retrieval-conditioned, not endpoint-deduplicated; no genomic mapping-quality filter. Positive phase requires independent alignment and molecular QC.'}
    (ROOT/'kmer_phase_audit.json').write_text(json.dumps(output,indent=2)+'\n')
    print(json.dumps(output,indent=2))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('lanes',type=int,nargs='+');run(p.parse_args().lanes)
