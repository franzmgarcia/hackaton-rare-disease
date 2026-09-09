"""Annotate canonical CDS SNVs locally. Indels require a separate normalizer."""
import json
from pathlib import Path
from itertools import product
from collections import Counter

BASES = 'TCAG'
AMINO = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
CODE = dict(zip((''.join(x) for x in product(BASES, repeat=3)), AMINO))
COMPLEMENT = str.maketrans('ACGT', 'TGCA')


def annotate(row, gene, sequence):
    r = dict(row)
    t = next(t for t in gene['Transcript'] if t.get('is_canonical'))
    r['transcript'] = t['id'] + '.' + str(t['version'])
    if len(r['ref']) != 1 or len(r['alt']) != 1:
        r['effect'] = 'indel_or_complex_requires_normalization'
        coding = [(max(e['start'], t['Translation']['start']),
                   min(e['end'], t['Translation']['end'])) for e in t['Exon']]
        coding = [(a,b) for a,b in coding if a <= b]
        if not any(r['pos'] <= b and r['pos'] + len(r['ref']) >= a for a,b in coding):
            r['effect'] = 'outside_canonical_CDS'
        return r
    coding = [(max(e['start'], t['Translation']['start']),
               min(e['end'], t['Translation']['end'])) for e in t['Exon']]
    coding = sorted(((a,b) for a,b in coding if a <= b), reverse=t['strand']==-1)
    offset = 0
    for a,b in coding:
        if a <= r['pos'] <= b:
            index = offset + (r['pos']-a if t['strand']==1 else b-r['pos'])
            ref = r['ref'] if t['strand']==1 else r['ref'].translate(COMPLEMENT)
            alt = r['alt'] if t['strand']==1 else r['alt'].translate(COMPLEMENT)
            if sequence[index] != ref:
                raise ValueError('Reference mismatch: ' + r['gene'])
            start = index//3*3
            codon = sequence[start:start+3]
            altered = codon[:index%3]+alt+codon[index%3+1:]
            before, after = CODE[codon], CODE[altered]
            effect = 'synonymous' if before==after else 'stop_gained' if after=='*' else 'stop_lost' if before=='*' else 'missense'
            if index//3==0 and before=='M' and after!='M': effect='start_lost'
            r.update(effect=effect, coding_position=index+1, protein_change=f'{before}{index//3+1}{after}')
            return r
        offset += b-a+1
    r['effect']='outside_canonical_CDS'
    return r


if __name__ == '__main__':
    refs=Path('data/reference')
    rows=[]
    for line in Path('data/candidates/panel.jsonl').read_text().splitlines():
        r=json.loads(line)
        gene=json.loads((refs/(r['gene']+'.json')).read_text())
        seq=json.loads((refs/(r['gene']+'_cds.json')).read_text())['seq']
        rows.append(annotate(r,gene,seq))
    Path('data/candidates/panel_cds.json').write_text(json.dumps(rows,indent=2)+'\n')
    print(json.dumps(dict(Counter(r['effect'] for r in rows))))
