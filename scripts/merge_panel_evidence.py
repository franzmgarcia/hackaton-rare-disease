"""Join exact local candidate alleles to downloaded public evidence."""
import csv
import json
from pathlib import Path

refs=Path('data/reference')
rows=json.loads(Path('data/candidates/panel_cds.json').read_text())
with (refs/'BUB1B_alphamissense_hg38.csv').open() as f:
    am={(r['CHROM'].removeprefix('chr'),int(r['POS']),r['REF'],r['ALT']):r
        for r in csv.DictReader(f)}
out=[]
for row in rows:
    if row['effect'] not in ('missense','stop_gained'):
        continue
    key=(row['chrom'],row['pos'],row['ref'],row['alt'])
    vid='-'.join(map(str,key))
    response=json.loads((refs/f'gnomad_{vid}.json').read_text())
    if response.get('errors'):
        raise ValueError('gnomAD returned an error for '+vid)
    v=response['data']['variant']
    if v is not None and v['variant_id']!=vid:
        raise ValueError('gnomAD allele mismatch')
    row['gnomad_r4']=v
    row['alphamissense']=am.get(key)
    row['population_priority_note']='no usable population record'
    if v:
        freqs=[v[k]['af'] for k in ('exome','genome') if v.get(k) and not v[k]['filters']]
        if freqs:
            top=max(freqs)
            row['population_priority_note']=('common alternate allele (>=1%); low priority for a highly penetrant rare cause' if top>=.01
                else 'population frequency >=0.1%; review disease model and homozygotes' if top>=.001
                else 'rare observed alternate allele; rarity alone does not establish pathogenicity')
    out.append(row)
Path('data/candidates/panel_population_evidence.json').write_text(json.dumps(out,indent=2)+'\n')
print('Joined evidence for',len(out),'coding SNVs. No clinical classifications inferred.')
