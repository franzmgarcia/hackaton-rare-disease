"""Match SNVs by exact SPDI and check their CDS codons against public references."""
import json
from pathlib import Path

reference = Path('data/reference')
index = {}
seen = set()
for path in reference.glob('clinvar_bub1b_summary_*.json'):
    result = json.loads(path.read_text())['result']
    for uid in result['uids']:
        seen.add(uid)
        record = result[uid]
        for allele in record['variation_set']:
            index.setdefault(allele.get('canonical_spdi', ''), []).append(record)
expected = set(json.loads((reference / 'clinvar_bub1b_search.json').read_text())['esearchresult']['idlist'])
assert seen == expected, 'ClinVar download is incomplete'
gene = json.loads((reference / 'BUB1B.json').read_text())
transcript = next(t for t in gene['Transcript'] if t.get('is_canonical'))
assert gene['assembly_name'] == 'GRCh38' and transcript['strand'] == 1
sequence = json.loads((reference / 'BUB1B_cds.json').read_text())['seq']
coding = sorted((max(e['start'], transcript['Translation']['start']),
                 min(e['end'], transcript['Translation']['end'])) for e in transcript['Exon'])
coding = [(a, b) for a, b in coding if a <= b]
rows = []
for line in Path('data/candidates/panel.jsonl').read_text().splitlines():
    row = json.loads(line)
    if row['gene'] != 'BUB1B':
        continue
    assert len(row['ref']) == len(row['alt']) == 1, 'This checker supports SNVs only'
    spdi = f"NC_000015.10:{row['pos'] - 1}:{row['ref']}:{row['alt']}"
    row['clinvar_exact_matches'] = [dict(accession=h['accession_version'], title=h['title'],
        classification=h['germline_classification']) for h in index.get(spdi, [])]
    offset = 0
    for start, end in coding:
        if start <= row['pos'] <= end:
            position = offset + row['pos'] - start + 1
            assert sequence[position - 1] == row['ref'], 'Reference base mismatch'
            codon_start = ((position - 1) // 3) * 3
            codon = sequence[codon_start:codon_start + 3]
            within = (position - 1) % 3
            alternate = codon[:within] + row['alt'] + codon[within + 1:]
            row['local_cds_check'] = dict(transcript=transcript['id'] + '.' + str(transcript['version']),
                coding_position=position, protein_position=codon_start // 3 + 1,
                reference_codon=codon, alternate_codon=alternate, reference_base_matches=True)
            break
        offset += end - start + 1
    else:
        row['local_cds_check'] = {'status': 'outside canonical CDS'}
    rows.append(row)
Path('data/candidates/bub1b_evidence.json').write_text(json.dumps(rows, indent=2) + '\n')
print(f'Annotated {len(rows)} BUB1B alleles with exact matching and local CDS checks.')
