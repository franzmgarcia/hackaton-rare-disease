"""Local exploratory panel extraction; no remote variant queries or clinical ranking."""
import argparse
from collections import Counter
import gzip
import json
from pathlib import Path


def load_panel(directory):
    panel = []
    for path in sorted(directory.glob('*.json')):
        gene = json.loads(path.read_text())
        if 'Transcript' not in gene:
            continue
        if gene.get('assembly_name') != 'GRCh38':
            raise ValueError(f'Unexpected assembly in {path.name}')
        exons = sorted({(e['start'], e['end']) for t in gene['Transcript']
                        if t.get('biotype') == 'protein_coding'
                        for e in t.get('Exon', [])})
        panel.append(dict(gene=gene['display_name'], chrom=gene['seq_region_name'],
                          start=gene['start'], end=gene['end'], exons=exons))
    if not panel:
        raise ValueError('No panel genes loaded')
    return panel


def extract(vcf, panel, output):
    counts = Counter()
    with gzip.open(vcf, 'rt') as source, output.open('w') as target:
        for line in source:
            if line.startswith('#'):
                continue
            f = line.rstrip('\n').split('\t')
            chrom, pos, ref = f[0].removeprefix('chr'), int(f[1]), f[3]
            end = pos + len(ref) - 1
            genes = [g for g in panel if g['chrom'] == chrom
                     and pos <= g['end'] + 10 and end >= g['start'] - 10]
            if not genes:
                continue
            if len(f) != 10:
                raise ValueError('Expected exactly one sample')
            sample = dict(zip(f[8].split(':'), f[9].split(':')))
            gt = sample.get('GT', '.')
            called = set(gt.replace('|', '/').split('/')) - {'.', '0'}
            for gene in genes:
                counts[gene['gene'] + ':all_records'] += 1
                nearby = any(pos <= b + 10 and end >= a - 10 for a, b in gene['exons'])
                if not nearby:
                    continue
                for index, alt in enumerate(f[4].split(','), 1):
                    if str(index) not in called:
                        continue
                    counts[gene['gene'] + ':exon_or_10bp_flank_alleles'] += 1
                    record = dict(gene=gene['gene'], chrom=chrom, pos=pos, ref=ref,
                                  alt=alt, alt_index=index, filter=f[6], qual=f[5],
                                  gt=gt, dp=sample.get('DP'), gq=sample.get('GQ'),
                                  ad=sample.get('AD'), pid=sample.get('PID'),
                                  pgt=sample.get('PGT'),
                                  region='exon_or_10bp_flank_any_protein_coding_transcript')
                    target.write(json.dumps(record) + '\n')
    return dict(counts)


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('vcf', type=Path)
    p.add_argument('--reference-dir', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    a = p.parse_args()
    a.output.parent.mkdir(parents=True, exist_ok=True)
    counts = extract(a.vcf, load_panel(a.reference_dir), a.output)
    a.output.with_suffix('.counts.json').write_text(json.dumps(counts, indent=2)+'\n')
    print(json.dumps(counts, indent=2))
