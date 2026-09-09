"""Streaming VCF quality summary. Never emits sample names or variant rows."""
import argparse
from collections import Counter
import gzip
import hashlib
import json
from pathlib import Path


def summarize(path):
    counts = Counter()
    filters = Counter()
    info_fields = set()
    format_fields = set()
    versions = []
    samples = None
    with gzip.open(path, 'rt') as stream:
        for line in stream:
            if line.startswith('##fileformat='):
                versions.append(line.strip().split('=', 1)[1])
            elif line.startswith('##INFO=<ID='):
                info_fields.add(line.split('ID=', 1)[1].split(',', 1)[0])
            elif line.startswith('##FORMAT=<ID='):
                format_fields.add(line.split('ID=', 1)[1].split(',', 1)[0])
            elif line.startswith('#CHROM\t'):
                samples = max(0, len(line.rstrip('\n').split('\t')) - 9)
            elif not line.startswith('#'):
                fields = line.rstrip('\n').split('\t')
                if samples is None or len(fields) < 8:
                    raise ValueError('Missing VCF header or malformed record')
                counts['records'] += 1
                ref, alt = fields[3:5]
                alleles = alt.split(',')
                counts['multiallelic_records'] += len(alleles) > 1
                for allele in alleles:
                    if allele == '.':
                        counts['missing_alt'] += 1
                    elif allele == '*' or allele.startswith('<') or '[' in allele or ']' in allele:
                        counts['symbolic_or_breakend_alleles'] += 1
                    elif len(ref) == len(allele) == 1:
                        counts['snv_alleles'] += 1
                    elif len(ref) != len(allele):
                        counts['indel_alleles'] += 1
                    else:
                        counts['other_sequence_alleles'] += 1
                filters[fields[6]] += 1
    if samples is None:
        raise ValueError('Missing VCF column header')
    with open(path, 'rb') as stream:
        digest = hashlib.file_digest(stream, 'sha256').hexdigest()
    return dict(sha256=digest, compressed_bytes=path.stat().st_size,
                vcf_versions=versions, sample_count=samples, counts=dict(counts),
                filters=dict(filters), info_fields=sorted(info_fields),
                format_fields=sorted(format_fields))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('vcf', type=Path)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    result = summarize(args.vcf)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + '\n')
    print('QC completed; aggregate summary saved locally.')
