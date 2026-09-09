"""Mask centromeric neighborhoods and thin SNPs; stress-test exploratory model."""
import gzip
import json
from pathlib import Path
import numpy as np
from mosaic_models import histogram, fit, logprob

ROOT = Path('data/mosaic')

def sensitivity():
    with np.load(ROOT / 'heterozygous_readouts.npz') as archive:
        data = {key: archive[key] for key in archive.files}
    centromeres = {}
    with gzip.open('data/reference/hg38_cytoBandIdeo.txt.gz', 'rt') as stream:
        for line in stream:
            chrom, start, end, _, stain = line.rstrip('\n').split('\t')
            if stain != 'acen' or not chrom[3:].isdigit():
                continue
            chrom = int(chrom[3:])
            previous = centromeres.get(chrom, (int(start), int(end)))
            centromeres[chrom] = (min(previous[0], int(start)), max(previous[1], int(end)))
    keep = np.zeros(len(data['pos']), dtype=bool)
    excluded = 0
    for chrom in range(1, 23):
        start, end = centromeres[chrom]
        ids = np.flatnonzero(data['chrom'] == chrom)
        last = -1000000
        for i in ids:
            pos = int(data['pos'][i])
            if start - 5000000 <= pos - 1 < end + 5000000:
                excluded += 1
                continue
            if pos - last >= 400:
                keep[i] = True
                last = pos
    target = ROOT / 'sensitivity'
    target.mkdir(exist_ok=True)
    np.savez_compressed(target / 'heterozygous_readouts.npz', **{key: values[keep] for key, values in data.items()})
    record = dict(input_sites=len(keep), retained_sites=int(keep.sum()), centromere_excluded=excluded,
                  centromere_intervals_0based_halfopen=centromeres, flank_bp=5000000, min_distance_bp=400,
                  limitations='Not a comprehensive mappability/repeat/segmental-duplication mask; thinning reduces but does not remove dependence.')
    (target / 'mask_audit.json').write_text(json.dumps(record, indent=2) + '\n')
    return record

def simulations():
    # Model-conditional checks, not calibrated detection limits for this patient's data.
    data = np.load(ROOT / 'heterozygous_readouts.npz')
    empirical_n = data['alt_reads'].astype(int) + data['ref_reads'].astype(int)
    results = []
    for seed in (11, 22, 33):
        for fraction in (0., .1, .2, .4):
            rng = np.random.default_rng(seed)
            delta = fraction / (2 * (2 + fraction))
            n = rng.choice(empirical_n, 20000)
            p = .49 + rng.choice([-1, 1], len(n)) * delta
            rho = .01
            concentration = (1 - rho) / rho
            k = rng.binomial(n, rng.beta(p * concentration, (1 - p) * concentration))
            # Apply the same minimum-AD rule as patient extraction; n is already 20..100.
            train = histogram(k[:10000], n[:10000]); hold = histogram(k[10000:], n[10000:])
            null = fit(train); mix = fit(train, True)
            gain = float(np.dot(hold[2], logprob(hold, mix.x, True) - logprob(hold, null.x, False)))
            results.append(dict(seed=seed, simulated_gain_clone_fraction=fraction, true_delta=delta,
                                fitted_delta=float(mix.x[2]), holdout_loglik_gain=gain))
    q = .2
    mean_a = (1 - 2*q) + 2*q
    mean_b = (1 - 2*q) + q + q
    assert abs(mean_a - 1) < 1e-12 and abs(mean_b - 1) < 1e-12
    record = dict(training_sites=10000, test_sites=10000, true_rho=.01, true_center=.49,
                  limitations='Generates from fitted model family. No alignment, genotype ascertainment, regional correlation or biological selection simulation. Not patient power or a clinical threshold.',
                  results=results, cancellation_example=dict(normal_AB=.6, gain_AAB=.2, loss_B=.2,
                  aneuploid_cell_fraction=.4, mean_A=mean_a, mean_B=mean_b, bulk_copy_number=mean_a+mean_b,
                  bulk_B_fraction=mean_b/(mean_a+mean_b)))
    (ROOT / 'simulation_checks.json').write_text(json.dumps(record, indent=2) + '\n')
    return record

if __name__ == '__main__':
    print(json.dumps(dict(sensitivity=sensitivity(), simulations=simulations()), indent=2))
