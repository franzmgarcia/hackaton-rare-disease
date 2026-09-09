# MVA Hackathon 2026 - Track 1 methodology

Participant: franzmg. This source-only repository documents an evidence-weighted,
manually curated causal comparison. It contains no patient VCF, reads, alignments,
genotype tables, credentials or analysis outputs. The prediction and methods report
are supplied separately through the official submission UI.

## Approach

Local exact-allele annotation, population evidence, targeted fragment validation,
read-backed phase falsification, bulk allele-balance sensitivity and competing
causal models. Unresolved phase is never treated as evidence for trans. The final
submission is a conditional pair hypothesis, not a proven molecular diagnosis.

## Setup

Use Python 3.14 and the pinned requirements. Build the two C++ retrieval programs
with C++17 and zlib. Install bcftools 1.22. Details and ordered commands are in
REPRODUCE.md. Run synthetic alignment tests without patient data:

```bash
python -m unittest discover -s scripts -p 'test_*.py'
```

## Authorized data only

Obtain the gated organizer dataset independently after accepting its terms. All
inputs and derived outputs must remain in ignored data/. Local reproduction
configuration and historical seed snapshots are not public. The public code has
fixed target intervals as method parameters, not observed genotypes. No source
file should be interpreted as an assertion of a patient's haplotypes.

## Repository policy

Commit only this reviewed source tree. Do not copy the parent workspace or its
data/, logs, private replay inputs, old exploratory reports or credentials into it.
No automatic upload or git push is part of the workflow. Respect the organizers'
30-day post-close deletion requirement for genomic data under your control.

The current scientific interpretation is in the separately supplied final report;
historical report-generating templates are intentionally not part of this tree.
See PRIVACY.md for the official distinction between findings and genomic datasets.
