# Reproduction instructions

## Boundaries

Run from this repository root using an independently authorized dataset. Outputs
always go to ignored data/. The source export contains no genomic data or saved
patient results. Do not run git add on the original working directory.

Two modes are distinguished: computational replay of the recorded run, and a fresh
analysis from raw data. Exact replay requires the historical reference/seed snapshots
and local preparation JSON as well as the authorized dataset. Those private inputs
are intentionally not public. Reference APIs change; a fresh download is not claimed
to reproduce old snapshots byte for byte. Original preparation included interactive
steps, so this is not represented as a single-command, raw-data-to-final-decision pipeline.
The final causal judgment and EPCR remain manually curated.

## Environment and synthetic verification

```bash
python3.14 -m venv .venv
.venv/bin/python -m pip install -r requirements.txt
.venv/bin/python -m unittest discover -s scripts -p 'test_*.py'
```

Use bcftools 1.22, C++17 and zlib. The original macOS run required rebuilding mappy
from source because its downloaded wheel linked to an unavailable zlib path:

```bash
.venv/bin/python -m pip install --force-reinstall --no-binary=mappy mappy==2.31
```

## Inputs

For exact replay from the original authorized workspace, in a fresh checkout:

```bash
.venv/bin/python scripts/stage_authorized_inputs.py --authorized-workspace /path/to/authorized/workspace
```

This stages a generic input.vcf.gz alias, reference files, lane metadata and historical
seed/preparation files, including the original public-panel query configuration. It does not copy patient outputs into the source tree. Reference
links are read-only inputs by convention; do not edit them. The staging command is an
access convenience, not a transfer of authorization. Otherwise obtain the dataset from:
https://huggingface.co/datasets/SageBio/mva-hackathon-2026-data

Pin revision 59e322d27f399006b398d366d33e703e48a29914. Authenticate interactively with
`hf auth login`; never write tokens in source or shell commands. The VCF alias and TBI
must refer to the one supplied sample. FASTQ filenames are read from the authorized
inventory; they are not copied into public configuration.

Reference sources needed under data/reference/: Ensembl expanded JSON and CDS for
the panel genes, exact public ClinVar/gnomAD response snapshots, UniProt/AlphaMissense,
UCSC hg38.fa.gz and hg38 cytoBandIdeo. Consult the final report's reference list and
the authorized workspace reference manifests for exact hashes. The hidden-hit length
screen also needs UCSC hg38 ncbiRefSeq.txt.gz in data/causal_models/. The large public
reference and external software binaries are downloaded independently, not redistributed.

## Ordered analysis commands

```bash
mkdir -p data/qc data/candidates data/mosaic data/causal_models
.venv/bin/python scripts/qc_vcf.py data/input.vcf.gz --output data/qc/vcf_summary.json
.venv/bin/python scripts/panel_candidates.py data/input.vcf.gz --reference-dir data/reference --output data/candidates/panel.jsonl
.venv/bin/python scripts/annotate_bub1b.py
.venv/bin/python scripts/panel_cds.py
.venv/bin/python scripts/merge_panel_evidence.py
.venv/bin/python scripts/mosaic_readouts.py
.venv/bin/python scripts/mosaic_models.py
.venv/bin/python scripts/mosaic_sensitivity.py
.venv/bin/python scripts/mosaic_models.py --directory data/mosaic/sensitivity
```

Public population/ClinVar annotation snapshots must already be available before their
join commands. Their original retrieval was interactive and rate-limited, not an
automated end-to-end annotation service. The code fails on incomplete/mismatched inputs.

## Targeted read workflow

```bash
c++ -O3 -std=c++17 scripts/phase_extract.cpp -lz -o data/phasing/phase_extract
c++ -O3 -std=c++17 scripts/phase_seed_specificity.cpp -lz -o data/phasing/phase_seed_specificity
data/phasing/phase_seed_specificity data/phasing/seeds.txt data/reference/hg38.fa.gz data/phasing/specific_seeds.txt data/phasing/seed_counts.txt
.venv/bin/python scripts/phase_acquire.py --lanes 1 2 3 4 --workers 1
.venv/bin/python scripts/phase_specific_extract.py 1 2 3 4
.venv/bin/python scripts/phase_align.py --lanes 1 2 3 4
.venv/bin/python scripts/phase_lower_quality.py
.venv/bin/python scripts/phase_kmer_pairs.py 1 2 3 4
.venv/bin/python scripts/phase_panel_audit.py
```

The acquisition script deletes only its newly downloaded successful scratch FASTQs
after extraction. It retains selected pairs and completion metadata. Its disk check
is per concurrent lane; provision additional space for index and retained pairs. Run
with one worker when storage is constrained. No runtime/cost benchmark is promised.
Targeted seeds are a replay input, not a claim of exhaustive error-tolerant WGS mapping.

Place the bcftools 1.22 executable at data/phasing/bcftools-1.22/bcftools. Avoid pysam's
bcftools subprocess wrapper for this workflow; the original run used the standalone
executable. Do not overwrite an index while simultaneously opening it for input.

## Hidden-hit and alternative analysis

```bash
.venv/bin/python scripts/causal_model_audit.py
.venv/bin/python scripts/causal_model_followup.py
.venv/bin/python scripts/causal_model_mosaic_screen.py
mkdir -p data/submission_audit
.venv/bin/python scripts/audit_tail_candidates.py
```

These generate private JSON inventories, not validated SV/CNV diagnoses. Their output
must remain in data/. The local report names original source files; public scripts use
input.vcf.gz as a privacy-neutral alias. These are path-only differences, not new analyses.

## Interpretation and final decision

The quantitative stages generate evidence. Final model ranking, source-ablation
interpretation and CSV choice are curated from that evidence; no automatic code claims
to establish trans or causality. Historical report templates with superseded text are
not distributed. The final submitted methods report is the authoritative interpretation.

## Historical public scoring audit (optional)

The synthetic-only script `audit_track1_scoring.py` expects the public evaluator at
`data/submission_audit/public/evaluation.py`. Obtain that file from the official Space,
record its hash, and create the parent directories before running it. Never download
or import an answer key. This is a historical format audit, not a required biological
reproduction stage, and it was not rerun for the final package.
