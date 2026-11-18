# OmegaPulldown

Structure-based protein–protein interaction screening with [Boltz-2](https://github.com/jwohlwend/boltz) and [AlphaFold 3](https://github.com/google-deepmind/alphafold3).

Loosely analogous to [AlphaPulldown](https://github.com/KosinskiLab/AlphaPulldown) (Yu et al., *Bioinformatics* 2023) but built on a Nextflow pipeline that runs two independently trained next-generation structure predictors under a single ranking schema.

<picture>
  <source media="(prefers-color-scheme: dark)" srcset="extras/pipeline_subway_dark.png">
  <img alt="OmegaPulldown pipeline" src="extras/pipeline_subway.png">
</picture>

## Features

- **Explicit pair lists** — every run takes `--sequences_fasta` (protein pool) and `--pairs_tsv` (the pairs to predict), allowing for flexible screens that can be shaped like baits x candidates, all-vs-all or something else.
- **Dual backend** — Boltz-2 and AlphaFold 3, selectable with `--predictor boltz2 | alphafold3`
- **Post-translational modifications** — specified via a single ID-keyed TSV/CSV using CCD codes (e.g. SEP, TPO, NAG)
- **Per-predictor MSA strategies** — Boltz-2 uses ColabFold-lineage search with taxonomy tagging; AlphaFold 3 uses its own native 4-database search (UniRef90, MGnify, Small BFD for unpaired MSA; UniProt for paired MSA)
- **MSA caching** — cached alignments are reused across pairs and across runs; proteins with existing cache entries skip MSA search entirely
- **MSA server mode** — Boltz-2 can fetch MSAs from the public ColabFold API (`--use_msa_server true`, the default), avoiding local database requirements for small screens
- **MSA-only mode** — `--msa_only` runs the MSA stages without submitting any GPU work, for validating infrastructure or pre-warming the cache
- **Unified ranking** — ipTM, pTM, pDockQ, and pLDDT reported in a single TSV; shared metrics use bare column names, while predictor-specific scores carry a `boltz2_` or `af3_` prefix
- **Chain-order contract** — column 1 of the pair list is the bait, which becomes chain A in every output structure

## Input model

Every run takes two inputs:

| Flag | Description |
|------|-------------|
| `--sequences_fasta` | A multi-sequence FASTA holding every protein in the study |
| `--pairs_tsv` | The pairs to predict: first two columns are protein IDs (header auto-detected; further columns are ignored) |

The pair list is the screen definition. If you want a cross-product or all-vs-all rather than an explicit list, generate it first by hand or with the `bin/make_pairs.py` helper:

```bash
# Classic pulldown: one bait against a candidate panel
bin/make_pairs.py --baits baits.fasta --candidates candidates.fasta > pairs/pulldown.tsv

# All-vs-all
bin/make_pairs.py --all-vs-all pool.fasta > pairs/all_vs_all.tsv
```

## Pipeline stages

1. **GENERATE_PAIRS** — resolve each ID in `--pairs_tsv` against `--sequences_fasta`, inline both sequences, collapse unordered duplicates, and emit the materialized `pairs.tsv`. This file also defines the MSA universe.
2. **MSA search** — compute MSAs per unique protein, with predictor-specific strategies:
   - *Boltz-2 local mode* (`use_msa_server=false`): MSA_SEARCH_DB_BOLTZ2 (UniRef30) → optionally MSA_SEARCH_ENV_BOLTZ2 (env/metagenomic) → MERGE_TAXONOMY_BOLTZ2 (taxonomy-ID tagging via redis + CSV output)
   - *Boltz-2 server mode* (`use_msa_server=true`, default): skipped entirely; Boltz-2 fetches MSAs from the ColabFold API at prediction time
   - *AlphaFold 3*: MSA_SEARCH_DB (UniRef90, MGnify, Small BFD independently) → MERGE_AF3_NATIVE_MSA (merge with deduplication) + MSA_SEARCH_UNIPROT_PAIRED (UniProt paired MSA)
3. **Input assembly** — MAKE_BOLTZ_YAML or MAKE_AF3_JSON per pair, with optional PTM injection
4. **Structure prediction** — BOLTZ_PREDICT or AF3_PREDICT (GPU)
5. **PARSE_METRICS / PARSE_METRICS_AF3** — extract confidence scores and compute pDockQ per pair
6. **RANK_PREDICTIONS** — aggregate all pairs into a single ranked TSV

Stages 3–6 are skipped under `--msa_only`.

## Quick start

```bash
# Boltz-2 with MSA server (default, simplest — no local databases needed)
nextflow run main.nf -profile slurm \
  --sequences_fasta inputs/sequences/proteins.fasta \
  --pairs_tsv inputs/pairs/actin_pulldown.tsv

# Boltz-2 with local MSA (no public API dependency; for large screens)
nextflow run main.nf -profile slurm \
  --sequences_fasta inputs/sequences/proteins.fasta \
  --pairs_tsv inputs/pairs/all_vs_all.tsv \
  --use_msa_server false \
  --boltz2_uniref30_db /path/to/uniref30_2302_db \
  --taxonomy_rdb /path/to/taxonomy.rdb

# AlphaFold 3 (always uses local AF3-native MSA search)
nextflow run main.nf -profile slurm \
  --predictor alphafold3 \
  --use_msa_server false \
  --sequences_fasta inputs/sequences/proteins.fasta \
  --pairs_tsv inputs/pairs/actin_pulldown.tsv

# Proteomics edge list: predict only reported pairs
nextflow run main.nf -profile slurm \
  --sequences_fasta all_proteins.fasta \
  --pairs_tsv ms_hits.tsv

# MSA-only mode (validate infrastructure, pre-warm cache, no GPU work)
nextflow run main.nf -profile slurm \
  --sequences_fasta all_proteins.fasta \
  --pairs_tsv ms_hits.tsv \
  --use_msa_server false \
  --msa_only true

# Reuse pre-computed alignments (skips MSA search for cached proteins)
nextflow run main.nf -profile slurm \
  --predictor alphafold3 \
  --use_msa_server false \
  --sequences_fasta all_proteins.fasta \
  --pairs_tsv ms_hits.tsv \
  --alignments_dir_af3 results/alignments_af3
```

The `-profile slurm` flag is required for GPU-accelerated prediction on HPC clusters. Without it, Nextflow runs locally and predictions will fail if no GPU is available.

## Profiles

| Profile | Description |
|---------|-------------|
| `standard` | Local execution (default) |
| `slurm` | Slurm scheduler with GPU partitions and escalate-on-retry |
| `sge` | Sun Grid Engine |
| `test` | Quick test with the bundled actin toy set (reduced sampling) |

## Post-translational modifications

Place a file named `ptms.tsv` (or `ptms.csv`) in the same directory as your `--sequences_fasta`, or specify it explicitly with `--ptms`. One table covers the entire run — modifications are resolved by protein ID, not by bait/candidate role:

```tsv
protein_id	position	ccd
MyProtein	84	SEP
MyProtein	142	NAG
```

Codes follow the [RCSB PDB Chemical Component Dictionary](https://www.rcsb.org/ligand).

## MSA strategies

The two predictors were trained on different MSA lineages, so OmegaPulldown maintains **separate alignment directories and search strategies** for each.

### Boltz-2

- **Server mode** (`use_msa_server=true`, default): Boltz-2 fetches MSAs from the ColabFold API at prediction time. Simplest setup: no local databases required. Recommended for small-to-moderate screens.
- **Local mode** (`use_msa_server=false`): pre-computes MSAs via mmseqs2 against UniRef30 (+ optionally an environmental database with `--use_env true`), performs taxonomy-ID tagging via a local redis instance loaded with Boltz-2's `taxonomy.rdb`, and writes per-protein CSV files. Recommended for large all-vs-all screens to avoid API rate limits and redundant remote queries.
- **Cache**: per-protein CSV files at `alignments_dir_boltz2/csv/<protein_id>.csv`.

### AlphaFold 3

- Always uses local MSA search (the ColabFold API produces Boltz-2-lineage MSAs, not AF3-native ones).
- Unpaired MSA: three independent mmseqs2 searches (UniRef90, MGnify, Small BFD) merged with exact-duplicate removal.
- Paired MSA: UniProt search, kept separate so AF3's own species-pairing code can do cross-chain matching at inference time.
- **Cache**: unpaired at `alignments_dir_af3/a3m/<protein_id>.a3m`, paired at `alignments_dir_af3/paired_a3m/<protein_id>.paired.a3m`.

## Sampling strategy

- **Boltz-2**: 5 diffusion samples per prediction (configurable via `diffusion_samples`). Boltz-2 uses conditional flow matching, which preserves structural diversity across samples. Physics-based steering potentials are enabled by default (`use_potentials=true`).
- **AlphaFold 3**: 5 seeds × 1 diffusion sample (configurable via `af3_seeds` and `af3_samples`). Within-seed samples show narrow variance in practice, so seed-to-seed variation (independent MSA subsampling and trunk computation) better captures prediction uncertainty.

## Ranking and metrics

All predictions are ranked by `protein_iptm` (mean across samples) by default. This can be changed via `--rank_by`. The `--rank_by` value is validated at startup.

The two rankings tables share a **comparable** core schema and nothing else. Column naming reflects this:

| Column | Comparable? | Notes |
|--------|-------------|-------|
| `ptm`, `iptm`, `protein_iptm` | Yes | Same definition in both predictors |
| `pdockq`, `n_interface_contacts`, `avg_interface_plddt` | Yes | Computed from CIF by `compute_pdockq.py` |
| `boltz2_confidence_score` | No | 0.8 × complex_plddt + 0.2 × iptm |
| `af3_ranking_score` | No | 0.8 × iptm + 0.2 × ptm + 0.5 × fraction_disordered − 100 × has_clash |
| `boltz2_complex_plddt`, `boltz2_complex_iplddt` | No | Token (per-residue) mean, 0–1 scale |
| `af3_complex_plddt`, `af3_complex_iplddt` | No | Atom-weighted mean from CIF B-factors, 0–100 scale |

The output TSV (`results/rankings_<predictor>.tsv`) reports mean and standard deviation for every metric across the samples for each pair.

## Repository structure

```
OmegaPulldown/
├── main.nf                        # Pipeline entry point and workflow logic
├── nextflow.config                # Parameters, resource labels, and profiles
├── modules/
│   ├── prepare_input.nf           # GENERATE_PAIRS, MAKE_BOLTZ_YAML(_NO_MSA)
│   ├── search_msa_boltz2.nf       # MSA_SEARCH_DB_BOLTZ2, MSA_SEARCH_ENV_BOLTZ2,
│   │                              #   MERGE_TAXONOMY_BOLTZ2
│   ├── search_msa_af3.nf          # MSA_SEARCH_DB, MERGE_AF3_NATIVE_MSA,
│   │                              #   MSA_SEARCH_UNIPROT_PAIRED
│   ├── predict_boltz2.nf          # BOLTZ_PREDICT
│   ├── predict_af3.nf             # MAKE_AF3_JSON, AF3_PREDICT
│   ├── rank.nf                    # PARSE_METRICS, RANK_PREDICTIONS
│   └── rank_af3.nf                # PARSE_METRICS_AF3
├── bin/
│   ├── generate_pairs.py          # Pair-list resolver (IDs → sequences)
│   ├── make_pairs.py              # Pair-list generator (cross-product, all-vs-all)
│   ├── make_boltz_yaml.py         # Boltz-2 YAML input builder with PTM support
│   ├── make_af3_json.py           # AF3 JSON builder with A3M cleaning + PTMs
│   ├── a3m_to_boltz_csv.py        # A3M → Boltz-2 CSV with taxonomy tagging
│   ├── aligned_fasta_to_a3m.py    # Aligned FASTA → A3M (UniProt paired MSA)
│   ├── merge_af3_native_msa.py    # Merge multi-DB A3Ms with deduplication
│   └── compute_pdockq.py          # pDockQ from predicted CIF structures
├── inputs/
│   ├── sequences/
│   │   └── proteins.fasta         # Example protein pool (actin toy set)
│   └── pairs/
│       ├── actin_pulldown.tsv     # Example: actin vs. four candidates
│       └── all_vs_all.tsv         # Example: all pairwise combinations
├── extras/
│   └── pipeline_subway.png        # Pipeline schematic
└── results/                       # Output directory (not tracked in git)
```

## Requirements

- [Nextflow](https://www.nextflow.io/) ≥ 22.10
- [mmseqs2](https://github.com/soedinglab/MMseqs2) (for local MSA computation; not needed if using `use_msa_server=true` with Boltz-2)
- [Boltz-2](https://github.com/jwohlwend/boltz) (conda environment)
- [AlphaFold 3](https://github.com/google-deepmind/alphafold3) (Singularity container). Model weights require a [license agreement with Google DeepMind](https://github.com/google-deepmind/alphafold3/blob/main/WEIGHTS_TERMS_OF_USE.md).
- Python 3.9+ with numpy ≥ 1.20, PyYAML ≥ 5.1, redis ≥ 4.0
- Redis (for Boltz-2 local MSA mode only — used to serve `taxonomy.rdb` during taxonomy tagging)

## Local setup

The shipped `nextflow.config` contains placeholder paths (e.g. `/path/to/...`) for the conda environment, mmseqs2 databases, AlphaFold 3 container, and model weights. Edit these to match your local installation before running.

Key parameters to configure:

| Parameter | Description |
|-----------|-------------|
| `conda_init` | Path to your conda executable |
| `boltz_env` | Name of the conda environment with Boltz-2 installed |
| `boltz2_uniref30_db` | mmseqs2 UniRef30 database (Boltz-2 local MSA mode) |
| `taxonomy_rdb` | Boltz-2 taxonomy.rdb file (Boltz-2 local MSA mode) |
| `af3_sif` | Singularity image for AlphaFold 3 |
| `af3_db_dir` | AlphaFold 3 genetic databases directory |
| `af3_model_dir` | AlphaFold 3 model weights directory |
| `gpu_partition` | Slurm GPU partition name |
| `gpu_partition_large` | Slurm large-GPU partition (fallback on retry) |

## References

- Yu D, Chojnowski G, Rosenthal M, Kosinski J. AlphaPulldown — a python package for protein–protein interaction screens using AlphaFold-Multimer. *Bioinformatics* 2023;39(1):btac749.
- Passaro S et al. Boltz-2: Exploring the frontiers of biomolecular structure prediction. *bioRxiv* 2025. DOI:10.1101/2025.06.14.659707.
- Abramson J et al. Accurate structure prediction of biomolecular interactions with AlphaFold 3. *Nature* 2024;630:493–500.
- Bryant P, Pozzati G, Elofsson A. Improved prediction of protein-protein interactions using AlphaFold2. *Nat Commun* 2022;13:6028.

## License

MIT
