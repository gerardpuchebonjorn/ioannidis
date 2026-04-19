# LAI Pipeline

A modular, maintainable refactor of the Local Ancestry Inference (LAI) preprocessing pipeline. This README serves as a summary of what the pipeline and CLI looks like now and which changes have been made.

---

## Project structure

```
ioannidis/
│
├── lai_pipeline/              ← Python package
│   ├── __init__.py
│   ├── models.py              ← shared dataclasses (ToolConfig, Templates, ChromStats, AlleleConcordanceStats)
│   ├── utils.py               ← subprocess helpers and centralized logging
│   ├── io.py                  ← VCF reading, writing, indexing via bcftools + bundle manifest readers
│   ├── harmonize.py           ← contig name detection and renaming (e.g. chr1 vs 1)
│   ├── qc.py                  ← allele concordance check against SNP manifest and VCF normalization
│   ├── phasing.py             ← phasing status detection
│   ├── impute.py              ← phasing and imputation via Beagle
│   ├── assembly.py            ← final VCF reconstruction in exact SNP manifest order
│   └── pipeline.py            ← LAIPipeline class: per-chromosome orchestration
│
├── cli.py                     ← command-line entry point (outside the package by design)
└── pipeline_original.py       ← original pipeline (kept for reference)
```

`cli.py` lives outside `lai_pipeline/` intentionally — it is the entry point that wires user arguments to the pipeline, not part of the library itself.

---

## Module dependency graph

```
cli.py
  └── pipeline.py
        ├── models.py
        ├── utils.py
        ├── io.py        ← depends on utils.py
        ├── harmonize.py ← depends on utils.py, io.py
        ├── qc.py        ← depends on utils.py, io.py, models.py
        ├── phasing.py   ← depends on utils.py, models.py
        ├── impute.py    ← depends on utils.py, io.py, models.py
        └── assembly.py  ← depends on utils.py, io.py, models.py
```

`models.py` was introduced during refactoring to break a circular import between `pipeline.py` and `qc.py` — all shared dataclasses now live there.

---

## Centralized logging

Rather than creating a separate logger in each module, all modules share a single `LOG` instance defined in `utils.py`:

```python
from lai_pipeline.utils import LOG
```

This means `import logging` is absent from most modules to keep logging consistent and centralized across the whole pipeline.

---

## Function ordering in `io.py`

The functions in `io.py` are ordered by logical level rather than following the order in `pipeline_original.py`:

1. Index management (`ensure_index`)
2. Metadata readers (`bcftools_count_records`, `get_vcf_contigs`, `read_samples_from_vcf_header`)
3. Data extraction (`extract_chrom_variant_vcf`)
4. Bundle readers (`load_bundle_manifest`, `bundle_entries_for_chrom`, `available_bundle_chroms`, `combined_snp_manifest_for_chrom`)
5. Iterators and builders (`_iter_pos_ref_alt`, `_iter_vcf_data_lines`, `build_key_to_tail_list`)

In the original file these were spread out as they were added over time.

---

## CLI usage

```bash
python cli.py \
  --input-vcf path/to/target.vcf.gz \
  --workdir path/to/workdir \
  --bundle-dir path/to/pclai_bundle \
  [--reference-vcf-template "path/to/ref/chr{chrom}.vcf.gz"] \
  [--reference-fasta path/to/ref.fa] \
  [--genetic-map-template "path/to/maps/chr{chrom}.map"] \
  [--impute-engine {beagle,none}] \
  [--beagle-jar path/to/beagle.jar] \
  [--bcftools bcftools] \
  [--java java] \
  [--threads 4] \
  [--qc-strict] \
  [--min-exact-match-pct 50.0] \
  [--allow-inversions] \
  [--allow-other-mismatch] \
  [--auto-normalize-on-qc-fail] \
  [--log-level {DEBUG,INFO,WARNING,ERROR}]
```

The CLI validates inputs before starting: it checks that the input VCF and bundle directory exist, creates the workdir if needed, and errors early if `--beagle-jar` is missing when `--impute-engine beagle` is selected.

---

## Test run results

Tested with `chr21.admx.vcf.gz` (654 samples, chr21) against the `pclai_1kg_bundle_cpu` bundle:

```
=== SUMMARY ===
  chr21 | qc=PASS | exact=96.386% | inv=0 | other=84 | records=366615

  filled_from_beagle     = 326432   (imputed by Beagle)
  filled_from_target     = 32407    (present in input)
  missing_unfilled       = 7776     (./.)
  total_manifest_records = 366615   (exact match with bundle)
```

Beagle reduced missing positions from 78,467 (without imputation) to 7,776.

---

## Changes from `pipeline_original.py`

### Bundle-based architecture (major change)

The pipeline now reads the SNP manifest from a PCLAI bundle directory instead of per-chromosome VCF files. This affects:

- **`--model-vcf-template`** is replaced by **`--bundle-dir`** — a single path to the bundle folder containing `manifest.json` and `snp_manifests/*.snps.tsv`
- **`io.py`** — four new functions to read the bundle: `load_bundle_manifest`, `bundle_entries_for_chrom`, `available_bundle_chroms`, `combined_snp_manifest_for_chrom`
- **`qc.py`** — `allele_concordance_check_streaming` replaced by `allele_concordance_check_streaming_vs_manifest`, which compares against a pandas DataFrame instead of streaming a VCF
- **`assembly.py`** — `write_final_vcf_in_model_order` replaced by `write_final_vcf_in_manifest_order`, which iterates over DataFrame rows instead of VCF lines
- **`pipeline.py`** — `LAIPipeline` now receives `bundle_dir` and uses `combined_snp_manifest_for_chrom` to load the SNP manifest per chromosome

### CLI argument renames

| Original | New | Reason |
|---|---|---|
| `--model-sites-template` | `--bundle-dir` | Pipeline now reads from bundle, not per-chrom VCFs |
| `--reference-split-template` | `--reference-vcf-template` | Consistent naming |

### CLI arguments removed

- **`--qc-nonstrict`** — was already marked as deprecated in the original. Removed entirely.
- **`--no-split-beagle-multiallelics`** — removed as multiallelic support is not needed.
- **`minimac4`** as an option for `--impute-engine` — not yet implemented. Reserved as a future option in the codebase but removed from the CLI choices for now. The method itself has also been commented in `pipeline.py`.

### Bug fix

- The error message inside `reference_split_vcf_for` referenced the old argument name `--reference-split-template`. Updated to `--reference-vcf-template` for consistency with the new CLI.

### Refactoring

- `_iter_pos_ref_alt` was implicitly defined inside `allele_concordance_check_streaming` in the original. It is now a proper named function in `io.py`, imported by `qc.py`.
- The four shared dataclasses (`ToolConfig`, `Templates`, `ChromStats`, `AlleleConcordanceStats`) were extracted into `models.py` to resolve a circular import between `pipeline.py` and `qc.py`.
- `import logging` removed from most modules — all use the centralized `LOG` from `utils.py`.

---

## Known limitations

The pipeline assumes the input VCF is already on the correct genome build (GRCh38). Build mismatch detection and LiftOver are not implemented — this is noted as work in progress in the broader project context.
