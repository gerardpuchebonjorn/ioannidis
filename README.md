# LAI Harmonization Pipeline

A modular, maintainable refactor of the Local Ancestry Inference (LAI) preprocessing pipeline. This tool prepares a target VCF file for inference by harmonizing it against the model's expected sites, handling contig renaming, optional normalization, phasing, imputation, and final VCF assembly in exact model record order.

---

## Project structure

```
ioannidis/
│
├── lai_pipeline/              ← Python package
│   ├── __init__.py
│   ├── models.py              ← shared dataclasses (ToolConfig, Templates, ChromStats, AlleleConcordanceStats)
│   ├── utils.py               ← subprocess helpers and centralized logging
│   ├── io.py                  ← VCF reading, writing, indexing via bcftools
│   ├── harmonize.py           ← contig name detection and renaming (e.g. chr1 vs 1)
│   ├── qc.py                  ← allele concordance check and VCF normalization
│   ├── phasing.py             ← phasing status detection
│   ├── impute.py              ← phasing and imputation via Beagle
│   ├── assembly.py            ← final VCF reconstruction in exact model record order
│   └── pipeline.py            ← LAIPipeline class: per-chromosome orchestration
│
├── cli.py                     ← command-line entry point (outside the package by design)
└── pipeline_original.py       ← original monolithic pipeline (kept for reference)
```

`cli.py` lives outside `lai_pipeline/` intentionally — it is the entry point that wires user arguments to the pipeline, not part of the library itself.

---

## Module dependency graph

```
cli.py
  └── pipeline.py
        ├── models.py
        ├── utils.py
        ├── io.py       ← depends on utils.py
        ├── harmonize.py← depends on utils.py, io.py
        ├── qc.py       ← depends on utils.py, io.py, models.py
        ├── phasing.py  ← depends on utils.py, models.py
        ├── impute.py   ← depends on utils.py, io.py, models.py
        └── assembly.py ← depends on utils.py, io.py, models.py
```

`models.py` was introduced during refactoring to break a circular import between `pipeline.py` and `qc.py` — all shared dataclasses now live there.

---

## Centralized logging

Rather than creating a separate logger in each module, all modules share a single `LOG` instance defined in `utils.py`:

```python
from lai_pipeline.utils import LOG
```

This means `import logging` is absent from most modules — a deliberate design decision to keep logging consistent and centralized across the whole pipeline.

---

## Function ordering in `io.py`

The functions in `io.py` are ordered by logical level rather than following the order in `pipeline_original.py`:

1. Index management (`ensure_index`)
2. Metadata readers (`bcftools_count_records`, `get_vcf_contigs`, `read_samples_from_vcf_header`)
3. Data extraction (`extract_chrom_variant_vcf`)
4. Iterators and builders (`_iter_pos_ref_alt`, `_iter_vcf_data_lines`, `build_key_to_tail_list`)

In the original file these were spread out as they were added over time. The new order reflects the natural dependency flow.

---

## CLI usage

```bash
python cli.py \
  --input-vcf path/to/target.vcf.gz \
  --workdir path/to/workdir \
  --model-vcf-template "path/to/model/chr{chrom}.vcf.gz" \
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
  [--no-split-beagle-multiallelics] \
  [--log-level {DEBUG,INFO,WARNING,ERROR}]
```

The CLI validates inputs before starting: it checks that the input VCF exists, creates the workdir if needed, and errors early if `--beagle-jar` is missing when `--impute-engine beagle` is selected.

---

## Changes from `pipeline_original.py`

### CLI argument renames
| Original | New | Reason |
|---|---|---|
| `--model-sites-template` | `--model-vcf-template` | Clearer — these are VCF files, not just "sites" |
| `--reference-split-template` | `--reference-vcf-template` | Clearer — consistent naming with model argument |

### CLI arguments removed
- **`--qc-nonstrict`** — was already marked as deprecated in the original (`"Deprecated alias; default behavior is nonstrict"`). Removed entirely.
- **`minimac4`** as an option for `--impute-engine` — was never implemented in the original (raised `RuntimeError` immediately). Removed from both the CLI and the `run` method in `pipeline.py`.

### Bug fix
- The error message inside `reference_split_vcf_for` referenced the old argument name `--reference-split-template`. Updated to `--reference-vcf-template` for consistency with the new CLI.

### Refactoring
- `_iter_pos_ref_alt` was implicitly defined inside `allele_concordance_check_streaming` in the original. It is now a proper named function in `io.py`, imported by `qc.py`.
- The four shared dataclasses (`ToolConfig`, `Templates`, `ChromStats`, `AlleleConcordanceStats`) were extracted into `models.py` to resolve a circular import between `pipeline.py` and `qc.py`.

---

## Known limitations

The pipeline assumes the input VCF is already on the correct genome build (GRCh38). Build mismatch detection and LiftOver are not implemented — this is noted as work in progress in the broader project context.
