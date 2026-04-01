"""
qc.py

Allele concordance check and VCF normalization.
"""

from __future__ import annotations

from pathlib import Path
from typing import List, Tuple, Iterable

from lai_pipeline.utils import run, popen_lines, LOG
from lai_pipeline.io import ensure_index
from lai_pipeline.models import AlleleConcordanceStats


def clean_snps_biallelic(cfg, in_vcf: Path, out_vcf: Path) -> Path:
    """
    Filter VCF to keep only biallelic SNPs.
    Removes indels, multiallelic sites, and missing ALT.
    """
    out_vcf.parent.mkdir(parents=True, exist_ok=True)
    run([
        cfg.bcftools, "view",
        "-v", "snps",
        "-m2", "-M2",
        "-e", 'ALT="."',
        "-Oz", "-o", str(out_vcf),
        str(in_vcf),
    ])
    ensure_index(cfg, out_vcf, prefer="tbi", force=True)
    return out_vcf


def normalize_vcf(cfg, in_vcf: Path, fasta: Path, out_vcf: Path) -> Path:
    """
    Normalize a VCF against a reference FASTA using bcftools norm.
    Splits multiallelic sites and left-aligns indels.
    """
    out_vcf.parent.mkdir(parents=True, exist_ok=True)
    run([
        cfg.bcftools, "norm",
        "-m", "-any",
        "-f", str(fasta),
        "--check-ref", "w",
        "-Oz", "-o", str(out_vcf),
        str(in_vcf),
    ])
    ensure_index(cfg, out_vcf, prefer="tbi", force=True)
    return out_vcf


def allele_concordance_check_streaming(
    cfg: ToolConfig,
    *,
    chrom: str,
    target_vcf: Path,
    model_vcf: Path,
    max_examples: int = 12,
) -> AlleleConcordanceStats:
    """
    For each target variant POS, check whether its (REF,ALT) is present among ANY
    model alleles at that POS.
    Classify target record as exact / inverted / other mismatch.
    """
    target_it = iter(_iter_pos_ref_alt(cfg, target_vcf))
    model_it = iter(_iter_pos_ref_alt(cfg, model_vcf))

    def next_model():
        try:
            return next(model_it)
        except StopIteration:
            return None

    mcur = next_model()

    shared = 0
    exact = 0
    inv = 0
    other = 0
    missing_in_model = 0
    examples: List[str] = []

    for tpos, tref, talt in target_it:
        shared += 1

        while mcur is not None and mcur[0] < tpos:
            mcur = next_model()

        if mcur is None or mcur[0] > tpos:
            missing_in_model += 1
            if len(examples) < max_examples:
                examples.append(f"POS={tpos} target={tref}>{talt} (NO_MODEL_AT_POS)")
            continue

        # mcur.pos == tpos: gather all model (ref,alt) pairs for this POS
        model_pairs = set()
        mp = mcur[0]
        while mcur is not None and mcur[0] == mp:
            _p, mref, malt = mcur
            for a in malt.split(","):
                model_pairs.add((mref, a))
            mcur = next_model()

        if (tref, talt) in model_pairs:
            exact += 1
        elif (talt, tref) in model_pairs:
            inv += 1
            if len(examples) < max_examples:
                examples.append(f"POS={tpos} target={tref}>{talt} (INVERTED_vs_model)")
        else:
            other += 1
            if len(examples) < max_examples:
                sample = ", ".join([f"{r}>{a}" for (r, a) in list(model_pairs)[:3]])
                examples.append(f"POS={tpos} target={tref}>{talt} model_candidates=[{sample}] (MISMATCH)")

    exact_pct = (100.0 * exact / shared) if shared else 0.0
    inv_pct = (100.0 * inv / shared) if shared else 0.0
    other_pct = (100.0 * other / shared) if shared else 0.0

    LOG.info("Allele concordance chr%s over TARGET records:", chrom)
    LOG.info("  target_records_checked = %d", shared)
    LOG.info("  exact_match            = %d (%.6f%%)", exact, exact_pct)
    LOG.info("  inverted               = %d (%.6f%%)", inv, inv_pct)
    LOG.info("  other_mismatch         = %d (%.6f%%)", other, other_pct)
    LOG.info("  missing_in_model_at_pos= %d", missing_in_model)

    if examples:
        LOG.warning("chr%s mismatch examples (first %d):", chrom, len(examples))
        for e in examples:
            LOG.warning("  %s", e)

    return AlleleConcordanceStats(
        shared_pos=shared,
        exact_match=exact,
        inverted_ref_alt=inv,
        other_mismatch=other,
        missing_in_model=missing_in_model,
        exact_match_pct=exact_pct,
        inverted_pct=inv_pct,
        other_mismatch_pct=other_pct,
        examples=examples,
    )