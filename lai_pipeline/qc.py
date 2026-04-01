"""
qc.py

Allele concordance check and VCF normalization.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import List

from lai_pipeline.utils import LOG


def clean_snps_biallelic(cfg, in_vcf: Path, out_vcf: Path) -> Path:
    """
    Filter VCF to keep only biallelic SNPs.
    Removes indels, multiallelic sites, and missing ALT.
    """
    raise NotImplementedError


def normalize_vcf(cfg, in_vcf: Path, fasta: Path, out_vcf: Path) -> Path:
    """
    Normalize a VCF against a reference FASTA using bcftools norm.
    Splits multiallelic sites and left-aligns indels.
    """
    raise NotImplementedError


def allele_concordance_check_streaming(
    cfg,
    *,
    chrom: str,
    target_vcf: Path,
    model_vcf: Path,
    max_examples: int = 12,
):
    """
    For each variant in target_vcf, check whether its (REF, ALT) matches
    the model VCF at the same position. Classifies each record as:
      - exact match
      - inverted (REF/ALT swapped)
      - other mismatch
      - missing in model
    Returns an AlleleConcordanceStats dataclass.
    """
    raise NotImplementedError