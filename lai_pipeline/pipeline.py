"""
pipeline.py

LAIPipeline class: orchestrates per-chromosome processing.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional

from lai_pipeline.utils import LOG


@dataclass
class ToolConfig:
    """External tool paths and runtime settings."""
    bcftools: str
    java: str
    beagle_jar: Optional[Path]
    minimac4: str
    threads: int


@dataclass
class Templates:
    """Path templates for per-chromosome model and reference files."""
    model_sites_vcf_template: str
    genetic_map_template: Optional[str]
    reference_split_template: Optional[str] = None


@dataclass
class ChromStats:
    """Per-chromosome processing results."""
    chrom: str
    input_contig: str
    total_model_records: int
    present_exact_in_target: int
    qc_passed: bool
    allele_exact_match_pct: float
    allele_inverted: int
    allele_other_mismatch: int
    target_is_phased: bool
    phased_or_imputed_vcf: Optional[str]
    final_model_order_vcf: str


@dataclass
class AlleleConcordanceStats:
    """Results of allele concordance QC check."""
    shared_pos: int
    exact_match: int
    inverted_ref_alt: int
    other_mismatch: int
    missing_in_model: int
    exact_match_pct: float
    inverted_pct: float
    other_mismatch_pct: float
    examples: List[str]


class LAIPipeline:
    def __init__(
        self,
        cfg: ToolConfig,
        templates: Templates,
        workdir: Path,
        *,
        impute_engine: str,
        window_size: int,
        low_cov_threshold: int,
        qc_strict: bool,
        min_exact_match_pct: float,
        require_zero_inversions: bool,
        require_zero_other_mismatch: bool,
        reference_fasta: Optional[Path],
        auto_normalize_on_qc_fail: bool,
        split_beagle_multiallelics: bool = True,
    ):
        raise NotImplementedError

    def run(self, input_vcf: Path) -> List[ChromStats]:
        """Run the full pipeline on input_vcf, one chromosome at a time."""
        raise NotImplementedError