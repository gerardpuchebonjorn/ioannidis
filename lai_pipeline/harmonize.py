"""
harmonize.py

Contig name detection and renaming (e.g. chr1 vs 1).
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Dict, List, Optional

from lai_pipeline.utils import LOG


def detect_canonical_chrom_mapping(contigs: List[str]) -> Dict[str, str]:
    """
    Given a list of contig names from a VCF, return a mapping
    of canonical chrom (e.g. '1', 'X') to the actual contig name used
    in that VCF (e.g. 'chr1', 'chrX').
    """
    raise NotImplementedError


def contig_for_canonical_chrom(cfg, vcf_gz: Path, chrom: str) -> str:
    """
    Return the contig name used in vcf_gz for a given canonical chrom.
    E.g. canonical '1' -> 'chr1' if that VCF uses chr-prefixed names.
    """
    raise NotImplementedError


def rename_chrom_if_needed(cfg, in_vcf: Path, old_contig: str, new_contig: str, out_vcf: Path) -> Path:
    """
    Rename a contig in a VCF using bcftools annotate --rename-chrs.
    If old_contig == new_contig, returns in_vcf unchanged.
    """
    raise NotImplementedError