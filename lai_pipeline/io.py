"""
io.py

Reading, writing, and indexing VCF files via bcftools.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

from lai_pipeline.utils import LOG


def ensure_index(cfg, vcf_gz: Path, *, prefer: str = "tbi", force: bool = False) -> None:
    """Create or rebuild a .tbi or .csi index for a bgzipped VCF."""
    raise NotImplementedError


def bcftools_count_records(cfg, vcf_gz: Path) -> int:
    """Count the number of variant records in a VCF using bcftools."""
    raise NotImplementedError


def read_samples_from_vcf_header(cfg, vcf: Path) -> List[str]:
    """Return the list of sample names from the VCF header."""
    raise NotImplementedError


def get_vcf_contigs(cfg, vcf_gz: Path) -> List[str]:
    """Return all contig IDs declared in the VCF header."""
    raise NotImplementedError


def extract_chrom_variant_vcf(cfg, input_vcf: Path, contig: str, out_vcf: Path) -> None:
    """Extract biallelic SNPs for a single contig into a new VCF."""
    raise NotImplementedError


def _iter_vcf_data_lines(cfg, vcf: Path) -> Iterable[str]:
    """Stream raw VCF data lines (no header) via bcftools view -H."""
    raise NotImplementedError


def build_key_to_tail_list(cfg, vcf: Path) -> Dict[Tuple[int, str, str], List[List[str]]]:
    """
    Build a mapping of (POS, REF, ALT) -> list of tails (columns[5:]).
    Used to fill genotypes in the final VCF assembly.
    """
    raise NotImplementedError