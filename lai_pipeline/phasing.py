"""
phasing.py

Detection of phasing status in a VCF file.
"""

from __future__ import annotations

import logging
from pathlib import Path

from lai_pipeline.utils import LOG


def is_vcf_phased(cfg, vcf_gz: Path, contig: str, max_lines: int = 2000) -> bool:
    """
    Detect whether a VCF is phased by inspecting genotype separators.
    Phased genotypes use '|' (e.g. 0|1), unphased use '/' (e.g. 0/1).
    Checks up to max_lines records for efficiency.
    """
    raise NotImplementedError