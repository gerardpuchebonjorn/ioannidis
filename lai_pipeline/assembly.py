"""
assembly.py

Final VCF reconstruction in exact model record order.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Optional

from lai_pipeline.utils import LOG


def write_final_vcf_in_model_order(
    cfg,
    *,
    chrom: str,
    model_vcf: Path,
    header_source_vcf: Path,
    beagle_fill_vcf: Optional[Path],
    target_fill_vcf: Path,
    out_vcf_gz: Path,
) -> None:
    """
    Emit the final VCF in the exact order of records in model_vcf.
    For each model record (POS, REF, ALT), fill FORMAT + samples from:
      1. beagle_fill_vcf  (preferred, if available and allele matches)
      2. target_fill_vcf  (fallback)
      3. missing (./.     if neither source has the record)
    Output is bgzipped and indexed.
    """
    raise NotImplementedError