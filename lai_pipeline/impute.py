"""
impute.py

Phasing and imputation using Beagle.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Optional

from lai_pipeline.utils import LOG


def run_beagle_phasing(cfg, in_vcf: Path, out_prefix: Path, genetic_map: Optional[Path]) -> Path:
    """
    Phase an unphased VCF using Beagle.
    Returns the path to the phased output VCF.
    """
    raise NotImplementedError


def run_beagle_imputation(cfg, gt_vcf: Path, ref_vcf: Path, out_prefix: Path, genetic_map: Optional[Path]) -> Path:
    """
    Impute missing variants using Beagle with a reference panel.
    Returns the path to the imputed output VCF.
    """
    raise NotImplementedError