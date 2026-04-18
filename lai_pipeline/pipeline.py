"""
pipeline.py

LAIPipeline class: orchestrates per-chromosome processing.
"""

from __future__ import annotations

from pathlib import Path
from typing import List, Optional

from lai_pipeline.utils import run, LOG
from lai_pipeline.models import ToolConfig, Templates, ChromStats, AlleleConcordanceStats
from lai_pipeline.io import (
    ensure_index, get_vcf_contigs, extract_chrom_variant_vcf,
    bcftools_count_records
)
from lai_pipeline.harmonize import (
    detect_canonical_chrom_mapping, contig_for_canonical_chrom,
    rename_chrom_if_needed
)
from lai_pipeline.qc import (
    allele_concordance_check_streaming, clean_snps_biallelic, normalize_vcf
)
from lai_pipeline.phasing import is_vcf_phased
from lai_pipeline.impute import run_beagle_phasing, run_beagle_imputation
from lai_pipeline.assembly import write_final_vcf_in_model_order


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
    ):
        self.cfg = cfg
        self.templates = templates
        self.workdir = workdir
        self.impute_engine = impute_engine

        self.window_size = window_size
        self.low_cov_threshold = low_cov_threshold

        self.qc_strict = qc_strict
        self.min_exact_match_pct = min_exact_match_pct
        self.require_zero_inversions = require_zero_inversions
        self.require_zero_other_mismatch = require_zero_other_mismatch

        self.reference_fasta = reference_fasta
        self.auto_normalize_on_qc_fail = auto_normalize_on_qc_fail

    def model_sites_vcf_for(self, chrom: str) -> Path:
        return Path(self.templates.model_sites_vcf_template.format(chrom=chrom))

    def reference_split_vcf_for(self, chrom: str) -> Path:
        if not self.templates.reference_split_template:
            raise RuntimeError("--reference-vcf-template is required for imputation.")
        return Path(self.templates.reference_split_template.format(chrom=chrom))

    def map_for(self, chrom: str) -> Optional[Path]:
        if not self.templates.genetic_map_template:
            return None
        return Path(self.templates.genetic_map_template.format(chrom=chrom))

    def _prepare_model_for_target_contig(self, chrom: str, target_contig: str, chrom_dir: Path) -> Path:
        ms_raw = self.model_sites_vcf_for(chrom)
        if not ms_raw.exists():
            raise RuntimeError(f"Missing model-sites VCF for chr{chrom}: {ms_raw}")

        ms_contig = contig_for_canonical_chrom(self.cfg, ms_raw, chrom)
        ms_renamed = rename_chrom_if_needed(
            self.cfg,
            ms_raw,
            old_contig=ms_contig,
            new_contig=target_contig,
            out_vcf=chrom_dir / f"modelsites.chr{chrom}.renamed.vcf.gz",
        )

        # IMPORTANT: DO NOT “clean”/dedup the model here.
        # Your network expects exact model record list (including duplicate POS).
        ensure_index(self.cfg, ms_renamed, prefer="tbi", force=True)
        return ms_renamed

    def _prepare_ref_for_target_contig(self, chrom: str, target_contig: str, chrom_dir: Path) -> Path:
        ref_raw = self.reference_split_vcf_for(chrom)
        if not ref_raw.exists():
            raise RuntimeError(f"Missing reference VCF for chr{chrom}: {ref_raw}")

        ref_contig = contig_for_canonical_chrom(self.cfg, ref_raw, chrom)
        ref_renamed = rename_chrom_if_needed(
            self.cfg,
            ref_raw,
            old_contig=ref_contig,
            new_contig=target_contig,
            out_vcf=chrom_dir / f"reference.chr{chrom}.renamed.vcf.gz",
        )
        ensure_index(self.cfg, ref_renamed, prefer="tbi", force=True)
        return ref_renamed

    def _maybe_norm_target(self, chrom: str, target_chr_vcf: Path, chrom_dir: Path, tag: str) -> Path:
        """
        Normalize + biallelic-clean target against FASTA (if provided).
        This changes target representation only (safer for model order requirements).
        """
        if self.reference_fasta is None:
            return target_chr_vcf

        norm = normalize_vcf(
            self.cfg,
            target_chr_vcf,
            self.reference_fasta,
            chrom_dir / f"target.chr{chrom}.{tag}.norm.vcf.gz",
        )
        clean = clean_snps_biallelic(
            self.cfg,
            norm,
            chrom_dir / f"target.chr{chrom}.{tag}.clean.vcf.gz",
        )
        return clean

    def _qc_gate(self, chrom: str, allele_stats: AlleleConcordanceStats) -> bool:
        failures: List[str] = []
        if allele_stats.exact_match_pct < self.min_exact_match_pct:
            failures.append(f"allele_exact_match_pct={allele_stats.exact_match_pct:.6f}% < min_exact_match_pct={self.min_exact_match_pct:.6f}%")
        if self.require_zero_inversions and allele_stats.inverted_ref_alt > 0:
            failures.append(f"inversions={allele_stats.inverted_ref_alt}")
        if self.require_zero_other_mismatch and allele_stats.other_mismatch > 0:
            failures.append(f"other_mismatch={allele_stats.other_mismatch}")

        if failures:
            msg = "QC FAILED for chr{}:\n  - ".format(chrom) + "\n  - ".join(failures)
            LOG.error(msg)
            if self.qc_strict:
                raise RuntimeError(msg)
            return False

        LOG.info("QC PASSED for chr%s", chrom)
        return True

    def run(self, input_vcf: Path) -> List[ChromStats]:
        LOG.info("=== PIPELINE START ===")
        LOG.info("Input VCF: %s", input_vcf)
        self.workdir.mkdir(parents=True, exist_ok=True)

        ensure_index(self.cfg, input_vcf, prefer="tbi", force=False)
        contigs = get_vcf_contigs(self.cfg, input_vcf)
        chrom_to_target_contig = detect_canonical_chrom_mapping(contigs)
        if not chrom_to_target_contig:
            raise RuntimeError("Could not map any canonical chromosomes (1-22,X,Y,MT) in input VCF header.")

        chroms_to_process: List[str] = []
        for chrom in chrom_to_target_contig:
            if self.model_sites_vcf_for(chrom).exists():
                chroms_to_process.append(chrom)
            else:
                LOG.warning("No model-sites VCF for chr%s; skipping", chrom)

        stats: List[ChromStats] = []

        for chrom in chroms_to_process:
            target_contig = chrom_to_target_contig[chrom]
            chrom_dir = self.workdir / f"chr{chrom}"
            chrom_dir.mkdir(parents=True, exist_ok=True)

            LOG.info("------------------------------------------------------------------")
            LOG.info("Processing chr%s (target_contig=%s)", chrom, target_contig)

            # 1) Extract target chrom
            target_chr_vcf = chrom_dir / f"target.chr{chrom}.vcf.gz"
            extract_chrom_variant_vcf(self.cfg, input_vcf, target_contig, target_chr_vcf)

            # 2) Prepare model + reference contig naming to match target contig
            model_vcf = self._prepare_model_for_target_contig(chrom, target_contig, chrom_dir)
            
            if self.impute_engine != "none":
                ref_vcf = self._prepare_ref_for_target_contig(chrom, target_contig, chrom_dir)

            # 3) Optional: normalize target if QC fails (or always, if you want by setting auto flag + fasta)
            target_pre = target_chr_vcf

            # 4) QC pass 1: compare TARGET to MODEL at POS (grouping model alleles by POS)
            allele_stats = allele_concordance_check_streaming(
                self.cfg,
                chrom=chrom,
                target_vcf=target_pre,
                model_vcf=model_vcf,
            )
            qc_passed = self._qc_gate(chrom, allele_stats)

            if (not qc_passed) and self.reference_fasta is not None and self.auto_normalize_on_qc_fail:
                LOG.warning("QC failed for chr%s; retrying QC after target normalization vs FASTA.", chrom)
                target_pre = self._maybe_norm_target(chrom, target_chr_vcf, chrom_dir, tag="QC_RETRY")
                allele_stats = allele_concordance_check_streaming(
                    self.cfg,
                    chrom=chrom,
                    target_vcf=target_pre,
                    model_vcf=model_vcf,
                )
                qc_passed = self._qc_gate(chrom, allele_stats)

            # 5) Decide phasing/imputation input
            target_is_phased = is_vcf_phased(self.cfg, target_pre, target_contig)

            phased_or_imputed: Optional[Path] = None
            if self.impute_engine == "beagle":
                # If target unphased, Beagle imputation will still output phased haplotypes.
                beagle_prefix = chrom_dir / f"target.imputed.by_beagle.chr{chrom}"
                imputed_vcf = run_beagle_imputation(
                    self.cfg,
                    gt_vcf=target_pre,
                    ref_vcf=ref_vcf,
                    out_prefix=beagle_prefix,
                    genetic_map=self.map_for(chrom),
                )
                # Split any multiallelics so we can match biallelic model records
                imputed_vcf = self._split_multiallelic_if_needed(imputed_vcf, chrom, chrom_dir, tag="beagle_out")
                phased_or_imputed = imputed_vcf

            elif self.impute_engine == "none":
                # If you requested none, we still ensure phased output by phasing if needed
                if target_is_phased:
                    phased_or_imputed = target_pre
                else:
                    phased_prefix = chrom_dir / f"target.phased.by_beagle.chr{chrom}"
                    phased_or_imputed = run_beagle_phasing(self.cfg, target_pre, phased_prefix, self.map_for(chrom))
            #elif self.impute_engine == "minimac4":
            #   raise NotImplementedError("minimac4 support is planned but not yet implemented.")
            else:
                raise RuntimeError(f"Unknown impute engine: {self.impute_engine}")

            # 6) FINAL: emit exactly MODEL record list in MODEL order, fill GT from beagle/target
            final_vcf = chrom_dir / f"final.for_model.chr{chrom}.vcf.gz"
            write_final_vcf_in_model_order(
                self.cfg,
                chrom=chrom,
                model_vcf=model_vcf,
                header_source_vcf=phased_or_imputed if phased_or_imputed else target_pre,
                beagle_fill_vcf=phased_or_imputed if phased_or_imputed else None,
                target_fill_vcf=target_pre,
                out_vcf_gz=final_vcf,
            )

            # Optional: log counts that actually matter
            model_n = bcftools_count_records(self.cfg, model_vcf)
            final_n = bcftools_count_records(self.cfg, final_vcf)
            if model_n != final_n:
                LOG.error("chr%s: FINAL != MODEL record count (this should not happen). model=%d final=%d", chrom, model_n, final_n)
            else:
                LOG.info("chr%s: FINAL matches MODEL record count exactly: %d", chrom, final_n)

            stats.append(
                ChromStats(
                    chrom=chrom,
                    input_contig=target_contig,
                    total_model_records=model_n,
                    present_exact_in_target=-1, 
                    qc_passed=qc_passed,
                    allele_exact_match_pct=allele_stats.exact_match_pct,
                    allele_inverted=allele_stats.inverted_ref_alt,
                    allele_other_mismatch=allele_stats.other_mismatch,
                    target_is_phased=target_is_phased,
                    phased_or_imputed_vcf=str(phased_or_imputed) if phased_or_imputed else None,
                    final_model_order_vcf=str(final_vcf),
                )
            )

            LOG.info("chr%s DONE. Final (model-order) VCF: %s", chrom, final_vcf)

        LOG.info("=== PIPELINE END ===")
        return stats