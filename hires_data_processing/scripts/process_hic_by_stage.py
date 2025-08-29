#!/usr/bin/env python3
"""
Unified schicluster pipeline runner:
- Stage loop: prepare-impute (hicluster) -> snakemake
- Supports relative or absolute base_dir
- All parameters configurable via CLI
"""

import os
import subprocess
import argparse
import logging
from datetime import datetime

# ---------------------------
# Logging
# ---------------------------
def setup_logging(log_dir: str) -> logging.Logger:
    os.makedirs(log_dir, exist_ok=True)
    ts = datetime.now().strftime("%Y%m%d_%H%M%S")
    log_file = os.path.join(log_dir, f"schicluster_processing_{ts}.log")
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s",
        handlers=[logging.FileHandler(log_file), logging.StreamHandler()],
    )
    return logging.getLogger(__name__)

# ---------------------------
# Core runners
# ---------------------------
def run_schicluster_prepare(
    stage: str,
    contact_table: str,
    output_dir: str,
    chrom_size_path: str,
    *,
    resolution: int = 100_000,
    batch_size: int = 1536,
    cpu_per_job: int = 30,
    pad: int = 1,
    chr1: int = 1,
    pos1: int = 2,
    chr2: int = 3,
    pos2: int = 4,
    output_dist: int = 500_000_000,
    window_size: int = 500_000_000,
    step_size: int = 500_000_000,
    timeout_sec: int = 3600,
    dry_run: bool = False,
    logger: logging.Logger | None = None,
) -> bool:
    """Run hicluster prepare-impute for a stage."""
    logger = logger or logging.getLogger(__name__)

    cmd = [
        "hicluster", "prepare-impute",
        "--cell_table", contact_table,
        "--batch_size", str(batch_size),
        "--pad", str(pad),
        "--cpu_per_job", str(cpu_per_job),
        "--chr1", str(chr1),
        "--pos1", str(pos1),
        "--chr2", str(chr2),
        "--pos2", str(pos2),
        "--output_dir", output_dir,
        "--chrom_size_path", chrom_size_path,
        "--output_dist", str(output_dist),
        "--window_size", str(window_size),
        "--step_size", str(step_size),
        "--resolution", str(resolution),
    ]

    logger.info(f"[Step 1] prepare-impute @ {stage}")
    logger.info(f"Command: {' '.join(cmd)}")
    logger.info(f"Output directory: {output_dir}")
    logger.info(f"Contact table: {contact_table}")

    if dry_run:
        logger.info("Dry-run: skipping execution.")
        return True

    try:
        result = subprocess.run(
            cmd,
            cwd=os.path.dirname(output_dir) or None,  # None -> current
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            timeout=timeout_sec,
        )
        if result.returncode == 0:
            logger.info(f"prepare-impute OK @ {stage}")
            if result.stdout:
                logger.info(f"STDOUT:\n{result.stdout}")
            return True
        else:
            logger.error(f"prepare-impute FAILED @ {stage} (code {result.returncode})")
            if result.stderr:
                logger.error(f"STDERR:\n{result.stderr}")
            if result.stdout:
                logger.error(f"STDOUT:\n{result.stdout}")
            return False
    except subprocess.TimeoutExpired:
        logger.error(f"prepare-impute TIMEOUT @ {stage} ({timeout_sec}s)")
        return False
    except Exception as e:
        logger.error(f"prepare-impute EXCEPTION @ {stage}: {e}")
        return False

def run_snakemake_imputation(
    stage: str,
    output_dir: str,
    *,
    snakemake_cmd_file: str = "snakemake_cmd.txt",
    timeout_sec: int = 86_400,
    dry_run: bool = False,
    logger: logging.Logger | None = None,
) -> bool:
    """Run snakemake imputation using command recorded by hicluster."""
    logger = logger or logging.getLogger(__name__)

    cmd_path = os.path.join(output_dir, snakemake_cmd_file)
    if not os.path.exists(cmd_path):
        logger.error(f"Snakemake command file not found: {cmd_path}")
        return False

    with open(cmd_path, "r") as f:
        snakemake_cmd = f.read().strip()

    logger.info(f"[Step 2] snakemake imputation @ {stage}")
    logger.info(f"Command: {snakemake_cmd}")

    if dry_run:
        logger.info("Dry-run: skipping execution.")
        return True

    try:
        result = subprocess.run(
            snakemake_cmd.split(),
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            timeout=timeout_sec,
        )
        if result.returncode == 0:
            logger.info(f"snakemake imputation OK @ {stage}")
            if result.stdout:
                logger.info(f"STDOUT:\n{result.stdout}")
            return True
        else:
            logger.error(f"snakemake FAILED @ {stage} (code {result.returncode})")
            if result.stderr:
                logger.error(f"STDERR:\n{result.stderr}")
            if result.stdout:
                logger.error(f"STDOUT:\n{result.stdout}")
            return False
    except subprocess.TimeoutExpired:
        logger.error(f"snakemake TIMEOUT @ {stage} ({timeout_sec}s)")
        return False
    except Exception as e:
        logger.error(f"snakemake EXCEPTION @ {stage}: {e}")
        return False

def run_full_pipeline(
    stage: str,
    contact_table: str,
    output_dir: str,
    chrom_size_path: str,
    *,
    resolution: int,
    batch_size: int,
    cpu_per_job: int,
    pad: int,
    chr1: int,
    pos1: int,
    chr2: int,
    pos2: int,
    output_dist: int,
    window_size: int,
    step_size: int,
    prep_timeout: int,
    snakemake_timeout: int,
    snakemake_cmd_file: str,
    dry_run: bool,
    logger: logging.Logger,
) -> bool:
    logger.info(f"Start pipeline @ {stage}")

    ok = run_schicluster_prepare(
        stage, contact_table, output_dir, chrom_size_path,
        resolution=resolution,
        batch_size=batch_size,
        cpu_per_job=cpu_per_job,
        pad=pad,
        chr1=chr1, pos1=pos1, chr2=chr2, pos2=pos2,
        output_dist=output_dist, window_size=window_size, step_size=step_size,
        timeout_sec=prep_timeout,
        dry_run=dry_run,
        logger=logger,
    )
    if not ok:
        logger.error(f"Abort: prepare failed @ {stage}")
        return False

    ok = run_snakemake_imputation(
        stage, output_dir,
        snakemake_cmd_file=snakemake_cmd_file,
        timeout_sec=snakemake_timeout,
        dry_run=dry_run,
        logger=logger,
    )
    if not ok:
        logger.error(f"Abort: snakemake failed @ {stage}")
        return False

    logger.info(f"Pipeline OK @ {stage}")
    return True

# ---------------------------
# Utils
# ---------------------------
def get_cell_count(contact_table: str) -> int:
    try:
        with open(contact_table, "r") as f:
            return sum(1 for line in f if line.strip())
    except Exception:
        return 0

# ---------------------------
# CLI
# ---------------------------
def main():
    p = argparse.ArgumentParser(
        description="Process Hi-C data using schicluster by developmental stage"
    )
    # Workflow selection
    p.add_argument("--stages", nargs="+",
                   default=["E70", "E75", "E80", "E85", "E95", "EX05", "EX15"],
                   help="Stages to process")
    p.add_argument("--specific_stage", type=str, default=None,
                   help="If set, only process this stage")

    # Params
    p.add_argument("--resolution", type=int, default=100_000,
                   help="Resolution (e.g., 100000 for 100K)")
    p.add_argument("--batch_size", type=int, default=1536,
                   help="Batch size for prepare-impute")
    p.add_argument("--cpu_per_job", type=int, default=30,
                   help="CPUs per job for prepare-impute")
    p.add_argument("--pad", type=int, default=1, help="pad")
    p.add_argument("--chr1", type=int, default=1, help="Column index for chr1")
    p.add_argument("--pos1", type=int, default=2, help="Column index for pos1")
    p.add_argument("--chr2", type=int, default=3, help="Column index for chr2")
    p.add_argument("--pos2", type=int, default=4, help="Column index for pos2")
    p.add_argument("--output_dist", type=int, default=500_000_000)
    p.add_argument("--window_size", type=int, default=500_000_000)
    p.add_argument("--step_size", type=int, default=500_000_000)

    # Timeouts
    p.add_argument("--prep_timeout", type=int, default=3600,
                   help="Timeout (s) for prepare step")
    p.add_argument("--snakemake_timeout", type=int, default=86_400,
                   help="Timeout (s) for snakemake step")
    p.add_argument("--snakemake_cmd_file", type=str, default="snakemake_cmd.txt",
                   help="Filename containing snakemake command")

    # Paths
    path_group = p.add_mutually_exclusive_group()
    path_group.add_argument("--base_dir", type=str, default=None,
                            help="Absolute base dir (overrides relative mode)")
    path_group.add_argument("--relative_base", action="store_true",
                            help="Use script-relative base dir (default)")

    p.add_argument("--contact_table_dir", type=str, default=None,
                   help="Override contact_table dir (else base_dir/contact_tables)")
    p.add_argument("--output_base_dir", type=str, default=None,
                   help="Override outputs dir (else base_dir/outputs)")
    p.add_argument("--chrom_size_path", type=str, default=None,
                   help="Override chrom sizes file (else base_dir/mm10_chrom_sizes_with_chrY.txt)")
    p.add_argument("--log_dir", type=str, default=None,
                   help="Override log dir (else base_dir/logs)")

    # Misc
    p.add_argument("--dry_run", action="store_true",
                   help="Print commands only; do not execute")

    args = p.parse_args()

    # Resolve base_dir
    if args.base_dir:
        base_dir = args.base_dir
    else:
        # default: relative to this script
        script_dir = os.path.dirname(os.path.abspath(__file__))
        base_dir = os.path.dirname(script_dir)

    # Resolve subdirs/files
    contact_table_dir = args.contact_table_dir or os.path.join(base_dir, "contact_tables")
    output_base_dir = args.output_base_dir or os.path.join(base_dir, "outputs")
    chrom_size_path = args.chrom_size_path or os.path.join(base_dir, "mm10_chrom_sizes_with_chrY.txt")
    log_dir = args.log_dir or os.path.join(base_dir, "logs")

    logger = setup_logging(log_dir)

    # Input checks
    if not os.path.exists(chrom_size_path):
        logger.error(f"Chromosome size file not found: {chrom_size_path}")
        return

    stages = [args.specific_stage] if args.specific_stage else args.stages
    logger.info(f"Stages: {stages}")
    logger.info(f"Params: resolution={args.resolution}, batch_size={args.batch_size}, cpu_per_job={args.cpu_per_job}")
    logger.info(f"Base: {base_dir}")
    logger.info(f"contact_table_dir: {contact_table_dir}")
    logger.info(f"output_base_dir: {output_base_dir}")
    logger.info(f"chrom_size_path: {chrom_size_path}")
    logger.info(f"logs: {log_dir}")
    if args.dry_run:
        logger.info("DRY-RUN mode is ON")

    results: dict[str, str] = {}

    for stage in stages:
        logger.info("\n" + "=" * 60)
        logger.info(f"Processing stage: {stage}")
        logger.info("=" * 60)

        contact_table = os.path.join(contact_table_dir, f"contact_table_{stage}.tsv")
        if not os.path.exists(contact_table):
            logger.error(f"Contact table not found: {contact_table}")
            results[stage] = "Contact table missing"
            continue

        # Output dir: .../<stage>/impute/<resolution K>
        resK = f"{args.resolution // 1000}K"
        output_dir = os.path.join(output_base_dir, stage, "impute", resK)
        os.makedirs(output_dir, exist_ok=True)

        cell_count = get_cell_count(contact_table)
        logger.info(f"Cell count (rough): {cell_count}")

        ok = run_full_pipeline(
            stage=stage,
            contact_table=contact_table,
            output_dir=output_dir,
            chrom_size_path=chrom_size_path,
            resolution=args.resolution,
            batch_size=args.batch_size,
            cpu_per_job=args.cpu_per_job,
            pad=args.pad,
            chr1=args.chr1, pos1=args.pos1, chr2=args.chr2, pos2=args.pos2,
            output_dist=args.output_dist,
            window_size=args.window_size,
            step_size=args.step_size,
            prep_timeout=args.prep_timeout,
            snakemake_timeout=args.snakemake_timeout,
            snakemake_cmd_file=args.snakemake_cmd_file,
            dry_run=args.dry_run,
            logger=logger,
        )
        results[stage] = "Success" if ok else "Failed"
        logger.info(f"Result @ {stage}: {results[stage]}")

    # Summary
    logger.info("\n" + "=" * 60)
    logger.info("PROCESSING SUMMARY")
    logger.info("=" * 60)
    for s, r in results.items():
        logger.info(f"{s}: {r}")
    ok_n = sum(1 for r in results.values() if r == "Success")
    logger.info(f"\nSuccessfully processed {ok_n}/{len(results)} stages")

if __name__ == "__main__":
    main()