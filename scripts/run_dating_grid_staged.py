#!/usr/bin/env python3
"""Staged PCR dating workflow.

Runs scripts/run_dating_grid.R five times in ONE output folder:

  1. RelTime (MEGA-CC) point estimate, Tao CI, and bootstrap CI
  2. treePL: prime -> cross-validation -> final dating (no CI yet)
  3. chronos grid for all four clock models (no CI yet), seeded by treePL
  4. treePL bootstrap CI for the stage-2 tree
  5. chronos bootstrap CI for the stage-3 winners

Each stage reuses the checkpointed trees of earlier stages, so the CI stages
add intervals to the same trees instead of refitting them. The R engine
refreshes MAIN_OUTPUT_TREES/ after every stage with the best version of each
candidate (CI tree when available, plain tree otherwise).
"""
import argparse
import csv
import subprocess
from pathlib import Path


SCRIPT_DIR = Path(__file__).resolve().parent
RUNNER = SCRIPT_DIR / "run_dating_grid.R"


def run_cmd(args, log_path: Path):
    log_path.parent.mkdir(parents=True, exist_ok=True)
    print(f"[stage] {log_path.name}", flush=True)
    with log_path.open("w") as handle:
        rc = subprocess.call(args, stdout=handle, stderr=subprocess.STDOUT)
    if rc != 0:
        raise RuntimeError(f"command failed ({rc}): {' '.join(map(str, args))}; see {log_path}")


def read_csv_rows(path: Path):
    if not path.exists():
        return []
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def load_tip_label_map(csv_path: Path | None):
    if csv_path is None:
        return {}
    rows = read_csv_rows(csv_path)
    if not rows:
        return {}
    first = rows[0]
    if "Old Name" in first and "NewName" in first:
        old_col, new_col = "Old Name", "NewName"
    else:
        keys = list(first.keys())
        if len(keys) < 2:
            return {}
        old_col, new_col = keys[0], keys[1]
    return {row[old_col]: row[new_col] for row in rows if row.get(old_col) and row.get(new_col)}


def relabel_main_trees(outdir: Path, tip_map: dict[str, str]):
    if not tip_map:
        return
    for path in (outdir / "MAIN_OUTPUT_TREES").glob("*.tre"):
        text = path.read_text()
        for old, new in tip_map.items():
            text = text.replace(old, new)
        path.write_text(text)


def get_chronos_seed(outdir: Path, prefix: str):
    """Seed chronos starting ages with the treePL tree, else the RelTime tree."""
    for row in read_csv_rows(outdir / "treepl" / f"{prefix}_treepl_runs.csv"):
        tree = row.get("tree_file", "")
        if tree and Path(tree).exists():
            return Path(tree), "treePL"
    reltime = outdir / "reltime_mega" / f"{prefix}_RelTime_MEGA.tre"
    if reltime.exists():
        return reltime, "RelTime_MEGA"
    raise RuntimeError(f"no chronos seed tree found in {outdir}")


def build_args(args, methods: str, chronos_ci_reps: int, treepl_bootstrap_reps: int,
               reltime_bootstrap_reps: int, chronos_age_start_tree: Path | None = None):
    cmd = [
        "Rscript", str(RUNNER),
        f"--phylogram={args.phylogram}",
        f"--outdir={args.outdir}",
        f"--out-prefix={args.out_prefix}",
        f"--methods={methods}",
        f"--chronos-ci-models={args.chronos_ci_models}",
        f"--chronos-ci-reps={chronos_ci_reps}",
        f"--chronos-allow-cal-drop={args.chronos_allow_cal_drop}",
        f"--treepl-bootstrap-reps={treepl_bootstrap_reps}",
        f"--treepl-bootstrap-jobs={args.treepl_bootstrap_jobs}",
        f"--treepl-smoothing={args.treepl_smoothing}",
        f"--treepl-cv={args.treepl_cv}",
        f"--treepl-numsites={args.treepl_numsites}",
        f"--reltime-bootstrap-reps={reltime_bootstrap_reps}",
        f"--ci-sites={args.ci_sites}",
        f"--seed={args.seed}",
    ]
    if args.calibrations_csv:
        cmd.append(f"--calibrations-csv={args.calibrations_csv}")
    else:
        cmd.append(f"--reference-time-tree={args.reference_time_tree}")
    for flag, value in (("additional-calibrations-csv", args.additional_calibrations_csv),
                        ("treepl-bin", args.treepl_bin), ("megacc-bin", args.megacc_bin),
                        ("treepl-opt", args.treepl_opt), ("root-age", args.root_age),
                        ("treepl-prime-reps", args.treepl_prime_reps),
                        ("treepl-cv-extend", args.treepl_cv_extend),
                        ("treepl-opt-check", args.treepl_opt_check),
                        ("treepl-cv-reps", args.treepl_cv_reps),
                        ("treepl-opt-screen", args.treepl_opt_screen),
                        ("treepl-opt-seeds", args.treepl_opt_seeds),
                        ("treepl-deep-seeds", args.treepl_deep_seeds),
                        ("treepl-screen-pliter", args.treepl_screen_pliter),
                        ("treepl-deep-pliter", args.treepl_deep_pliter),
                        ("treepl-jobs", args.treepl_jobs),
                        ("subset-tips", args.subset_tips), ("subset-min", args.subset_min),
                        ("subset-max", args.subset_max)):
        if value:
            cmd.append(f"--{flag}={value}")
    if chronos_age_start_tree is not None:
        cmd.append(f"--chronos-age-start-tree={chronos_age_start_tree}")
    return cmd


def parse_args():
    p = argparse.ArgumentParser(description="Run the staged PCR dating workflow in one output folder.")
    p.add_argument("--phylogram", required=True)
    cal = p.add_mutually_exclusive_group(required=True)
    cal.add_argument("--calibrations-csv")
    cal.add_argument("--reference-time-tree")
    p.add_argument("--additional-calibrations-csv", default="",
                   help="extra calibrations (e.g. fossils) merged with the primary set; "
                        "they win over congruified ages when the two disagree")
    p.add_argument("--outdir", required=True)
    p.add_argument("--out-prefix", required=True)
    p.add_argument("--root-age", default="")
    p.add_argument("--chronos-ci-models", default="clock,correlated,discrete")
    p.add_argument("--chronos-allow-cal-drop", default="FALSE")
    p.add_argument("--chronos-ci-reps", type=int, default=100)
    p.add_argument("--treepl-ci-reps", type=int, default=100)
    p.add_argument("--treepl-bootstrap-jobs", type=int, default=1)
    p.add_argument("--treepl-smoothing", default="0.01,0.1,1,10,100")
    p.add_argument("--treepl-cv", choices=["loo", "random"], default="loo")
    p.add_argument("--treepl-numsites", type=int, default=1000,
                   help="alignment length passed to treePL (numsites)")
    p.add_argument("--treepl-prime-reps", default="")
    p.add_argument("--treepl-cv-extend", default="")
    p.add_argument("--treepl-opt-check", default="")
    p.add_argument("--treepl-cv-reps", default="")
    p.add_argument("--treepl-opt-screen", choices=["", "prime", "full"], default="",
                   help="full adds A. Pyron's 144-setting optimizer screen and deep seed search")
    p.add_argument("--treepl-opt-seeds", default="")
    p.add_argument("--treepl-deep-seeds", default="")
    p.add_argument("--treepl-screen-pliter", default="")
    p.add_argument("--treepl-deep-pliter", default="")
    p.add_argument("--treepl-jobs", default="")
    p.add_argument("--subset-tips", default="")
    p.add_argument("--subset-min", default="")
    p.add_argument("--subset-max", default="")
    p.add_argument("--treepl-opt", default="")
    p.add_argument("--treepl-bin", default="")
    p.add_argument("--megacc-bin", default="")
    p.add_argument("--reltime-bootstrap-reps", type=int, default=100)
    p.add_argument("--ci-sites", type=int, default=1000)
    p.add_argument("--seed", type=int, default=1)
    p.add_argument("--tip-label-map-csv", default="")
    p.add_argument("--start-stage", type=int, default=1, choices=range(1, 6))
    return p.parse_args()


def main():
    args = parse_args()
    args.phylogram = Path(args.phylogram).resolve()
    if args.calibrations_csv:
        args.calibrations_csv = Path(args.calibrations_csv).resolve()
    else:
        args.reference_time_tree = Path(args.reference_time_tree).resolve()
    if args.additional_calibrations_csv:
        args.additional_calibrations_csv = Path(args.additional_calibrations_csv).resolve()
    args.outdir = Path(args.outdir).resolve()
    args.outdir.mkdir(parents=True, exist_ok=True)
    tip_map = load_tip_label_map(Path(args.tip_label_map_csv).resolve() if args.tip_label_map_csv else None)
    out = args.outdir

    if args.start_stage <= 1:
        run_cmd(build_args(args, "reltime_mega", 0, 0, args.reltime_bootstrap_reps),
                out / "stage1_reltime_mega.log")
        relabel_main_trees(out, tip_map)

    if args.start_stage <= 2:
        run_cmd(build_args(args, "treepl", 0, 0, 0), out / "stage2_treepl_plain.log")
        relabel_main_trees(out, tip_map)

    seed, seed_method = get_chronos_seed(out, args.out_prefix)
    print(f"Using {seed_method} tree to seed chronos: {seed}", flush=True)

    if args.start_stage <= 3:
        run_cmd(build_args(args, "chronos", 0, 0, 0, chronos_age_start_tree=seed),
                out / "stage3_chronos_plain.log")
        relabel_main_trees(out, tip_map)

    if args.start_stage <= 4:
        run_cmd(build_args(args, "treepl", 0, args.treepl_ci_reps, 0), out / "stage4_treepl_ci.log")
        relabel_main_trees(out, tip_map)

    if args.start_stage <= 5:
        run_cmd(build_args(args, "chronos", args.chronos_ci_reps, 0, 0, chronos_age_start_tree=seed),
                out / "stage5_chronos_ci.log")
        relabel_main_trees(out, tip_map)


if __name__ == "__main__":
    main()
