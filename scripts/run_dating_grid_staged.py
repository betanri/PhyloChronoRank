#!/usr/bin/env python3
import argparse
import csv
import shutil
import subprocess
from pathlib import Path


SCRIPT_DIR = Path(__file__).resolve().parent
RUNNER = SCRIPT_DIR / "run_dating_grid.R"


def run_cmd(args, log_path: Path):
    log_path.parent.mkdir(parents=True, exist_ok=True)
    with log_path.open("w") as handle:
        proc = subprocess.Popen(args, stdout=handle, stderr=subprocess.STDOUT)
        rc = proc.wait()
    if rc != 0:
        raise RuntimeError(f"command failed ({rc}): {' '.join(args)}; see {log_path}")


def read_candidates(outdir: Path):
    path = outdir / "candidates.csv"
    if not path.exists():
        return []
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def choose_first(rows, pred):
    for row in rows:
        if pred(row):
            return row
    return None


def choose_chronos(rows, model):
    return choose_first(rows, lambda r: r.get("method") == "chronos" and r.get("model") == model)


def copy_if_exists(src: Path, dst: Path):
    if src.exists():
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(src, dst)
        return True
    return False


def rm_glob(directory: Path, pattern: str):
    if not directory.exists():
        return
    for path in directory.glob(pattern):
        try:
            path.unlink()
        except FileNotFoundError:
            pass


def resolve_tree_source(source_outdir: Path, prefix: str, row: dict):
    candidate = row["candidate"]
    method = row["method"]
    model = row.get("model", "")
    if method == "treePL":
        ci = source_outdir / "selected_ci" / f"{prefix}_{candidate}_with_CI.tre"
        if ci.exists():
            return ci
        plain_full = source_outdir / "treepl" / "trees" / f"{prefix}_{candidate}_fulltree.tre"
        if plain_full.exists():
            return plain_full
        plain = source_outdir / "treepl" / "trees" / f"{prefix}_{candidate}.tre"
        if plain.exists():
            return plain
    elif method == "chronos":
        if model != "relaxed":
            ci = source_outdir / "selected_ci" / f"{prefix}_{candidate}_with_CI.tre"
            if ci.exists():
                return ci
        plain = source_outdir / "chronos" / "trees" / f"{prefix}_{candidate}.tre"
        if plain.exists():
            return plain
    elif method == "RelTime_MEGA":
        ci = source_outdir / "selected_ci" / f"{prefix}_RelTime_MEGA_with_bootstrap_CI.tre"
        if ci.exists():
            return ci
        plain = source_outdir / "reltime_mega" / f"{prefix}_RelTime_MEGA.tre"
        if plain.exists():
            return plain
    fallback = Path(row.get("tree_file", ""))
    return fallback if fallback.exists() else None


def sync_main(main_outdir: Path, prefix: str, treepl_source: Path | None = None, chronos_source: Path | None = None):
    main_dir = main_outdir / "MAIN_OUTPUT_TREES"
    ci_dir = main_dir / "ci"
    main_dir.mkdir(parents=True, exist_ok=True)
    ci_dir.mkdir(parents=True, exist_ok=True)

    copy_if_exists(main_outdir / "selected_ci" / f"{prefix}_RelTime_MEGA_with_bootstrap_CI.tre",
                   main_dir / "RelTime_MEGA_with_bootstrap_CI.tre")
    copy_if_exists(main_outdir / "selected_ci" / f"{prefix}_RelTime_MEGA_bootstrap_ci.csv",
                   ci_dir / "RelTime_MEGA_bootstrap_ci.csv")
    copy_if_exists(main_outdir / "reltime_mega" / f"{prefix}_RelTime_MEGA_with_tao_CI.tre",
                   main_dir / "RelTime_MEGA_with_tao_CI.tre")
    copy_if_exists(main_outdir / "reltime_mega" / f"{prefix}_RelTime_MEGA_tao_ci.csv",
                   ci_dir / "RelTime_MEGA_tao_ci.csv")

    main_candidates = read_candidates(main_outdir)

    treepl_src_dir = treepl_source or main_outdir
    treepl_rows = read_candidates(treepl_src_dir)
    treepl_row = choose_first(treepl_rows, lambda r: r.get("method") == "treePL")
    rm_glob(main_dir, "treePL*.tre")
    rm_glob(ci_dir, "treePL*_ci.csv")
    if treepl_row:
        src = resolve_tree_source(treepl_src_dir, prefix, treepl_row)
        if src:
            dest_name = f"{treepl_row['candidate']}_with_CI.tre" if src.name.endswith("_with_CI.tre") else f"{treepl_row['candidate']}.tre"
            copy_if_exists(src, main_dir / dest_name)
        copy_if_exists(treepl_src_dir / "selected_ci" / f"{prefix}_{treepl_row['candidate']}_ci.csv",
                       ci_dir / f"{treepl_row['candidate']}_ci.csv")

    chronos_src_dir = chronos_source or main_outdir
    chronos_rows = read_candidates(chronos_src_dir)
    models = ["clock", "correlated", "discrete", "relaxed"]
    selected = {m: choose_chronos(chronos_rows, m) for m in models}
    for model in models:
        rm_glob(main_dir, f"chronos_{model}_*.tre")
        rm_glob(ci_dir, f"chronos_{model}_*_ci.csv")
        row = selected[model]
        if not row:
            continue
        src = resolve_tree_source(chronos_src_dir, prefix, row)
        if not src:
            continue
        dest_name = f"{row['candidate']}_with_CI.tre" if src.name.endswith("_with_CI.tre") else f"{row['candidate']}.tre"
        copy_if_exists(src, main_dir / dest_name)
        copy_if_exists(chronos_src_dir / "selected_ci" / f"{prefix}_{row['candidate']}_ci.csv",
                       ci_dir / f"{row['candidate']}_ci.csv")

    main_candidates = [r for r in main_candidates if r.get("method") == "RelTime_MEGA"]
    if treepl_row:
        row = dict(treepl_row)
        row["tree_file"] = str((main_dir / (f"{treepl_row['candidate']}_with_CI.tre" if (main_dir / f"{treepl_row['candidate']}_with_CI.tre").exists() else f"{treepl_row['candidate']}.tre")).resolve())
        main_candidates.append(row)
    for model in models:
        row = selected[model]
        if not row:
            continue
        row = dict(row)
        dst = main_dir / (f"{row['candidate']}_with_CI.tre" if (main_dir / f"{row['candidate']}_with_CI.tre").exists() else f"{row['candidate']}.tre")
        row["tree_file"] = str(dst.resolve())
        main_candidates.append(row)
    for row in main_candidates:
        if row.get("method") == "RelTime_MEGA":
            row["tree_file"] = str((main_dir / "RelTime_MEGA_with_bootstrap_CI.tre").resolve())
    if main_candidates:
        fieldnames = ["candidate", "tree_file", "method", "model", "lambda", "nb_rate_cat", "smooth", "selected_from_candidate", "selection_score", "selection_rule"]
        with (main_outdir / "candidates.csv").open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fieldnames)
            writer.writeheader()
            for row in main_candidates:
                writer.writerow({k: row.get(k, "") for k in fieldnames})


def get_treepl_seed(main_outdir: Path, prefix: str):
    for row in read_candidates(main_outdir):
        if row.get("method") == "treePL":
            cand = row["candidate"]
            full = main_outdir / "treepl" / "trees" / f"{prefix}_{cand}_fulltree.tre"
            plain = main_outdir / "treepl" / "trees" / f"{prefix}_{cand}.tre"
            if full.exists():
                return full
            if plain.exists():
                return plain
    raise RuntimeError(f"treePL seed not found in {main_outdir}")


def build_args(args, outdir: Path, methods: str, chronos_ci_reps: int,
               treepl_bootstrap_reps: int, reltime_bootstrap_reps: int,
               chronos_age_start_tree: Path | None = None):
    cmd = [
        "Rscript", str(RUNNER),
        f"--phylogram={args.phylogram}",
        f"--calibrations-csv={args.calibrations_csv}",
        f"--outdir={outdir}",
        f"--out-prefix={args.out_prefix}",
        f"--methods={methods}",
        f"--chronos-ci-models={args.chronos_ci_models}",
        f"--chronos-ci-reps={chronos_ci_reps}",
        f"--treepl-bootstrap-reps={treepl_bootstrap_reps}",
        "--treepl-bootstrap-jobs=1",
        f"--reltime-bootstrap-reps={reltime_bootstrap_reps}",
    ]
    if "chronos" in methods and chronos_age_start_tree is not None:
        cmd.append(f"--chronos-age-start-tree={chronos_age_start_tree}")
    return cmd


def parse_args():
    p = argparse.ArgumentParser(description="Run staged dating workflow with early MAIN_OUTPUT_TREES syncing.")
    p.add_argument("--phylogram", required=True)
    p.add_argument("--calibrations-csv", required=True)
    p.add_argument("--outdir", required=True)
    p.add_argument("--out-prefix", required=True)
    p.add_argument("--chronos-ci-models", default="clock,correlated,discrete")
    p.add_argument("--reltime-bootstrap-reps", type=int, default=100)
    p.add_argument("--treepl-ci-reps", type=int, default=100)
    p.add_argument("--chronos-ci-reps", type=int, default=100)
    p.add_argument("--tmp-root", default="/tmp")
    return p.parse_args()


def main():
    args = parse_args()
    args.phylogram = Path(args.phylogram).resolve()
    args.calibrations_csv = Path(args.calibrations_csv).resolve()
    outdir = Path(args.outdir).resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    run_cmd(
        build_args(args, outdir, "reltime_mega", 0, 0, args.reltime_bootstrap_reps),
        outdir / "stage1_reltime_mega.log"
    )
    sync_main(outdir, args.out_prefix)

    run_cmd(
        build_args(args, outdir, "treepl", 0, 0, 0),
        outdir / "stage2_treepl_plain.log"
    )
    sync_main(outdir, args.out_prefix)

    seed = get_treepl_seed(outdir, args.out_prefix)
    run_cmd(
        build_args(args, outdir, "chronos", 0, 0, 0, chronos_age_start_tree=seed),
        outdir / "stage3_chronos_plain.log"
    )
    sync_main(outdir, args.out_prefix)

    treepl_ci_out = Path(args.tmp_root) / f"{args.out_prefix}_treepl_ci"
    if treepl_ci_out.exists():
        shutil.rmtree(treepl_ci_out)
    run_cmd(
        build_args(args, treepl_ci_out, "treepl", 0, args.treepl_ci_reps, 0),
        outdir / "stage4_treepl_ci.log"
    )
    sync_main(outdir, args.out_prefix, treepl_source=treepl_ci_out)

    chronos_ci_out = Path(args.tmp_root) / f"{args.out_prefix}_chronos_ci"
    if chronos_ci_out.exists():
        shutil.rmtree(chronos_ci_out)
    run_cmd(
        build_args(args, chronos_ci_out, "chronos", args.chronos_ci_reps, 0, 0, chronos_age_start_tree=seed),
        outdir / "stage5_chronos_ci.log"
    )
    sync_main(outdir, args.out_prefix, chronos_source=chronos_ci_out)


if __name__ == "__main__":
    main()
