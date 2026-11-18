"""
Compare the pipeline's placements with the simulated truth.

Usage: python tests/check_results.py SIMDIR
For every result folder results/<Link>_<sb_index>/<section>/ :
  * normal nuclei: most placed, placement error small, high score
  * nuclei with ~no spatial UMIs: not placed
  * doublet-like nuclei (two clusters): may be placed at their dominant cluster (as in the R
    pipeline), but get a clearly lower score
  * droplets shared by a nucleus of each section: both nuclei placed on their own section
  * droplet shared by two nuclei of the same section: left out
"""
import glob
import json
import os
import sys

import numpy as np
import pandas as pd


def main(sim):
    truth = pd.read_csv(os.path.join(sim, "truth.csv"))
    sheet = pd.read_csv(os.path.join(sim, "sample_sheet.csv"))
    ok = True
    for _, line in sheet.dropna(subset=["sb_index"]).iterrows():
        for res in sorted(glob.glob(os.path.join(sim, "results", f"{line.Link}_{line.sb_index}", "*", ""))):
            sec = os.path.basename(os.path.dirname(res))
            t = truth[(truth.gem_well == line.Link) & (truth.section_id == sec)]  # noqa
            c = pd.read_csv(os.path.join(res, "coords.csv")).drop(columns="gem_barcode")
            m = t.merge(c, on="cellranger_barcode", how="left", suffixes=("_true", ""))
            placed = m["placed"].fillna(False).astype(bool)
            err = np.hypot(m["x"] - m["x_true"], m["y"] - m["y_true"])
            normal, doublet, low = (m["kind"] == k for k in ("normal", "doublet", "low"))
            other = truth.loc[(truth.gem_well == line.Link) & (truth.section_id != sec), "gem_barcode"]
            shared_cross = m["gem_barcode"].isin(other)
            shared_same = m["gem_barcode"].duplicated(keep=False)
            log = json.load(open(os.path.join(res, "log.json")))
            r = dict(
                normal_placed=f"{(placed & normal).sum()}/{normal.sum()}",
                median_error=round(float(err[placed & normal].median()), 2),
                max_error=round(float(err[placed & normal].max()), 2),
                doublet_score=round(float(m.loc[doublet, "score"].median()), 2),
                normal_score=round(float(m.loc[normal, "score"].median()), 2),
                low_placed=int((placed & low).sum()),
                cross_shared_ok=f"{int((placed & shared_cross & (err < 30)).sum())}/{int(shared_cross.sum())}",
                same_shared_left_out=bool(m.loc[shared_same, "placed"].isna().all()),
                tag=log["tag"],
                reads_wrong_puck=log["STEP 5 reads (this section)"]["removed: bead not on this section's puck"],
            )
            checks = [(placed & normal).sum() / normal.sum() > 0.95, r["median_error"] < 10, r["max_error"] < 40,
                      r["doublet_score"] < 0.8, r["normal_score"] > 0.95, (m.loc[normal, "score"] < 0.5).sum() == 0,
                      r["low_placed"] == 0, r["same_shared_left_out"],
                      int((placed & shared_cross & (err < 30)).sum()) == int(shared_cross.sum())]
            ok &= all(checks)
            print(f"{'PASS' if all(checks) else 'FAIL'} {line.Link}_{line.sb_index} {sec}: {r}")
    ok &= check_summaries(sim, truth, sheet)
    print("ALL CHECKS PASSED" if ok else "SOME CHECKS FAILED")
    return ok


def check_summaries(sim, truth, sheet):
    """results/summary/<section>/: pooled positions vs truth, RNA-only GEM well, files.csv, PDF."""
    ok = True
    for sec, t in truth.groupby("section_id"):
        folder = os.path.join(sim, "results", "summary", sec)
        c = pd.read_csv(os.path.join(folder, "coords.csv"))
        files = pd.read_csv(os.path.join(folder, "files.csv"))
        m = t.merge(c, on=["cellranger_barcode"], how="left", suffixes=("_true", ""))
        normal = m["kind"] == "normal"
        err = np.hypot(m["x"] - m["x_true"], m["y"] - m["y_true"])
        rna_only = m["kind"] == "rna_only"
        pdf = os.path.join(folder, f"{sec}_summary.pdf")
        r = dict(nuclei=len(c), placed=int(c["placed"].sum()), median_error=round(float(err[normal & m["placed"]].median()), 2),
                 rna_only_marked=bool((m.loc[rna_only, "spatial_status"] == "RNA only (no spatial library)").all()),
                 files_rows=len(files), libraries=sorted(files["spatial_library"].fillna("").tolist()),
                 pdf_kb=os.path.getsize(pdf) // 1024)
        checks = [len(c) == len(t), r["median_error"] < 10, r["rna_only_marked"], rna_only.sum() > 0,
                  len(files) == len(sheet), (m.loc[normal, "placed"]).mean() > 0.9, r["pdf_kb"] > 100]
        ok &= all(checks)
        print(f"{'PASS' if all(checks) else 'FAIL'} summary {sec}: {r}")
    return ok


if __name__ == "__main__":
    sys.exit(0 if main(sys.argv[1]) else 1)
