"""
Two checks of STEP 6 of slidetags_flex.py:
  1. our DBSCAN (dbscan) gives the same clusters as scikit-learn's DBSCAN with sample weights
  2. our directly computed i1 / i2 (one_cluster_range) equal what the R search loops in
     positioning.R (ported line by line below) return, on simulated nuclei

Usage: python tests/test_dbscan.py [umis.csv.gz ...]
"""
import sys

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from sklearn.cluster import DBSCAN

sys.path.insert(0, ".")
import slidetags_flex as S  # noqa: E402


class Cell:
    """The beads of one nucleus, most UMIs first (what the R code works on)."""
    def __init__(self, x, y, umis, eps):
        order = np.argsort(-np.asarray(umis), kind="stable")
        self.x, self.y = np.asarray(x, float)[order], np.asarray(y, float)[order]
        self.umis = np.asarray(umis, int)[order]
        self.n, self.eps = len(self.umis), eps
        self.beads = np.column_stack([self.x, self.y, self.umis])
        nbs = cKDTree(self.beads[:, :2]).query_ball_point(self.beads[:, :2], r=eps) if self.n else []
        self.neighbourhood_umis = np.array([self.umis[nb].sum() for nb in nbs])


def dbscan(cell, min_pts):
    return S.dbscan(cell.beads, cell.eps, min_pts)


# --- literal port of positioning.R lines 248-334 ----------------------------------------
def r_iestimate(cell, idx):
    """sum of UMIs within eps of the max-UMI bead(s), searched within the subset idx"""
    pts = np.column_stack([cell.x[idx], cell.y[idx]])
    umi = cell.umis[idx]
    tree = cKDTree(pts)
    q = pts[umi == umi.max()]
    return max(umi[nb].sum() for nb in tree.query_ball_point(q, r=cell.eps))


def r_ipercluster(cell, v):
    out = []
    for c in sorted(set(v[v > 0])):
        idx = np.flatnonzero(v == c)
        pts = np.column_stack([cell.x[idx], cell.y[idx]])
        umi = cell.umis[idx]
        nbs = cKDTree(pts).query_ball_point(pts, r=cell.eps)
        out.append(max(umi[nb].sum() for nb in nbs))  # neighbours (excl. self) + self
    return out


def r_i1(cell):
    if cell.n == 0:
        return 0
    if cell.n == 1:
        return int(cell.umis.sum())
    i = r_iestimate(cell, np.arange(cell.n))
    v = dbscan(cell, i)
    i = max(r_ipercluster(cell, v))
    v = dbscan(cell, i)
    while v.max() > 0:
        i += 1
        v = dbscan(cell, i)
    return i - 1


def r_i2(cell):
    if cell.n <= 1:
        return 0
    i = r_iestimate(cell, np.arange(cell.n))
    v = dbscan(cell, i)
    seen = []
    while True:
        seen.append(i)
        if (v != 1).any():
            i = r_iestimate(cell, np.flatnonzero(v != 1))
            v = dbscan(cell, i)
        else:
            return 0
        if v.max() >= 2:
            i = sorted(r_ipercluster(cell, v))[-2]
            v = dbscan(cell, i)
        if v.max() >= 2 or i in seen:
            break
    if v.max() < 2:
        while v.max() < 2 and i >= 1:
            i -= 1
            v = dbscan(cell, i)
        return i
    while v.max() >= 2:
        i += 1
        v = dbscan(cell, i)
    return i - 1


# --- tests ---------------------------------------------------------------------------------
def same_partition(a, b):
    """labels a (0 = noise) vs b (-1 = noise): same noise set and same grouping"""
    if not np.array_equal(a == 0, b == -1):
        return False
    pairs = set(zip(a[a > 0], b[b >= 0]))
    return len(pairs) == len(set(a[a > 0])) == len(set(b[b >= 0]))


def main(matrices):
    rng = np.random.default_rng(0)
    cells = []
    for path in matrices:
        m = pd.read_csv(path)
        for _, g in m.groupby("cellranger_barcode"):
            cells.append((g["x"].to_numpy(), g["y"].to_numpy(), g["umis"].to_numpy()))
    for _ in range(300):  # plus random synthetic cells with 1-3 blobs
        k = rng.integers(1, 4)
        pts = np.concatenate([rng.normal(rng.uniform(0, 500, 2), 15, (rng.integers(1, 60), 2)) for _ in range(k)])
        cells.append((pts[:, 0], pts[:, 1], rng.integers(1, 6, len(pts))))
    eps = 42.0

    n_ok_sk = n_sk = 0
    agree_i1 = agree_i2 = 0
    diffs = []
    for x, y, umi in cells:
        cell = Cell(x, y, umi, eps)
        # 1. compare with scikit-learn at a few min_pts values
        for mp in sorted(set([1, 5, 20, int(cell.neighbourhood_umis.max()) if cell.n else 1])):
            ours = dbscan(cell, mp)
            sk = DBSCAN(eps=eps, min_samples=mp).fit(np.column_stack([cell.x, cell.y]),
                                                     sample_weight=cell.umis).labels_
            n_sk += 1
            n_ok_sk += same_partition(ours, sk) or _border_only_difference(cell, ours, sk, mp)
        # 2. compare exact boundaries with the R search
        a1, a2 = S.one_cluster_range(cell.beads, eps, 0)  # cmes = 0: the pure search
        b1, b2 = r_i1(cell), r_i2(cell)
        agree_i1 += a1 == b1
        agree_i2 += a2 == b2
        if a2 != b2:
            diffs.append((cell.n, a2, b2, dbscan(cell, a2).max() if a2 else None, dbscan(cell, b2).max()))
    n = len(cells)
    print(f"DBSCAN vs scikit-learn : {n_ok_sk}/{n_sk} identical clusterings")
    print(f"i1 exact vs R search   : {agree_i1}/{n} equal")
    print(f"i2 exact vs R search   : {agree_i2}/{n} equal")
    for d in diffs[:10]:
        print(f"   differs: beads={d[0]} exact i2={d[1]} (clusters {d[3]}) | R i2={d[2]} (clusters {d[4]})")
    return n_ok_sk == n_sk and agree_i1 == n


def _border_only_difference(cell, ours, sk, mp):
    """A border bead within eps of two clusters may join either one; that is not an error."""
    core = cell.neighbourhood_umis >= mp
    if not np.array_equal(ours[core] > 0, sk[core] >= 0):
        return False
    return same_partition(np.where(core, ours, 0), np.where(core, sk, -1))


if __name__ == "__main__":
    ok = main(sys.argv[1:])
    sys.exit(0 if ok else 1)
