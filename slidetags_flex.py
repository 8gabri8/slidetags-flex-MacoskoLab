#!/usr/bin/env python3
"""
Slide-tags for 10x Flex: give every nucleus that Cell Ranger called a position on its puck.
A Python version of the MacoskoLab slide-tags pipeline (spatial-count.jl + positioning.R).

    python3 slidetags_flex.py SAMPLE_SHEET.csv            every line of the sheet
    python3 slidetags_flex.py SAMPLE_SHEET.csv --line 2   only line 2

The sample sheet is the one you run Cumulus (Cell Ranger) with, plus 4 columns:
    sb_index               FASTQ prefix of the spatial library, e.g. SI-TT-A6 (empty = skip the line)
    pucks                  each section's puck: brain_section1=gs://.../Puck.csv;brain_section2=gs://...
    cellranger_output_dir  Cumulus output_directory (Cell Ranger results are in <this>/<Link>)
    slidetags_output_dir   results are written to <this>/<Link>_<sb_index>/
It also uses Cumulus' Link, Flowcell (folder with the FASTQs, spatial ones included) and AuxFile
(its sample_id column lists the sections).

Every line is processed on its own:
    STEP 1  download the AuxFile, the pucks, the spatial FASTQs
    STEP 2  load the pucks
    STEP 3  read the spatial FASTQs, keep the good reads, count them per (cell barcode, UMI, bead)
    for each section:
    STEP 4  load the nuclei Cell Ranger called in this section
    STEP 5  give each nucleus its spatial UMIs
    STEP 6  place each nucleus (DBSCAN)
    STEP 7  save the results
"""
import argparse
import csv
import glob
import gzip
import itertools
import json
import os
import shutil
import subprocess
from collections import Counter
from multiprocessing import Pool

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy.spatial import cKDTree  # noqa: E402

CELL_BARCODE_LENGTH = 16   # read 1: 16 bases of 10x GEM barcode ...
UMI_LENGTH = 12            # ... then 12 bases of UMI

# Read 2 holds the bead's spatial barcode in two parts plus a constant "UP" sequence.
# Where they are depends on the bead type (0-based positions).
UP_LONG, UP_SHORT = "TCTTCAGCGTTCCCGAGA", "CTGTTTCCTG"
BEAD_TYPES = {
    "V10": dict(part1=0, part2=26, up=(8, 26), up_seq=UP_LONG, min_length=33),
    "V17": dict(part1=0, part2=27, up=(9, 27), up_seq=UP_LONG, min_length=35),
    "V19": dict(part1=0, part2=19, up=(9, 19), up_seq=UP_LONG[:10], min_length=27),
    "V15": dict(part1=0, part2=8, up=(15, 25), up_seq=UP_SHORT, min_length=25),
    "V16": dict(part1=0, part2=9, up=(17, 27), up_seq=UP_SHORT, min_length=27),
}
PART1_LENGTH = {"V10": 8, "V17": 9, "V19": 9, "V15": 8, "V16": 9}


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("sample_sheet")
    parser.add_argument("--line", type=int, help="process only this line of the sheet (1 = first line after the header)")
    parser.add_argument("--count-lines", action="store_true", help="print the line numbers to process and stop (for the WDL)")
    parser.add_argument("--barcode-length", default="auto", help="Cell Ranger barcode length without -1: auto, 16, 24, ...")
    parser.add_argument("--knn", type=int, default=36, help="eps = median distance from a bead to its knn-th nearest bead")
    parser.add_argument("--cmes", type=float, default=10.0, help="clusters closer than eps*cmes count as one (0 = never)")
    parser.add_argument("--cores", type=int, default=os.cpu_count())
    parser.add_argument("--workdir", default="work")
    args = parser.parse_args()

    with open(args.sample_sheet, newline="") as f:
        lines = [{k.strip(): (v or "").strip() for k, v in row.items()} for row in csv.DictReader(f)]
    to_do = [n for n, line in enumerate(lines, start=1) if line.get("sb_index") and args.line in (None, n)]
    if args.count_lines:
        print("\n".join(map(str, to_do)))
        return
    for n in to_do:
        run_line(lines[n - 1], args)


def run_line(line, args):
    name = f"{line['Link']}_{line['sb_index']}"
    work = os.path.join(args.workdir, name)
    output_dir = f"{line['slidetags_output_dir'].rstrip('/')}/{name}"
    print(f"\n======== {name} ========")

    # STEP 1 - download the inputs of this line
    print("STEP 1: downloading")
    sections = read_auxfile(download(line["AuxFile"], work)[0])
    pucks = dict(pair.strip().split("=", 1) for pair in line["pucks"].split(";") if pair.strip())
    missing = [s for s in sections if s not in pucks]
    if missing:
        raise SystemExit(f"{name}: no puck given for section(s) {missing} (sections come from the AuxFile)")
    puck_files = {s: download(pucks[s], f"{work}/pucks/{s}")[0] for s in sections}  # own folder: all are Puck.csv
    fastq = line["Flowcell"].rstrip("/")
    read1_files = download(f"{fastq}/{line['sb_index']}_*_R1_001.fastq.gz", f"{work}/fastqs")
    read2_files = download(f"{fastq}/{line['sb_index']}_*_R2_001.fastq.gz", f"{work}/fastqs")

    # STEP 2 - load the pucks
    print("STEP 2: loading the pucks")
    beads = load_pucks(puck_files)

    # STEP 3 - read the spatial FASTQs, filter, count
    print("STEP 3: reading the spatial reads")
    counts, cell_barcodes, bead_barcodes, read_log = count_spatial_reads(list(zip(read1_files, read2_files)), beads)

    summary = []
    for section in sections:
        print(f"\n---- section {section} ----")
        folder = os.path.join(work, "results", section)
        os.makedirs(folder, exist_ok=True)

        # STEP 4 - the nuclei Cell Ranger called in this section
        print("STEP 4: loading the nuclei")
        barcodes_file = f"{line['cellranger_output_dir'].rstrip('/')}/{line['Link']}/per_sample_outs/{section}" \
                        f"/sample_filtered_feature_bc_matrix/barcodes.tsv.gz"
        nuclei, nuclei_log = load_nuclei(download(barcodes_file, f"{work}/cellranger/{section}")[0], args.barcode_length)

        # STEP 5 - their spatial UMIs, from this section's puck only
        print("STEP 5: giving each nucleus its spatial UMIs")
        puck = beads[beads["section"] == section]
        umis, umi_log = umis_per_nucleus(counts, cell_barcodes, bead_barcodes, puck, nuclei)

        # STEP 6 - their positions
        print("STEP 6: placing the nuclei")
        coords, place_log = place_nuclei(umis, nuclei, puck, args.knn, args.cmes, args.cores)

        # STEP 7 - save
        print("STEP 7: saving")
        log = {"line": name, "section": section, "STEP 3 reads": read_log, "STEP 4 nuclei": nuclei_log,
               "STEP 5 reads": umi_log, "STEP 6 placement": place_log}
        save_results(folder, coords, umis, nuclei, puck, log)
        upload(folder, output_dir)
        summary.append({"line": name, "section": section, **nuclei_log, **place_log})

    pd.DataFrame(summary).to_csv(f"{work}/results/summary.csv", index=False)
    upload(f"{work}/results/summary.csv", output_dir)
    print(pd.DataFrame(summary).to_string(index=False))


# =====================================================================================
# STEP 1 - downloading (gs:// with gcloud; plain paths are used as they are)
# =====================================================================================
def download(source, folder):
    """Copy a file, or all files matching a * pattern, into folder. Returns the local paths."""
    os.makedirs(folder, exist_ok=True)
    if source.startswith("gs://"):
        listing = subprocess.run(["gcloud", "storage", "ls", source], capture_output=True, text=True).stdout
        files = sorted(listing.split())
        if files:
            subprocess.run(["gcloud", "storage", "cp", *files, folder + "/"], check=True)
    else:
        files = sorted(glob.glob(source))
        for f in files:
            shutil.copy(f, folder)
    if not files:
        raise SystemExit(f"Not found: {source}")
    return [os.path.join(folder, os.path.basename(f)) for f in files]


def upload(local, destination):
    if destination.startswith("gs://"):
        subprocess.run(["gcloud", "storage", "cp", "-r", local, destination + "/"], check=True)
    elif os.path.isdir(local):
        shutil.copytree(local, os.path.join(destination, os.path.basename(local)), dirs_exist_ok=True)
    else:
        os.makedirs(destination, exist_ok=True)
        shutil.copy(local, destination)


def read_auxfile(path):
    """The sections = the sample_id column of the Flex AuxFile."""
    with open(path, newline="") as f:
        return [row["sample_id"].strip() for row in csv.DictReader(f) if row["sample_id"].strip()]


# =====================================================================================
# STEP 2 - pucks: Puck.csv = spatial barcode, x, y (no header), one line per bead
# =====================================================================================
def load_pucks(puck_files):
    """One table with the beads of all pucks: barcode, x, y, section.
    Removed (as in the R pipeline): barcodes with an N, and degenerate barcodes (one base > 4/5 of
    the barcode, or a run of one base > 5/7 of it). A barcode that is on two pucks is kept: we know
    each nucleus' section, so we know which puck its reads came from."""
    pucks = []
    for section, path in puck_files.items():
        puck = pd.read_csv(path, header=None, usecols=[0, 1, 2], names=["barcode", "x", "y"], dtype={0: str})
        puck["section"] = section
        pucks.append(puck)
    beads = pd.concat(pucks, ignore_index=True)
    bad = beads["barcode"].str.contains("N") | beads["barcode"].apply(is_degenerate)
    print(f"  {len(beads):,} beads on {len(pucks)} puck(s); removed {bad.sum()} with N or degenerate barcodes")
    return beads[~bad].reset_index(drop=True)


def is_degenerate(barcode):
    longest_run = max(len(list(run)) for _, run in itertools.groupby(barcode))
    most_common = max(barcode.count(base) for base in "ACGT")
    return most_common > round(len(barcode) * 4 / 5) or longest_run > round(len(barcode) * 5 / 7)


# =====================================================================================
# STEP 3 - spatial reads
# =====================================================================================
def count_spatial_reads(fastq_pairs, beads):
    """Read every read pair, keep the good ones, count reads per (cell barcode, UMI, bead).
    Returns counts (columns cell, umi, bead, reads: cell and bead are row numbers in the two
    returned barcode lists), the cell barcode list, the bead barcode list and the filter log."""

    # 3a. Which read holds the bead and which bead type: count UP-site hits in the first 100,000 reads
    bead_in_read1, bead_type = detect_layout(fastq_pairs)
    layout = BEAD_TYPES[bead_type]
    part1 = PART1_LENGTH[bead_type]
    part2 = len(beads["barcode"].iat[0]) - part1  # old 14-base V10 pucks: last base not used
    print(f"  bead type {bead_type}, bead in read {1 if bead_in_read1 else 2}")

    # 3b. What counts as a good UP site / bad UMI (1 mismatch allowed per 6 bases of UP)
    up_ok = within_mismatches(layout["up_seq"], round(len(layout["up_seq"]) / 6))
    up_all_g = within_mismatches("G" * len(layout["up_seq"]), round(len(layout["up_seq"]) / 6))
    umi_bad = set().union(*(within_mismatches(base * UMI_LENGTH, 2) for base in "ACGT"))
    bead_number = {barcode: i for i, barcode in enumerate(sorted(set(beads["barcode"])))}

    # 3c. Go through all reads
    log = Counter()
    cell_number = {}
    kept, tables = [], []
    for read1_file, read2_file in fastq_pairs:
        if bead_in_read1:
            read1_file, read2_file = read2_file, read1_file
        print(f"  reading {os.path.basename(read1_file)}", flush=True)
        for read1, read2 in zip(fastq_sequences(read1_file), fastq_sequences(read2_file)):
            log["reads"] += 1
            # filter 1: long enough?
            if len(read1) < CELL_BARCODE_LENGTH + UMI_LENGTH or len(read2) < layout["min_length"]:
                log["removed: too short"] += 1
                continue
            # filter 2: UMI without N and not one base repeated
            umi = read1[CELL_BARCODE_LENGTH:CELL_BARCODE_LENGTH + UMI_LENGTH]
            if "N" in umi or umi in umi_bad:
                log["removed: bad UMI"] += 1
                continue
            # filter 3: UP site present
            up = read2[layout["up"][0]:layout["up"][1]]
            if up not in up_ok:
                log["removed: all-G read" if up in up_all_g else "removed: no UP site"] += 1
                continue
            # filter 4: spatial barcode = a bead on the pucks (exact or 1 mismatch to one bead)
            spatial_barcode = read2[layout["part1"]:layout["part1"] + part1] + read2[layout["part2"]:layout["part2"] + part2]
            bead = find_bead(spatial_barcode, bead_number)
            if bead is None:
                log["removed: no bead / ambiguous bead"] += 1
                continue
            # keep it
            cell = cell_number.setdefault(read1[:CELL_BARCODE_LENGTH], len(cell_number))
            kept.append((cell, umi_to_number(umi), bead))
            if len(kept) == 5_000_000:  # add up duplicates now and then to save memory
                tables.append(pd.DataFrame(kept, columns=["cell", "umi", "bead"]).value_counts().rename("reads"))
                kept = []
    tables.append(pd.DataFrame(kept, columns=["cell", "umi", "bead"]).value_counts().rename("reads"))
    counts = pd.concat(tables).groupby(level=[0, 1, 2]).sum().reset_index()
    log["reads kept"] = int(counts["reads"].sum())
    print(f"  {log['reads']:,} reads, {log['reads kept']:,} kept")
    return counts, list(cell_number), list(bead_number), dict(log)


def detect_layout(fastq_pairs):
    def up_hits(path):
        hits = Counter()
        for seq in fastq_sequences(path, max_reads=100_000):
            for bead_type, layout in BEAD_TYPES.items():
                if len(seq) >= layout["min_length"] and seq[layout["up"][0]:layout["up"][1]] == layout["up_seq"]:
                    hits[bead_type] += 1
        hits["V19"] -= hits["V17"]  # V17 reads also contain V19's shorter UP site
        return hits
    read1_hits, read2_hits = up_hits(fastq_pairs[0][0]), up_hits(fastq_pairs[0][1])
    bead_in_read1 = max(read1_hits.values()) > max(read2_hits.values())
    hits = read1_hits if bead_in_read1 else read2_hits
    return bead_in_read1, max(BEAD_TYPES, key=lambda t: hits[t])


def find_bead(barcode, bead_number):
    """Exact match, or the only bead 1 mismatch away. None if no bead or 2+ beads are 1 mismatch away."""
    if barcode in bead_number:
        return bead_number[barcode]
    hits = {bead_number[b] for b in one_mismatch(barcode) if b in bead_number}
    return hits.pop() if len(hits) == 1 else None


def one_mismatch(seq):
    return [seq[:i] + base + seq[i + 1:] for i in range(len(seq)) for base in "ACGTN" if base != seq[i]]


def within_mismatches(seq, n):
    result = {seq}
    for _ in range(n):
        result |= {variant for s in result for variant in one_mismatch(s)}
    return result


def umi_to_number(umi):
    """12-base UMI as one number (2 bits per base, same numbering as spatial-count.jl)."""
    return int(umi.translate(str.maketrans("ACTG", "0123"))[::-1], 4)


def fastq_sequences(path, max_reads=None):
    """The sequence of every read in a .fastq.gz (line 2 of every 4)."""
    process = subprocess.Popen([shutil.which("pigz") or "gzip", "-dc", path], stdout=subprocess.PIPE, text=True)
    for n, row in enumerate(process.stdout):
        if n % 4 == 1:
            yield row.rstrip()
            if max_reads and n // 4 + 1 >= max_reads:
                break
    process.kill()
    process.wait()


# =====================================================================================
# STEP 4 - nuclei called by Cell Ranger in one section
# =====================================================================================
def load_nuclei(barcodes_file, barcode_length="auto"):
    """Cell Ranger barcode (e.g. 24 bases + -1 in multiplex Flex) -> GEM barcode (first 16 bases),
    which is what the spatial reads carry. Two nuclei of this section with the same GEM barcode
    (same droplet) cannot be told apart: both are left out."""
    with gzip.open(barcodes_file, "rt") as f:
        cellranger = [row.strip() for row in f if row.strip()]
    length = len(cellranger[0].split("-")[0])
    if barcode_length != "auto" and int(barcode_length) != length:
        raise SystemExit(f"Cell Ranger barcodes have {length} bases, not {barcode_length} ({cellranger[0]})")
    nuclei = pd.DataFrame({"cellranger_barcode": cellranger, "gem_barcode": [b[:CELL_BARCODE_LENGTH] for b in cellranger]})
    shared = nuclei["gem_barcode"].duplicated(keep=False)
    print(f"  {len(nuclei):,} nuclei, barcode length {length} (e.g. {cellranger[0]} -> {nuclei['gem_barcode'][0]}); "
          f"{shared.sum()} left out (same droplet)")
    return nuclei[~shared].reset_index(drop=True), {"nuclei called": len(nuclei), "barcode length": length,
                                                    "nuclei left out (same droplet)": int(shared.sum())}


# =====================================================================================
# STEP 5 - spatial UMIs of every nucleus of one section
# =====================================================================================
def umis_per_nucleus(counts, cell_barcodes, bead_barcodes, puck, nuclei):
    log = {"reads": int(counts["reads"].sum())}

    # 5a. keep only reads whose bead is on this section's puck (the rest belong to other sections)
    puck_barcodes = set(puck["barcode"])
    on_puck = np.array([b in puck_barcodes for b in bead_barcodes])
    reads = counts[on_puck[counts["bead"]]].copy()
    log["removed: bead on the other section's puck"] = log["reads"] - int(reads["reads"].sum())

    # 5b. match each read's cell barcode to a nucleus: exact, or the only nucleus 1 mismatch away.
    #     If swapping bases 8-9 (as 10x does between some library types) explains more reads, use that.
    nucleus_of = {b: i for i, b in enumerate(nuclei["gem_barcode"])}
    swap = {"TC": "AG", "AG": "TC", "GT": "CA", "CA": "GT"}
    swapped = {b[:7] + swap[b[7:9]] + b[9:]: i for b, i in nucleus_of.items() if b[7:9] in swap}
    reads_per_cell = reads.groupby("cell")["reads"].sum()
    barcodes_seen = [cell_barcodes[c] for c in reads_per_cell.index]
    exact_plain = sum(n for b, n in zip(barcodes_seen, reads_per_cell) if b in nucleus_of)
    exact_swapped = sum(n for b, n in zip(barcodes_seen, reads_per_cell) if b in swapped)
    if len(swapped) == len(nucleus_of) and exact_swapped > exact_plain:
        print("  using the 10x bases 8-9 swap for cell barcodes")
        nucleus_of = swapped
    nucleus = {}
    for cell, barcode in zip(reads_per_cell.index, barcodes_seen):
        if barcode in nucleus_of:
            nucleus[cell] = nucleus_of[barcode]
        else:
            hits = {nucleus_of[b] for b in one_mismatch(barcode) if b in nucleus_of}
            if len(hits) == 1:
                nucleus[cell] = hits.pop()
    reads["nucleus"] = reads["cell"].map(nucleus)

    # 5c. remove reads whose cell barcode is not a nucleus of this section (empty droplets, ambient)
    log["removed: not a nucleus of this section"] = int(reads.loc[reads["nucleus"].isna(), "reads"].sum())
    reads = reads.dropna(subset=["nucleus"]).astype({"nucleus": int})
    reads = reads.groupby(["nucleus", "umi", "bead"], as_index=False)["reads"].sum()

    # 5d. remove chimeras: a UMI of a nucleus comes from ONE bead. Keep the bead with strictly most reads.
    most = reads.groupby(["nucleus", "umi"])["reads"].transform("max")
    top = reads["reads"] == most
    single_top = top.groupby([reads["nucleus"], reads["umi"]]).transform("sum") == 1
    log["removed: chimeric UMI"] = int(reads.loc[~(top & single_top), "reads"].sum())
    reads = reads[top & single_top]
    log["reads used"] = int(reads["reads"].sum())

    # 5e. UMIs per (nucleus, bead), with the bead's position on this puck
    umis = reads.groupby(["nucleus", "bead"]).size().rename("umis").reset_index()
    umis["barcode"] = np.array(bead_barcodes)[umis["bead"]]
    umis = umis.merge(puck[["barcode", "x", "y"]], on="barcode")
    log["UMIs used"] = int(umis["umis"].sum())
    for what, n in log.items():
        print(f"  {what}: {n:,}")
    return umis, log


# =====================================================================================
# STEP 6 - place every nucleus with DBSCAN
#
# A nucleus collected UMIs mostly from the beads under it. DBSCAN looks for that dense spot:
# a bead is "core" when the UMIs within distance eps of it (itself included) add up to min_pts or
# more; core beads closer than eps form one cluster; other beads near a core bead join it.
#   min_pts too high -> no cluster (not placed)     min_pts too low -> 2+ clusters (ambiguous)
# Every nucleus has a range of min_pts giving exactly one cluster: (i2, i1]. One min_pts is chosen
# for the whole section: the one inside the most ranges. A nucleus is placed at the UMI-weighted
# centre of its cluster when it has exactly one cluster at that min_pts.
# =====================================================================================
def place_nuclei(umis, nuclei, puck, knn, cmes, cores):
    # 6a. eps = median distance from a bead to its knn-th nearest bead (10,000 beads sampled)
    xy = puck[["x", "y"]].to_numpy()
    sample = np.random.default_rng(1).choice(len(xy), min(10_000, len(xy)), replace=False)
    eps = float(np.median(cKDTree(xy).query(xy[sample], k=knn)[0][:, knn - 1]))

    # 6b. the one-cluster range (i2, i1] of every nucleus
    beads_of = {n: g for n, g in umis.groupby("nucleus")}
    nuclei_beads = [beads_of.get(n, umis.iloc[:0])[["x", "y", "umis"]].to_numpy() for n in range(len(nuclei))]
    with Pool(cores) as pool:
        ranges = pool.starmap(one_cluster_range, [(b, eps, cmes) for b in nuclei_beads])
        i1 = np.array([r[0] for r in ranges])
        i2 = np.array([r[1] for r in ranges])

        # 6c. one min_pts for the section: the value inside the most ranges (ties: the highest)
        values = np.arange(1, max(i1.max(initial=0), 1) + 1)
        nuclei_ok = np.array([np.sum((i2 < v) & (i1 >= v)) for v in values])
        min_pts = int(values[np.flatnonzero(nuclei_ok == nuclei_ok.max())[-1]])

        # 6d. final DBSCAN of every nucleus at that min_pts
        rows = pool.starmap(final_position, [(b, eps, cmes, min_pts, i2_n) for b, i2_n in zip(nuclei_beads, i2)])

    coords = pd.concat([nuclei, pd.DataFrame(rows)], axis=1)
    coords["score"] = placement_score(coords["umi_main"].to_numpy(float), coords["umi_second"].to_numpy(float))
    placed = int(coords["placed"].sum())
    print(f"  eps {eps:.1f}, min_pts {min_pts}: placed {placed:,} of {len(coords):,} nuclei")
    return coords, {"eps": round(eps, 3), "min_pts": min_pts, "nuclei positioned": len(coords),
                    "nuclei placed": placed, "% placed": round(100 * placed / max(len(coords), 1), 2)}


def dbscan(beads, eps, min_pts):
    """Cluster number of each bead (0 = background). Beads are sorted by UMIs, so cluster 1 is the
    one found first starting from the bead with most UMIs (as R's dbscan)."""
    order = np.argsort(-beads[:, 2], kind="stable")
    beads = beads[order]
    neighbours = cKDTree(beads[:, :2]).query_ball_point(beads[:, :2], r=eps)
    core = np.array([beads[nb, 2].sum() for nb in neighbours]) >= min_pts
    cluster = np.zeros(len(beads), int)
    for start in np.flatnonzero(core):
        if cluster[start]:
            continue
        cluster[start] = cluster.max() + 1
        to_visit = [start]
        while to_visit:
            for nb in neighbours[to_visit.pop()]:
                if cluster[nb] == 0:
                    cluster[nb] = cluster[start]
                    if core[nb]:
                        to_visit.append(nb)
    result = np.zeros(len(beads), int)
    result[order] = cluster  # back to the input order
    return result


def merge_close(beads, cluster, max_distance):
    """Clusters whose UMI-weighted centre is within max_distance of cluster 1's centre join cluster 1."""
    if cluster.max() < 2 or max_distance <= 0:
        return cluster
    centres = [np.average(beads[cluster == c, :2], axis=0, weights=beads[cluster == c, 2]) for c in range(1, cluster.max() + 1)]
    close = [c for c in range(2, cluster.max() + 1) if np.linalg.norm(centres[c - 1] - centres[0]) < max_distance]
    cluster = np.where(np.isin(cluster, close), 1, cluster)
    numbers = {old: new for new, old in enumerate(sorted(set(cluster[cluster > 0])), start=1)}
    return np.array([numbers.get(c, 0) for c in cluster])


def one_cluster_range(beads, eps, cmes):
    """i1 = highest min_pts with a cluster; i2 = highest min_pts with 2+ clusters (0 if none).
    Clusters only change at the beads' neighbourhood UMI sums, so only those values are tried.
    With cmes, i2 is lowered while all clusters at i2 lie within eps*cmes of cluster 1."""
    if len(beads) == 0:
        return 0, 0
    neighbours = cKDTree(beads[:, :2]).query_ball_point(beads[:, :2], r=eps)
    candidates = np.unique([beads[nb, 2].sum() for nb in neighbours])[::-1]  # high to low
    i1 = int(candidates[0])
    i2 = next((int(v) for v in candidates if dbscan(beads, eps, v).max() >= 2), 0)
    while cmes > 0 and i2 >= 1 and merge_close(beads, dbscan(beads, eps, i2), eps * cmes).max() < 2:
        lower = candidates[candidates < i2]
        i2 = int(lower[0]) if len(lower) else 0
    return i1, i2


def final_position(beads, eps, cmes, min_pts, i2):
    """Position at the section's min_pts + the two largest clusters at the nucleus' own i2 (for the score)."""
    row = dict(x=np.nan, y=np.nan, placed=False, clusters=0, umi=0, beads=0, spread=np.nan,
               umi_main=np.nan, umi_second=np.nan)
    if len(beads) == 0:
        return row
    cluster = merge_close(beads, dbscan(beads, eps, min_pts), eps * cmes)
    row.update(clusters=int(cluster.max()), umi=int(beads[:, 2].sum()), beads=len(beads))
    if row["clusters"] == 1:
        inside = beads[cluster == 1]
        centre = np.average(inside[:, :2], axis=0, weights=inside[:, 2])
        row.update(x=centre[0], y=centre[1], placed=True,
                   spread=float(np.sqrt(np.average(((inside[:, :2] - centre) ** 2).sum(1), weights=inside[:, 2]))))
    at_i2 = merge_close(beads, dbscan(beads, eps, i2), eps * cmes)
    sizes = sorted((beads[at_i2 == c, 2].sum() for c in range(1, at_i2.max() + 1)), reverse=True)
    row["umi_main"] = sizes[0] if sizes else np.nan
    row["umi_second"] = sizes[1] if len(sizes) > 1 else np.nan
    return row


def placement_score(umi_main, umi_second):
    """How much the main cluster dominates the second one (positioning.R):
    (1 - F_main(umi_second) / F_main(umi_main)) * F_second(umi_main), with F = cumulative
    distribution over all nuclei; no second cluster -> F_second(umi_main). 1 = clear single cluster."""
    def cdf(values):
        v = np.sort(values[~np.isnan(values)])
        return lambda q: np.where(np.isnan(q), np.nan, np.searchsorted(v, q, side="right") / max(len(v), 1)) \
            if len(v) else np.ones_like(q)
    F_main, F_second = cdf(umi_main), cdf(umi_second)
    with np.errstate(invalid="ignore", divide="ignore"):
        score = (1 - F_main(umi_second) / F_main(umi_main)) * F_second(umi_main)
    return np.where(np.isnan(umi_second), F_second(umi_main), score)


# =====================================================================================
# STEP 7 - save
# =====================================================================================
def save_results(folder, coords, umis, nuclei, puck, log):
    """coords.csv (one row per nucleus, join with RNA on cellranger_barcode), umis.csv.gz (UMIs per
    nucleus and bead), log.json (reads removed at every filter), summary.png (QC)."""
    coords.to_csv(f"{folder}/coords.csv", index=False)
    umis.assign(cellranger_barcode=nuclei["cellranger_barcode"].to_numpy()[umis["nucleus"]])[
        ["cellranger_barcode", "barcode", "umis", "x", "y"]].to_csv(f"{folder}/umis.csv.gz", index=False)
    with open(f"{folder}/log.json", "w") as f:
        json.dump(log, f, indent=2, default=int)

    fig, (a, b, c) = plt.subplots(1, 3, figsize=(18, 5.5))
    funnel = {**{k: v for k, v in log["STEP 3 reads"].items()}, **log["STEP 5 reads"]}
    a.barh(list(funnel)[::-1], list(funnel.values())[::-1])
    a.set_title("spatial reads at each filter")
    b.hist(np.log10(coords["umi"].clip(lower=1)), bins=40)
    b.set_xlabel("log10 spatial UMIs per nucleus")
    b.set_title(f"{log['STEP 6 placement']['% placed']}% of nuclei placed")
    c.scatter(puck["x"], puck["y"], s=0.1, c="0.85", rasterized=True)
    c.scatter(coords["x"], coords["y"], s=2, c=coords["score"], vmin=0, vmax=1)
    c.set_aspect("equal")
    c.set_title(f"{log['section']}: placed nuclei (colour = score)")
    fig.tight_layout()
    fig.savefig(f"{folder}/summary.png", dpi=120)
    plt.close(fig)


if __name__ == "__main__":
    main()
