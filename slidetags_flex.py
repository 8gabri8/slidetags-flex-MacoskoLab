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
(the Flex sample table: sample_id = the sections; optional column spatial_tag, see below).

Section tag: the oligo that links a section's spatial beads to the 10x gel bead carries a tag, so
every spatial read says which section it comes from (read 2 of the spatial library). By default a
section's tag is its Flex probe barcode, taken from its Cell Ranger barcodes (GEM barcode + probe
barcode). If you used another sequence, put it in the AuxFile column spatial_tag.

Every line is processed on its own (--line N, one Terra job per line):
    STEP 1  download the AuxFile, the pucks, the Cell Ranger barcodes, the spatial FASTQs
    STEP 2  load the pucks
    STEP 3  load the nuclei Cell Ranger called, and each section's tag
    STEP 4  read the spatial FASTQs, keep the good reads, count them per (section, cell barcode, UMI, bead)
    for each section:
    STEP 5  give each nucleus its spatial UMIs
    STEP 6  place each nucleus (DBSCAN)
    STEP 7  save the results           -> <slidetags_output_dir>/<Link>_<sb_index>/

Then every section is summarised on its own (--summarize SECTION, one Terra job per section),
over all lines (GEM wells x spatial libraries) of the sheet:
    STEP 8  merge_sections.py: pool the spatial libraries of each GEM well, place every nucleus
            once, RNA UMAP/clusters, index of all files, merged positions
    STEP 9  report_section.py: the PDF summary
                                       -> <slidetags_output_dir>/summary/<section>/

    python3 slidetags_flex.py SHEET.csv                     all lines, then all sections
    python3 slidetags_flex.py SHEET.csv --summarize brain_section1
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
import time
from collections import Counter
from multiprocessing import Pool, get_context

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
    parser.add_argument("--list-sections", action="store_true", help="print the sections of the sheet and stop (for the WDL)")
    parser.add_argument("--summarize", metavar="SECTION", help="only make the summary of this section (steps 8-9)")
    parser.add_argument("--barcode-length", default="auto", help="Cell Ranger barcode length without -1: auto, 16, 26, ...")
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
    if args.list_sections:
        print("\n".join(sections_of_sheet(lines, args.workdir)))
        return
    if args.summarize:
        summarize(args.summarize, lines, args)
        return
    for n in to_do:
        run_line(lines[n - 1], args)
    if args.line is None:  # whole sheet in one go (local runs): also the summaries
        for section in sections_of_sheet(lines, args.workdir):
            summarize(section, lines, args)


def sections_of_sheet(lines, workdir):
    """All sections of the sheet (the sample_ids of every line's AuxFile), in order."""
    sections = []
    for line in lines:
        for section in read_auxfile(download(line["AuxFile"], os.path.join(workdir, "auxfiles", line["Link"]))[0]):
            if section not in sections:
                sections.append(section)
    return sections


def summarize(section, lines, args):
    """STEPS 8-9 for one section: see merge_sections.py and report_section.py."""
    import merge_sections
    import report_section
    print(f"\n======== summary of {section} ========")
    folder, merged = merge_sections.merge_section(section, lines, args)
    pdf = report_section.make_report(folder, merged)
    for path in [os.path.join(folder, f) for f in ("files.csv", "coords.csv", "qc_per_library.csv")] + [pdf]:
        upload(path, merged["output_dir"])


def step(text):
    """Progress line with the time, so the log shows how long each step took."""
    print(f"[{time.strftime('%H:%M:%S')}] {text}", flush=True)


def run_line(line, args):
    name = f"{line['Link']}_{line['sb_index']}"
    work = os.path.join(args.workdir, name)
    output_dir = f"{line['slidetags_output_dir'].rstrip('/')}/{name}"
    print(f"\n======== {name} ========")

    # STEP 1 - download the inputs of this line
    step("STEP 1: downloading")
    sections = read_auxfile(download(line["AuxFile"], work)[0])  # {section: spatial_tag or ""}
    pucks = dict(pair.strip().split("=", 1) for pair in line["pucks"].split(";") if pair.strip())
    missing = [s for s in sections if s not in pucks]
    if missing:
        raise SystemExit(f"{name}: no puck given for section(s) {missing} (sections come from the AuxFile)")
    puck_files = {s: download(pucks[s], f"{work}/pucks/{s}")[0] for s in sections}  # own folder: all are Puck.csv
    cellranger = f"{line['cellranger_output_dir'].rstrip('/')}/{line['Link']}/per_sample_outs"
    barcode_files = {s: download(f"{cellranger}/{s}/sample_filtered_feature_bc_matrix/barcodes.tsv.gz",
                                 f"{work}/cellranger/{s}")[0] for s in sections}
    fastq = line["Flowcell"].rstrip("/")
    read1_files = download(f"{fastq}/{line['sb_index']}_*_R1_001.fastq.gz", f"{work}/fastqs")
    read2_files = download(f"{fastq}/{line['sb_index']}_*_R2_001.fastq.gz", f"{work}/fastqs")

    # STEP 2 - load the pucks
    step("STEP 2: loading the pucks")
    beads = load_pucks(puck_files)

    # STEP 3 - the nuclei of every section, and every section's tag
    step("STEP 3: loading the nuclei")
    nuclei, nuclei_log, tags = {}, {}, {}
    for section, spatial_tag in sections.items():
        nuclei[section], nuclei_log[section], probe_barcode = load_nuclei(section, barcode_files[section], args.barcode_length)
        tags[section] = spatial_tag or probe_barcode
        print(f"    tag {tags[section]} ({'spatial_tag in AuxFile' if spatial_tag else 'Flex probe barcode'})")

    # STEP 4 - read the spatial FASTQs, filter, count
    step("STEP 4: reading the spatial reads")
    counts, cell_barcodes, bead_barcodes, read_log, library = count_spatial_reads(
        list(zip(read1_files, read2_files)), beads, tags, args.cores)
    library.update(line=name, gem_well=line["Link"], sb_index=line["sb_index"], tags=tags, read_filters=read_log,
                   read1_files=[os.path.basename(f) for f in read1_files],
                   read2_files=[os.path.basename(f) for f in read2_files])
    os.makedirs(os.path.join(work, "results"), exist_ok=True)
    with open(os.path.join(work, "results", "library.json"), "w") as f:
        json.dump(library, f, indent=1, default=int)
    upload(os.path.join(work, "results", "library.json"), output_dir)

    summary = []
    for section in sections:
        print(f"\n---- section {section} ----")
        folder = os.path.join(work, "results", section)
        os.makedirs(folder, exist_ok=True)

        # STEP 5 - the section's reads -> spatial UMIs of its nuclei, on its own puck
        step("STEP 5: giving each nucleus its spatial UMIs")
        puck = beads[beads["section"] == section].reset_index(drop=True)
        reads = counts[counts["section"] == section]
        nucleus_reads, umi_log, spatial_qc = reads_of_nuclei(reads, cell_barcodes, bead_barcodes, puck, nuclei[section])
        umis, chimera_log = umis_per_nucleus(nucleus_reads, puck)
        umi_log.update(chimera_log)
        spatial_qc["bead_umis_nuclei"] = np.bincount(umis["bead_row"], weights=umis["umis"], minlength=len(puck)).astype(np.int32)
        for what, n in umi_log.items():
            print(f"  {what}: {n:,}" if isinstance(n, (int, np.integer)) else f"  {what}: {n}")

        # STEP 6 - their positions
        step("STEP 6: placing the nuclei")
        coords, place_log = place_nuclei(umis, nuclei[section], puck, args.knn, args.cmes, args.cores)

        # STEP 7 - save
        step("STEP 7: saving")
        log = {"line": name, "section": section, "tag": tags[section], "STEP 3 nuclei": nuclei_log[section],
               "STEP 4 reads (whole library)": read_log, "STEP 5 reads (this section)": umi_log,
               "STEP 6 placement": place_log}
        save_results(folder, coords, umis, nuclei[section], puck, log, nucleus_reads, spatial_qc)
        upload(folder, output_dir)
        summary.append({"line": name, "section": section, "tag": tags[section],
                        "spatial reads": read_log["reads per section"][section], **nuclei_log[section], **place_log})

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
    """The Flex AuxFile: sample_id = the sections; spatial_tag (optional column, may be empty)."""
    with open(path, newline="") as f:
        return {row["sample_id"].strip(): (row.get("spatial_tag") or "").strip().upper()
                for row in csv.DictReader(f) if row["sample_id"].strip()}


# =====================================================================================
# STEP 2 - pucks: Puck.csv = spatial barcode, x, y (no header), one line per bead
# =====================================================================================
def load_pucks(puck_files):
    """One table with the beads of all pucks: barcode, x, y, section.
    Removed (as in the R pipeline): barcodes with an N, and degenerate barcodes (one base > 4/5 of
    the barcode, or a run of one base > 5/7 of it). A barcode that is on two pucks is kept: the
    read's tag says which section, so which puck, it came from."""
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
# STEP 3 - nuclei called by Cell Ranger in one section
# =====================================================================================
def load_nuclei(section, barcodes_file, barcode_length="auto"):
    """Cell Ranger barcode = GEM barcode (16 bases) + Flex probe barcode (multiplex Flex) + "-1".
    The spatial reads carry the GEM barcode in read 1 (and the section tag in read 2).
    Two nuclei of this section with the same GEM barcode (same droplet) cannot be told apart:
    both are left out. Also returns the section's probe barcode (empty if there is none)."""
    with gzip.open(barcodes_file, "rt") as f:
        cellranger = [row.strip() for row in f if row.strip()]
    length = len(cellranger[0].split("-")[0])
    if barcode_length != "auto" and int(barcode_length) != length:
        raise SystemExit(f"Cell Ranger barcodes have {length} bases, not {barcode_length} ({cellranger[0]})")
    probe_barcode = Counter(b.split("-")[0][CELL_BARCODE_LENGTH:] for b in cellranger).most_common(1)[0][0]
    nuclei = pd.DataFrame({"cellranger_barcode": cellranger, "gem_barcode": [b[:CELL_BARCODE_LENGTH] for b in cellranger]})
    shared = nuclei["gem_barcode"].duplicated(keep=False)
    print(f"  {section}: {len(nuclei):,} nuclei, barcode {cellranger[0]} = GEM {cellranger[0][:CELL_BARCODE_LENGTH]}"
          f" + probe barcode {probe_barcode or '(none)'}; {shared.sum()} left out (same droplet)")
    return nuclei[~shared].reset_index(drop=True), {"nuclei called": len(nuclei), "barcode length": length,
                                                    "nuclei left out (same droplet)": int(shared.sum())}, probe_barcode


# =====================================================================================
# STEP 4 - spatial reads
# =====================================================================================
def count_spatial_reads(fastq_pairs, beads, tags, cores=1):
    """Read every read pair, keep the good ones, count reads per (section, cell barcode, UMI, bead).
    Returns counts (columns section, cell, umi, bead, reads: cell and bead are row numbers in the
    two returned barcode lists), the cell barcode list, the bead barcode list and the filter log."""

    # 4a. Which read holds the bead and which bead type: count UP-site hits in the first 100,000 reads
    bead_in_read1, bead_type = detect_layout(fastq_pairs)
    layout = BEAD_TYPES[bead_type]
    part1 = PART1_LENGTH[bead_type]
    part2 = len(beads["barcode"].iat[0]) - part1  # old 14-base V10 pucks: last base not used
    print(f"  bead type {bead_type}, bead in read {1 if bead_in_read1 else 2}")

    # 4b. Where the section tag is in the bead read: the position where the tags are found most
    #     (as written, or reverse-complemented), in the first 100,000 reads
    tag_start, tag_length, section_of_tag, tag_how, tag_fraction = find_tags(fastq_pairs[0][0 if bead_in_read1 else 1], tags)
    min_read2_length = max(layout["min_length"], tag_start + tag_length)

    # 4c. What counts as a good UP site / bad UMI (1 mismatch allowed per 6 bases of UP)
    up_ok = within_mismatches(layout["up_seq"], round(len(layout["up_seq"]) / 6))
    up_all_g = within_mismatches("G" * len(layout["up_seq"]), round(len(layout["up_seq"]) / 6))
    umi_bad = set().union(*(within_mismatches(base * UMI_LENGTH, 2) for base in "ACGT"))
    bead_number = {barcode: i for i, barcode in enumerate(sorted(set(beads["barcode"])))}

    # 4d. Go through all reads. Each lane (R1/R2 file pair) is read by its own process, in parallel.
    global READ_SETUP  # what every process needs; shared with them by fork (not copied up front)
    READ_SETUP = dict(bead_in_read1=bead_in_read1, layout=layout, part1=part1, part2=part2, up_ok=up_ok,
                      up_all_g=up_all_g, umi_bad=umi_bad, bead_number=bead_number, tag_start=tag_start,
                      tag_length=tag_length, section_of_tag=section_of_tag, min_read2_length=min_read2_length,
                      section_code={s: i for i, s in enumerate(tags)})
    processes = max(1, min(cores, len(fastq_pairs)))
    print(f"  reading {len(fastq_pairs)} lane(s) with {processes} process(es)", flush=True)
    if processes == 1:
        lanes = [read_one_lane(pair) for pair in fastq_pairs]
    else:
        with get_context("fork").Pool(processes) as pool:
            lanes = pool.map(read_one_lane, fastq_pairs, chunksize=1)

    # 4e. Put the lanes together: one cell numbering, add up the reads of the same (section, cell, UMI, bead)
    log, match, fuzzy_position = Counter(), {"UP site": Counter(), "bead barcode": Counter(), "section tag": Counter()}, Counter()
    cell_number, tables = {}, []
    for lane_log, lane_match, lane_fuzzy, lane_cells, lane_counts in lanes:
        log.update(lane_log)
        fuzzy_position.update(lane_fuzzy)
        for part in match:
            match[part].update(lane_match[part])
        to_global = np.array([cell_number.setdefault(c, len(cell_number)) for c in lane_cells], dtype=np.int64)
        tables.append(lane_counts.assign(cell=to_global[lane_counts["cell"].to_numpy()] if len(lane_cells) else lane_counts["cell"]))
    counts = pd.concat(tables).groupby(["section", "cell", "umi", "bead"], as_index=False)["reads"].sum()
    counts["section"] = pd.Categorical.from_codes(counts["section"].to_numpy(), categories=list(tags))
    log["reads kept"] = int(counts["reads"].sum())
    log = dict(log)
    log["reads per section"] = {s: int(n) for s, n in counts.groupby("section", observed=True)["reads"].sum().items()}
    log["reads per section"] = {s: log["reads per section"].get(s, 0) for s in tags}
    print(f"  {log['reads']:,} reads, {log['reads kept']:,} kept: " +
          ", ".join(f"{s} {n:,}" for s, n in log["reads per section"].items()))

    # 4f. Library numbers for the summary PDF
    reads_per_umi = counts["reads"].value_counts().sort_index()  # how many (cell, UMI, bead) have 1, 2, ... reads
    fractions = np.linspace(0, 1, 21)
    library = {
        "bead type": bead_type, "bead in read 1 (R1 <-> R2)": bool(bead_in_read1),
        "tag position (1-based)": f"{tag_start + 1}-{tag_start + tag_length}",
        "tag orientation": tag_how, "reads with a tag (first 100,000)": round(tag_fraction, 4),
        "matching (reads)": {k: dict(v) for k, v in match.items()},
        "bead barcode mismatch position (reads)": {int(k) + 1: v for k, v in sorted(fuzzy_position.items())},
        "reads per UMI": {int(k): int(v) for k, v in reads_per_umi.items()},
        "sequencing saturation": round(1 - len(counts) / max(int(counts["reads"].sum()), 1), 4),
        # expected UMIs if only a fraction of the reads had been sequenced: each (cell, UMI, bead)
        # with r reads is still seen with probability 1 - (1 - f)^r
        "downsampling": {"fraction of reads": fractions.round(2).tolist(),
                         "reads": (fractions * log["reads"]).round().tolist(),
                         "UMIs": [float((reads_per_umi * (1 - (1 - f) ** reads_per_umi.index)).sum()) for f in fractions]},
    }
    return counts, list(cell_number), list(bead_number), log, library


def read_one_lane(fastq_pair):
    """STEP 4d for one lane: every read pair through the 5 filters. Returns the filter log, the
    matching counts, the bead-barcode mismatch positions, the lane's cell barcodes, and the reads
    kept, counted per (section, cell, UMI, bead) (section and cell as numbers)."""
    r = READ_SETUP
    layout, up_start, up_end = r["layout"], r["layout"]["up"][0], r["layout"]["up"][1]
    read1_file, read2_file = fastq_pair[::-1] if r["bead_in_read1"] else fastq_pair
    print(f"  reading {os.path.basename(read1_file)}", flush=True)
    log = Counter()
    match = {"UP site": Counter(), "bead barcode": Counter(), "section tag": Counter()}
    fuzzy_position = Counter()  # where the 1 mismatch is, for beads found with 1 mismatch
    cell_number = {}
    kept, tables = [], []
    columns = ["section", "cell", "umi", "bead"]
    for read1, read2 in zip(fastq_sequences(read1_file), fastq_sequences(read2_file)):
        log["reads"] += 1
        # filter 1: long enough?
        if len(read1) < CELL_BARCODE_LENGTH + UMI_LENGTH or len(read2) < r["min_read2_length"]:
            log["removed: too short"] += 1
            continue
        # filter 2: UMI without N and not one base repeated
        umi = read1[CELL_BARCODE_LENGTH:CELL_BARCODE_LENGTH + UMI_LENGTH]
        if "N" in umi or umi in r["umi_bad"]:
            log["removed: bad UMI"] += 1
            continue
        # filter 3: UP site present
        up = read2[up_start:up_end]
        if up not in r["up_ok"]:
            match["UP site"]["GG" if up in r["up_all_g"] else "none"] += 1
            log["removed: all-G read" if up in r["up_all_g"] else "removed: no UP site"] += 1
            continue
        match["UP site"]["exact" if up == layout["up_seq"] else "fuzzy"] += 1
        # filter 4: spatial barcode = a bead on the pucks (exact or 1 mismatch to one bead)
        spatial_barcode = read2[layout["part1"]:layout["part1"] + r["part1"]] + read2[layout["part2"]:layout["part2"] + r["part2"]]
        bead, how, position = find_match(spatial_barcode, r["bead_number"])
        match["bead barcode"][how] += 1
        if how == "fuzzy":
            fuzzy_position[position] += 1
        if bead is None:
            log["removed: no bead / ambiguous bead"] += 1
            continue
        # filter 5: section tag (exact or 1 mismatch to one section's tag)
        section, how, _ = find_match(read2[r["tag_start"]:r["tag_start"] + r["tag_length"]], r["section_of_tag"])
        match["section tag"][how] += 1
        if section is None:
            log["removed: no section tag"] += 1
            continue
        # keep it
        cell = cell_number.setdefault(read1[:CELL_BARCODE_LENGTH], len(cell_number))
        kept.append((r["section_code"][section], cell, umi_to_number(umi), bead))
        if len(kept) == 2_000_000:  # add up duplicates now and then to save memory
            tables.append(pd.DataFrame(kept, columns=columns).value_counts().rename("reads").reset_index())
            kept = []
    tables.append(pd.DataFrame(kept, columns=columns).value_counts().rename("reads").reset_index())
    counts = pd.concat(tables).groupby(columns, as_index=False)["reads"].sum()
    return log, match, fuzzy_position, list(cell_number), counts


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


def find_tags(bead_read_file, tags):
    """Find where the section tags sit in the bead read. Returns (start, length, {tag as read: section})."""
    if len({len(t) for t in tags.values()}) != 1 or not all(tags.values()):
        raise SystemExit(f"Every section needs a tag, all of the same length: {tags} (AuxFile column spatial_tag)")
    length = len(next(iter(tags.values())))
    reverse = str.maketrans("ACGT", "TGCA")
    options = {"as written": {t: s for s, t in tags.items()},
               "reverse complement": {t.translate(reverse)[::-1]: s for s, t in tags.items()}}
    reads = list(fastq_sequences(bead_read_file, max_reads=100_000))
    best = max(((sum(r[p:p + length] in found for r in reads), p, how)
                for how, found in options.items() for p in range(max(map(len, reads)) - length + 1)))
    hits, start, how = best
    print(f"  section tags found at read bases {start + 1}-{start + length} ({how}) in {100 * hits / len(reads):.1f}% of reads")
    if hits < 0.1 * len(reads):
        raise SystemExit(f"The section tags {tags} are in less than 10% of the spatial reads: check spatial_tag in the AuxFile")
    return start, length, options[how], how, hits / len(reads)


def find_match(barcode, table):
    """Look a barcode up in table {barcode: value}: exact, or the only entry 1 mismatch away.
    Returns (value or None, "exact" / "fuzzy" / "ambiguous" / "none", 0-based mismatch position or -1).
    Used for the bead barcodes and the section tags."""
    if barcode in table:
        return table[barcode], "exact", -1
    hits = {}
    for position, variant in one_mismatch_with_position(barcode):
        if variant in table:
            hits[table[variant]] = position
    if len(hits) == 1:
        value, position = hits.popitem()
        return value, "fuzzy", position
    return None, ("ambiguous" if hits else "none"), -1


def find_bead(barcode, bead_number):
    """Exact match, or the only bead 1 mismatch away. None if no bead or 2+ beads are 1 mismatch away."""
    return find_match(barcode, bead_number)[0]


def one_mismatch_with_position(seq):
    return [(i, seq[:i] + base + seq[i + 1:]) for i in range(len(seq)) for base in "ACGTN" if base != seq[i]]


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
# STEP 5 - spatial UMIs of every nucleus of one section
# =====================================================================================
def reads_of_nuclei(counts, cell_barcodes, bead_barcodes, puck, nuclei):
    """counts = the reads with this section's tag. Returns the reads of this section's nuclei
    (nucleus, umi, bead_row, reads; bead_row = row of the bead in puck), the log, and QC numbers."""
    log = {"reads with this section's tag": int(counts["reads"].sum())}

    # 5a. keep only reads whose bead is on this section's puck (a check: should be almost all)
    row_of = {b: i for i, b in reversed(list(enumerate(puck["barcode"])))}  # first row if a barcode is twice
    bead_row = np.array([row_of.get(b, -1) for b in bead_barcodes], dtype=np.int64)
    reads = counts.assign(bead_row=bead_row[counts["bead"].to_numpy()])
    reads = reads[reads["bead_row"] >= 0]
    log["removed: bead not on this section's puck"] = log["reads with this section's tag"] - int(reads["reads"].sum())

    # 5b. match each read's cell barcode to a nucleus: exact, or the only nucleus 1 mismatch away.
    #     If swapping bases 8-9 (as 10x does between some library types) explains more reads, use that.
    nucleus_of = {b: i for i, b in enumerate(nuclei["gem_barcode"])}
    swap = {"TC": "AG", "AG": "TC", "GT": "CA", "CA": "GT"}
    swapped = {b[:7] + swap[b[7:9]] + b[9:]: i for b, i in nucleus_of.items() if b[7:9] in swap}
    reads_per_cell = reads.groupby("cell")["reads"].sum()
    barcodes_seen = [cell_barcodes[c] for c in reads_per_cell.index]
    exact_plain = sum(n for b, n in zip(barcodes_seen, reads_per_cell) if b in nucleus_of)
    exact_swapped = sum(n for b, n in zip(barcodes_seen, reads_per_cell) if b in swapped)
    remap = len(swapped) == len(nucleus_of) and exact_swapped > exact_plain
    if remap:
        print("  using the 10x bases 8-9 swap for cell barcodes")
        nucleus_of = swapped
    nucleus, cell_matching = {}, Counter()
    for cell, barcode, n in zip(reads_per_cell.index, barcodes_seen, reads_per_cell):
        found, how, _ = find_match(barcode, nucleus_of)
        cell_matching[how] += int(n)
        if found is not None:
            nucleus[cell] = found
    reads["nucleus"] = reads["cell"].map(nucleus)

    # QC numbers: spatial UMIs per cell barcode (all barcodes / nuclei) and per bead (all barcodes)
    umis_per_cell = reads.groupby("cell").size()
    called = umis_per_cell.index.isin(list(nucleus))
    qc = {"cell_umis_all": np.sort(umis_per_cell.to_numpy())[::-1].astype(np.int32),
          "cell_umis_called": np.sort(umis_per_cell.to_numpy()[called])[::-1].astype(np.int32),
          "bead_umis_all": np.bincount(reads["bead_row"], minlength=len(puck)).astype(np.int32)}
    log["cell barcode matching (reads)"] = dict(cell_matching)
    log["10x bases 8-9 swap used"] = bool(remap)

    # 5c. remove reads whose cell barcode is not a nucleus of this section (empty droplets, ambient)
    log["removed: not a nucleus of this section"] = int(reads.loc[reads["nucleus"].isna(), "reads"].sum())
    reads = reads.dropna(subset=["nucleus"]).astype({"nucleus": int})
    reads = reads.groupby(["nucleus", "umi", "bead_row"], as_index=False)["reads"].sum()
    return reads, log, qc


def umis_per_nucleus(reads, puck):
    """reads = nucleus, umi, bead_row, reads (one section; can be several libraries pooled).
    Returns UMIs per (nucleus, bead) with the bead's position, and the log."""
    log = {}
    # 5d. remove chimeras: a UMI of a nucleus comes from ONE bead. Keep the bead with strictly most reads.
    most = reads.groupby(["nucleus", "umi"])["reads"].transform("max")
    top = reads["reads"] == most
    single_top = top.groupby([reads["nucleus"], reads["umi"]]).transform("sum") == 1
    log["removed: chimeric UMI"] = int(reads.loc[~(top & single_top), "reads"].sum())
    reads = reads[top & single_top]
    log["reads used"] = int(reads["reads"].sum())

    # 5e. UMIs per (nucleus, bead), with the bead's position on this puck
    umis = reads.groupby(["nucleus", "bead_row"]).size().rename("umis").reset_index()
    umis["barcode"] = puck["barcode"].to_numpy()[umis["bead_row"]]
    umis["x"] = puck["x"].to_numpy()[umis["bead_row"]]
    umis["y"] = puck["y"].to_numpy()[umis["bead_row"]]
    log["UMIs used"] = int(umis["umis"].sum())
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
    beads_of = {n: g for n, g in umis.groupby("nucleus")}  # nucleus = row of the nuclei table
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
        rows = pool.starmap(final_position, [(b, eps, cmes, min_pts, i2_n, i1_n)
                                             for b, i2_n, i1_n in zip(nuclei_beads, i2, i1)])

    coords = pd.concat([nuclei, pd.DataFrame(rows)], axis=1)
    coords["score"] = placement_score(coords["umi_main"].to_numpy(float), coords["umi_second"].to_numpy(float))
    coords["i1"], coords["i2"], coords["min_pts"], coords["eps"] = i1, i2, min_pts, eps
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


def final_position(beads, eps, cmes, min_pts, i2, i1=None):
    """One nucleus.
    Global DBSCAN (the section's min_pts): placed at cluster 1 if there is exactly 1 cluster.
      x, y, clusters, umi/beads (all its beads), x1, y1, umi1, beads1, rmsd1 (spread), h1 (h-index), snr
      (= umi1 / umi) of cluster 1, and x2, y2, umi2 of cluster 2 when there are 2+ clusters.
    Dynamic DBSCAN (the nucleus' own min_pts): x_dynamic, y_dynamic = its cluster at i1 (the highest
      min_pts that still gives a cluster), and umi_main / umi_second = UMIs of the two largest
      clusters at i2 (the highest min_pts giving 2+ clusters), used for the score."""
    row = dict(x=np.nan, y=np.nan, placed=False, clusters=0, umi=0, beads=0, x1=np.nan, y1=np.nan, umi1=np.nan, beads1=np.nan,
               rmsd1=np.nan, h1=np.nan, snr=np.nan, x2=np.nan, y2=np.nan, umi2=np.nan,
               x_dynamic=np.nan, y_dynamic=np.nan, umi_main=np.nan, umi_second=np.nan)
    if len(beads) == 0:
        return row
    cluster = merge_close(beads, dbscan(beads, eps, min_pts), eps * cmes)
    row.update(clusters=int(cluster.max()), umi=int(beads[:, 2].sum()), beads=len(beads))
    for c in (1, 2):
        if cluster.max() >= c:
            inside = beads[cluster == c]
            centre = np.average(inside[:, :2], axis=0, weights=inside[:, 2])
            if c == 1:
                ranked = np.sort(inside[:, 2])[::-1]
                row.update(x1=centre[0], y1=centre[1], umi1=int(inside[:, 2].sum()), beads1=len(inside),
                           rmsd1=float(np.sqrt(np.average(((inside[:, :2] - centre) ** 2).sum(1), weights=inside[:, 2]))),
                           h1=int(np.sum(ranked >= np.arange(1, len(ranked) + 1))))
                row["snr"] = row["umi1"] / row["umi"]
            else:
                row.update(x2=centre[0], y2=centre[1], umi2=int(inside[:, 2].sum()))
    if row["clusters"] == 1:
        row.update(x=row["x1"], y=row["y1"], placed=True)
    if i1:
        at_i1 = dbscan(beads, eps, i1)
        if at_i1.max() == 1:
            inside = beads[at_i1 == 1]
            row["x_dynamic"], row["y_dynamic"] = np.average(inside[:, :2], axis=0, weights=inside[:, 2])
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
def save_results(folder, coords, umis, nuclei, puck, log, nucleus_reads, spatial_qc):
    """coords.csv (one row per nucleus, join with RNA on cellranger_barcode), umis.csv.gz (UMIs per
    nucleus and bead), log.json (reads removed at every filter), summary.png (quick QC), and for the
    section summary (step 8-9): nucleus_reads.npz (reads of the nuclei, to pool libraries of one
    GEM well) and spatial_qc.npz (UMIs per cell barcode and per bead)."""
    coords.to_csv(f"{folder}/coords.csv", index=False)
    np.savez_compressed(f"{folder}/nucleus_reads.npz", nucleus=nucleus_reads["nucleus"].to_numpy(np.int32),
                        umi=nucleus_reads["umi"].to_numpy(np.int64), bead_row=nucleus_reads["bead_row"].to_numpy(np.int32),
                        reads=nucleus_reads["reads"].to_numpy(np.int32),
                        cellranger_barcode=nuclei["cellranger_barcode"].to_numpy(str), puck_beads=len(puck))
    np.savez_compressed(f"{folder}/spatial_qc.npz", **spatial_qc)
    umis.assign(cellranger_barcode=nuclei["cellranger_barcode"].to_numpy()[umis["nucleus"]])[
        ["cellranger_barcode", "barcode", "umis", "x", "y"]].to_csv(f"{folder}/umis.csv.gz", index=False)
    with open(f"{folder}/log.json", "w") as f:
        json.dump(log, f, indent=2, default=int)

    fig, (a, b, c) = plt.subplots(1, 3, figsize=(18, 5.5))
    library = {k: v for k, v in log["STEP 4 reads (whole library)"].items() if k != "reads per section"}
    funnel = {**library, **{k: v for k, v in log["STEP 5 reads (this section)"].items() if not isinstance(v, (dict, bool))}}
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
