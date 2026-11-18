"""
Simulate a small multiplexed-Flex Slide-tags experiment with known answers.

Design (like the real pilot): 2 sections, each on its own puck, Flex-barcoded, pooled and run in
2 GEM wells. Per GEM well we write:
  * SB FASTQs (R1 = cell barcode + UMI, R2 = bead) with sequencing errors, chimeras, ambient
    barcodes and broken reads. GEM well 1: two lane files, plus a second spatial library (SI-TT-A8,
    so 2 sample-sheet lines for GEM well 1). GEM well 2: R1/R2 swapped.
  * Cell Ranger multi-like outputs, laid out like Cumulus (<output_directory>/<Link>/per_sample_outs/
    <section>/sample_filtered_feature_bc_matrix/) with 26-base barcodes (16 GEM + 10 probe barcode) + "-1".
  * the section tag in spatial read 2, as in the real reads: section 1 = its probe barcode,
    section 2 = a custom tag given in the AuxFile column spatial_tag
  * the pucks, both named Puck.csv in different folders (as on the lab bucket)
plus truth.csv (true position of every nucleus), the Flex AuxFile, and sample_sheet.csv in Cumulus
format (one line per spatial library) with the slide-tags columns added. A third GEM well has Cell
Ranger results but no spatial library (empty sb_index): its nuclei are "RNA only" in the summaries.

Usage: python tests/simulate.py OUTDIR
"""
import gzip
import os
import sys

import numpy as np
import pandas as pd

rng = np.random.default_rng(42)
B = np.array(list("ACGT"))
UP2 = "CTGTTTCCTG"  # V15 bead: [15-base barcode][UP2][9 bases]
PROBE_BARCODE = {"brain_section1": "ACCTGGTTGT", "brain_section2": "AGCGAACCCT"}  # Flex probe barcodes (A-D04, A-E04)
# Section tag in the spatial reads: section 1 uses its probe barcode (AuxFile spatial_tag empty),
# section 2 a custom sequence (given in the AuxFile spatial_tag column)
SPATIAL_TAG = {"brain_section1": "ACCTGGTTGT", "brain_section2": "TTGCAGGACA"}
ADAPTER_LEFT, ADAPTER_RIGHT = "CCCATATAAGAAAACCTGAATACGCGGTT", "CGGTCCTAGCAA"  # around the tag (real reads)


def rand_seq(n, k):
    return ["".join(s) for s in B[rng.integers(0, 4, size=(n, k))]]


def mutate(s, pos=None):
    pos = rng.integers(len(s)) if pos is None else pos
    return s[:pos] + rng.choice([c for c in "ACGT" if c != s[pos]]) + s[pos + 1:]


GENES = [f"Gene{g}" for g in range(197)] + ["mt-Co1", "mt-Nd1", "mt-Atp6"]


def write_cellranger(out, gw, nuc, sections):
    """Cell Ranger multi-like per-section output: filtered matrix (MEX) + metrics_summary.csv.
    RNA: 6 cell types that follow the tissue (bands across the puck), 20 marker genes each, ~2% mito."""
    for s in sections:
        cells = nuc[nuc.section_id == s]
        folder = os.path.join(out, "cellranger_outputs", gw, "per_sample_outs", s)
        mdir = os.path.join(folder, "sample_filtered_feature_bc_matrix")
        os.makedirs(mdir, exist_ok=True)
        with gzip.open(os.path.join(mdir, "barcodes.tsv.gz"), "wt") as f:
            f.write("\n".join(cells["cellranger_barcode"]) + "\n")
        with gzip.open(os.path.join(mdir, "features.tsv.gz"), "wt") as f:
            f.write("".join(f"ENSMUSG{g:05d}\t{name}\tGene Expression\n" for g, name in enumerate(GENES)))
        cell_type = np.minimum(((cells["x"] + cells["y"]) / 4000 * 6).astype(int).to_numpy(), 5)
        rate = np.full((len(cells), len(GENES)), 0.5)
        for t in range(6):
            rate[cell_type == t, t * 20:(t + 1) * 20] = 6
        rate[:, -3:] = 0.4
        depth = rng.lognormal(0, 0.4, len(cells))[:, None]
        counts = rng.poisson(rate * depth)
        g, c = np.nonzero(counts.T)
        with gzip.open(os.path.join(mdir, "matrix.mtx.gz"), "wt") as f:
            f.write("%%MatrixMarket matrix coordinate integer general\n")
            f.write(f"{len(GENES)} {len(cells)} {len(g)}\n")
            f.write("".join(f"{a + 1} {b + 1} {counts[b, a]}\n" for a, b in zip(g, c)))
        umis = counts.sum(1)
        metrics = [("Cells", "Gene Expression", "Cells", f"{len(cells):,}"),
                   ("Cells", "Gene Expression", "Median UMI counts per cell", f"{int(np.median(umis)):,}"),
                   ("Cells", "Gene Expression", "Median genes per cell", f"{int(np.median((counts > 0).sum(1))):,}"),
                   ("Cells", "Gene Expression", "Total genes detected", f"{int((counts.sum(0) > 0).sum()):,}"),
                   ("Cells", "Gene Expression", "Confidently mapped reads in cells", "64.12%"),
                   ("Library", "Gene Expression", "Number of reads", "2,176,185,260"),
                   ("Library", "Gene Expression", "Sequencing saturation", "90.12%"),
                   ("Library", "Gene Expression", "Valid barcodes", "93.80%"),
                   ("Library", "Gene Expression", "Reads mapped confidently to intronic regions", "0.1%")]
        pd.DataFrame([dict(Category=a, **{"Library Type": b, "Grouped By": "", "Group Name": "",
                                          "Metric Name": m, "Metric Value": v}) for a, b, m, v in metrics]
                     ).to_csv(os.path.join(folder, "metrics_summary.csv"), index=False)


def make_puck(n=25_000, size=2000.0):
    sb = list(dict.fromkeys(rand_seq(int(n * 1.05), 15)))[:n]
    df = pd.DataFrame({"sb": sb, "x": rng.uniform(0, size, n), "y": rng.uniform(0, size, n)})
    return df


def main(out):
    os.makedirs(out, exist_ok=True)
    sections = list(PROBE_BARCODE)

    # --- pucks (both files named Puck.csv) + some bad beads -------------------------------
    pucks = {s: make_puck() for s in sections}
    p1, p2 = pucks[sections[0]], pucks[sections[1]]
    p2.loc[:19, "sb"] = p1.loc[:19, "sb"].to_numpy()             # 20 barcodes on both pucks
    p1.loc[20:29, "sb"] = [s[:5] + "N" + s[6:] for s in p1.loc[20:29, "sb"]]  # one N
    p1.loc[30:34, "sb"] = [s[:3] + "NN" + s[5:] for s in p1.loc[30:34, "sb"]]  # two N
    p1.loc[35:39, "sb"] = ["AAAAAAAAAAAAACG", "CCCCCCCCCCCCGTA", "GGGGGGGGGGGGTAC",
                           "TTTTTTTTTTTTACG", "ACGTAAAAAAAAAAC"]  # degenerate
    puck_paths = {}
    for i, s in enumerate(sections, start=1):
        d = os.path.join(out, "recon", f"puck_run{i}", "UMAP2_params")
        os.makedirs(d, exist_ok=True)
        puck_paths[s] = os.path.join(d, "Puck.csv")
        pucks[s].to_csv(puck_paths[s], header=False, index=False)
    good = {s: pucks[s][~pucks[s]["sb"].str.contains("N")].reset_index(drop=True) for s in sections}

    truth, sheet = [], []
    for gw_no, gw in enumerate(["brain_GEMwell1", "brain_GEMwell2"], start=1):
        # --- nuclei: 200 per section, 6 droplets shared across sections, 1 within section 1 ---
        n = 200
        gem = {s: rand_seq(n, 16) for s in sections}
        gem[sections[1]][:6] = gem[sections[0]][:6]          # cross-section shared droplets
        gem[sections[0]][199] = gem[sections[0]][198]         # same-section shared droplet
        records = []
        for s in sections:
            for i in range(n):
                kind = "normal"
                if i in (100, 101, 102, 103, 104, 105, 106, 107, 108, 109):
                    kind = "doublet"   # second cluster far away
                elif i in (110, 111, 112):
                    kind = "low"       # almost no spatial UMIs
                records.append(dict(gem_well=gw, section_id=s, gem_barcode=gem[s][i], kind=kind,
                                    cellranger_barcode=gem[s][i] + PROBE_BARCODE[s] + "-1",
                                    x=rng.uniform(200, 1800), y=rng.uniform(200, 1800)))
        nuc = pd.DataFrame(records)
        truth.append(nuc)

        # --- Cell Ranger multi-like output ----------------------------------------------------
        write_cellranger(out, gw, nuc, sections)

        # --- spatial reads ----------------------------------------------------------------------
        reads = []  # (cell barcode, umi, bead barcode, section tag)
        def add_umi(cb, bead_sb, tag):
            umi = rand_seq(1, 12)[0]
            for _ in range(1 + rng.geometric(0.4)):
                reads.append((cb, umi, bead_sb, tag))
            if rng.random() < 0.03:  # chimera: same UMI, other bead
                reads.append((cb, umi, good_beads.loc[rng.integers(len(good_beads)), "sb"], tag))

        for _, r in nuc.iterrows():
            good_beads = good[r.section_id]
            xy = good_beads[["x", "y"]].to_numpy()
            n_umi = 3 if r.kind == "low" else int(rng.lognormal(np.log(150), 0.4))
            centres = [(r.x, r.y)]
            if r.kind == "doublet":
                centres.append((2000 - r.x, 2000 - r.y))
            for u in range(n_umi):
                if rng.random() < 0.15:  # background bead anywhere on the puck
                    b = rng.integers(len(good_beads))
                else:
                    cx, cy = centres[1] if (len(centres) == 2 and u % 5 < 2) else centres[0]
                    d2 = (xy[:, 0] - cx) ** 2 + (xy[:, 1] - cy) ** 2
                    near = np.flatnonzero(d2 < 60 ** 2)
                    w = np.exp(-d2[near] / (2 * 15 ** 2))
                    b = near[rng.choice(len(near), p=w / w.sum())]
                add_umi(r.gem_barcode, good_beads.loc[b, "sb"], SPATIAL_TAG[r.section_id])
        for cb in rand_seq(2000, 16):  # ambient barcodes (empty droplets)
            s = sections[rng.integers(2)]
            good_beads = good[s]
            for _ in range(rng.integers(1, 4)):
                add_umi(cb, good_beads.loc[rng.integers(len(good_beads)), "sb"], SPATIAL_TAG[s])

        # sequencing: errors and broken reads, then write FASTQs
        order = rng.permutation(len(reads))
        r1s, r2s = [], []
        for k in order:
            cb, umi, sb, tag = reads[k]
            if rng.random() < 0.02:
                cb = mutate(cb)
            if rng.random() < 0.03:
                sb = mutate(sb)
            if rng.random() < 0.005:
                umi = umi[:4] + "N" + umi[5:]
            if rng.random() < 0.02:
                tag = mutate(tag)  # 1 error: still assigned
            elif rng.random() < 0.01:
                tag = rand_seq(1, 10)[0]  # unreadable tag: removed
            up = UP2 if rng.random() > 0.02 else "ACGTACGTAC"
            r1s.append(cb + umi)
            r2s.append(sb + up + rand_seq(1, 9)[0] + ADAPTER_LEFT + tag + ADAPTER_RIGHT)
        fq_dir = os.path.join(out, "fastqs", "run_lane1-2")
        os.makedirs(fq_dir, exist_ok=True)
        prefix = f"SI-TT-A{5 + gw_no}"
        lanes = [("L001", slice(0, len(r1s) // 2)), ("L002", slice(len(r1s) // 2, None))] if gw_no == 1 \
            else [("L001", slice(None))]
        for lane, part in lanes:
            first, second = (r1s, r2s) if gw_no == 1 else (r2s, r1s)  # GEM well 2: SB read in R1
            for name, seqs in (("R1", first[part]), ("R2", second[part]), ("I1", ["ACGTACGTAC"] * len(first[part]))):
                with gzip.open(os.path.join(fq_dir, f"{prefix}_S{gw_no}_{lane}_{name}_001.fastq.gz"), "wt") as f:
                    for i, s in enumerate(seqs):
                        f.write(f"@r{i}\n{s}\n+\n{'F' * len(s)}\n")
        print(f"{gw}: {len(reads):,} spatial reads, {len(nuc)} nuclei")
        libraries = [prefix]
        if gw_no == 1:  # a second spatial library made from the same GEM well (another 60% of its reads)
            libraries.append("SI-TT-A8")
            pick = np.sort(rng.choice(len(r1s), int(len(r1s) * 0.6), replace=False))
            for name, seqs in (("R1", r1s), ("R2", r2s)):
                with gzip.open(os.path.join(fq_dir, f"SI-TT-A8_S3_L001_{name}_001.fastq.gz"), "wt") as f:
                    for i in pick:
                        f.write(f"@r{i}\n{seqs[i]}\n+\n{'F' * len(seqs[i])}\n")

        for library in libraries:  # one sample-sheet line per spatial library
            sheet.append(dict(Link=gw, Sample=f"SI-TS-G{gw_no}", Reference="GRCm39-2024-A", Flowcell=fq_dir,
                              AuxFile=os.path.join(out, "flex_samples.csv"), DataType="flex-v2", Chemistry="auto",
                              sb_index=library, pucks=";".join(f"{s}={puck_paths[s]}" for s in sections),
                              cellranger_output_dir=os.path.join(out, "cellranger_outputs"),
                              slidetags_output_dir=os.path.join(out, "results")))

    # a third GEM well with Cell Ranger results but no spatial library (RNA only)
    nuc3 = pd.DataFrame([dict(gem_well="brain_GEMwell3", section_id=s, gem_barcode=g, kind="rna_only",
                              cellranger_barcode=g + PROBE_BARCODE[s] + "-1", x=rng.uniform(200, 1800), y=rng.uniform(200, 1800))
                         for s in sections for g in rand_seq(150, 16)])
    write_cellranger(out, "brain_GEMwell3", nuc3, sections)
    truth.append(nuc3)
    sheet.append(dict(Link="brain_GEMwell3", Sample="SI-TS-G3", Reference="GRCm39-2024-A", Flowcell=fq_dir,
                      AuxFile=os.path.join(out, "flex_samples.csv"), DataType="flex-v2", Chemistry="auto",
                      sb_index="", pucks="", cellranger_output_dir=os.path.join(out, "cellranger_outputs"),
                      slidetags_output_dir=os.path.join(out, "results")))

    with open(os.path.join(out, "flex_samples.csv"), "w") as f:  # the Flex AuxFile (as for Cumulus)
        f.write("sample_id,probe_barcode_ids,description,emptydrops_minimum_umis,spatial_tag\n"
                "brain_section1,A-D04,coronal_section_1 - regular (regular puck),100,\n"
                f"brain_section2,A-E04,coronal_section_2 - preloaded_puck,100,{SPATIAL_TAG['brain_section2']}\n")
    pd.concat(truth).to_csv(os.path.join(out, "truth.csv"), index=False)
    pd.DataFrame(sheet).to_csv(os.path.join(out, "sample_sheet.csv"), index=False)
    print(f"Wrote {out}/sample_sheet.csv and truth.csv")


if __name__ == "__main__":
    main(sys.argv[1])
