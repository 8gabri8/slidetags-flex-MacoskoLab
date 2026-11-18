"""
STEP 8 - Merge everything about ONE section (= one sample = one puck).

A section's nuclei can be spread over several GEM wells (lines with different Link), and the
spatial reads of one GEM well over several spatial libraries (lines with the same Link and
different sb_index). This step collects all of them:

    8a. find the lines of the sheet that contain this section (its sample_id is in their AuxFile)
    8b. download, per GEM well, the Cell Ranger results of this section (RNA counts, metrics,
        unfiltered matrix) and, per spatial library, the results of steps 1-7
    8c. POOL the spatial libraries of each GEM well and place every nucleus ONCE:
          - libraries of the same GEM well come from the same cDNA, so they contain the same
            nuclei (same GEM barcodes). Their reads of each (nucleus, UMI, bead) are added up
            (the same molecule seen in two libraries counts once), chimeras are removed again,
            UMIs are counted per (nucleus, bead), and DBSCAN places each nucleus once.
          - min_pts is chosen per GEM well (as the R pipeline does per library): GEM wells can
            have very different spatial depth, and one min_pts for all would under-place the
            shallow ones. eps comes from the puck, so it is the same for all.
          - a GEM well with one library gives the same positions as its line (steps 1-7).
    8d. RNA: counts of all nuclei of the section (all GEM wells together, no batch correction):
        UMIs, genes and %mito per nucleus, then scanpy: normalize, log, 2,000 variable genes,
        PCA, neighbours, UMAP, Leiden clusters (resolution 1.0)
    8e. write, in <slidetags_output_dir>/summary/<section>/:
          files.csv            where every file of this section is (one row per GEM well x library)
          coords.csv           every nucleus of the section: GEM well, barcode, position, score,
                               spatial and RNA numbers, RNA cluster, UMAP
          qc_per_library.csv   QC numbers per GEM well x library (also page 2 of the PDF)
"""
import json
import os
import sys

import numpy as np
import pandas as pd
import scipy.io
import scipy.sparse

import slidetags_flex as core

CELLRANGER_FILES = {  # Cumulus layout: <cellranger_output_dir>/<Link>/per_sample_outs/<section>/<file>
    "filtered_matrix_h5": "sample_filtered_feature_bc_matrix.h5",
    "filtered_matrix_mex_dir": "sample_filtered_feature_bc_matrix/",
    "raw_matrix_h5": "sample_raw_feature_bc_matrix.h5",
    "molecule_info_h5": "sample_molecule_info.h5",
    "metrics_summary_csv": "metrics_summary.csv",
    "web_summary_html": "web_summary.html",
    "cloupe": "sample_cloupe.cloupe",
}
SPATIAL_FILES = {  # in <slidetags_output_dir>/<Link>_<sb_index>/<section>/
    "spatial_coords_csv": "coords.csv", "spatial_umis_csv": "umis.csv.gz", "spatial_log_json": "log.json",
    "nucleus_reads_npz": "nucleus_reads.npz", "spatial_qc_npz": "spatial_qc.npz",
}


def merge_section(section, lines, args):
    work = os.path.join(args.workdir, "summary", section)
    os.makedirs(work, exist_ok=True)

    # 8a. the lines of this section, grouped by GEM well
    mine = [line for line in lines
            if section in core.read_auxfile(core.download(line["AuxFile"], os.path.join(args.workdir, "auxfiles", line["Link"]))[0])]
    if not mine:
        raise SystemExit(f"Section {section} is in no AuxFile of the sample sheet")
    gem_wells = list(dict.fromkeys(line["Link"] for line in mine))
    spatial = [line for line in mine if line.get("sb_index")]
    output_dir = f"{(spatial or mine)[0]['slidetags_output_dir'].rstrip('/')}/summary/{section}"
    core.step(f"STEP 8: {section}: {len(gem_wells)} GEM well(s), {len(spatial)} spatial library(ies)")

    # 8b. Cell Ranger results of this section, per GEM well
    rna, metrics, raw = {}, {}, {}
    for gem_well in gem_wells:
        line = next(line for line in mine if line["Link"] == gem_well)
        source = cellranger_dir(line, section)
        folder = os.path.join(work, "cellranger", gem_well)
        mex = [core.download(f"{source}/sample_filtered_feature_bc_matrix/{f}", folder)[0]
               for f in ("barcodes.tsv.gz", "features.tsv.gz", "matrix.mtx.gz")]
        rna[gem_well] = read_mex(*mex)
        metrics[gem_well] = optional(lambda: read_metrics(core.download(f"{source}/metrics_summary.csv", folder)[0]))
        raw[gem_well] = optional(lambda: raw_barcode_totals(core.download(f"{source}/sample_raw_feature_bc_matrix.h5", folder)[0]))

    # ... the puck (the same for every line of the section)
    puck_paths = {dict(p.strip().split("=", 1) for p in line["pucks"].split(";") if p.strip()).get(section) for line in spatial}
    if len(puck_paths) > 1:
        raise SystemExit(f"{section}: different pucks on different lines of the sheet: {puck_paths}")
    puck_path = puck_paths.pop() if puck_paths else ""
    puck = core.load_pucks({section: core.download(puck_path, os.path.join(work, "puck"))[0]}) if puck_path else None
    beads_in_file = len(pd.read_csv(os.path.join(work, "puck", os.path.basename(puck_path)), header=None, usecols=[0])) if puck_path else 0

    # ... and the results of steps 1-7, per spatial library
    libraries = []
    for line in spatial:
        name = f"{line['Link']}_{line['sb_index']}"
        source = f"{line['slidetags_output_dir'].rstrip('/')}/{name}"
        folder = os.path.join(work, "lines", name)
        library = {"name": name, "gem_well": line["Link"], "sb_index": line["sb_index"], "source": source}
        with open(core.download(f"{source}/library.json", folder)[0]) as f:
            library["library"] = json.load(f)
        with open(core.download(f"{source}/{section}/log.json", folder)[0]) as f:
            library["log"] = json.load(f)
        library["coords"] = pd.read_csv(core.download(f"{source}/{section}/coords.csv", folder)[0])
        library["reads"] = dict(np.load(core.download(f"{source}/{section}/nucleus_reads.npz", folder)[0]))
        library["qc"] = dict(np.load(core.download(f"{source}/{section}/spatial_qc.npz", folder)[0]))
        if int(library["reads"]["puck_beads"]) != len(puck):
            raise SystemExit(f"{name}: its puck had {int(library['reads']['puck_beads'])} beads, the section's puck has {len(puck)}")
        libraries.append(library)

    # 8c. pool the libraries of each GEM well and place every nucleus once
    print("  placing the nuclei (spatial libraries of one GEM well pooled)")
    placements = []
    for gem_well in gem_wells:
        mine_libraries = [lib for lib in libraries if lib["gem_well"] == gem_well]
        if not mine_libraries:
            continue
        nuclei, _, _ = core.load_nuclei(section, os.path.join(work, "cellranger", gem_well, "barcodes.tsv.gz"), args.barcode_length)
        row_of = {b: i for i, b in enumerate(nuclei["cellranger_barcode"])}
        pooled = []
        for lib in mine_libraries:
            r = lib["reads"]
            nucleus = np.array([row_of.get(b, -1) for b in r["cellranger_barcode"]])[r["nucleus"]]
            pooled.append(pd.DataFrame({"nucleus": nucleus, "umi": r["umi"], "bead_row": r["bead_row"], "reads": r["reads"]}))
        pooled = pd.concat(pooled)
        pooled = pooled[pooled["nucleus"] >= 0].groupby(["nucleus", "umi", "bead_row"], as_index=False)["reads"].sum()
        umis, _ = core.umis_per_nucleus(pooled, puck)
        coords, place_log = core.place_nuclei(umis, nuclei, puck, args.knn, args.cmes, args.cores)
        coords.insert(0, "gem_well", gem_well)
        coords["spatial_libraries"] = "+".join(lib["sb_index"] for lib in mine_libraries)
        placements.append(coords)
    placements = pd.concat(placements, ignore_index=True) if placements else pd.DataFrame()

    # 8d. RNA of all nuclei of the section
    print("  RNA: UMIs, genes, %mito, UMAP, clusters")
    cells = []
    for gem_well, (matrix, barcodes, genes) in rna.items():
        mito = genes["name"].str.lower().str.startswith("mt-").to_numpy()
        umi = np.asarray(matrix.sum(1)).ravel()
        cells.append(pd.DataFrame({
            "gem_well": gem_well, "cellranger_barcode": barcodes,
            "gem_barcode": [b[:core.CELL_BARCODE_LENGTH] for b in barcodes],
            "rna_umis": umi, "rna_genes": np.diff(matrix.indptr),
            "pct_mito": 100 * np.asarray(matrix[:, mito].sum(1)).ravel() / np.maximum(umi, 1) if mito.any() else np.nan}))
    cells = pd.concat(cells, ignore_index=True)
    all_genes = rna[gem_wells[0]][2]
    counts = scipy.sparse.vstack([align_genes(m, g, all_genes) for m, _, g in rna.values()]).tocsr()
    umap, clusters = rna_analysis(counts)
    cells["rna_cluster"], cells["umap_1"], cells["umap_2"] = clusters, umap[:, 0], umap[:, 1]

    # 8e. one table with every nucleus: RNA + pooled placement
    merged = cells.merge(placements.drop(columns=["gem_barcode"], errors="ignore"),
                         on=["gem_well", "cellranger_barcode"], how="left") if len(placements) else cells.copy()
    merged["spatial_status"] = spatial_status(merged, cells, set(lib["gem_well"] for lib in libraries))
    merged["placed"] = merged.get("placed", pd.Series(False, index=merged.index)).fillna(False).astype(bool)
    merged.to_csv(os.path.join(work, "coords.csv"), index=False)

    files = files_table(section, mine, output_dir)
    files.to_csv(os.path.join(work, "files.csv"), index=False)
    qc = qc_table(merged, libraries, gem_wells)
    qc.to_csv(os.path.join(work, "qc_per_library.csv"), index=False)
    print(f"  {len(merged):,} nuclei, {merged['placed'].sum():,} placed -> {output_dir}")

    return work, {"section": section, "output_dir": output_dir, "coords": merged, "files": files, "qc": qc,
                  "gem_wells": gem_wells, "libraries": libraries, "puck": puck, "puck_path": puck_path,
                  "beads_in_file": beads_in_file,
                  "metrics": metrics, "raw": raw, "cmes": args.cmes, "knn": args.knn,
                  "mito_measured": bool(cells["pct_mito"].notna().any())}


def cellranger_dir(line, section):
    return f"{line['cellranger_output_dir'].rstrip('/')}/{line['Link']}/per_sample_outs/{section}"


def optional(make):
    """Files that only some pages need: a missing or unreadable file gives None (the page says so)."""
    try:
        return make()
    except (SystemExit, Exception) as error:  # noqa: BLE001
        print(f"  (skipped: {str(error)[:200]})")
        return None


def spatial_status(merged, cells, gem_wells_with_spatial):
    """Why a nucleus is or is not placed."""
    status = np.where(merged["gem_well"].isin(gem_wells_with_spatial), "not placed (no cluster)", "RNA only (no spatial library)")
    shared = cells.duplicated(["gem_well", "gem_barcode"], keep=False).to_numpy()
    if "clusters" in merged:
        status = np.where(merged["clusters"] >= 2, "not placed (2+ clusters)", status)
        status = np.where(merged["clusters"] == 1, "placed", status)
    status = np.where(shared & merged["gem_well"].isin(gem_wells_with_spatial), "left out (droplet shared within the section)", status)
    return status


# ---------------------------------------------------------------------------------------------
# RNA
# ---------------------------------------------------------------------------------------------
def read_mex(barcodes_file, features_file, matrix_file):
    """Cell Ranger matrix folder -> (cells x genes CSR, barcodes, genes table)."""
    barcodes = pd.read_csv(barcodes_file, header=None)[0].astype(str).to_numpy()
    genes = pd.read_csv(features_file, header=None, sep="\t").iloc[:, :2]
    genes.columns = ["id", "name"]
    matrix = scipy.sparse.csr_matrix(scipy.io.mmread(matrix_file).T)
    return matrix, barcodes, genes


def align_genes(matrix, genes, all_genes):
    """Put a GEM well's matrix in the gene order of the first GEM well (normally identical)."""
    if genes["id"].equals(all_genes["id"]):
        return matrix
    column = pd.Series(range(len(genes)), index=genes["id"]).reindex(all_genes["id"]).to_numpy()
    present = ~np.isnan(column)
    out = scipy.sparse.lil_matrix((matrix.shape[0], len(all_genes)), dtype=matrix.dtype)
    out[:, np.flatnonzero(present)] = matrix[:, column[present].astype(int)]
    return out.tocsr()


def rna_analysis(counts, seed=0):
    """Standard scanpy clustering, no batch correction: normalize to 10,000 per nucleus, log,
    2,000 most variable genes, PCA (30), 15 neighbours, UMAP, Leiden (resolution 1.0)."""
    import anndata
    import scanpy as sc
    data = anndata.AnnData(scipy.sparse.csr_matrix(counts, dtype=np.float32))
    sc.pp.normalize_total(data, target_sum=1e4)
    sc.pp.log1p(data)
    sc.pp.highly_variable_genes(data, n_top_genes=min(2000, data.n_vars - 1), flavor="seurat")
    data = data[:, data.var["highly_variable"].to_numpy()].copy()
    sc.pp.pca(data, n_comps=max(2, min(30, data.n_obs - 1, data.n_vars - 1)), random_state=seed)
    sc.pp.neighbors(data, n_neighbors=min(15, data.n_obs - 1), random_state=seed)
    sc.tl.umap(data, random_state=seed)
    sc.tl.leiden(data, resolution=1.0, random_state=seed)
    return np.asarray(data.obsm["X_umap"]), data.obs["leiden"].astype(int).to_numpy()


def read_metrics(path):
    """Cell Ranger metrics_summary.csv -> (metric, value) table. Metrics about genome / intron /
    exon mapping are dropped: Flex reads come from probes, not from the genome."""
    table = pd.read_csv(path, dtype=str)
    if "Metric Name" in table:
        grouped = table["Grouped By"].fillna("") if "Grouped By" in table else ""
        table = table[grouped == ""]
        table = pd.DataFrame({"metric": table.get("Library Type", "").fillna("").str.replace("Gene Expression", "GEX") + ": "
                              + table["Metric Name"], "value": table["Metric Value"]})
        table["metric"] = table["metric"].str.lstrip(": ")
    else:
        table = pd.DataFrame({"metric": table.columns, "value": table.iloc[0].to_numpy()})
    drop = table["metric"].str.contains("ntron|enom|ntergen|xonic|ntisense", regex=True)
    return table[~drop].drop_duplicates("metric").reset_index(drop=True)


def raw_barcode_totals(path):
    """UMIs and genes of every barcode of the unfiltered matrix (10x .h5), for the barcode rank plot."""
    import h5py
    with h5py.File(path, "r") as f:
        indptr = f["matrix/indptr"][:]
        data = f["matrix/data"][:]
        barcodes = f["matrix/barcodes"][:].astype(str)
    cumulative = np.concatenate([[0], np.cumsum(data, dtype=np.int64)])
    return pd.DataFrame({"barcode": barcodes, "umis": cumulative[indptr[1:]] - cumulative[indptr[:-1]],
                         "genes": np.diff(indptr)})


# ---------------------------------------------------------------------------------------------
# Tables
# ---------------------------------------------------------------------------------------------
def files_table(section, mine, output_dir):
    """files.csv: everything needed to build the section's object, one row per GEM well x spatial
    library (a GEM well without spatial library gets one row with empty spatial columns)."""
    rows = []
    for gem_well in dict.fromkeys(line["Link"] for line in mine):
        gem_lines = [line for line in mine if line["Link"] == gem_well]
        spatial = [line for line in gem_lines if line.get("sb_index")] or [None]
        source = cellranger_dir(gem_lines[0], section)
        for line in spatial:
            row = {"section": section, "gem_well": gem_well, "spatial_library": line["sb_index"] if line else "",
                   "cellranger_sample_dir": source,
                   **{key: f"{source}/{name}" for key, name in CELLRANGER_FILES.items()},
                   "cellranger_probe_set_csv": f"{gem_lines[0]['cellranger_output_dir'].rstrip('/')}/{gem_well}/probe_set.csv"}
            if line:
                results = f"{line['slidetags_output_dir'].rstrip('/')}/{gem_well}_{line['sb_index']}"
                pucks = dict(p.strip().split("=", 1) for p in line["pucks"].split(";") if p.strip())
                row.update({"spatial_results_dir": results, "library_json": f"{results}/library.json",
                            **{key: f"{results}/{section}/{name}" for key, name in SPATIAL_FILES.items()},
                            "puck_csv": pucks.get(section, ""),
                            "spatial_fastqs": f"{line['Flowcell'].rstrip('/')}/{line['sb_index']}_*_R[12]_001.fastq.gz"})
            row.update({"section_coords_csv": f"{output_dir}/coords.csv",
                        "section_summary_pdf": f"{output_dir}/{section}_summary.pdf"})
            rows.append(row)
    return pd.DataFrame(rows)


def qc_table(merged, libraries, gem_wells):
    """One row per GEM well x spatial library, plus one 'pooled' row per GEM well with 2+ libraries."""
    rows = []
    for gem_well in gem_wells:
        nuclei = merged[merged["gem_well"] == gem_well]
        base = {"GEM well": gem_well, "nuclei (Cell Ranger)": len(nuclei),
                "median RNA UMIs": nuclei["rna_umis"].median(), "median genes": nuclei["rna_genes"].median(),
                "median %mito": round(nuclei["pct_mito"].median(), 2) if nuclei["pct_mito"].notna().any() else np.nan}
        mine = [lib for lib in libraries if lib["gem_well"] == gem_well]
        if not mine:
            rows.append({**base, "spatial library": "none (RNA only)"})
        for lib in mine:
            log, coords = lib["log"], lib["coords"]
            step5 = log["STEP 5 reads (this section)"]
            rows.append({**base, "spatial library": lib["sb_index"],
                         "spatial reads (library)": lib["library"]["read_filters"]["reads"],
                         "reads with this section's tag": step5["reads with this section's tag"],
                         "% of tagged reads in nuclei": round(100 * step5["reads used"] / max(step5["reads with this section's tag"], 1), 1),
                         "median spatial UMIs per nucleus": coords["umi"].median(),
                         "eps": round(coords["eps"].iloc[0], 1) if len(coords) else np.nan,
                         "min_pts": coords["min_pts"].iloc[0] if len(coords) else np.nan,
                         "% placed": round(100 * coords["placed"].mean(), 1) if len(coords) else np.nan})
        if len(mine) > 1:  # the pooled placement of this GEM well (what coords.csv contains)
            pooled = nuclei[nuclei["clusters"].notna()]
            rows.append({**base, "spatial library": "pooled: " + "+".join(lib["sb_index"] for lib in mine),
                         "median spatial UMIs per nucleus": pooled["umi"].median(),
                         "eps": round(pooled["eps"].iloc[0], 1), "min_pts": pooled["min_pts"].iloc[0],
                         "% placed": round(100 * pooled["placed"].mean(), 1)})
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------------------------
# Build-time check of the scanpy / h5py parts (Dockerfile: python3 merge_sections.py --self-test)
# ---------------------------------------------------------------------------------------------
def self_test():
    import h5py
    rng = np.random.default_rng(0)
    profiles = rng.gamma(0.3, 1, size=(3, 300))  # 3 cell types, 300 genes
    counts = scipy.sparse.csr_matrix(rng.poisson(profiles[rng.integers(0, 3, 400)] * 3))
    umap, clusters = rna_analysis(counts)
    assert umap.shape == (400, 2) and clusters.min() >= 0 and len(set(clusters)) >= 2, (umap.shape, set(clusters))
    path = "/tmp/self_test_raw.h5"
    csc = counts.T.tocsc()  # 10x .h5: genes x barcodes, compressed by barcode
    with h5py.File(path, "w") as f:
        g = f.create_group("matrix")
        g["data"], g["indices"], g["indptr"] = csc.data, csc.indices, csc.indptr
        g["barcodes"] = np.array([f"BC{i}-1".encode() for i in range(400)])
    totals = raw_barcode_totals(path)
    assert (totals["umis"].to_numpy() == np.asarray(counts.sum(1)).ravel()).all()
    print(f"self-test ok: {len(set(clusters))} clusters, UMAP {umap.shape}, .h5 reader ok")


if __name__ == "__main__" and "--self-test" in sys.argv:
    self_test()
