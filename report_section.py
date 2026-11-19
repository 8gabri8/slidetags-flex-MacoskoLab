"""
STEP 9 - The PDF summary of ONE section (= sample = puck): <section>_summary.pdf

Pages (each page is drawn on its own: if one fails, it shows the error and the others are still made)
     1  Overview                 what went in, headline numbers, where every nucleus ended up
     2  QC per library           one row per GEM well x spatial library (+ pooled)
     3  GEX metrics              Cell Ranger metrics of this section, per GEM well
     4  Cell calling             barcode rank plot and genes vs UMIs (unfiltered matrix)
     5  RNA                      UMAP: clusters, GEM well, UMIs per cluster, %mito, genes
     6  Spatial metrics          per spatial library: reads, puck, matching, filters, FASTQs
     7  Spatial counts           spatial UMIs per cell / per bead, downsampling, reads per UMI
     8  Puck                     spatial UMIs on the puck (all barcodes / nuclei only)
     9  Global DBSCAN            choice of min_pts
    10  DBSCAN = 1               nuclei with one cluster: positions, SNR, spread, h-index
    11  DBSCAN = 2               nuclei with two clusters
    12  Dynamic DBSCAN           each nucleus with its own min_pts, and the score
    13  Placement by cluster     placed nuclei coloured by RNA cluster
    14  Placement by UMIs        placed nuclei coloured by RNA / spatial UMIs
    15  RNA vs spatial           RNA UMIs vs spatial UMIs, % placed per cluster
    16  Pipeline checks          section tags, pooled vs single library, why nuclei are not placed
"""
import os
import traceback

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.backends.backend_pdf import PdfPages  # noqa: E402
from scipy.stats import gaussian_kde  # noqa: E402

PAGE = (8.5, 11)
SCORE = matplotlib.colors.LinearSegmentedColormap.from_list("score", ["#132B43", "#56B1F7"])  # ggplot2 default
COLOURS = list(plt.get_cmap("tab20").colors) + list(plt.get_cmap("tab20b").colors) + list(plt.get_cmap("tab20c").colors)
plt.rcParams.update({"font.size": 8, "axes.titlesize": 10, "axes.spines.top": False, "axes.spines.right": False})


def make_report(folder, m):
    path = os.path.join(folder, f"{m['section']}_summary.pdf")
    import slidetags_flex as core
    core.step(f"STEP 9: {path}")
    pages = [("overview", overview), ("QC per library", qc_per_library), ("GEX metrics", gex_metrics),
             ("cell calling", cell_calling), ("RNA", rna)]
    pages += [(f"spatial metrics {lib['name']}", lambda m, fig, lib=lib: spatial_metrics(m, fig, lib)) for lib in m["libraries"]]
    pages += [("spatial counts", spatial_counts), ("puck", puck_umis), ("global DBSCAN", dbscan_optimisation),
              ("DBSCAN=1", dbscan_one), ("DBSCAN=2", dbscan_two), ("dynamic DBSCAN", dynamic_dbscan),
              ("placement by cluster", placement_by_cluster), ("placement by UMIs", placement_by_umis),
              ("RNA vs spatial", rna_vs_spatial), ("pipeline checks", pipeline_checks)]
    with PdfPages(path) as pdf:
        for name, page in pages:
            fig = plt.figure(figsize=PAGE)
            try:
                page(m, fig)
            except Exception:  # noqa: BLE001 - one broken page must not stop the report
                fig.clf()
                fig.text(0.05, 0.95, f"Page '{name}' could not be drawn:", va="top", fontsize=11)
                fig.text(0.05, 0.9, traceback.format_exc()[-3000:], va="top", fontsize=6, family="monospace")
                print(f"  page '{name}' failed: {traceback.format_exc().splitlines()[-1]}")
            pdf.savefig(fig)
            plt.close("all")
    return path


# ---------------------------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------------------------
def spatial(m):
    """Nuclei that went through the spatial placement (a GEM well with a spatial library, not left out)."""
    c = m["coords"]
    return c[c["clusters"].notna()] if "clusters" in c else c.iloc[:0]


def title(fig, text, sub=None):
    fig.suptitle(text, fontsize=14, y=0.985)
    if sub:
        fig.text(0.5, 0.955, sub, ha="center", fontsize=8, color="0.3")


def message(ax, text):
    ax.axis("off")
    ax.text(0.5, 0.5, text, ha="center", va="center", fontsize=9, color="0.4", wrap=True, transform=ax.transAxes)


def table(ax, rows, columns=None, title_text=None, font=7, col_widths=None):
    ax.axis("off")
    if title_text:
        ax.set_title(title_text)
    rows = [[fmt(v) for v in r] for r in rows]
    if not rows:
        return message(ax, "(nothing)")
    t = ax.table(cellText=rows, colLabels=columns, loc="upper center", cellLoc="center", colWidths=col_widths)
    t.auto_set_font_size(False)
    t.set_fontsize(font)
    t.scale(1, 1.25)
    for (r, _), cell in t.get_celld().items():
        cell.set_edgecolor("white")
        cell.set_facecolor("#d9d9d9" if r == 0 and columns else ("#f2f2f2" if r % 2 else "#e6e6e6"))


def fmt(v):
    if isinstance(v, (bool, np.bool_)):
        return str(v)
    if isinstance(v, (int, np.integer)):
        return f"{v:,}"
    if isinstance(v, (float, np.floating)):
        return "" if np.isnan(v) else (f"{v:,.0f}" if abs(v) >= 1000 or v == int(v) else f"{v:,.2f}")
    return str(v)


def puck_axes(ax, m, label=True):
    p = m["puck"]
    if p is not None and len(p):
        pad = 0.02 * (p["x"].max() - p["x"].min())
        ax.set_xlim(p["x"].min() - pad, p["x"].max() + pad)
        ax.set_ylim(p["y"].min() - pad, p["y"].max() + pad)
    ax.set_aspect("equal")
    if not label:
        ax.set_xticks([]), ax.set_yticks([])


def outline(ax, m):
    """The puck in light grey (all beads, binned) under the nuclei."""
    p = m["puck"]
    if p is not None and len(p):
        h, xe, ye = np.histogram2d(p["x"], p["y"], bins=100)
        ax.imshow(np.where(h.T > 0, 1, np.nan), origin="lower", extent=(xe[0], xe[-1], ye[0], ye[-1]),
                  cmap="Greys", vmin=0, vmax=10, alpha=0.5, interpolation="nearest")


def density(ax, values, name, percent=False):
    v = np.asarray(values, float)
    v = v[np.isfinite(v)]
    if len(v) < 3 or np.ptp(v) == 0:
        return message(ax, f"{name}: too few values")
    grid = np.linspace(v.min(), v.max(), 200)
    ax.plot(grid, gaussian_kde(v)(grid), color="k")
    ax.axvline(v.mean(), color="red", ls="--")
    ax.text(v.mean(), ax.get_ylim()[1], f" Mean: {100 * v.mean():.1f}%" if percent else f" Mean: {v.mean():.2f}",
            color="red", va="top")
    ax.set(xlabel=name, ylabel="Density", title=f"{name} distribution")


def cluster_colours(clusters):
    return [COLOURS[int(c) % len(COLOURS)] for c in clusters]


def violin(ax, groups, labels):
    keep = [(g, l) for g, l in zip(groups, labels) if len(g) > 1]
    if not keep:
        return message(ax, "too few values")
    ax.violinplot([g for g, _ in keep], showmedians=True)
    ax.set_xticks(range(1, len(keep) + 1), [str(l) for _, l in keep])


def lines_of(m):
    return [f"{lib['gem_well']} / {lib['sb_index']}" for lib in m["libraries"]]


# ---------------------------------------------------------------------------------------------
# pages
# ---------------------------------------------------------------------------------------------
def overview(m, fig):
    c, s = m["coords"], spatial(m)
    title(fig, f"Slide-tags summary: {m['section']}")
    placed = int(c["placed"].sum())
    tags = sorted({lib["library"]["tags"][m["section"]] for lib in m["libraries"]})
    min_pts = s.groupby("gem_well")["min_pts"].first().to_dict() if len(s) else {}
    rows = [["Section (sample)", m["section"]],
            ["Puck", os.path.basename(os.path.dirname(os.path.dirname(m["puck_path"]))) + "/…/" + os.path.basename(m["puck_path"]) if m["puck_path"] else "none"],
            ["Section tag", ", ".join(tags) or "-"],
            ["GEM wells", ", ".join(m["gem_wells"])],
            ["Spatial libraries", ", ".join(lines_of(m)) or "none"],
            ["Nuclei (Cell Ranger)", len(c)],
            ["Nuclei with a spatial library", int((c["spatial_status"] != "RNA only (no spatial library)").sum())],
            ["Nuclei placed", f"{placed:,} ({100 * placed / max(len(s), 1):.1f}% of nuclei with a spatial library)"],
            ["Median RNA UMIs / genes per nucleus", f"{c['rna_umis'].median():,.0f} / {c['rna_genes'].median():,.0f}"],
            ["Median spatial UMIs per nucleus", fmt(s["umi"].median()) if len(s) else "-"],
            ["eps (µm) / knn / cmes", f"{s['eps'].iloc[0]:.1f} / {m['knn']} / {m['cmes']}" if len(s) else "-"],
            ["min_pts per GEM well", ", ".join(f"{k}: {int(v)}" for k, v in min_pts.items()) or "-"],
            ["RNA clusters (Leiden)", c["rna_cluster"].nunique()]]
    table(fig.add_axes([0.05, 0.58, 0.9, 0.35]), rows, None, font=9, col_widths=[0.35, 0.65])
    ax = fig.add_axes([0.33, 0.08, 0.18, 0.42])
    counts = c["spatial_status"].value_counts()
    ax.barh(counts.index[::-1], counts.to_numpy()[::-1], color="0.4")
    for i, v in enumerate(counts.to_numpy()[::-1]):
        ax.text(v, i, f" {v:,}", va="center")
    ax.set_title("What happened to each nucleus")
    ax = fig.add_axes([0.58, 0.08, 0.4, 0.42])
    outline(ax, m)
    p = c[c["placed"]]
    ax.scatter(p["x"], p["y"], s=2, c=cluster_colours(p["rna_cluster"]), rasterized=True)
    puck_axes(ax, m)
    ax.set_title(f"Placed nuclei by RNA cluster ({placed:,})")


def qc_per_library(m, fig):
    title(fig, "QC per library", "one row per GEM well x spatial library; 'pooled' = the libraries of one GEM well together (what coords.csv contains)")
    qc = m["qc"]
    table(fig.add_axes([0.02, 0.45, 0.96, 0.47]), [[k] + list(qc[k]) for k in qc.columns if k not in ("GEM well", "spatial library")],
          ["metric"] + [f"{g}\n{lib}" for g, lib in zip(qc["GEM well"], qc["spatial library"])], font=6.5,
          col_widths=[0.24] + [0.76 / len(qc)] * len(qc))
    s = qc.dropna(subset=["% placed"])
    if len(s):
        names = [f"{g}\n{lib}" for g, lib in zip(s["GEM well"], s["spatial library"])]
        ax = fig.add_axes([0.08, 0.06, 0.4, 0.3])
        ax.bar(names, s["% placed"], color="0.4")
        ax.set(ylabel="% placed", ylim=(0, 100), title="% of nuclei placed")
        ax.tick_params(axis="x", labelsize=6, rotation=30)
        ax = fig.add_axes([0.57, 0.06, 0.4, 0.3])
        ax.bar(names, s["median spatial UMIs per nucleus"], color="0.4")
        ax.set(ylabel="spatial UMIs", title="Median spatial UMIs per nucleus")
        ax.tick_params(axis="x", labelsize=6, rotation=30)


def gex_metrics(m, fig):
    title(fig, "GEX Metrics Summary", "Cell Ranger, this section only; genome/intron/exon mapping metrics left out (Flex reads come from probes)")
    tables = {g: t for g, t in m["metrics"].items() if t is not None}
    if not tables:
        return message(fig.add_subplot(), "No metrics_summary.csv found for this section")
    merged = None
    for g, t in tables.items():
        t = t.rename(columns={"value": g})
        merged = t if merged is None else merged.merge(t, on="metric", how="outer")
    table(fig.add_axes([0.03, 0.03, 0.94, 0.9]), merged.fillna("").values.tolist()[:45], ["metric"] + list(tables),
          font=7, col_widths=[0.5] + [0.5 / len(tables)] * len(tables))


def cell_calling(m, fig):
    c = m["coords"]
    title(fig, f"Cell calling (Cell Ranger): {len(c):,} nuclei", "unfiltered matrix of this section (sample_raw_feature_bc_matrix.h5)")
    ax1, ax2 = fig.add_subplot(2, 1, 1), fig.add_subplot(2, 1, 2)
    drawn = False
    for i, (g, raw) in enumerate(m["raw"].items()):
        if raw is None:
            continue
        drawn = True
        cells = set(c.loc[c["gem_well"] == g, "cellranger_barcode"])
        raw = raw[raw["umis"] > 0].sort_values("umis", ascending=False).reset_index(drop=True)
        is_cell = raw["barcode"].isin(cells).to_numpy()
        rank = np.arange(1, len(raw) + 1)
        ax1.plot(rank, raw["umis"], color=COLOURS[2 * i + 1], lw=1, label=f"{g} background")
        ax1.plot(rank[is_cell], raw["umis"][is_cell], ".", ms=1, color=COLOURS[2 * i], label=f"{g} nuclei")
        top = raw.iloc[:100_000]
        ax2.scatter(np.log10(top["umis"]), top["genes"], s=1, alpha=0.3, rasterized=True,
                    c=np.where(top["barcode"].isin(cells), "tab:cyan", "tab:red"))
    if not drawn:
        message(ax1, "Unfiltered matrix not available\n(sample_raw_feature_bc_matrix.h5 missing or unreadable)")
        nuclei = c.sort_values("rna_umis", ascending=False)
        ax2.plot(np.arange(1, len(nuclei) + 1), nuclei["rna_umis"], color="tab:cyan")
        ax2.set(xscale="log", yscale="log", xlabel="nucleus rank", ylabel="RNA UMIs", title="RNA UMIs of the called nuclei")
        return
    ax1.set(xscale="log", yscale="log", xlabel="Barcode rank", ylabel="UMI counts", title="Cell barcode rank plot")
    ax1.text(0.02, 0.05, f"Median UMI counts per nucleus: {c['rna_umis'].median():,.0f}", transform=ax1.transAxes)
    ax1.legend(fontsize=6, markerscale=5)
    ax2.set(xlabel="log10(UMI)", ylabel="genes", title="Genes vs UMIs (top 100,000 barcodes; cyan = nuclei, red = background)")


def rna(m, fig):
    c = m["coords"]
    if m.get("rna_error"):
        return message(fig.add_subplot(), f"RNA clustering (UMAP, Leiden) failed:\n{m['rna_error'][:300]}\n\n"
                                          "Everything else in this summary is fine.")
    fig.subplots_adjust(hspace=0.35, wspace=0.25)
    title(fig, "RNA", f"{len(c):,} nuclei, all GEM wells together (no batch correction); scanpy: 2,000 variable genes, PCA, UMAP, Leiden 1.0")
    ax = fig.add_subplot(3, 2, 1)
    ax.scatter(c["umap_1"], c["umap_2"], s=1, c=cluster_colours(c["rna_cluster"]), rasterized=True)
    for k, g in c.groupby("rna_cluster"):
        ax.text(g["umap_1"].median(), g["umap_2"].median(), str(k), fontsize=8, weight="bold")
    ax.set(title="UMAP (Leiden clusters)", xticks=[], yticks=[])
    ax = fig.add_subplot(3, 2, 2)
    order = sorted(c["rna_cluster"].unique())
    violin(ax, [np.log10(c.loc[c["rna_cluster"] == k, "rna_umis"].clip(lower=1)) for k in order], order)
    ax.set(title="log10 RNA UMIs per cluster")
    ax = fig.add_subplot(3, 2, 3)
    if m["mito_measured"]:
        sc = ax.scatter(c["umap_1"], c["umap_2"], s=1, c=c["pct_mito"], cmap="Purples", rasterized=True)
        fig.colorbar(sc, ax=ax, shrink=0.7)
        ax.set(title=f"%Mito (median {c['pct_mito'].median():.2f}%, mean {c['pct_mito'].mean():.2f}%)", xticks=[], yticks=[])
    else:
        message(ax, "%Mito: not measured\n(no mt- probes in the probe set)")
    ax = fig.add_subplot(3, 2, 4)
    sc = ax.scatter(c["umap_1"], c["umap_2"], s=1, c=np.log10(c["rna_genes"].clip(lower=1)), cmap="viridis", rasterized=True)
    fig.colorbar(sc, ax=ax, shrink=0.7)
    ax.set(title="log10 genes per nucleus", xticks=[], yticks=[])
    ax = fig.add_subplot(3, 2, 5)
    for i, (g, d) in enumerate(c.groupby("gem_well")):
        ax.scatter(d["umap_1"], d["umap_2"], s=1, color=COLOURS[2 * i], label=g, rasterized=True)
    ax.legend(fontsize=6, markerscale=5)
    ax.set(title="GEM well (batch check)", xticks=[], yticks=[])
    ax = fig.add_subplot(3, 2, 6)
    n = c.groupby(["rna_cluster", "gem_well"]).size().unstack(fill_value=0)
    (n / n.sum()).mul(100).plot.bar(ax=ax, width=0.8, legend=False, color=COLOURS[::2])
    ax.set(title="% of each GEM well's nuclei per cluster", xlabel="cluster")


def spatial_metrics(m, fig, lib):
    L, log = lib["library"], lib["log"]
    title(fig, f"Spatial Metrics Summary: {lib['gem_well']} / {lib['sb_index']}")
    reads = L["read_filters"]["reads"]
    step5 = log["STEP 5 reads (this section)"]
    table(fig.add_axes([0.03, 0.72, 0.44, 0.2]), [
        ["Total reads", reads], ["Reads kept (all sections)", L["read_filters"]["reads kept"]],
        ["Reads with this section's tag", step5["reads with this section's tag"]],
        ["Final UMIs (this section)", step5["UMIs used"]],
        ["R1 <-> R2", L["bead in read 1 (R1 <-> R2)"]], ["Remap CB (10x bases 8-9)", step5.get("10x bases 8-9 swap used", "")],
        ["Bead type", L["bead type"]], ["Sequencing saturation", f"{100 * L['sequencing saturation']:.1f}%"]],
        None, "Library information", col_widths=[0.6, 0.4])
    p = m["puck"]
    table(fig.add_axes([0.53, 0.72, 0.44, 0.2]), [
        ["Puck file", os.path.basename(m["puck_path"])], ["Beads total", m["beads_in_file"]],
        ["Beads removed (N / degenerate)", m["beads_in_file"] - len(p)],
        ["Diameter", f"{max(np.ptp(p['x']), np.ptp(p['y'])):,.0f}"],
        ["Section tag", f"{L['tags'][m['section']]} (read bases {L['tag position (1-based)']}, {L['tag orientation']})"]],
        None, "Puck information", col_widths=[0.45, 0.55])

    def pct(d):
        total = max(sum(d.values()), 1)
        return [[k, f"{100 * d.get(k, 0) / total:.2f}%"] for k in ("exact", "fuzzy", "GG", "ambiguous", "none") if k in d]
    matching = L["matching (reads)"]
    for i, (name, d) in enumerate([("UP matching", matching["UP site"]), ("SB matching", matching["bead barcode"]),
                                   ("Tag matching", matching["section tag"]),
                                   ("CB matching", step5.get("cell barcode matching (reads)", {}))]):
        table(fig.add_axes([0.03 + i * 0.245, 0.55, 0.22, 0.13]), pct(d), None, name, col_widths=[0.55, 0.45])
    filters = {k.replace("removed: ", ""): v for k, v in L["read_filters"].items() if k.startswith("removed")}
    filters.update({k.replace("removed: ", "") + " (this section)": v for k, v in step5.items()
                    if k.startswith("removed") and isinstance(v, (int, float))})
    table(fig.add_axes([0.03, 0.27, 0.4, 0.24]), [[k, f"{100 * v / reads:.2f}%"] for k, v in filters.items()],
          ["Filter", "% of reads"], "SB filtering", col_widths=[0.75, 0.25])
    ax = fig.add_axes([0.52, 0.31, 0.44, 0.18])
    pos = {int(k): v for k, v in L["bead barcode mismatch position (reads)"].items()}
    length = len(m["puck"]["barcode"].iat[0])
    ax.bar(list(pos), list(pos.values()), color="0.35")
    ax.set(xlabel="Spatial barcode base position", ylabel="Fuzzy matches", title="Location of spatial barcode fuzzy match",
           xlim=(0.5, length + 0.5), xticks=range(1, length + 1))
    files = list(zip(L["read1_files"], L["read2_files"]))
    rows = [[a, b] for a, b in files[:6]] + ([[f"({len(files) - 6} more)", ""]] if len(files) > 6 else [])
    table(fig.add_axes([0.15, 0.03, 0.7, 0.2]), rows, ["R1s", "R2s"], "Spatial library FASTQs")


def spatial_counts(m, fig):
    title(fig, "Spatial counts", "one line per spatial library")
    ax1, ax2, ax3, ax4 = (fig.add_subplot(2, 2, i) for i in range(1, 5))
    width = 0.8 / max(len(m["libraries"]), 1)
    for i, lib in enumerate(m["libraries"]):
        name, q, L = f"{lib['gem_well']}/{lib['sb_index']}", lib["qc"], lib["library"]
        ax1.plot(np.arange(1, len(q["cell_umis_all"]) + 1), q["cell_umis_all"], color=COLOURS[2 * i + 1], label=f"{name} all")
        ax1.plot(np.arange(1, len(q["cell_umis_called"]) + 1), q["cell_umis_called"], color=COLOURS[2 * i], label=f"{name} nuclei")
        share = q["cell_umis_called"].sum() / max(q["cell_umis_all"].sum(), 1)
        ax1.text(0.98, 0.95 - 0.06 * i, f"{name}: {100 * share:.1f}% of UMIs in nuclei", transform=ax1.transAxes, ha="right", fontsize=7)
        for arr, shade, lab in ((q["bead_umis_all"], 1, "all"), (q["bead_umis_nuclei"], 0, "nuclei")):
            arr = np.sort(arr[arr > 0])[::-1]
            ax2.plot(np.arange(1, len(arr) + 1), arr, color=COLOURS[2 * i + shade], label=f"{name} {lab}")
        d = L["downsampling"]
        ax3.plot(np.array(d["reads"]) / 1e6, np.array(d["UMIs"]) / 1e6, "o-", ms=3, color=COLOURS[2 * i], label=name)
        rpu = {int(k): v for k, v in L["reads per UMI"].items()}
        x = np.arange(1, 11)
        y = [rpu.get(k, 0) for k in range(1, 10)] + [sum(v for k, v in rpu.items() if k >= 10)]
        ax4.bar(x + i * width - 0.4 + width / 2, np.array(y) / 1e6, width=width, color=COLOURS[2 * i],
                label=f"{name}: saturation {100 * L['sequencing saturation']:.1f}%")
    ax1.set(xscale="log", yscale="log", xlabel="Cell barcode rank", ylabel="SB UMI counts", title="SB UMI per cell (this section's tag)")
    ax2.set(xscale="log", yscale="log", xlabel="Beads", ylabel="SB UMI counts", title="SB UMI per bead")
    ax3.set(xlabel="Millions of SB reads", ylabel="Millions of filtered SB UMIs", title="SB downsampling curve (whole library)")
    ax4.set(xlabel="Reads per UMI", ylabel="Millions of filtered SB UMIs", title="Reads per UMI (whole library)",
            xticks=range(1, 11), xticklabels=[str(k) for k in range(1, 10)] + ["10+"])
    for ax in (ax1, ax2, ax3, ax4):
        ax.legend(fontsize=6)


def puck_umis(m, fig):
    title(fig, "SB UMI per bead on the puck", "all libraries of this section added up; beads binned (sum of UMIs per bin)")
    p = m["puck"]
    total = {k: sum(lib["qc"][k].astype(np.int64) for lib in m["libraries"]) for k in ("bead_umis_all", "bead_umis_nuclei")}
    for i, (k, name) in enumerate([("bead_umis_all", "SB UMI per bead (all cell barcodes)"),
                                   ("bead_umis_nuclei", "SB UMI per bead (called nuclei only)")]):
        ax = fig.add_subplot(2, 1, i + 1)
        h, xe, ye = np.histogram2d(p["x"], p["y"], bins=150, weights=total[k])
        present, _, _ = np.histogram2d(p["x"], p["y"], bins=[xe, ye])
        im = ax.imshow(np.where(present.T > 0, np.log10(h.T + 1), np.nan), origin="lower", extent=(xe[0], xe[-1], ye[0], ye[-1]),
                       cmap="inferno", interpolation="nearest")
        fig.colorbar(im, ax=ax, shrink=0.8, label="log10(UMI + 1)")
        ax.set(title=name, aspect="equal")


def dbscan_optimisation(m, fig):
    s = spatial(m)
    title(fig, "Global DBSCAN Optimization", "min_pts is chosen per GEM well (libraries of a GEM well pooled)")
    ax = fig.add_subplot(2, 2, 1)
    k = s["clusters"].clip(upper=5).astype(int).value_counts().reindex(range(6), fill_value=0)
    pct = 100 * k / max(len(s), 1)
    ax.bar(["0", "1", "2", "3", "4", "5+"], pct, color="0.35")
    for i, v in enumerate(pct):
        ax.text(i, v, f"{v:.0f}%", ha="center", va="bottom")
    ax.set(ylim=(0, 105), xlabel="DBSCAN clusters", ylabel="Percent", title="Cluster distribution")
    ax = fig.add_subplot(2, 2, 2)
    note = []
    for i, (g, d) in enumerate(s.groupby("gem_well")):
        top = int(max(d["i1"].max(), 1))
        values = np.arange(1, top + 1)
        placed = np.array([np.mean((d["i2"] < v) & (d["i1"] >= v)) for v in values]) * 100
        mp = int(d["min_pts"].iloc[0])
        show = values <= max(20, 2 * mp)
        ax.plot(values[show], placed[show], color=COLOURS[2 * i], label=g)
        ax.axvline(mp, color=COLOURS[2 * i], ls="--")
        note.append(f"{g}: minPts {mp}, placed {placed[mp - 1]:.2f}%")
    ax.text(0.98, 0.04, "\n".join(note + [f"eps: {s['eps'].iloc[0]:.2f} | k: {m['knn']}"]), transform=ax.transAxes,
            ha="right", va="bottom", color="red", fontsize=7)
    ax.set(xlabel="minPts", ylabel="%Placed", title="Parameter optimization")
    ax.legend(fontsize=6, loc="upper right")
    ax = fig.add_subplot(2, 2, 3)
    groups = [np.log10(s.loc[s["clusters"].clip(upper=5) == n, "umi"].clip(lower=1)) for n in range(6)]
    violin(ax, groups, ["0", "1", "2", "3", "4", "5+"])
    ax.set(xlabel="DBSCAN clusters", ylabel="log10 SB UMI", title="SB UMI distribution")
    ax.text(0.98, 0.02, f"Mean SB UMI per nucleus: {s['umi'].mean():,.1f}", transform=ax.transAxes, ha="right")
    ax = fig.add_subplot(2, 2, 4)
    colour = np.select([s["clusters"] == 0, s["clusters"] == 1], ["#619CFF", "#F8766D"], "#00BA38")
    ax.scatter(s["i2"] + 1, (s["i1"] - s["i2"]).clip(lower=0.8), s=1, c=colour, rasterized=True)
    for i, (g, d) in enumerate(s.groupby("gem_well")):
        mp = int(d["min_pts"].iloc[0])
        ax.axvline(mp, color="k", lw=0.8)
        x = np.arange(1, mp + 1)
        ax.plot(x, mp - x + 1, color="k", lw=0.8)
    ax.set(xscale="log", yscale="log", xlabel="Minimum minPts", ylabel="Placement Interval Width",
           title="Placement interval (blue 0, red 1, green 2+ clusters)")


def dbscan_one(m, fig):
    s = spatial(m)
    one = s[s["clusters"] == 1]
    title(fig, "Global DBSCAN=1 Results", f"{len(one):,} nuclei with exactly one cluster")
    ax = fig.add_subplot(2, 2, 1)
    ax.scatter(one["x"], one["y"], s=0.5, c="k", rasterized=True)
    puck_axes(ax, m)
    ax.set_title("DBSCAN=1 Placements")
    density(fig.add_subplot(2, 2, 2), one["snr"], "SNR", percent=True)
    density(fig.add_subplot(2, 2, 3), one["rmsd1"], "RMSD")
    ax = fig.add_subplot(2, 2, 4)
    h = one["h1"].dropna().astype(int)
    if len(h):
        ax.bar(*np.unique(h, return_counts=True), color="0.35")
        ax.axvline(h.mean(), color="red", ls="--")
        ax.text(h.mean(), ax.get_ylim()[1], f" Mean: {h.mean():.1f}", color="red", va="top")
    ax.set(xlabel="h-index", ylabel="Count", title="h-index distribution")


def dbscan_two(m, fig):
    s = spatial(m)
    two = s[s["clusters"] == 2]
    title(fig, "Global DBSCAN=2 Results", f"{len(two):,} nuclei with two clusters (not placed)")
    for i, (x, y, name) in enumerate([("x1", "y1", "DBSCAN=1/2 Positions"), ("x2", "y2", "DBSCAN=2/2 Positions")]):
        ax = fig.add_subplot(2, 2, i + 1)
        ax.scatter(two[x], two[y], s=1, c="k", rasterized=True)
        puck_axes(ax, m)
        ax.set_title(name)
    ax = fig.add_subplot(2, 2, 3)
    ax.scatter(np.log10(two["umi1"]), np.log10(two["umi2"]), s=2, c="k", rasterized=True)
    ax.set(xlabel="log10(UMI1)", ylabel="log10(UMI2)", title="SB UMI counts")
    ax = fig.add_subplot(2, 2, 4)
    dist = (np.hypot(two["x1"] - two["x2"], two["y1"] - two["y2"])).to_numpy()
    if len(dist):
        ax.hist(dist[np.isfinite(dist)], bins=30, color="0.35")
    if len(s):
        ax.axvline(m["cmes"] * s["eps"].iloc[0], color="red", ls="--")
        ax.text(m["cmes"] * s["eps"].iloc[0], ax.get_ylim()[1], f" CMES: {m['cmes']}", color="red", va="top")
    ax.set(xlabel="Distance", ylabel="Nuclei", title="DBSCAN 1-2 Distance")


def dynamic_dbscan(m, fig):
    s = spatial(m)
    title(fig, "Dynamic DBSCAN Results",
          "each nucleus clustered with its own min_pts; score = how much its main cluster dominates the second (1 = clear)")
    ax = fig.add_subplot(2, 2, 1)
    d = s.dropna(subset=["x_dynamic"])
    sc = ax.scatter(d["x_dynamic"], d["y_dynamic"], s=1, c=d["score"], cmap=SCORE, vmin=0, vmax=1, rasterized=True)
    fig.colorbar(sc, ax=ax, shrink=0.7)
    puck_axes(ax, m)
    ax.set_title(f"Scored DBSCAN Positions ({len(d):,})")
    ax = fig.add_subplot(2, 2, 2)
    score = s["score"].dropna()
    ax.hist(score, bins=50, color="0.35")
    if len(score):
        ax.axvline(score.median(), color="red", ls="--")
        ax.text(score.median(), ax.get_ylim()[1], f"Median: {score.median():.2f} ", color="red", va="top", ha="right")
    ax.set(xlabel="Score", ylabel="Count", title="Score distribution")
    ax = fig.add_subplot(2, 2, 3)
    for col, colour, name in (("umi_main", "red", "DBSCAN=1"), ("umi_second", "green", "DBSCAN=2")):
        v = np.log10(s[col].dropna().clip(lower=1))
        if len(v) > 2 and np.ptp(v) > 0:
            grid = np.linspace(0, v.max(), 200)
            ax.plot(grid, gaussian_kde(v)(grid), color=colour, label=name)
    ax.legend()
    ax.set(xlabel="log10(UMI)", ylabel="Density", title="DBSCAN UMI Density")
    ax = fig.add_subplot(2, 2, 4)
    both = s.dropna(subset=["umi_main", "umi_second"])
    if len(both):
        h = ax.hist2d(np.log10(both["umi_main"]), np.log10(both["umi_second"].clip(lower=1)), bins=30, cmap=SCORE, cmin=1)
        fig.colorbar(h[3], ax=ax, shrink=0.7, label="count")
    ax.set(xlabel="log10(UMI1)", ylabel="log10(UMI2)", title="DBSCAN UMI Distribution")


def placement_by_cluster(m, fig):
    c, s = m["coords"], spatial(m)
    p = c[c["placed"]]
    title(fig, f"%placed: {100 * len(p) / max(len(s), 1):.2f} ({len(p):,}/{len(s):,})")
    ax = fig.add_axes([0.08, 0.52, 0.84, 0.4])
    outline(ax, m)
    ax.scatter(p["x"], p["y"], s=2, c=cluster_colours(p["rna_cluster"]), rasterized=True)
    puck_axes(ax, m)
    ax.set(xlabel="x-position", ylabel="y-position")
    clusters = sorted(c["rna_cluster"].unique())
    columns = 6
    rows = int(np.ceil(len(clusters) / columns))
    for i, k in enumerate(clusters):
        ax = fig.add_axes([0.02 + (i % columns) * 0.163, 0.45 - (i // columns + 1) * (0.43 / max(rows, 1)),
                           0.15, 0.4 / max(rows, 1)])
        d = p[p["rna_cluster"] == k]
        ax.scatter(d["x"], d["y"], s=1, color=COLOURS[int(k) % len(COLOURS)], rasterized=True)
        puck_axes(ax, m, label=False)
        ax.set_title(f"{k} ({len(d)})", fontsize=7)
        for spine in ax.spines.values():
            spine.set_visible(False)


def placement_by_umis(m, fig):
    c = m["coords"]
    p = c[c["placed"]]
    title(fig, "Placed nuclei coloured by UMIs")
    for i, (values, name, cmap) in enumerate([
            (p["rna_umis"].clip(upper=p["rna_umis"].quantile(0.99)), "RNA UMIs (capped at 99th percentile)", "viridis"),
            (np.log10(p["rna_umis"].clip(lower=1)), "log10 RNA UMIs", "viridis"),
            (np.log10(p["umi"].clip(lower=1)), "log10 spatial (SB) UMIs", "magma"),
            (p["score"], "Dynamic DBSCAN score", SCORE)]):
        ax = fig.add_subplot(2, 2, i + 1)
        outline(ax, m)
        sc = ax.scatter(p["x"], p["y"], s=2, c=values, cmap=cmap, rasterized=True)
        fig.colorbar(sc, ax=ax, shrink=0.6)
        puck_axes(ax, m)
        ax.set_title(name)


def rna_vs_spatial(m, fig):
    s = spatial(m)
    title(fig, "RNA vs. SB metrics")
    ax = fig.add_axes([0.08, 0.55, 0.4, 0.37])
    for placed, colour in ((False, "#F8766D"), (True, "#00BFC4")):
        d = s[s["placed"] == placed]
        ax.scatter(np.log10(d["rna_umis"].clip(lower=1)), np.log10(d["umi"].clip(lower=1)), s=1, color=colour,
                   label=f"placed {placed}", rasterized=True)
    ax.legend(markerscale=5)
    ax.set(xlabel="log10 RNA UMI", ylabel="log10 SB UMI", title="SB UMI vs. RNA UMI")
    ax = fig.add_axes([0.08, 0.08, 0.4, 0.37])
    violin(ax, [np.log10(s.loc[s["clusters"].clip(upper=5) == n, "rna_umis"].clip(lower=1)) for n in range(6)],
           ["0", "1", "2", "3", "4", "5+"])
    ax.set(xlabel="DBSCAN clusters", ylabel="log10 RNA UMI", title="RNA UMI vs. DBSCAN cluster")
    per = s.groupby("rna_cluster").agg(nuclei=("placed", "size"), placed=("placed", "mean")).reset_index()
    rows = [[int(r.rna_cluster), int(r.nuclei), f"{100 * r.placed:.1f}%"] for r in per.itertuples()]
    half = (len(rows) + 1) // 2
    table(fig.add_axes([0.53, 0.1, 0.21, 0.8]), rows[:half], ["cluster", "nuclei", "placed"], font=7)
    table(fig.add_axes([0.76, 0.1, 0.21, 0.8]), rows[half:], ["cluster", "nuclei", "placed"], font=7)


def pipeline_checks(m, fig):
    c, s = m["coords"], spatial(m)
    title(fig, "Pipeline checks", "Flex-specific: section tags, pooling of libraries, depth per GEM well")
    ax = fig.add_subplot(2, 2, 1)
    names, this, other, none = [], [], [], []
    for lib in m["libraries"]:
        f = lib["library"]["read_filters"]
        per = f["reads per section"]
        names.append(f"{lib['gem_well']}\n{lib['sb_index']}")
        this.append(per.get(m["section"], 0) / f["reads"] * 100)
        other.append(sum(v for k, v in per.items() if k != m["section"]) / f["reads"] * 100)
        none.append(f.get("removed: no section tag", 0) / f["reads"] * 100)
    if names:
        ax.bar(names, this, color="tab:blue", label="this section's tag")
        ax.bar(names, other, bottom=this, color="tab:orange", label="other sections' tags")
        ax.bar(names, none, bottom=np.add(this, other), color="tab:red", label="no tag found")
        ax.legend(fontsize=6)
        ax.tick_params(axis="x", labelsize=6)
    ax.set(ylabel="% of library reads", title="Reads by section tag (per spatial library)")
    ax = fig.add_subplot(2, 2, 2)
    distances = []
    for lib in m["libraries"]:
        pooled = s[s["gem_well"] == lib["gem_well"]]
        if sum(other["gem_well"] == lib["gem_well"] for other in m["libraries"]) < 2:
            continue  # this GEM well has one library: pooled = this library
        one = lib["coords"][["cellranger_barcode", "x", "y"]].merge(pooled[["cellranger_barcode", "x", "y"]],
                                                                   on="cellranger_barcode", suffixes=("_lib", "_pooled")).dropna()
        distances.append((lib["sb_index"], np.hypot(one["x_lib"] - one["x_pooled"], one["y_lib"] - one["y_pooled"])))
    if distances:
        for i, (name, d) in enumerate(distances):
            ax.hist(d, bins=40, histtype="step", color=COLOURS[2 * i], label=f"{name}: median {d.median():.1f}")
        ax.legend(fontsize=6)
        ax.set(xlabel="distance single library vs pooled position", ylabel="nuclei", title="Pooled vs single-library positions")
    else:
        message(ax, "Pooled vs single-library positions:\nevery GEM well has one spatial library")
    ax = fig.add_subplot(2, 2, 3)
    wells = list(s.groupby("gem_well"))
    violin(ax, [np.log10(d["umi"].clip(lower=1)) for _, d in wells], [g for g, _ in wells])
    ax.set(ylabel="log10 SB UMI per nucleus", title="Spatial depth per GEM well")
    ax.tick_params(axis="x", labelsize=6)
    ax = fig.add_subplot(2, 2, 4)
    rows = []
    for lib in m["libraries"]:
        step5 = lib["log"]["STEP 5 reads (this section)"]
        nuclei = lib["log"]["STEP 3 nuclei"]
        rows.append([f"{lib['gem_well']}/{lib['sb_index']}", nuclei["nuclei left out (same droplet)"],
                     step5["removed: bead not on this section's puck"], step5["removed: chimeric UMI"]])
    table(ax, rows, ["library", "nuclei left out\n(shared droplet)", "reads: bead not\non this puck", "reads:\nchimeric"],
          "Should be small", font=6.5, col_widths=[0.4, 0.2, 0.2, 0.2])
