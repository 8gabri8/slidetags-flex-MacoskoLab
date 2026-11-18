# slidetags_flex: Slide-tags positioning for 10x Flex

Cell Ranger (Cumulus) runs first, separately. This workflow uses Cell Ranger's results to give every
nucleus a position on its puck, and makes one summary per section (= sample = puck). It is a Python
version of the MacoskoLab slide-tags pipeline (`spatial-count.jl` + `positioning.R`), adapted to Flex.

## Overview

```
Terra workflow slidetags_flex (one launch)
├── list_lines                     which sample-sheet lines have a spatial library, which sections exist
├── run_line   x each line         STEPS 1-7, in parallel     -> <slidetags_output_dir>/<Link>_<sb_index>/
└── summarize  x each section      STEPS 8-9, in parallel,    -> <slidetags_output_dir>/summary/<section>/
                                   after ALL run_line jobs
```

| File | Role |
|---|---|
| `slidetags_flex.py` | **main file**: steps 1-7 per line; `--summarize SECTION` calls the two files below |
| `merge_sections.py` | **step 8**: collects a section across GEM wells and libraries, pools libraries, places once, RNA UMAP/clusters, `files.csv`, `coords.csv` |
| `report_section.py` | **step 9**: the PDF summary |
| `slidetags_flex.wdl` | the Terra workflow |
| `Dockerfile` | the image (Python, scanpy, h5py, gcloud). Building it also tests scanpy and the .h5 reader |
| `.github/workflows/docker.yml` | GitHub builds the image and pushes it to Docker Hub |
| `.dockstore.yml` | tells Dockstore where the WDL is, so it can send it to Terra |
| `examples/Pilot_FFPE_092026/` | the pilot's 3 input files, named as in the bucket (`gs://…/flex/cumulus_cellranger_inputs/Pilot_FFPE_092026/`): sample sheet, Flex sample file (AuxFile), and one Terra inputs JSON for both workflows (Cumulus `cellranger_workflow` and `slidetags_flex`) |
| `tests/` | simulator and checks |

## Steps 1-7: each line of the sample sheet

Each line is one spatial library (`sb_index`) and is processed on its own.

| Step | What it does |
|---|---|
| **1. Download** | The AuxFile (its `sample_id`s are the sections), each section's `Puck.csv`, each section's Cell Ranger `barcodes.tsv.gz`, and the spatial FASTQs `<Flowcell>/<sb_index>_*_R1/R2_001.fastq.gz`. |
| **2. Load the pucks** | `Puck.csv` = barcode, x, y. Removes beads whose barcode has an N or is degenerate (one base > 4/5 of the barcode, or a run > 5/7). A barcode found on two pucks is **kept**: the read's section tag says which puck it came from. |
| **3. Load the nuclei and the section tags** | A Cell Ranger barcode is the GEM barcode (16 bases) + the Flex probe barcode (10 bases) + `-1`, e.g. `AAACCAATCATGGAAC` + `ACCTGGTTGT` + `-1`. Spatial read 1 carries the GEM barcode. Two nuclei of the same section in the same droplet can't be told apart, so both are left out. **Section tag:** the AuxFile's `spatial_tag` if filled in, otherwise the section's probe barcode. |
| **4. Read the spatial reads** | From the first 100,000 reads it finds the bead type, which read holds the bead, and where the section tags sit (pilot: read 2, bases 66–75, as written). Every read pair then goes through 5 filters: **(1)** too short; **(2)** UMI with an N, or one base repeated; **(3)** no UP site (1 mismatch per 6 bases allowed); **(4)** spatial barcode matches no bead, or 2+ beads (1 mismatch allowed); **(5)** section tag matches no section (1 mismatch allowed). Each lane (R1/R2 file pair) is read by its own process, in parallel (up to `num_cpu`). The reads kept are counted per (section, cell barcode, UMI, bead). It also records the numbers for the summary: how each part matched (exact / 1 mismatch / …), where the bead-barcode mismatches are, reads per UMI, and a downsampling curve. |
| *then for each section:* | |
| **5. Spatial UMIs per nucleus** | Uses only the reads with **this section's tag**. **(a)** Check: removes reads whose bead isn't on this section's puck (should be ~0; many means a wrong tag). **(b)** Matches the cell barcode to a nucleus of this section: exact, or the only nucleus 1 mismatch away. The 10x bases 8-9 swap is used if it explains more reads. **(c)** Removes reads that aren't from a nucleus (empty droplets, ambient). **(d)** Removes chimeras: a UMI comes from one bead, so it keeps the bead with strictly the most reads. **(e)** Counts UMIs per (nucleus, bead) and adds the bead's x, y. |
| **6. Place the nuclei (DBSCAN)** | **(a)** `eps` = median distance from a bead to its 36th nearest bead. **(b)** For each nucleus, the range of `min_pts` giving exactly one cluster (clusters closer than eps × cmes count as one). **(c)** One `min_pts`: the value inside the most nuclei's ranges. **(d)** Final DBSCAN: a nucleus is placed at its cluster's UMI-weighted centre if it has exactly one cluster. **(e)** Dynamic DBSCAN: each nucleus also gets a position with its *own* `min_pts` and a score (1 = one clear cluster, 0 = two equal ones). |
| **7. Save** | To `<slidetags_output_dir>/<Link>_<sb_index>/` (see Results). |

### The section tag
The spatial oligo and the 10x gel-bead oligo are joined by an adapter that carries a tag for the
section, so every spatial read says which section it comes from. A read from the pilot (GEM well 1, read 2):

```
GNCCGCCCGACCTGGAT CTGTTTCCTG CCAATTTCC CCCATATAAGAAAACCTGAATACGCGGTT ACCTGGTTGT CGGTCCTAGCAA...
bead barcode 1-17 UP 18-27   random    adapter 37-65                 tag 66-75  adapter
```

Because of the tag:
- **One nucleus from each section in a droplet:** resolved, because they have the same GEM barcode but different tags.
- **A bead barcode found on two pucks:** resolved too.
- **Several sections on one puck:** list the same `Puck.csv` for each of them in `pucks`; each section is placed using only its own reads.
- **Two nuclei of the *same* section in one droplet:** these still can't be resolved and are left out.

**Custom tags:** add a column `spatial_tag` to the Flex AuxFile, filled only for the sections whose
tag is not their probe barcode. All tags must have the same length.

```
sample_id,probe_barcode_ids,description,emptydrops_minimum_umis,spatial_tag
brain_section1,A-D04,coronal_section_1,100,
brain_section2,A-E04,coronal_section_2,100,TTGCAGGACA
```

Cumulus copies the AuxFile into Cell Ranger's `[samples]` section as it is, and Cell Ranger may not
accept an extra column. If it complains, give Cell Ranger a copy without `spatial_tag`. If every tag
is its probe barcode (as in the pilot), leave the column out.

## Steps 8-9: one summary per section

A section's nuclei can be in several GEM wells (lines with different `Link`), and one GEM well's
spatial reads can be in several libraries (lines with the same `Link` and different `sb_index`).

**8. Merge** (`merge_sections.py`)
- **8a.** It finds the lines of the section: those whose AuxFile lists the section's `sample_id`.
- **8b.** It downloads two things:
  - **per GEM well, from Cell Ranger:** the RNA matrix, `metrics_summary.csv` and the unfiltered matrix;
  - **per spatial library, from steps 1-7:** the results.
- **8c. Pool the libraries of each GEM well, then place each nucleus once.**
  - **Pooling:** libraries made from the same GEM well come from the same cDNA, so they contain the same nuclei. Their reads of each (nucleus, UMI, bead) are added up, so a molecule seen in two libraries counts once. Chimeras are then removed again and UMIs counted.
  - **Placement:** DBSCAN places each nucleus once, with more UMIs than any single library gives. A GEM well with one library gets exactly the positions of its line.
  - **`min_pts`:** chosen **per GEM well**, as the R pipeline does per library. GEM wells can differ a lot in spatial depth, and one `min_pts` for all would place fewer nuclei in the shallower ones. `eps` comes from the puck, so it is the same for all.
- **8d. RNA.**
  - **Per nucleus:** UMIs, genes and %mito. %mito is computed only if the probe set has `mt-` probes; otherwise the PDF says "not measured".
  - **Clustering:** scanpy on all GEM wells together, with no batch correction. Normalize to 10,000 per nucleus, log, take the 2,000 most variable genes, then PCA (30), 15 neighbours, UMAP and Leiden (resolution 1.0).
- **8e.** It writes the three tables below.

**9. PDF** (`report_section.py`): 16 pages, plus one per extra spatial library. If one page fails,
it shows the error and the rest are still made.

| Page | Content |
|---|---|
| 1 Overview | Section, puck, tag, GEM wells and libraries, headline numbers, what happened to each nucleus, map of placed nuclei |
| 2 QC per library | One column per GEM well × library (+ pooled): nuclei, median UMIs and genes, %mito, spatial reads, tagged reads, % in nuclei, spatial UMIs, eps, min_pts, % placed |
| 3 GEX metrics | Cell Ranger `metrics_summary.csv`, one column per GEM well. Genome, intron and exon mapping metrics are dropped (Flex reads come from probes) |
| 4 Cell calling | Barcode rank plot and genes vs UMIs from the unfiltered matrix, nuclei vs background |
| 5 RNA | UMAP by cluster and by GEM well, UMIs per cluster, %mito, genes, share of each GEM well in each cluster |
| 6 Spatial metrics (one per library) | Library and puck tables; UP, bead, tag and cell-barcode matching; % of reads per filter; where the bead-barcode mismatches are; FASTQs |
| 7 Spatial counts | Spatial UMIs per cell barcode and per bead (rank plots), downsampling curve, reads per UMI and saturation |
| 8 Puck | Spatial UMIs on the puck: all barcodes / nuclei only |
| 9 Global DBSCAN | Clusters per nucleus, % placed vs min_pts (one curve per GEM well), UMIs vs clusters, placement interval |
| 10 DBSCAN=1 | Positions, SNR (UMIs in the cluster / all), spread (RMSD), h-index |
| 11 DBSCAN=2 | Positions of both clusters, their UMIs, their distance vs eps × cmes |
| 12 Dynamic DBSCAN | Positions with each nucleus' own min_pts coloured by score, score histogram, UMIs of cluster 1 vs 2 |
| 13 Placement by cluster | % placed (n/N), map by RNA cluster, one small map per cluster |
| 14 Placement by UMIs | Maps coloured by RNA UMIs, log10 RNA UMIs, log10 spatial UMIs, score |
| 15 RNA vs spatial | RNA vs spatial UMIs (placed / not), RNA UMIs vs DBSCAN clusters, % placed per RNA cluster |
| 16 Pipeline checks | Reads by section tag per library, pooled vs single-library positions, spatial depth per GEM well, numbers that should be small |

## Results

```
<slidetags_output_dir>/
├── <Link>_<sb_index>/                 one folder per sample-sheet line (steps 1-7)
│   ├── summary.csv                    per section: tag, spatial reads, nuclei, eps, min_pts, % placed
│   ├── library.json                   whole-library numbers (filters, matching, reads per UMI, downsampling)
│   └── <section>/
│       ├── coords.csv                 one row per nucleus: position from THIS library alone
│       ├── umis.csv.gz                spatial UMIs per (nucleus, bead) with the bead's x, y
│       ├── log.json                   reads removed at every filter, eps, min_pts, % placed
│       ├── nucleus_reads.npz          reads of the nuclei (used to pool libraries in step 8)
│       ├── spatial_qc.npz             UMIs per cell barcode and per bead (for the PDF)
│       └── summary.png                quick look
└── summary/<section>/                 one folder per section (steps 8-9): USE THESE
    ├── files.csv                      where every file of the section is (see below)
    ├── coords.csv                     every nucleus of the section, all GEM wells (see below)
    ├── qc_per_library.csv             page 2 of the PDF as a table
    └── <section>_summary.pdf
```

**`summary/<section>/files.csv`**: one row per GEM well × spatial library. A GEM well without a
spatial library gets one row with the spatial columns empty. The columns hold every path needed to
build the section's object:
- **Cell Ranger:** `cellranger_sample_dir`, `filtered_matrix_h5`, `filtered_matrix_mex_dir`, `raw_matrix_h5`, `molecule_info_h5`, `metrics_summary_csv`, `web_summary_html`, `cloupe`, `cellranger_probe_set_csv`
- **Slide-tags:** `spatial_results_dir`, `library_json`, `spatial_coords_csv`, `spatial_umis_csv`, `spatial_log_json`, `nucleus_reads_npz`, `spatial_qc_npz`, `puck_csv`, `spatial_fastqs`
- **Section summary:** `section_coords_csv`, `section_summary_pdf`

**`summary/<section>/coords.csv`**: one row per nucleus Cell Ranger called in the section, across all GEM wells. Join it to the RNA on `gem_well` + `cellranger_barcode`.
- `spatial_status`: placed / not placed (no cluster) / not placed (2+ clusters) / left out (droplet shared within the section) / RNA only (no spatial library)
- position (pooled libraries): `x`, `y`, `placed`, `clusters`, `spatial_libraries`
- spatial numbers: `umi`, `beads`; for cluster 1: `x1`, `y1`, `umi1`, `beads1`, `rmsd1`, `h1`, `snr`; for cluster 2: `x2`, `y2`, `umi2`
- Dynamic DBSCAN: `x_dynamic`, `y_dynamic`, `score`, `umi_main`, `umi_second`, `i1`, `i2`
- placement settings: `eps`, `min_pts`
- RNA: `rna_umis`, `rna_genes`, `pct_mito`, `rna_cluster`, `umap_1`, `umap_2`

## The sample sheet = your Cumulus sheet + 4 columns

| Column | Example | Used for |
|---|---|---|
| `Link` | `brain_GEMwell1` | the GEM well; Cell Ranger results are in `<cellranger_output_dir>/<Link>` |
| `Flowcell` | `gs://…/fastqs/260903_SL-EXC_0924_A2535KKLT4_lane1-6` | the spatial FASTQs are taken from here (directly in this folder) |
| `AuxFile` | `gs://…/Pilot_FFPE_092026-flex_samples.csv` | its `sample_id`s are the sections; optional `spatial_tag` column |
| **`sb_index`** | `SI-TT-A6` | FASTQ prefix of the spatial library (empty = RNA only, no spatial library) |
| **`pucks`** | `brain_section1=gs://…/Puck.csv;brain_section2=gs://…/Puck.csv` | each section's puck |
| **`cellranger_output_dir`** | `gs://…/flex/outputs/Pilot_FFPE_092026` | the Cumulus `output_directory` |
| **`slidetags_output_dir`** | `gs://…/flex/outputs/Pilot_FFPE_092026/flex_slidetag_out` | where all results go |

A second spatial library from the same GEM well is one more line, with the same `Link` and the new
`sb_index`. Cumulus would read that line as a duplicate gene-expression run, so keep it out of the
copy of the sheet you give Cumulus.

## Settings

| Setting | Default | Meaning |
|---|---|---|
| `only_summarize` | false | true: skip `run_line` and rebuild only the summaries from the results already in the bucket (e.g. to redo the PDFs) |
| `barcode_length` | `auto` | Cell Ranger barcode length without `-1` (26 for multiplexed Flex v2). The run stops if it differs |
| `knn` | 36 | eps = median distance to the knn-th nearest bead |
| `cmes` | 10.0 | clusters within eps × cmes of the main one count as one (lab default) |
| `num_cpu`, `memory`, `disk_space`, `preemptible`, `zones` | 8, 32G, 100, 1, us-central1-a/b/c/f | `run_line` VM (same names as in Cumulus). Lanes are read in parallel, so `num_cpu` = number of lanes is the fastest. The disk must hold one line's spatial FASTQs |
| `summarize_num_cpu`, `summarize_memory`, `summarize_disk_space` | 8, 64G, 100 | `summarize` VM. It reads all RNA matrices (filtered and unfiltered) of one section |

The sample sheet input is called `input_csv_file`, as in Cumulus: give both workflows the same file.

## Running

**On Terra: one-time setup, all in the browser.** The path is GitHub (code) → Docker Hub (image) + Dockstore (workflow) → Terra.

1. **GitHub repository.** Create a public repository and upload the contents of this folder, so that `slidetags_flex.wdl` sits at the top. Drag-and-drop often skips hidden files, so create `.dockstore.yml` and `.github/workflows/docker.yml` with *Add file → Create new file* and paste their contents.
2. **Docker image.**
   1. Create a Docker Hub access token.
   2. In GitHub, go to *Settings → Secrets and variables → Actions* and add `DOCKERHUB_USERNAME` and `DOCKERHUB_TOKEN`.
   3. Go to *Actions → Build Docker image → Run workflow* and type a tag, e.g. `v3.0.0`.
3. **Dockstore.**
   1. Register with GitHub.
   2. Go to *My Dockstore → + → Register using GitHub Apps* and install the app on this repository.
   3. Publish the workflow.
4. **Terra.**
   1. On Dockstore, click *Launch with → Terra*.
   2. In Terra, choose *inputs defined by file paths*.
   3. Upload `examples/Pilot_FFPE_092026/Pilot_FFPE_092026-terra_inputs.json`, set `docker`, and click *Save*.

**Each run.** Upload the sample sheet, set `input_csv_file`, and click *Launch*. *Job History* shows
`list_lines`, the `run_line` jobs, then the `summarize` jobs. The PDFs are also listed as workflow
outputs.

**Updating the code.**
1. Replace the changed files on GitHub.
2. If a `.py` file or the `Dockerfile` changed, run *Actions → Build Docker image* again: with the same tag to overwrite the image, or with a new tag (then set `docker` to it).
3. Dockstore syncs the WDL by itself. In Terra, use the `main` version.

**Files outside the workspace bucket** (e.g. pucks in another lab bucket): that bucket must give
read access to your Terra *proxy group*, which is shown on your Terra profile page.

**Anywhere else.**
- `python3 slidetags_flex.py SHEET.csv` runs every line, then every section.
- `--line N` runs one line; `--summarize SECTION` makes one section's summary.

## Tests

```
python3 tests/simulate.py SIM            # pilot-like experiment with known positions: 3 GEM wells (one RNA only),
                                         # 4 lines (2 spatial libraries from GEM well 1), 2 sections, a custom tag
python3 slidetags_flex.py SIM/sample_sheet.csv --workdir SIM/work
python3 tests/check_results.py SIM       # per-line and per-section placements vs truth, RNA-only, files.csv, PDF
python3 tests/test_dbscan.py SIM/results/brain_*/*/umis.csv.gz   # DBSCAN vs scikit-learn; i1/i2 vs the R search loops
python3 tests/test_matching.py           # bead matching + UMI numbers vs the Julia code
```

Without scanpy installed, add `PYTHONPATH=tests/stand_ins` to the commands. It uses simple stand-ins
(sklearn PCA + k-means) instead of UMAP and Leiden, for testing only. The Docker image uses the real
scanpy, and building the image runs `merge_sections.py --self-test`, which checks scanpy and the
.h5 reader.
