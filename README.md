# slidetags_flex: Slide-tags positioning for 10x Flex

This workflow gives every nucleus that Cell Ranger called a position on its puck. It is a Python
version of the MacoskoLab slide-tags pipeline (`spatial-count.jl` + `positioning.R`), adapted to Flex.

**All the code is in one file, `slidetags_flex.py`.** Read it top to bottom: `run_line()` is the
pipeline, and each function below it is one step.

## What happens, for each line of the sample sheet

Each line is one spatial library (`sb_index`) and is processed on its own. Several lines can share
the same GEM well (`Link`).

| Step | What it does |
|---|---|
| **1. Download** | Downloads the AuxFile (its `sample_id`s are the sections), each section's `Puck.csv`, and the spatial FASTQs `<Flowcell>/<sb_index>_*_R1/R2_001.fastq.gz`. |
| **2. Load the pucks** | `Puck.csv` = barcode, x, y. Removes beads whose barcode has an N or is degenerate (one base > 4/5 of the barcode, or a run > 5/7). A barcode that is on both pucks is **kept** (see below). |
| **3. Read the spatial reads** | Finds the bead type and which read holds the bead from the first 100,000 reads. Then every read pair goes through 4 filters: **(1)** read too short; **(2)** UMI with N or one repeated base; **(3)** no UP site (1 mismatch allowed per 6 bases); **(4)** spatial barcode matches no bead, or 2+ beads, allowing 1 mismatch. The reads kept are counted per (cell barcode, UMI, bead). |
| *then for each section:* | |
| **4. Load the nuclei** | Reads `barcodes.tsv.gz` from `<cellranger_output_dir>/<Link>/per_sample_outs/<section>/sample_filtered_feature_bc_matrix/`. Detects the barcode length (24 for multiplexed Flex), or checks the one you give. Keeps the first 16 bases, the GEM barcode the spatial reads carry. Two nuclei of the same section in the same droplet can't be told apart, so both are left out. |
| **5. Spatial UMIs per nucleus** | **(a)** Keeps reads whose bead is on *this section's* puck. **(b)** Matches the cell barcode to a nucleus: exact, or the only nucleus 1 mismatch away. The 10x bases 8-9 swap is used if it explains more reads. **(c)** Removes reads that aren't from a nucleus of this section. **(d)** Removes chimeras: a UMI comes from one bead, so it keeps the bead with strictly the most reads. **(e)** Counts UMIs per (nucleus, bead) and adds the bead's x, y. |
| **6. Place the nuclei (DBSCAN)** | **(a)** `eps` = median distance from a bead to its 36th nearest bead. **(b)** For each nucleus, finds the range of `min_pts` giving exactly one cluster (clusters closer than eps × cmes count as one). **(c)** One `min_pts` for the section: the value inside the most ranges. **(d)** Final DBSCAN. A nucleus is placed at its cluster's UMI-weighted centre if it has exactly one cluster. It also gets a score (1 = one clear cluster). |
| **7. Save** | Writes the result files and copies them to `<slidetags_output_dir>/<Link>_<sb_index>/<section>/`. |

### Why a barcode on two pucks is no longer removed
The R pipeline drops beads whose barcode is on two pucks, because it can't tell which puck a read
came from. With Flex we know each nucleus's section from its Flex probe barcode. In step 5 a nucleus
only looks at its own section's puck, so the duplicated barcode has only one possible position.

The same holds for a droplet that holds one nucleus from each section. Its two nuclei have the same
GEM barcode, but each is placed on its own puck. Only two nuclei of the *same* section in one droplet
can't be resolved.

## The sample sheet = your Cumulus sheet + 4 columns

Use the same CSV for Cumulus and for this workflow; Cumulus ignores the extra columns.

| Column | Example | Used for |
|---|---|---|
| `Link` | `brain_GEMwell1` | the Cell Ranger results folder `<cellranger_output_dir>/<Link>` |
| `Flowcell` | `gs://…/fastqs/260903_sl-exc_0924_a2535kklt4_lane1` | the spatial FASTQs are taken from here |
| `AuxFile` | `gs://…/Pilot_FFPE_092026-flex_samples.csv` | its `sample_id`s are the sections |
| **`sb_index`** | `SI-TT-A6` | FASTQ prefix of the spatial library (empty = skip the line) |
| **`pucks`** | `brain_section1=gs://…/Puck.csv;brain_section2=gs://…/Puck.csv` | each section's puck |
| **`cellranger_output_dir`** | `gs://…/flex/outputs/Pilot_FFPE_092026` | the Cumulus `output_directory` |
| **`slidetags_output_dir`** | `gs://…/slide-tags-flex/Pilot_FFPE_092026` | results go to `/<Link>_<sb_index>/` |

To add a second spatial library from the same GEM well, add a line with the same `Link` and the new
`sb_index`.

## Results: `<slidetags_output_dir>/<Link>_<sb_index>/`

```
summary.csv                 one line per section: nuclei called, left out, eps, min_pts, % placed
<section>/coords.csv        one row per nucleus: cellranger_barcode, x, y, placed, clusters, umi, beads,
                            spread, umi_main, umi_second, score   (join to the RNA on cellranger_barcode)
<section>/umis.csv.gz       spatial UMIs per (nucleus, bead) with the bead's x, y
<section>/log.json          reads removed by every filter, eps, min_pts, % placed
<section>/summary.png       read filters, UMIs per nucleus, map of placed nuclei
```

## Settings

| Setting | Default | Meaning |
|---|---|---|
| `barcode_length` | `auto` | Cell Ranger barcode length without `-1` (16, 24, …). The run stops if it differs |
| `knn` | 36 | eps = median distance to the knn-th nearest bead |
| `cmes` | 10.0 | clusters within eps × cmes of the main one count as one (lab default) |
| `cpu`, `memory`, `disk_gb` | 8, 32G, 200 | VM per line. The disk must hold that line's spatial FASTQs |

## Files

| File | Role |
|---|---|
| `slidetags_flex.py` | the whole pipeline |
| `slidetags_flex.wdl` | Terra: `list_lines` (which lines have an `sb_index`), then `run_line` on every line in parallel |
| `Dockerfile` | the image the workflow runs in |
| `.github/workflows/docker.yml` | GitHub builds the image and pushes it to Docker Hub when you publish a release |
| `.dockstore.yml` | tells Dockstore where the WDL is, so it can send it to Terra |
| `examples/` | your sample sheet and Terra inputs |
| `tests/` | simulator and checks |

## Running

**On Terra (one-time setup, all in the browser)**

The path is GitHub (code) → Docker Hub (image) + Dockstore (workflow) → Terra.

1. **GitHub repository.** Create a new public repository, e.g. `slidetags-flex`. Upload the contents
   of this folder so that `slidetags_flex.wdl` sits at the top. Hidden files are often skipped by
   drag-and-drop, so create `.dockstore.yml` and `.github/workflows/docker.yml` with
   *Add file → Create new file* and paste their contents.
2. **Docker image.**
   1. Make a free Docker Hub account, then create an access token (Read & Write).
   2. In GitHub, go to *Settings → Secrets and variables → Actions* and add `DOCKERHUB_USERNAME` and `DOCKERHUB_TOKEN`.
   3. Go to *Releases → Draft a new release*, create tag `v2.0.0`, and publish it.
   4. Under *Actions* you should see "Build Docker image" run. It pushes `<user>/slidetags-flex:v2.0.0` (takes ~5 min).
3. **Dockstore.**
   1. At dockstore.org, register with GitHub.
   2. Go to *My Dockstore → + → Register using GitHub Apps* and install the app on this repository only.
   3. Open the workflow, go to *Versions → set default version*, and click *Publish*.
4. **Terra.**
   1. On the Dockstore public page, under *Launch with*, click **Terra** and pick your workspace.
   2. In the workspace, open the workflow and choose *Run workflow with inputs defined by file paths*.
   3. Upload `examples/terra_inputs.json` and set `docker` to `<user>/slidetags-flex:v2.0.0`.
   4. Save.

**Each run**
1. Fill in the slide-tags columns of your sample sheet (`sb_index`, `pucks`, `cellranger_output_dir`,
   `slidetags_output_dir`) and upload it to the workspace bucket.
2. Set `sample_sheet` to its gs:// path, then *Run analysis → Launch*.
3. *Job History* shows one `run_line` job per line. Results appear in `slidetags_output_dir`.

**When you change the code:** push it to GitHub, publish a new release (e.g. `v2.0.1`), and set
`docker` to the new tag. Dockstore and Terra pick up the new WDL by themselves; choose the version in the
workflow's dropdown in Terra.

**Files outside the workspace bucket** (e.g. pucks in another lab bucket): the bucket must give read
access to your Terra *proxy group* (shown on your Terra profile page).

**Anywhere else**
`python3 slidetags_flex.py SHEET.csv [--line 2]`

## Tests

```
python3 tests/simulate.py SIM                     # a pilot-like experiment with known positions (3 lines, 2 from one GEM well)
python3 slidetags_flex.py SIM/sample_sheet.csv --workdir SIM/work
python3 tests/check_results.py SIM                # placements vs truth
python3 tests/test_dbscan.py SIM/results/*/*/umis.csv.gz   # DBSCAN vs scikit-learn; i1/i2 vs the R search loops
python3 tests/test_matching.py                    # bead matching + UMI numbers vs the Julia code
```
