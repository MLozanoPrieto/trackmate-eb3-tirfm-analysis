# TrackMate EB3 TIRFM Analysis Pipeline

**Purpose.** This repository provides a reproducible, post‑processing pipeline to summarize microtubule plus‑end dynamics from ImageJ/Fiji **TrackMate** exports of EB3‑GFP TIRF movies in primary human CD4 T cells. It:
1) parses TrackMate CSVs per cell,  
2) computes biologically meaningful **per‑cell metrics** (robust to outliers and track length),  
3) optionally runs **non‑parametric statistics** (Mann–Whitney; Van Elteren‑like stratified by Experiment) across conditions, and  
4) generates **publication‑quality plots**.

The pipeline assumes EB1/EB3 comet tracking (LoG detector + LAP tracker in TrackMate) and consistent spatial/temporal calibration. It does **not** perform detection or tracking itself—only analysis of TrackMate outputs.

---

## Table of contents

- [Quick start](#quick-start)
- [Inputs](#inputs)
  - [Expected folder layout](#expected-folder-layout)
  - [Required filenames](#required-filenames)
  - [CSV header parsing & multi-line headers](#csv-header-parsing--multi-line-headers)
  - [TrackMate exports must be calibrated](#trackmate-exports-must-be-calibrated)
- [Outputs](#outputs)
  - [`1_METRICS/cell_metrics.csv` (columns)](#1_metricscell_metricscsv-columns)
  - [`1_METRICS/track_metrics_from_edges.csv`](#1_metricstrack_metrics_from_edgescsv)
  - [`2_STATS/mannwhitney_all_metrics.csv` (if `-s`)](#2_statsmannwhitney_all_metricscsv-if--s)
  - [`3_PLOTS` (if `-p`)](#3_plots-if--p)
- [Biological meaning of metrics](#biological-meaning-of-metrics)
  - [Tracks (most informative for growth dynamics)](#tracks-most-informative-for-growth-dynamics)
  - [Edges (instantaneous, useful for distribution tails)](#edges-instantaneous-useful-for-distribution-tails)
- [Why medians and two-stage aggregation?](#why-medians-and-twostage-aggregation)
- [Installation](#installation)
- [Usage details](#usage-details)
  - [Command-line options](#commandline-options)
  - [What the script does not do](#what-the-script-does-not-do)
- [Reproducibility & logging](#reproducibility--logging)
- [Demo dataset (quick test)](#demo-dataset-quick-test)
- [Troubleshooting](#troubleshooting)
- [Citation](#citation)

---

## Quick start

```bash
# 1) Metrics only
python trackmate_eb3_tirfm_analysis.py -i /path/to/root

# 2) Metrics + statistics (exactly 2 conditions required)
python trackmate_eb3_tirfm_analysis.py -i /path/to/root -s

# 3) Metrics + plots (with/without stats)
python trackmate_eb3_tirfm_analysis.py -i /path/to/root -p

# 4) Full pipeline: metrics → stats → plots, writing outputs elsewhere
python trackmate_eb3_tirfm_analysis.py -i /path/to/input -o /path/to/output -s -p
```

A tee log is saved to `logs/trackmate_pipeline_<timestamp>.log` under the chosen output base.

---

## Inputs

### Expected folder layout
The script walks all subfolders under the input **root** and groups CSVs **per cell** using path metadata:
```
<root>/Experiment/Condition/[Plate]/cellName-edges.csv
<root>/Experiment/Condition/[Plate]/cellName-tracks.csv
```
*`Plate` is optional, but `Experiment` and `Condition` subfolders are mandatory, even if there is only one of each. Any extra nesting above the files is ignored for metadata. Make sure all files from one cell  have the same basename and they are unique to one cell/video*

### Required filenames
The file kind is inferred by suffix:
- `*-edges.csv` → **edges**
- `*-tracks.csv` → **tracks**
- Files that match `*all-spots*.csv` are **ignored**.

*Make sure suffixes match exactly*

### CSV header parsing & multi-line headers

This pipeline is tolerant to **multi-line headers** and varied CSV encodings. The reader:

- **Scans the first 25 lines** of each CSV to guess where the real header starts.
- **Auto-detects the delimiter** (comma, semicolon, tab, …).
- Tries common encodings (**UTF-8**, **UTF-8-SIG**, **Latin-1**).
- Matches required columns **case-insensitively** and ignores units in parentheses and extra symbols (e.g., `Track duration (s)` → recognized as `TRACK_DURATION`).

As long as the final header row includes the expected metrics, the file will be parsed correctly—even if there are introductory lines above.

### TrackMate exports must be calibrated
Ensure **pixel size** (µm/pixel) and **frame interval** (s) are correctly set before tracking. The pipeline trusts the numeric values written to CSV by TrackMate; it does not convert units.

---

## Outputs

All outputs are written under the chosen base (default = input root).

```
1_METRICS/
  ├─ cell_metrics.csv                   # one row per cell (final table)
  └─ track_metrics_from_edges.csv       # per‑track summaries derived from edges

2_STATS/                                 # only if -s
  └─ mannwhitney_all_metrics.csv

3_PLOTS/                                 # only if -p
  ├─ BY_CONDITION/
  │   ├─ UNCORRECTED/                    # p (pooled or stratified primary p)
  │   └─ FDR_CORRECTED/                  # q (BH‑FDR on primary p)
  └─ BY_EXPERIMENT/                      # per‑experiment box+strip
logs/
  └─ trackmate_pipeline_<timestamp>.log
```

### `1_METRICS/cell_metrics.csv` (columns)
Metadata:
- `Experiment`, `Condition`, `Plate`, `Cell`

Per‑cell metrics from **tracks** (medians unless noted):
- `N_TRACKS` — number of **unique tracks** in the cell (unique `TRACK_ID` / `TRACK_INDEX` in *tracks.csv*).
- `TRACK_MEDIAN_SPEED_median` — median of each track’s median speed, then median across tracks per cell.  
- `TRACK_DURATION_median`, `TRACK_DURATION_p90` — track lifetimes (continuous growth).  
- `CATASTROPHE_RATE` — **approximate** disappearance rate: `1 / mean(TRACK_DURATION)` (use only if few tracks end at the last video frame; see caveats).  
- `TOTAL_DISTANCE_TRAVELED_median` — run length of growth (sum of step lengths).  
- `TRACK_DISPLACEMENT_median` — net displacement (start→end, straight line).  
- `CONFINEMENT_RATIO_median` — rectilinearity = displacement / total distance (1 = straight).  
- `LINEARITY_OF_FORWARD_PROGRESSION_median` — ratio of straight‑line speed to mean speed.  
- `MEAN_STRAIGHT_LINE_SPEED_median` — straight‑line (vectorial) speed.  
- `MEAN_DIRECTIONAL_CHANGE_RATE_median`, `MEAN_DIRECTIONAL_CHANGE_RATE_p90` — tortuosity (angular change rate).

Per‑cell metrics **from edges** (via per‑track summaries; medians across tracks):
- `SPEED_TRACK_MEDIAN_median`, `SPEED_TRACK_P10_median`, `SPEED_TRACK_P90_median`, `SPEED_TRACK_IQR_median`.  
- `DIRECTIONAL_CHANGE_RATE_TRACK_MEDIAN_median`, `_P10_median`, `_P90_median`, `_IQR_median`.

Cells with all‑empty or all‑zero metric values are omitted with a warning.

### `1_METRICS/track_metrics_from_edges.csv`
One row per **track**, containing per‑track distribution summaries derived from edges:
- `SPEED_median`, `SPEED_p10`, `SPEED_p90`, `SPEED_IQR`  
- `DIRECTIONAL_CHANGE_RATE_median`, `_p10`, `_p90`, `_IQR`  
plus the metadata columns.

### `2_STATS/mannwhitney_all_metrics.csv` (if `-s`)
For **exactly two** conditions, the script computes:
- **Primary p‑value** per metric: **Van Elteren‑like stratified** Mann–Whitney across `Experiment` (if ≥2 usable strata), otherwise **pooled** Mann–Whitney.  
- Multiple testing: **BH‑FDR** → `q_value_fdr_bh`.  
- Effect sizes: **CLES** and **Cliff’s delta**.  
- Descriptives: means, medians, sample sizes, and median‑based direction.

If the dataset does not contain exactly two conditions, stats are skipped with a console warning (the pipeline continues).

### `3_PLOTS` (if `-p`)
- **BY_CONDITION**: box + jittered dots, sample sizes in x‑tick labels. If stats exist, two versions are saved: **UNCORRECTED** (annotates `*` at p<0.05) and **FDR_CORRECTED** (annotates `*` at q<0.05); otherwise a single unannotated plot is saved.  
- **BY_EXPERIMENT**: the same per experiment (no significance annotations).

You can freely adjust the plotting style. Tweak the constants below in the plotting section of the script and re-run with `-p` to regenerate all figures.

```python
# Style
FIGSIZE = (4, 4.8)          # width, height in inches
DPI = 600                   # export resolution (pixels per inch)
PALETTE_NAME = "Set2"       # seaborn palette (e.g., "Set2", "colorblind", "deep", "viridis")
DOT_ALPHA = 0.65            # dot transparency for jittered points
DOT_SIZE = 3.6              # dot size (points)
JITTER = 0.18               # horizontal jitter amplitude
BOX_LINEWIDTH = 1.4         # box/whisker line width
MEDIAN_LINEWIDTH = 2.0      # thickness of the median line in the boxplot
LABEL_FONTSIZE = 11.5       # axis-label font size
TICK_FONTSIZE = 10.5        # tick-label font size
STAR_FONTSIZE = 16          # significance star (“*”) font size
```

---

## Biological meaning of metrics

Below, “per cell” always means aggregation **across tracks in that cell**.

### Tracks (most informative for growth dynamics)

- **`N_TRACKS` (per cell):**  
  *What it is.* Count of EB3-GFP growth events detected as **unique tracks** within the acquisition window.  
  *Why it matters.* Approximates the **frequency of plus-end growth initiation**; informative about nucleation when interpreted together with **speed** and **duration**.  
  *Caveats.* Sensitive to **track fragmentation** (linking/gap settings), **quality thresholds**, **frame interval**, **video duration**, and **EB3-GFP expression**. Compare only across datasets with **matched** parameters and acquisition settings.

- **`TRACK_MEDIAN_SPEED` (per‑cell: median of medians; also see edges percentiles):**  
  *What it is.* The typical **growth speed** of an EB3 comet along each track.  
  *Why it matters.* Reflects the **tubulin polymerization rate** at the plus‑end.  
  *Why median.* Comet speeds have **skew** and occasional outliers (false links, pauses); the median is **more robust** than the mean and is less biased by **track length**.

- **`TRACK_STD_SPEED` (per‑cell: median)**  
  *What.* Within‑track **variability** in speed.  
  *Why.* Captures intra‑track heterogeneity (e.g., intermittent pauses/accelerations).

- **`TRACK_DURATION` (per‑cell: median and P90):**  
  *What.* Continuous **growth lifetime** of a comet (until disappearance).  
  *Why.* Longer lifetimes indicate sustained polymerization; P90 highlights the **longest runs**.  
  *Derived (use with caution):* **Catastrophe rate ≈ 1 / mean(TRACK_DURATION)** per cell. Valid only when **few tracks are right‑censored** at the video end.

- **`TOTAL_DISTANCE_TRAVELED` & `TRACK_DISPLACEMENT` (per‑cell: medians):**  
  *What.* Total path length vs. **net** start→end displacement.  
  *Why.* Distinguishes **run length** from straightness; large distance with small displacement indicates tortuous growth.

- **`CONFINEMENT_RATIO` (displacement/total distance; per‑cell: median):**  
  *What.* **Rectilinearity** (1 = straight, 0 = highly tortuous).  
  *Why.* Persistence of plus‑end growth direction.

- **`LINEARITY_OF_FORWARD_PROGRESSION (LFP)` and `MEAN_STRAIGHT_LINE_SPEED` (per‑cell: medians):**  
  *What.* LFP compares straight‑line speed vs. mean speed; mean straight‑line speed is the **vectorial** speed.  
  *Why.* Both summarize **directional persistence** at the track level.

- **`MEAN_DIRECTIONAL_CHANGE_RATE` (per‑cell: median and P90):**  
  *What.* Average **turning rate** of a track.  
  *Why.* Quantifies **tortuosity/flexion**. P90 highlights the **sharpest turns**.


### Edges (instantaneous, useful for distribution tails)

- **`SPEED` (edge‑to‑edge):**  
  *What.* Instantaneous speed between two consecutive detections.  
  *Why.* Gives the **full distribution**, revealing very slow/fast **tails** that are lost by per‑track averaging.  
  *How summarized here.* Two‑stage aggregation to avoid biasing long tracks:  
  (i) per‑track summaries (**median, P10, P90, IQR**) →  
  (ii) per‑cell **median** across tracks.

- **`DIRECTIONAL_CHANGE_RATE` (edge‑scale):**  
  *What.* Frame‑to‑frame angular change.  
  *Why.* Captures **tortuosity** at the finest time scale.  
  *Summary.* Same two‑stage approach as for speed.  

---

## Why medians and two‑stage aggregation?

- **Medians** are resistant to outliers (e.g., rare false links) and do not over‑emphasize unusually long tracks.  
- **Two‑stage aggregation (per‑track → per‑cell)** prevents long tracks from dominating per‑cell distributions, which would happen if you pooled all edges directly. This respects the **cell as the statistical unit**.

---

## Installation

- Python ≥ 3.10
- Packages: `numpy`, `pandas`, `scipy`, `seaborn`, `matplotlib`

```bash
pip install -r requirements.txt
# or
pip install numpy pandas scipy seaborn matplotlib
```

---

## Usage details

### Command‑line options
- `-i/--input` (required): root folder that contains your experiment/condition/[plate]/CSV tree.  
- `-o/--output` (optional): base folder for `1_METRICS`, `2_STATS`, `3_PLOTS`, `logs`. If omitted, outputs go under the input root.  
- `-s/--stats` (optional): compute Mann–Whitney statistics (**exactly 2 conditions** required). If more than two conditions are present, stats are **skipped** with a warning.  
- `-p/--plots` (optional): generate plots. Plots work with or without `-s`; when no stats are available, significance annotations are omitted.

### What the script does **not** do
- It does **not** detect/tracks comets (TrackMate must export the CSVs).  
- It ignores any `*all-spots*.csv` files and does not parse **spots** reports as metrics.  - It does **not** normalize `N_TRACKS` by time or area; it reports the **raw count per cell**.
- It does **not** compute appearance rates from `TRACK_START/STOP`, gap metrics, or normalize by **cell area**.  
- It does **not** convert/standardize units; it trusts the TrackMate CSV values.  
- It does **not** Δt‑weight edge summaries; it assumes that frame intervals (Δt) are the same across tracks.
- Statistical testing currently supports **two conditions**; >2 conditions are not compared here.

---

## Reproducibility & logging

- Every run tees stdout/stderr to `logs/trackmate_pipeline_<timestamp>.log`.  
- Save alongside: microscope settings, TrackMate detector/tracker parameters, and software versions.

---

## Demo dataset (quick test)

A minimal, self-contained dataset is provided under `examples/dummy_dataset/` so you can run the pipeline end-to-end in seconds.

**Run the demo**
```bash
# Metrics only
python trackmate_eb3_tirfm_analysis.py -i examples/dummy_dataset/input

# Metrics + stats + plots, writing results to a separate "run" folder
python trackmate_eb3_tirfm_analysis.py   -i examples/dummy_dataset/input   -o examples/dummy_dataset/run   -s -p
```

**What you should see**
- `examples/dummy_dataset/run/1_METRICS/cell_metrics.csv`
- `examples/dummy_dataset/run/1_METRICS/track_metrics_from_edges.csv`
- (if `-s`) `examples/dummy_dataset/run/2_STATS/mannwhitney_all_metrics.csv`
- (if `-p`) figures in `examples/dummy_dataset/run/3_PLOTS/`

The expected outputs are provided under `examples/dummy_dataset/expected_output/`

---

## Troubleshooting

- **“Cell omitted due to empty metrics (NaN/0)”**  
  Likely no valid numeric columns were parsed (header shift, encoding), or all selected metrics were zero/NaN. Check the log for the **detected header/encoding** and verify your TrackMate export columns.

- **No stats were generated**  
  The dataset didn’t have exactly two conditions. Plots can still be produced with `-p`.

- **Weird units or magnitudes**  
  Confirm pixel size and frame interval before tracking. The script does not rescale values.

- **Too many false links**  
  Reduce **linking distance** and/or **gap‑closing distance** in TrackMate; consider mild pre‑filtering and a higher quality threshold.
  
- **Unexpectedly high/low `N_TRACKS`**  
  Check TrackMate **linking** and **gap-closing** distances and **quality thresholds**. Over-fragmentation inflates counts; overly strict thresholds deflate them. Ensure comparable **video duration** across cells and conditions.

---

## Citation

If you use this pipeline, please cite **both** the software release and the Methods chapter:

**Software:**  
Lozano-Prieto M. *TrackMate EB3 TIRFM Analysis Pipeline* (v1.0.0). Zenodo. **DOI:** 10.5281/zenodo.XXXXXXX.

**Methods chapter:**  
Gómez-Morón Á, Lozano-Prieto M, Martín Cófreces NB. *: Studying actin and tubulin cytoskeleton dynamics at planar immunological synapses in primary human T lymphocytes through TIRF Microscopy*. Methods in Molecular Biology. 2025.
