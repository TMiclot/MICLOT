# MDQ: Molecular Dynamics DockQ

A Python toolkit to analyse protein–protein (or protein–chain) interfaces across molecular dynamics (MD) trajectories.

MDQ uses a symmetric, DockQ-inspired metric. **No native reference structure is needed**: all frames are compared to each other, pair by pair. The package calculates contact conservation, interface RMSD (iRMSD), ligand RMSD (lRMSD) and an MDQ similarity score. It then clusters and visualises the trajectory frames.

---

## Contents

- [Overview](#overview)
- [Methodology](#methodology)
  - [1. Fraction of Shared Contact (FSC)](#1-fraction-of-shared-contact-fsc)
  - [2. Weighted interface RMSD (iRMSD)](#2-weighted-interface-rmsd-irmsd)
  - [3. Weighted ligand RMSD (lRMSD)](#3-weighted-ligand-rmsd-lrmsd)
  - [4. MDQ score](#4-mdq-score)
  - [Score interpretation](#score-interpretation)
- [Analysis workflow](#analysis-workflow)
- [Input files and conventions](#input-files-and-conventions)
- [Step 1: Calculate interface metrics (`TrajectoryAnalyzer`)](#step-1-calculate-interface-metrics-trajectoryanalyzer)
- [Step 2: Cluster trajectory frames (`TrajectoryClusterVisualizer`)](#step-2-cluster-trajectory-frames-trajectoryclustervisualizer)
- [Step 3: Export and summarise clusters](#step-3-export-and-summarise-clusters)
- [Usage examples](#usage-examples)
- [Output files](#output-files)
- [Interpretation](#interpretation)
- [Known issues and limitations](#known-issues-and-limitations)
- [Troubleshooting](#troubleshooting)
- [Dependencies](#dependencies)
- [Reproducibility](#reproducibility)

---

## Overview

The package has two classes.

### `TrajectoryAnalyzer`

It loads an MD trajectory and calculates pairwise frame metrics for every pair of chains you define:

- **FSC (Fraction of Shared Contact):** similarity between the residue–residue contacts of two frames.
- **Weighted iRMSD:** interface backbone RMSD, penalised when few atoms are shared.
- **Weighted lRMSD:** RMSD of the shorter chain (the ligand) after receptor superposition, also penalised.
- **MDQ score:** a combined score from FSR, iRMSD and lRMSD (weighted and unweighted versions).

### `TrajectoryClusterVisualizer`

It reads a pairwise distance table (`distance = 1 - MDQ`) and provides:

- DBSCAN clustering;
- agglomerative hierarchical clustering with dendrograms;
- MDS, UMAP and t-SNE projections;
- distance heatmaps and clustermaps;
- representative frame (medoid) of each cluster;
- export of representative PDB files, cluster tables and cluster statistics.

The two classes are separate. The first one calculates the metrics. The second one uses a distance table for clustering and plotting.


---

## Methodology

### 1. Fraction of Shared Contact (FSC)

**Contact definition.** In DockQ, two residues are in contact if the distance between a heavy atom of one residue and a heavy atom of a residue in another chain is 5 Å or less.

**The asymmetry problem.** The classical DockQ parameter `Fnat` is the fraction of native contacts that are kept in the model:

```text
Fnat = number of preserved native contacts / total number of native contacts
```

`Fnat` changes depending on which structure is the reference. A contact in the native structure can be longer than 5 Å in the model, and the opposite is also true.

**The FSR solution.** FSR is symmetric. The method is:

1. Get the residue indices of each chain and their heavy atoms.
2. Calculate the distances for all residue pairs between two chains (closest heavy-atom distance).
3. Keep the contacts inside the distance interval (default 2–5 Å).
4. Compare the two contact sets:

```text
FSR(i,j) = 2 × |C_i ∩ C_j| / (|C_i| + |C_j|)
```

- `C_i`, `C_j`: sets of contacts in frame `i` and frame `j`.
- `|C_i ∩ C_j|`: number of contacts (same residue pair) shared by both frames.

If both frames have no contact, FSR is set to `0.0`.

### 2. Weighted interface RMSD (iRMSD)

**Background.** In DockQ, iRMSD is the backbone RMSD of the interface. The interface is defined with a cutoff of 10 Å, which is twice the `Fnat` cutoff.

**The asymmetry problem.** The number of atoms used for the superposition can be different between two frames. Then the RMSD values cannot be compared directly.

**The solution.** For each interface and each pair of frames:

1. Find the interface backbone atoms (`N`, `CA`, `C`, `O`) of each frame. The cutoff is 10 Å (1.0 nm) between backbone atoms of the two chains.
2. Keep only the atoms that are common to both frames.
3. Superimpose the frames on these common atoms and calculate the RMSD (in Å).
4. Calculate the fraction of common atoms:

```text
f(i,j) = 2 × |A_i ∩ A_j| / (|A_i| + |A_j|)
```

5. Apply a penalty. A low overlap gives a large penalty, so an RMSD calculated on very few atoms is not artificially low:

```text
loss = -log2(f² + 1e-15) + 1
irmsd_weighted = irmsd_unweighted × loss
```

The value `1e-15` avoids `log2(0)`.

### 3. Weighted ligand RMSD (lRMSD)

**Background.** lRMSD is the RMSD of the backbone of the shorter chain (the ligand) after superposition of the longer chain (the receptor).

> In DockQ version 2, the superposition uses only the interface of the receptor, not its full backbone.

**In MDQ:**

1. The longer chain is the receptor. The shorter chain is the ligand. If both chains have the same length, the alphabetical order of the labels decides.
2. Find the receptor backbone atoms at the interface that are common to the two frames.
3. Superimpose the frames on these atoms.
4. Calculate the RMSD of all ligand backbone atoms.
5. Apply the same penalty as for iRMSD. Here, `f` is the fraction of common receptor interface backbone atoms.

### 4. MDQ score

MDQ follows the CAPRI (DockQ) definition:

```text
MDQ = 1/3 × [ FSR + 1/(1+(iRMSD/1.5)²) + 1/(1+(lRMSD/8.5)²) ]
```

Compared with DockQ:

- `Fnat` is replaced by **FSR**;
- `iRMSD` is replaced by **weighted iRMSD**;
- `lRMSD` is replaced by **weighted lRMSD**.

An unweighted version (`mdq_unweighted`) is also calculated for comparison. The constants 1.5 and 8.5 are in Å.

### Score interpretation

MDQ goes from 0 to 1. A higher value means that the two frames are more similar.

| Score range | Interpretation |
|-------------|----------------|
| > 0.80 | Very high similarity (almost identical interfaces) |
| 0.60 – 0.80 | High similarity |
| 0.40 – 0.60 | Moderate similarity |
| 0.20 – 0.40 | Low similarity |
| < 0.20 | Very low similarity |

---

## Analysis workflow

1. Prepare a trajectory and a topology file for the same system.
2. Define the residue indices of each chain.
3. Create a `TrajectoryAnalyzer`.
4. Calculate FSR, iRMSD and lRMSD.
5. Calculate the MDQ score.
6. Save the result as `mdq.csv`.
7. Create a `TrajectoryClusterVisualizer`.
8. Run DBSCAN and/or hierarchical clustering.
9. Look at the projections and distance plots.
10. Export representative frames and cluster statistics.

All metrics are calculated for every pair of frames. The number of comparisons grows with the **square** of the number of frames.

---

## Input files and conventions

### Trajectory and topology

The trajectory is loaded with:

```python
md.load(traj_file, top=top_file, stride=stride)
```

- `traj_file`: trajectory file (`.xtc`, `.dcd`, or another MDTraj format).
- `top_file`: topology file (`.pdb`, `.gro`, `.prmtop`, ...).
- The topology must match the atom order and atom count of the trajectory.

### Chain dictionary

`chains` links each chain label to a list of **MDTraj residue indices**:

```python
chains = {
    "A": [0, 1, 2, 3, 4],
    "B": [100, 101, 102, 103, 104]
}
```

These are **not** PDB residue numbers (`resSeq`). MDTraj residue indices start at 0. Always check that they match the loaded topology. Use `range(a, b)` to create the lists.

### Stride

`stride` sets how many frames are loaded:

- `stride=1`: every frame;
- `stride=10`: every 10th frame.

A larger stride makes the calculation faster and uses less memory, but it lowers the time resolution. The frame indices in the outputs refer to the **loaded (strided) trajectory**, not to the original frame numbers.

---

## Step 1: Calculate interface metrics (`TrajectoryAnalyzer`)

### Initialisation

```python
analyzer = TrajectoryAnalyzer(
    traj_file="production.xtc",
    top_file="topology.pdb",
    chains={"A": list(range(0, 314)), "B": list(range(314, 634))},
    stride=1
)
```

| Argument | Type | Description |
|----------|------|-------------|
| `traj_file` | `str` | Path to the trajectory file |
| `top_file` | `str` | Path to the topology file |
| `chains` | `dict` | Chain labels and lists of MDTraj residue indices |
| `stride` | `int` | Frame-loading interval (default `1`) |

During initialisation, the class loads the trajectory (`analyzer.traj`), finds the backbone atoms (`N`, `CA`, `C`, `O`) of each chain (`self.backbone_indices`), and creates a one-frame slice for each frame (`self.slice_frame`).

### `compute_shared_fraction(min_A=2, max_A=5)`

Calculates FSR for all frame pairs (including self-comparisons).

| Parameter | Default | Description |
|-----------|---------|-------------|
| `min_A` | `2` | Minimum contact distance (Å) |
| `max_A` | `5` | Maximum contact distance (Å) |

- **Result:** `analyzer.df_shared`
- **File written:** `shared_fraction.csv`

### `compute_irmsd()`

Calculates unweighted and weighted iRMSD for each pair of frames. It has no parameter. The interface cutoff is fixed at 1.0 nm. It is different from the contact interval used for FSR.

- **Result:** `analyzer.df_irmsd` with the columns `interface`, `frame_i`, `frame_j`, `irmsd_unweighted`, `irmsd_weighted` (in Å).
- **Intermediate table:** `analyzer.df_interface_backbone` (interface backbone atoms of each frame). It is used by `compute_lrmsd()`.

### `compute_lrmsd()`

Calculates unweighted and weighted lRMSD. **Run `compute_irmsd()` first**, because this method needs `df_interface_backbone`.

- **Result:** `analyzer.df_lrmsd` with the columns `interface`, `frame_i`, `frame_j`, `lrmsd_unweighted`, `lrmsd_weighted` (in Å).

### `compute_mdq()`

Merges the three tables on `interface`, `frame_i` and `frame_j`, then calculates the MDQ scores. It also reads `shared_fraction.csv` from the working directory.

- **Result:** `analyzer.df_score` with all the metric columns plus `mdq_weighted` and `mdq_unweighted`.

**Important:** `compute_mdq()` does **not** write `mdq.csv`. Save it yourself:

```python
analyzer.df_score.to_csv("mdq.csv", index=False)
```

---

## Step 2: Cluster trajectory frames (`TrajectoryClusterVisualizer`)

### Distance table

Clustering needs distances, not similarities:

```text
distance(i,j) = 1 - mdq_weighted(i,j)
```

A high MDQ gives a low distance.

If `distance.csv` already exists:

```python
visualizer = TrajectoryClusterVisualizer(distance_file="distance.csv")
```

If you give no file, the constructor calls `create_distance_from_mdq()`. It reads `mdq.csv` (columns `interface`, `frame_i`, `frame_j`, `mdq_weighted`) and writes `distance.csv`:

```python
visualizer = TrajectoryClusterVisualizer()
```

### `run_dbscan_analysis(...)`

DBSCAN works on the precomputed distance matrix, one interface at a time. Frames that do not belong to a cluster get the label `-1` (noise).

```python
visualizer.run_dbscan_analysis(
    eps_dict={"A-B": 0.20},
    min_samples_pct=0.05,
    min_pct_dyn=0.05,
    total_frames_dyn=1000
)
```

| Parameter | Default | Description |
|-----------|---------|-------------|
| `eps_dict` | required | `eps` distance threshold for each interface. If an interface is missing, `0.2` is used |
| `min_samples_pct` | `0.05` | Fraction of the frames used for `min_samples` |
| `min_pct_dyn` | `0.05` | Minimum cluster size (fraction). Smaller clusters become noise (`-1`) |
| `total_frames_dyn` | required | Total number of frames used to calculate the cluster-size percentage |

For each interface, the method builds the distance matrix, runs DBSCAN with `metric="precomputed"` (`min_samples = max(1, int(min_samples_pct * n_frames))`), filters small clusters, stores the results, and plots the MDS, UMAP and t-SNE projections.

- **Result:** `visualizer.dbscan_labels_dict`. Each interface gives `(frames, distance_matrix, labels)`.

**About `eps`:** it is a *distance* threshold, not an MDQ threshold. A smaller `eps` means a stricter similarity requirement. Choose it for each system and interface, using the distribution of distances (see the heatmap).

### `run_hierarchical_analysis(...)`

Agglomerative clustering with a dendrogram. You do not need to choose the number of clusters.

```python
visualizer.run_hierarchical_analysis(
    method="average",
    threshold=0.3,
    min_pct=0.05,
    total_frames=1000
)
```

| Parameter | Default | Description |
|-----------|---------|-------------|
| `method` | `"average"` | Linkage method: `"single"`, `"complete"`, `"average"` or `"ward"` |
| `threshold` | `0.3` | Distance used to cut the dendrogram |
| `min_pct` | `0.05` | Minimum cluster size (fraction) |
| `total_frames` | required | Total number of frames used for the cluster-size percentage |

The method converts the distance matrix to condensed form (`squareform`), calculates the SciPy linkage, plots the dendrogram, cuts it with `fcluster(..., criterion="distance")`, relabels the clusters, filters small clusters, and plots the projections.

- **Result:** `visualizer.hierarchical_labels_dict` (same structure as for DBSCAN).

> `"ward"` is designed for Euclidean distances. Use it with care on an MDQ distance matrix.

### Projection plots

Both methods draw **MDS**, **UMAP** and **t-SNE** projections of the distance matrix. Points are coloured by cluster. Noise (`-1`) is light grey. The representative frame of each cluster has a larger marker with a black outline. Frame indices are written on the plots.

Projections are only a visual help. Distances in a 2D plot are not the real MDQ distances.

### `plot_distance_heatmap()`

Shows the distance matrix of each interface as a heatmap. The colour scale is 0 to 1. Red cells are missing values (frames that never interact).

```python
visualizer.plot_distance_heatmap()
```

### `plot_distance_clustermap()`

Shows a seaborn clustermap (heatmap with hierarchical clustering) for each interface. It helps to see natural groups.

```python
visualizer.plot_distance_clustermap()
```

---

## Step 3: Export and summarise clusters

### Get frames and representatives

```python
dbscan_table, hierarchical_table = visualizer.get_cluster_frames_table()
```

For each interface and each cluster, you get `frames` (frame indices) and `representative` (representative frame). Noise points are excluded.

The representative is the frame with the smallest mean distance to the other frames of its cluster (medoid-like).

You can also use:

```python
representatives = visualizer.get_cluster_representatives()
# {"DBSCAN": {"A-B": {0: {"frames": [...], "representative": 1}}}, "Hierarchical": {...}}
```

### Export representative structures (PDB)

```python
visualizer.export_representative_frames_pdb(
    traj_file="production.xtc",
    top_file="topology.pdb",
    output_dir="representative_structures",
    cluster_tables=(dbscan_table, hierarchical_table)
)
```

File names:

```text
DBSCAN_<interface>_cluster<label>_rep.pdb
Hierarchical_<interface>_cluster<label>_rep.pdb
```

**Careful with the frame index.** The trajectory used for the export must have the same frame order and sampling as the analysed trajectory. If you used a stride, use the same stride.

### Export representatives (CSV)

```python
visualizer.export_cluster_representatives_csv(output_file="cluster_representatives.csv")
```

| Column | Description |
|--------|-------------|
| `method` | `DBSCAN` or `Hierarchical` |
| `interface` | Interface name |
| `cluster` | Cluster label |
| `n_frames` | Number of frames in the cluster |
| `representative_frame` | Representative frame index |
| `frames` | Frame indices, separated by commas |

### Save cluster statistics

```python
visualizer.save_cluster_stats_csv(output_prefix="cluster_stats")
```

It writes `cluster_stats_DBSCAN.csv` and `cluster_stats_Hierarchical.csv`. For each interface and each non-noise cluster, the statistics use the mean within-cluster distance of each frame:

- number of frames;
- mean and standard deviation;
- Q1, median and Q3;
- lower and upper bounds of a 95% Student's t confidence interval.

For a cluster with one frame, the standard deviation and the confidence interval are not defined. Check these cases carefully.

---

## Usage examples

### Two chains

```python
from mdq import TrajectoryAnalyzer, TrajectoryClusterVisualizer

traj_file = "production.xtc"
top_file = "topology.pdb"
chains = {"A": range(0, 314), "B": range(314, 634)}

# 1. Metrics
analyzer = TrajectoryAnalyzer(traj_file, top_file, chains, stride=1)
analyzer.compute_shared_fraction(min_A=2, max_A=5)
analyzer.compute_irmsd()
analyzer.compute_lrmsd()
analyzer.compute_mdq()
analyzer.df_score.to_csv("mdq.csv", index=False)

# 2. Distance table (reads mdq.csv, writes distance.csv)
visualizer = TrajectoryClusterVisualizer()

# 3. Clustering
visualizer.run_dbscan_analysis(
    eps_dict={"A-B": 0.20},
    min_samples_pct=0.05,
    min_pct_dyn=0.05,
    total_frames_dyn=1000
)
visualizer.run_hierarchical_analysis(
    method="average", threshold=0.3, min_pct=0.05, total_frames=1000
)

# 4. Plots
visualizer.plot_distance_heatmap()
visualizer.plot_distance_clustermap()

# 5. Export
dbscan_table, hierarchical_table = visualizer.get_cluster_frames_table()
visualizer.export_representative_frames_pdb(
    traj_file=traj_file,
    top_file=top_file,
    output_dir="representative_structures",
    cluster_tables=(dbscan_table, hierarchical_table)
)
visualizer.export_cluster_representatives_csv("cluster_representatives.csv")
visualizer.save_cluster_stats_csv(output_prefix="cluster_stats")
```

### Several interfaces

```python
chains = {
    "A": range(0, 200),
    "B": range(200, 400),
    "C": range(400, 600)
}
analyzer = TrajectoryAnalyzer("traj.xtc", "top.pdb", chains)
analyzer.compute_shared_fraction()
analyzer.compute_irmsd()
analyzer.compute_lrmsd()
analyzer.compute_mdq()
analyzer.df_score.to_csv("mdq.csv", index=False)

visualizer = TrajectoryClusterVisualizer()
eps_dict = {"A-B": 0.12, "A-C": 0.18, "B-C": 0.15}

visualizer.run_dbscan_analysis(
    eps_dict, min_samples_pct=0.05, min_pct_dyn=0.05, total_frames_dyn=1000
)
visualizer.run_hierarchical_analysis(
    method="average", threshold=0.25, min_pct=0.05, total_frames=1000
)
```

---

## Output files

| File | Created by | Description |
|------|-----------|-------------|
| `shared_fraction.csv` | `compute_shared_fraction()` | Pairwise FSR values |
| `mdq.csv` | `analyzer.df_score.to_csv(...)` | Pairwise metrics and MDQ scores. Needed by the default visualiser |
| `distance.csv` | `create_distance_from_mdq()` | Pairwise distance (`1 - mdq_weighted`) |
| `representative_structures/` | `export_representative_frames_pdb()` | Representative PDB files |
| `cluster_representatives.csv` | `export_cluster_representatives_csv()` | Cluster members and representative frames |
| `cluster_stats_DBSCAN.csv` | `save_cluster_stats_csv()` | DBSCAN cluster statistics |
| `cluster_stats_Hierarchical.csv` | `save_cluster_stats_csv()` | Hierarchical cluster statistics |

The DataFrames `df_shared`, `df_irmsd`, `df_lrmsd` and `df_score` are also available in memory on the `TrajectoryAnalyzer` object.

### File formats

**`shared_fraction.csv`**

| interface | frame_i | frame_j | shared_fraction |
|-----------|---------|---------|-----------------|
| A-B | frame_0 | frame_0 | 1.0 |
| A-B | frame_0 | frame_1 | 0.42 |

**`mdq.csv`**

| interface | frame_i | frame_j | shared_fraction | irmsd_weighted | lrmsd_weighted | mdq_weighted |
|-----------|---------|---------|-----------------|----------------|----------------|--------------|
| A-B | 0 | 1 | 0.42 | 2.15 | 3.82 | 0.67 |
| A-B | 0 | 2 | 0.38 | 2.48 | 4.21 | 0.61 |

(The real file also contains the unweighted columns and `mdq_unweighted`.)

**`distance.csv`**

| interface | frame_i | frame_j | distance |
|-----------|---------|---------|----------|
| A-B | 0 | 1 | 0.33 |
| A-B | 0 | 2 | 0.39 |

---

## Interpretation

- **Similarity and distance.** MDQ is a similarity. Clustering uses `D = 1 - MDQ_weighted`. Low distance means high similarity.
- **Weighted and unweighted.** The weighted MDQ uses RMSD values multiplied by the overlap penalty. Compare both versions to see if the shared interface changes your conclusion.
- **Representative frame.** It is the frame with the smallest mean distance to the other frames of the cluster. It is not necessarily the most frequent frame, and it is not the frame closest to a reference structure.
- **Clustering parameters.** The results depend on `eps`, `min_samples_pct`, the cluster-size filter, the linkage method and the dendrogram threshold. Choose them from the distance distribution and your scientific question, and report them.

---

## Known issues and limitations

These points describe the current `mdq.py`. Check them before a production analysis.

1. **Method signatures.** In `run_dbscan_analysis()` and `run_hierarchical_analysis()`, required arguments must come before arguments with default values. Otherwise Python raises a `SyntaxError` when importing the file.
2. **No automatic `mdq.csv`.** `compute_mdq()` only fills `self.df_score`. Save the file yourself.
3. **Frame sampling.** Frame indices refer to the loaded (strided) trajectory. Keep the same order when you export structures.
4. **Cluster-size filter.** The minimum cluster size is calculated with the `total_frames` value you give. Use a value that matches your definition of the minimum fraction.
5. **Scaling.** The number of pair calculations grows like N², where N is the number of frames. Large trajectories produce large tables and long run times.
6. **Backbone definition.** Only atoms named `N`, `CA`, `C` and `O` are used.
7. **Receptor and ligand.** The longer chain is the receptor and the shorter is the ligand. For equal lengths, the alphabetical order decides. This is a convention of the code, not a biological annotation.

---

## Troubleshooting

**`SyntaxError: non-default argument follows default argument`**
Fix the argument order in `run_dbscan_analysis()` and `run_hierarchical_analysis()` (see [Known issues](#known-issues-and-limitations)).

**`FileNotFoundError: mdq.csv`**
Save the table first:

```python
analyzer.df_score.to_csv("mdq.csv", index=False)
```

Run the visualiser from the folder that contains `mdq.csv`, or give an explicit `distance_file`.

**Chains or residues look wrong**
Use zero-based MDTraj residue indices from the topology loaded with the trajectory. PDB residue numbers and MDTraj indices are different.

**Representative PDB files do not match the analysed frames**
Use the same frame order and the same `stride` for the export trajectory.

**Unexpected number of clusters**
Look at the distance heatmap and the dendrogram. Check `eps` and `min_samples_pct` (DBSCAN) or `threshold` and `method` (hierarchical). Also check `total_frames` or `total_frames_dyn`.

---

## Dependencies

| Package | Version | Purpose |
|---------|---------|---------|
| `mdtraj` | ≥ 1.9.0 | Trajectory reading, distances and RMSD |
| `numpy` | ≥ 1.19.0 | Numerical operations |
| `pandas` | ≥ 1.2.0 | DataFrames and CSV files |
| `scikit-learn` | ≥ 0.24.0 | DBSCAN, MDS and t-SNE |
| `scipy` | ≥ 1.6.0 | Hierarchical clustering, `squareform`, statistics |
| `umap-learn` | – | UMAP projection |
| `matplotlib` | ≥ 3.3.0 | Plots |
| `seaborn` | ≥ 0.11.0 | Heatmaps and clustermaps |
| `tqdm` | ≥ 4.50.0 | Progress bars |

`ast` and `warnings` are in the Python standard library.

---

## Reproducibility

When you report results from this code, write down:

- the trajectory and topology files (and how they were prepared);
- the chain definitions and the residue-index convention;
- the stride and the number of analysed frames;
- the FSR distance interval (`min_A`, `max_A`);
- the iRMSD interface cutoff (1.0 nm);
- the MDQ formula and constants (1.5 Å and 8.5 Å);
- whether you used weighted or unweighted MDQ;
- the distance definition (`1 - mdq_weighted`);
- the clustering algorithm and all its parameters;
- the software versions and any change made to `mdq.py`.
