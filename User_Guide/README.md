# MDQ: Molecular Dynamics DockQ

A Python toolkit for analysing protein--protein or protein--chain
interfaces across molecular dynamics (MD) trajectories. The package
calculates contact conservation, interface RMSD (iRMSD), ligand RMSD
(lRMSD), an MDQ similarity score, and clustering/visualisation of
trajectory frames.

------------------------------------------------------------------------

## Contents

-   [Overview](#overview)
-   [Analysis workflow](#analysis-workflow)
-   [Requirements](#requirements)
-   [Input files and topology
    conventions](#input-files-and-topology-conventions)
-   [Installation](#installation)
-   [1. Calculate interface metrics](#1-calculate-interface-metrics)
    -   [Initialise
        `TrajectoryAnalyzer`](#initialise-trajectoryanalyzer)
    -   [Shared fraction (SF)](#shared-fraction-sf)
    -   [Interface RMSD (iRMSD)](#interface-rmsd-irmsd)
    -   [Ligand RMSD (lRMSD)](#ligand-rmsd-lrmsd)
    -   [MDQ score](#mdq-score)
-   [2. Cluster trajectory frames](#2-cluster-trajectory-frames)
    -   [Create a distance matrix from
        MDQ](#create-a-distance-matrix-from-mdq)
    -   [DBSCAN clustering](#dbscan-clustering)
    -   [Hierarchical clustering](#hierarchical-clustering)
    -   [Projection plots](#projection-plots)
    -   [Distance heatmaps and
        clustermaps](#distance-heatmaps-and-clustermaps)
-   [3. Export and summarise clusters](#3-export-and-summarise-clusters)
-   [Output files](#output-files)
-   [Example end-to-end workflow](#example-end-to-end-workflow)
-   [Interpretation](#interpretation)
-   [Implementation notes and
    limitations](#implementation-notes-and-limitations)
-   [Troubleshooting](#troubleshooting)

------------------------------------------------------------------------

## Overview

The analysis is organised into two classes:

### `TrajectoryAnalyzer`

Loads an MD trajectory and calculates pairwise frame metrics for every
pair of chains defined by the user:

-   **Shared fraction (SF):** similarity between the sets of
    residue--residue contacts observed in two frames.
-   **Interface RMSD (iRMSD):** backbone RMSD calculated using atoms
    identified as belonging to the interface.
-   **Ligand RMSD (lRMSD):** RMSD of the shorter chain (treated as the
    ligand) after receptor-based alignment, according to the
    implementation in the source file.
-   **MDQ score:** combines SF, iRMSD and lRMSD into weighted and
    unweighted similarity values.

### `TrajectoryClusterVisualizer`

Reads a pairwise distance table and performs clustering and
visualisation:

-   DBSCAN clustering;
-   agglomerative hierarchical clustering;
-   MDS, UMAP and t-SNE projections;
-   distance heatmaps and clustermaps;
-   representative-frame identification;
-   export of representative structures and cluster tables/statistics.

The two classes are separate: the first calculates the metrics, while
the second consumes a distance table for clustering and visualisation.

------------------------------------------------------------------------

## Analysis workflow

The intended workflow is:

1.  Prepare a trajectory and a topology file that describe the same
    molecular system.
2.  Define the residue indices belonging to each chain.
3.  Instantiate `TrajectoryAnalyzer`.
4.  Calculate SF, iRMSD and lRMSD.
5.  Calculate the MDQ score.
6.  Save the MDQ results as `mdq.csv` if you want to use the
    visualisation class' automatic distance-file creation.
7.  Instantiate `TrajectoryClusterVisualizer`.
8.  Run DBSCAN and/or hierarchical clustering.
9.  Inspect projections and distance plots.
10. Export representative frames and cluster statistics if required.

The metrics are calculated for every pair of frames, so the number of
comparisons grows quadratically with the number of frames.

------------------------------------------------------------------------

## Requirements

The source file imports the following Python packages:

  -----------------------------------------------------------------------
  Package                             Main use
  ----------------------------------- -----------------------------------
  `mdtraj`                            Reading trajectories, calculating
                                      contacts/distances and RMSD

  `numpy`                             Numerical operations and distance
                                      matrices

  `pandas`                            DataFrames and CSV input/output

  `tqdm`                              Progress bars

  `matplotlib`                        Plotting

  `seaborn`                           Heatmaps and clustermaps

  `scikit-learn`                      DBSCAN, MDS and t-SNE

  `scikit-learn-extra`                KMedoids import (not used by the
                                      shown class methods)

  `scipy`                             Hierarchical clustering, square
                                      distance matrices and statistics

  `umap-learn`                        UMAP projection
  -----------------------------------------------------------------------

`ast` and `warnings` are part of the Python standard library.

### Installation

Install the dependencies in the Python environment used to run the
analysis:

``` bash
python -m pip install mdtraj numpy pandas tqdm matplotlib seaborn scikit-learn scikit-learn-extra scipy umap-learn
```

Depending on the operating system and Python version,
`scikit-learn-extra` may require additional build dependencies. If
installation fails, consult that package's installation instructions.

------------------------------------------------------------------------

## Input files and topology conventions

### Trajectory and topology

The constructor uses:

``` python
md.load(traj_file, top=top_file, stride=stride)
```

-   `traj_file`: trajectory file, for example `.xtc`, `.dcd` or another
    format supported by MDTraj.
-   `top_file`: topology file used to interpret the trajectory, for
    example `.pdb` or `.gro`.
-   The topology must match the atom ordering and atom count of the
    trajectory.

### Chain dictionary

The `chains` argument maps chain labels to **MDTraj residue indices**:

``` python
chains = {
    "A": [0, 1, 2, 3, 4],
    "B": [100, 101, 102, 103, 104]
}
```

In this example, the values are residue indices, not PDB residue numbers
(`resSeq`). MDTraj residue indices are zero-based. Make sure the indices
correspond to the topology loaded by MDTraj.

The implementation identifies backbone atoms named `N`, `CA`, `C` and
`O` for each chain. These atom indices are then used in the iRMSD and
lRMSD calculations.

### Choosing a stride

`stride` controls how many trajectory frames are read:

-   `stride=1`: load every frame;
-   `stride=10`: load every tenth frame;
-   `stride=100`: load every hundredth frame.

A larger stride reduces the number of pairwise comparisons and therefore
the computational and memory cost. It also reduces temporal resolution.
The frame indices used in the output correspond to the loaded, strided
trajectory, not necessarily the original trajectory's absolute frame
numbering.

------------------------------------------------------------------------

## Installation and import

Place `mdq.py` in your project directory, then import the classes:

``` python
from mdq import TrajectoryAnalyzer, TrajectoryClusterVisualizer
```

Before importing, fix the invalid parameter ordering noted at the top of
this guide.

------------------------------------------------------------------------

# 1. Calculate interface metrics

## Initialise `TrajectoryAnalyzer`

``` python
from mdq import TrajectoryAnalyzer

traj_file = "production.xtc"
top_file = "topology.pdb"

chains = {
    "A": list(range(0, 314)),
    "B": list(range(314, 634))
}

analyzer = TrajectoryAnalyzer(
    traj_file=traj_file,
    top_file=top_file,
    chains=chains,
    stride=1
)
```

### Constructor arguments

  -----------------------------------------------------------------------
  Argument                Type                    Description
  ----------------------- ----------------------- -----------------------
  `traj_file`             `str`                   Path to the trajectory
                                                  file

  `top_file`              `str`                   Path to the topology
                                                  file

  `chains`                `dict`                  Mapping of chain labels
                                                  to lists of MDTraj
                                                  residue indices

  `stride`                `int`                   Frame-loading interval;
                                                  defaults to `1`
  -----------------------------------------------------------------------

During initialisation, the class:

1.  loads the trajectory with MDTraj;
2.  builds a residue-index-to-chain mapping;
3.  identifies backbone atoms (`N`, `CA`, `C`, `O`) for each specified
    chain;
4.  stores the backbone atom indices in `self.backbone_indices`;
5.  creates a one-frame trajectory slice for every loaded frame in
    `self.slice_frame`.

The object stores the loaded trajectory in `analyzer.traj`.

------------------------------------------------------------------------

## Shared fraction (SF)

### Purpose

`compute_shared_fraction()` compares the contact sets of every pair of
frames for each chain interface.

The code builds all residue pairs between two chains and calculates
their closest-heavy-atom distances using MDTraj. Distances are converted
from nanometres to ångströms, then contacts are selected within the
inclusive interval defined by `min_A` and `max_A`.

The shared fraction is calculated as:

\[ SF(i,j) = `\frac{2|C_i \cap C_j|}{|C_i| + |C_j|}`{=tex} \]

where:

-   (C_i) is the set of contacts in frame (i);
-   (C_j) is the set of contacts in frame (j);
-   (\|C_i `\cap `{=tex}C_j\|) is the number of contacts shared by the
    two frames.

If both frames contain no contacts, the implementation assigns `0.0`.

### Usage

``` python
analyzer.compute_shared_fraction(
    min_A=2,
    max_A=5
)
```

### Parameters

  Parameter     Default Description
  ----------- --------- -------------------------------
  `min_A`           `2` Minimum contact distance in Å
  `max_A`           `5` Maximum contact distance in Å

### Output

The method stores the results in:

``` python
analyzer.df_shared
```

It also writes:

``` text
shared_fraction.csv
```

The CSV contains:

  Column              Description
  ------------------- ----------------------------------------------------
  `interface`         Chain pair, for example `A-B`
  `frame_i`           First frame identifier, stored as `frame_<index>`
  `frame_j`           Second frame identifier, stored as `frame_<index>`
  `shared_fraction`   Shared fraction between the two contact sets

The comparison includes all ordered frame pairs generated by the code,
including self-comparisons.

------------------------------------------------------------------------

## Interface RMSD (iRMSD)

### Purpose

`compute_irmsd()` calculates the interface RMSD for pairs of frames.

For each chain pair, the method:

1.  identifies backbone atom pairs between the two chains;
2.  calculates interatomic distances;
3.  selects atom pairs within a cutoff of `1.0 nm`;
4.  collects the backbone atoms involved in at least one selected pair
    for each frame;
5.  finds the atoms common to the interfaces of two frames;
6.  aligns the selected frame structures using those common atoms;
7.  calculates the RMSD over the common atoms.

The RMSD is converted from nanometres to ångströms.

### Usage

``` python
analyzer.compute_irmsd()
```

There are no user-facing parameters in the current method.

### Interface cutoff

The current implementation uses:

``` python
distance_cutoff = 1.0  # nm
```

This is a fixed 1.0 nm cutoff between backbone atoms of the two chains.
It is not the same cutoff as the contact interval used for SF.

### Weighted and unweighted iRMSD

The method stores two values:

-   `irmsd_unweighted`: the calculated iRMSD in Å;
-   `irmsd_weighted`: the iRMSD multiplied by a penalty derived from the
    fraction of common interface atoms.

The common-atom fraction is:

\[ f(i,j) = `\frac{2|A_i \cap A_j|}{|A_i| + |A_j|}`{=tex} \]

The penalty applied in the implementation is:

\[ L(i,j) = -`\log`{=tex}\_2(f(i,j)\^2 + 10\^{-15}) + 1 \]

and:

\[ iRMSD\_{`\mathrm{weighted}`{=tex}} =
iRMSD\_{`\mathrm{unweighted}`{=tex}} `\times `{=tex}L(i,j) \]

Here, (A_i) and (A_j) are the sets of interface backbone atom indices
for the two frames.

### Output

The method stores its results in:

``` python
analyzer.df_irmsd
```

Columns:

  Column               Description
  -------------------- --------------------------------
  `interface`          Chain pair
  `frame_i`            First frame index
  `frame_j`            Second frame index
  `irmsd_unweighted`   Unweighted interface RMSD in Å
  `irmsd_weighted`     Weighted interface RMSD in Å

It also creates:

``` python
analyzer.df_interface_backbone
```

This intermediate DataFrame contains the interface backbone atom indices
identified for each interface and frame. The lRMSD method uses this
information.

------------------------------------------------------------------------

## Ligand RMSD (lRMSD)

### Purpose

`compute_lrmsd()` calculates the ligand RMSD between trajectory frames
for each chain interface.

The method assigns the longer chain as the **receptor** and the shorter
chain as the **ligand**. If both chains have the same number of
residues, the implementation uses the alphabetical ordering of their
labels to decide which is the receptor and which is the ligand.

For each interface and frame, the method extracts the receptor backbone
atoms that belong to the interface. It then uses receptor-based
alignment and calculates the RMSD of the ligand backbone.

### Usage

``` python
analyzer.compute_lrmsd()
```

Run `compute_irmsd()` first: `compute_lrmsd()` relies on the
`df_interface_backbone` attribute created by that method.

### Output

The method stores the results in:

``` python
analyzer.df_lrmsd
```

Columns:

  Column               Description
  -------------------- -----------------------------
  `interface`          Chain pair
  `frame_i`            Reference frame index
  `frame_j`            Compared frame index
  `lrmsd_unweighted`   Unweighted ligand RMSD in Å
  `lrmsd_weighted`     Weighted ligand RMSD in Å

The weighted value uses the same logarithmic penalty form as the iRMSD,
with the common receptor-interface atoms used to calculate the fraction.

------------------------------------------------------------------------

## MDQ score

### Purpose

`compute_mdq()` combines the three metrics:

-   shared fraction;
-   iRMSD;
-   lRMSD.

It first merges the three result tables using `interface`, `frame_i` and
`frame_j`. It then calculates a weighted and an unweighted MDQ score.

### Usage

``` python
analyzer.compute_mdq()
```

The method expects the metric DataFrames to have been calculated first:

``` python
analyzer.compute_shared_fraction()
analyzer.compute_irmsd()
analyzer.compute_lrmsd()
analyzer.compute_mdq()
```

The current implementation also reads `shared_fraction.csv` from the
working directory during this step.

### Formula

The weighted score is:

\[ MDQ\_{`\mathrm{weighted}`{=tex}} = `\frac{
SF +
\frac{1}{1 + (iRMSD_{\mathrm{weighted}}/1.5)^2} +
\frac{1}{1 + (lRMSD_{\mathrm{weighted}}/8.5)^2}
}{3}`{=tex} \]

The unweighted score is:

\[ MDQ\_{`\mathrm{unweighted}`{=tex}} = `\frac{
SF +
\frac{1}{1 + (iRMSD_{\mathrm{unweighted}}/1.5)^2} +
\frac{1}{1 + (lRMSD_{\mathrm{unweighted}}/8.5)^2}
}{3}`{=tex} \]

The constants `1.5` and `8.5` are the scaling values used by the current
implementation for iRMSD and lRMSD, respectively.

Each component contributes one third of the final score. The RMSD
contributions are transformed so that smaller RMSD values produce larger
similarity contributions.

### Output

The merged results are stored in:

``` python
analyzer.df_score
```

The output contains the columns from the three metric tables and:

-   `mdq_weighted`;
-   `mdq_unweighted`.

**Important:** in the uploaded implementation, `compute_mdq()` does not
write `analyzer.df_score` to `mdq.csv` automatically. Save it explicitly
if you want to use the visualiser's default input:

``` python
analyzer.df_score.to_csv("mdq.csv", index=False)
```

------------------------------------------------------------------------

# 2. Cluster trajectory frames

## Create a distance matrix from MDQ

The clustering class uses a distance rather than a similarity. Its
helper method `create_distance_from_mdq()` converts the weighted MDQ
score into a distance:

\[ D(i,j) = 1 - MDQ\_{`\mathrm{weighted}`{=tex}}(i,j) \]

Thus, high MDQ similarity corresponds to low distance.

The method expects a file named:

``` text
mdq.csv
```

with at least these columns:

-   `interface`;
-   `frame_i`;
-   `frame_j`;
-   `mdq_weighted`.

It writes:

``` text
distance.csv
```

with the columns:

-   `interface`;
-   `frame_i`;
-   `frame_j`;
-   `distance`.

### Initialise the visualiser

If `distance.csv` already exists:

``` python
from mdq import TrajectoryClusterVisualizer

visualizer = TrajectoryClusterVisualizer(
    distance_file="distance.csv"
)
```

If no distance file is supplied, the constructor calls
`create_distance_from_mdq()` and reads the generated `distance.csv`:

``` python
visualizer = TrajectoryClusterVisualizer()
```

In this case, `mdq.csv` must already exist in the current working
directory.

------------------------------------------------------------------------

## DBSCAN clustering

### Purpose

DBSCAN groups frames based on the precomputed pairwise distance matrix.
Frames that do not belong to a cluster are labelled `-1` (noise).

The implementation processes each interface separately.

### Usage

The source currently defines the method with these intended arguments:

``` python
visualizer.run_dbscan_analysis(
    eps_dict=eps_dict,
    min_samples_pct=0.05,
    min_pct_dyn=0.05,
    total_frames_dyn=1000
)
```

This example shows the intended call after correcting the parameter
order in the method definition.

### Parameters

  ------------------------------------------------------------------------
  Parameter                        Default in source Description
  --------------------- ---------------------------- ---------------------
  `eps_dict`                                Required Dictionary assigning
                                                     an `eps` distance
                                                     threshold to each
                                                     interface

  `min_samples_pct`                           `0.05` Fraction of the
                                                     number of frames used
                                                     to calculate DBSCAN's
                                                     `min_samples`

  `min_pct_dyn`                               `0.05` Minimum cluster-size
                                                     fraction used by the
                                                     cluster filtering
                                                     step

  `total_frames_dyn`                        Required Denominator used by
                                                     the cluster filtering
                                                     step
  ------------------------------------------------------------------------

Example interface-specific thresholds:

``` python
eps_dict = {
    "A-B": 0.20,
    "A-C": 0.25,
    "B-C": 0.20
}
```

For each interface, the method:

1.  builds a distance matrix;
2.  selects `eps` from `eps_dict`, using `0.2` if the interface is not
    present;
3.  calculates `min_samples` as
    `max(1, int(min_samples_pct * n_frames))`;
4.  runs DBSCAN with `metric="precomputed"`;
5.  filters clusters smaller than the specified minimum fraction;
6.  stores the frames, distance matrix and labels;
7.  plots MDS, UMAP and t-SNE projections.

Results are stored in:

``` python
visualizer.dbscan_labels_dict
```

Each interface maps to a tuple:

``` python
(frames, distance_matrix, labels)
```

### Understanding `eps`

`eps` is a distance threshold, not an MDQ similarity threshold. Since
the distance is defined as `1 - MDQ_weighted`, smaller `eps` values
impose a stricter similarity requirement.

The appropriate value depends on the score distribution and should be
assessed for each system/interface.

------------------------------------------------------------------------

## Hierarchical clustering

### Purpose

Hierarchical clustering groups frames using an agglomerative clustering
procedure applied to the precomputed distance matrix.

### Usage

``` python
visualizer.run_hierarchical_analysis(
    method="average",
    threshold=0.3,
    min_pct=0.05,
    total_frames=1000
)
```

This is the intended call after correcting the parameter order in the
method definition.

### Parameters

  ------------------------------------------------------------------------
  Parameter                        Default in source Description
  --------------------- ---------------------------- ---------------------
  `method`                               `"average"` Linkage method passed
                                                     to SciPy

  `threshold`                                  `0.3` Distance threshold
                                                     used to cut the
                                                     dendrogram

  `min_pct`                                   `0.05` Minimum cluster-size
                                                     fraction

  `total_frames`                            Required Denominator used by
                                                     the cluster filtering
                                                     step
  ------------------------------------------------------------------------

The method:

1.  constructs the distance matrix for an interface;
2.  converts it to condensed form with `squareform`;
3.  calculates the linkage matrix with SciPy;
4.  plots a dendrogram;
5.  assigns clusters using `fcluster(..., criterion="distance")`;
6.  relabels clusters sequentially;
7.  filters clusters below the minimum size;
8.  stores the results;
9.  plots MDS, UMAP and t-SNE projections.

Results are stored in:

``` python
visualizer.hierarchical_labels_dict
```

As with DBSCAN, each interface maps to:

``` python
(frames, distance_matrix, labels)
```

------------------------------------------------------------------------

## Projection plots

Both clustering methods call `_plot_projections()` to visualise the
distance matrix in two dimensions using:

-   **MDS** (Multidimensional Scaling);
-   **UMAP**;
-   **t-SNE**.

The projections are generated from the precomputed distance matrix.
Points are coloured by cluster label. Noise points (`-1`) are shown in
light grey.

The representative frame of each cluster is highlighted with a larger
marker and a black outline. Frame indices are annotated on the plots.

These projections are visual aids. They do not replace the clustering
performed on the original distance matrix, and distances measured in a
2D projection should not be interpreted as the original MDQ distances.

------------------------------------------------------------------------

## Distance heatmaps and clustermaps

### `plot_distance_heatmap()`

Displays the pairwise distance matrix for each interface as a heatmap.

``` python
visualizer.plot_distance_heatmap()
```

The colour scale is set from 0 to 1. Missing matrix entries are
represented using the colormap's bad-value colour.

**Source note:** `mdq.py` contains two definitions named
`plot_distance_heatmap`. The second definition overrides the first in
Python. Consequently, the effective method does not accept `start_frame`
or `end_frame`, even though these arguments appear in the earlier
definition.

### `plot_distance_clustermap()`

``` python
visualizer.plot_distance_clustermap()
```

Displays a seaborn clustermap for each interface. Missing entries are
replaced with `-1` for plotting and the colormap is configured with a
lower-range colour.

------------------------------------------------------------------------

# 3. Export and summarise clusters

## Get cluster frames and representatives

``` python
dbscan_table, hierarchical_table = visualizer.get_cluster_frames_table()
```

Returns two nested dictionaries, one for each clustering method. For
each interface and cluster, the dictionary contains:

-   `frames`: frame indices assigned to the cluster;
-   `representative`: representative frame index.

Noise points labelled `-1` are excluded.

The representative is selected as the frame with the smallest mean
distance to the other frames in its cluster (a medoid-like
representative).

An alternative method returns both clustering methods in a named
dictionary:

``` python
representatives = visualizer.get_cluster_representatives()
```

Its structure is:

``` python
{
    "DBSCAN": {
        "A-B": {
            0: {
                "frames": [0, 1, 2],
                "representative": 1
            }
        }
    },
    "Hierarchical": {
        # ...
    }
}
```

The values above are illustrative only.

## Export representative structures as PDB

``` python
visualizer.export_representative_frames_pdb(
    traj_file="production.xtc",
    top_file="topology.pdb",
    output_dir="representative_structures",
    cluster_tables=(dbscan_table, hierarchical_table)
)
```

The method reloads the trajectory and exports one PDB file per
representative frame.

The output filenames follow this pattern:

``` text
DBSCAN_<interface>_cluster<cluster_label>_rep.pdb
Hierarchical_<interface>_cluster<cluster_label>_rep.pdb
```

**Frame-index caution:** the representative indices come from the
trajectory loaded by the visualiser. The trajectory passed to the export
method must have the same frame ordering and sampling. If the analysis
used a non-default stride, ensure that the export trajectory is sampled
consistently; otherwise the exported frame may not correspond to the
analysed frame.

## Export cluster representatives as CSV

``` python
visualizer.export_cluster_representatives_csv(
    output_file="cluster_representatives.csv"
)
```

The CSV contains:

  Column                   Description
  ------------------------ ---------------------------------------
  `method`                 `DBSCAN` or `Hierarchical`
  `interface`              Interface identifier
  `cluster`                Cluster label
  `n_frames`               Number of frames in the cluster
  `representative_frame`   Representative frame index
  `frames`                 Comma-separated list of frame indices

## Save cluster statistics

``` python
visualizer.save_cluster_stats_csv(
    output_prefix="cluster_stats"
)
```

The method writes one CSV per clustering method:

``` text
cluster_stats_DBSCAN.csv
cluster_stats_Hierarchical.csv
```

For each interface and non-noise cluster, the method calculates
statistics from the mean within-cluster distance of each member frame to
the other frames in that cluster:

-   number of frames;
-   mean distance;
-   standard deviation;
-   first quartile (Q1);
-   median;
-   third quartile (Q3);
-   lower and upper bounds of a 95% Student's t confidence interval.

For a cluster containing only one frame, the sample standard deviation
and confidence interval are not meaningfully defined; inspect such cases
when interpreting the exported statistics.

------------------------------------------------------------------------

# Example end-to-end workflow

The following example illustrates the intended sequence. Adjust the file
paths, chain indices, thresholds and frame counts to your system.

``` python
from mdq import TrajectoryAnalyzer, TrajectoryClusterVisualizer

# --------------------------------------------------
# 1. Input files and chain definitions
# --------------------------------------------------

traj_file = "production.xtc"
top_file = "topology.pdb"

chains = {
    "A": list(range(0, 314)),
    "B": list(range(314, 634))
}

# --------------------------------------------------
# 2. Calculate interface metrics
# --------------------------------------------------

analyzer = TrajectoryAnalyzer(
    traj_file=traj_file,
    top_file=top_file,
    chains=chains,
    stride=1
)

analyzer.compute_shared_fraction(min_A=2, max_A=5)
analyzer.compute_irmsd()
analyzer.compute_lrmsd()
analyzer.compute_mdq()

# Save the MDQ table for the visualisation workflow.
analyzer.df_score.to_csv("mdq.csv", index=False)

# --------------------------------------------------
# 3. Create the pairwise distance table
# --------------------------------------------------

visualizer = TrajectoryClusterVisualizer()
# The default constructor reads mdq.csv and creates distance.csv.

# --------------------------------------------------
# 4. Run clustering
# --------------------------------------------------

eps_dict = {
    "A-B": 0.20
}

# The method signatures in mdq.py must first be corrected
# so that required arguments precede arguments with defaults.
visualizer.run_dbscan_analysis(
    eps_dict=eps_dict,
    min_samples_pct=0.05,
    min_pct_dyn=0.05,
    total_frames_dyn=1000
)

visualizer.run_hierarchical_analysis(
    method="average",
    threshold=0.3,
    min_pct=0.05,
    total_frames=1000
)

# --------------------------------------------------
# 5. Export representative frames and cluster results
# --------------------------------------------------

dbscan_table, hierarchical_table = visualizer.get_cluster_frames_table()

visualizer.export_representative_frames_pdb(
    traj_file=traj_file,
    top_file=top_file,
    output_dir="representative_structures",
    cluster_tables=(dbscan_table, hierarchical_table)
)

visualizer.export_cluster_representatives_csv(
    "cluster_representatives.csv"
)

visualizer.save_cluster_stats_csv(
    output_prefix="cluster_stats"
)
```

------------------------------------------------------------------------

## Output files

Depending on the methods used, the workflow creates or expects the
following files:

  ---------------------------------------------------------------------------------------------------
  File                               Created by / used by                     Description
  ---------------------------------- ---------------------------------------- -----------------------
  `shared_fraction.csv`              `compute_shared_fraction()`              Pairwise
                                                                              shared-fraction results

  `mdq.csv`                          Save `analyzer.df_score` explicitly      Pairwise MDQ results;
                                                                              expected by the default
                                                                              visualiser workflow

  `distance.csv`                     `create_distance_from_mdq()`             Pairwise distance table
                                                                              calculated as
                                                                              `1 - mdq_weighted`

  `representative_structures/`       `export_representative_frames_pdb()`     Representative PDB
                                                                              structures

  `cluster_representatives.csv`      `export_cluster_representatives_csv()`   Cluster membership and
                                                                              representative-frame
                                                                              table

  `cluster_stats_DBSCAN.csv`         `save_cluster_stats_csv()`               DBSCAN cluster
                                                                              statistics

  `cluster_stats_Hierarchical.csv`   `save_cluster_stats_csv()`               Hierarchical cluster
                                                                              statistics
  ---------------------------------------------------------------------------------------------------

The metric DataFrames (`df_shared`, `df_irmsd`, `df_lrmsd` and
`df_score`) are also available in memory on the `TrajectoryAnalyzer`
instance.

------------------------------------------------------------------------

## Interpretation

### Similarity and distance

The MDQ score is a similarity measure in the current implementation:
higher values indicate greater similarity according to the three
combined components.

The clustering workflow converts the weighted score to a distance:

\[ D = 1 - MDQ\_{`\mathrm{weighted}`{=tex}} \]

Therefore:

-   low distance corresponds to high MDQ similarity;
-   high distance corresponds to low MDQ similarity.

### Weighted versus unweighted metrics

The unweighted MDQ uses the raw iRMSD and lRMSD values in the two RMSD
transformations. The weighted MDQ uses the corresponding values
multiplied by the logarithmic interface-overlap penalty.

Compare both versions when assessing whether interface conservation
changes the interpretation of structural similarity.

### Cluster representatives

The representative frame is selected using the smallest average distance
to the other frames in the cluster. It is a medoid-like representative,
not necessarily the most populated frame or the frame with the lowest
RMSD to a single reference structure.

### Clustering parameters

DBSCAN's `eps`, `min_samples_pct` and cluster-size filtering threshold
affect the resulting clusters. Hierarchical clustering depends on the
linkage method and dendrogram distance threshold. These parameters
should be selected based on the distance distributions and the
scientific question, and reported alongside the results.

------------------------------------------------------------------------

## Implementation notes and limitations

The following points describe the uploaded implementation and should be
checked before using it for production analyses:

1.  **Invalid method signatures:** as noted at the beginning, the
    current DBSCAN and hierarchical method definitions place required
    parameters after parameters with defaults. Python will raise a
    `SyntaxError` when importing the file until these signatures are
    corrected.
2.  **Duplicate heatmap method:** `plot_distance_heatmap()` is defined
    twice. The second definition overrides the first, so the version
    with `start_frame` and `end_frame` is not available at runtime.
3.  **MDQ CSV export:** `compute_mdq()` stores the result in
    `self.df_score` but does not export `mdq.csv`. Save the DataFrame
    explicitly before using `TrajectoryClusterVisualizer()` with its
    default argument.
4.  **Distance-column condition:** `create_distance_from_mdq()` checks
    for `mdq_weighted` twice. The second check is identical to the
    first; there is no separate fallback for another MDQ column in the
    current code.
5.  **Frame sampling:** frame indices refer to the trajectory as loaded,
    including the selected stride. Keep the same frame ordering when
    exporting representative structures.
6.  **Cluster-size filtering:** the filtering helper compares cluster
    size to the `total_frames` value supplied to the clustering method.
    Supply a denominator consistent with the intended definition of the
    minimum cluster fraction.
7.  **Pairwise scaling:** pairwise frame calculations scale
    approximately as (N\^2), where (N) is the number of loaded frames.
    Large trajectories can generate large DataFrames and long runtimes.
8.  **Backbone definition:** only atoms named `N`, `CA`, `C` and `O` are
    included in the backbone selections.
9.  **Receptor/ligand assignment:** the longer chain is designated as
    receptor and the shorter as ligand. Equal-length chains are resolved
    alphabetically. This is an implementation convention, not an
    automatic biological annotation.
10. **MDQ terminology:** MDQ is the name used for the combined score in
    this project. The code itself does not define an expanded form of
    the acronym; provide the intended definition in the associated
    manuscript or project documentation if needed.

------------------------------------------------------------------------

## Troubleshooting

### `SyntaxError: non-default argument follows default argument`

Correct the argument order in `run_dbscan_analysis()` and
`run_hierarchical_analysis()`. Required arguments must appear before
arguments with default values.

### `FileNotFoundError: mdq.csv`

`compute_mdq()` does not save this file automatically. Export it:

``` python
analyzer.df_score.to_csv("mdq.csv", index=False)
```

Run the visualiser from the directory containing `mdq.csv`, or provide
an explicit `distance_file` after creating `distance.csv`.

### Chain or residue selections do not match the expected structure

Check that `chains` contains zero-based MDTraj residue indices from the
topology actually loaded with the trajectory. Do not assume that PDB
residue numbers and MDTraj residue indices are interchangeable.

### Representative PDB frames do not match the analysed frames

Make sure the trajectory used for export has the same frame order and
sampling as the trajectory analysed by `TrajectoryAnalyzer`. In
particular, account for any non-default `stride`.

### Unexpected cluster counts

Inspect the distance heatmap and dendrogram, and review the DBSCAN `eps`
/ `min_samples_pct` or hierarchical `threshold` / `method`. Also check
the value passed as `total_frames` or `total_frames_dyn`, because it is
used by the minimum-cluster-size filter.

------------------------------------------------------------------------

## Citation and reproducibility

When reporting results generated with this code, document at least:

-   trajectory and topology files (including preparation details);
-   chain definitions and residue-index convention;
-   trajectory stride and number of analysed frames;
-   SF distance interval (`min_A`, `max_A`);
-   iRMSD interface cutoff (currently 1.0 nm);
-   MDQ formula and scaling constants (currently 1.5 Å and 8.5 Å);
-   whether weighted or unweighted MDQ was used;
-   distance definition (`1 - mdq_weighted`);
-   clustering algorithm and all its parameters;
-   software versions and any modifications made to `mdq.py`.

This information is necessary for another user to reproduce the
analysis.
