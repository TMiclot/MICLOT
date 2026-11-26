import mdtraj as md
import itertools
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns
import ast
import warnings

from tqdm import tqdm
from sklearn_extra.cluster import KMedoids
from sklearn.cluster import DBSCAN
from sklearn.manifold import MDS, TSNE
from scipy.cluster.hierarchy import linkage, dendrogram, fcluster
from scipy.spatial.distance import squareform
from scipy import stats
import umap


"""
Trajectory analysis module for molecular dynamics simulations.

This module provides the TrajectoryAnalyzer class to compute various metrics
on protein-protein or protein-chain interfaces from MD trajectories:
shared fraction, interface RMSD (iRMSD), ligand RMSD (lRMSD), and a MDQ score.
"""

class TrajectoryAnalyzer:
    """
    Analyze molecular dynamics trajectories focusing on inter-chain interfaces.

    The class loads a trajectory, extracts backbone atoms per chain, and provides
    methods to compute:
        - Shared fraction of contacts of residues between frames.
        - Interface RMSD (iRMSD) weighted/unweighted.
        - Ligand RMSD (lRMSD) weighted/unweighted.
        - MDQ score combining the above.

    Attributes:
        traj_file (str): Path to trajectory file (e.g., .xtc, .dcd).
        top_file (str): Path to topology file (e.g., .pdb, .gro).
        chains (dict): Mapping chain identifiers to lists of residue indices.
        stride (int): Frame stride for loading the trajectory.
        traj (md.Trajectory): Loaded trajectory object.
        backbone_indices (dict): Atom indices of backbone atoms per chain.
        slice_frame (dict): Dictionary mapping frame index to a trajectory slice.
        df_shared (pd.DataFrame): Shared fraction results.
        df_irmsd (pd.DataFrame): iRMSD results.
        df_lrmsd (pd.DataFrame): lRMSD results.
        df_score (pd.DataFrame): Final MDQ scores.
    """

    def __init__(self, traj_file, top_file, chains, stride=1):
        """
        Initialize the trajectory analyzer.

        Args:
            traj_file (str): Path to trajectory file.
            top_file (str): Path to topology file.
            chains (dict): Dictionary where keys are chain IDs (str) and values are
                           lists of residue indices belonging to that chain.
                           Example: {'A': [0,1,2], 'B': [100,101,102]}
            stride (int, optional): Read every `stride` frame. Defaults to 1.
        """
        self.traj_file = traj_file
        self.top_file = top_file
        self.chains = chains
        self.stride = stride

        print("Loading trajectory...")
        self.traj = md.load(traj_file, top=top_file, stride=stride)

        # -----------------------------
        # Backbone atoms per chain
        # -----------------------------
        backbone_atoms = {'N', 'CA', 'C', 'O'}
        top = self.traj.topology

        # Build residue -> chain mapping
        res_to_chain = {}
        for chain_id, res_range in chains.items():
            for r in res_range:
                res_to_chain[r] = chain_id

        self.backbone_indices = {c: [] for c in chains}

        for atom in top.atoms:
            if atom.name in backbone_atoms:
                chain_id = res_to_chain.get(atom.residue.index)
                if chain_id:
                    self.backbone_indices[chain_id].append(atom.index)

        # Pre-slice trajectory per frame for faster access
        self.slice_frame = {
            i: self.traj.slice(i) for i in range(self.traj.n_frames)
        }

    # -----------------------------
    # Shared fraction
    # -----------------------------
    def compute_shared_fraction(self, min_A=2, max_A=5):
        """
        Compute the shared fraction of contacts between all pairs of frames.

        For each interface (pair of chains) and each pair of frames, the fraction
        of common contacting residue pairs is calculated.

        Args:
            min_A (float, optional): Minimum contact distance in Angstroms.
                                    Defaults to 2.
            max_A (float, optional): Maximum contact distance in Angstroms.
                                    Defaults to 5.

        Returns:
            None. Results are stored in `self.df_shared` and saved to a CSV file.
        """
        print("Computing shared fraction...")

        residue_indices = {chain: list(ids) for chain, ids in self.chains.items()}
        interfaces = list(itertools.combinations(self.chains.keys(), 2))
        df_list = []

        for chain1, chain2 in tqdm(interfaces, desc="Interfaces"):
            pairs = list(itertools.product(
                residue_indices[chain1],
                residue_indices[chain2]
            ))

            if len(pairs) == 0:
                continue

            distances, _ = md.compute_contacts(
                self.traj,
                contacts=pairs,
                scheme='closest-heavy'
            )

            if distances.shape[0] < distances.shape[1]:
                distances = distances.T

            n_pairs, n_frames = distances.shape
            distances_A = distances * 10

            contacts = (distances_A >= min_A) & (distances_A <= max_A)

            for i, j in itertools.product(range(n_frames), repeat=2):

                idx_i = np.where(contacts[:, i])[0]
                idx_j = np.where(contacts[:, j])[0]

                if len(idx_i) + len(idx_j) == 0:
                    frac = 0.0
                else:
                    frac = 2 * len(np.intersect1d(idx_i, idx_j)) / (len(idx_i) + len(idx_j))

                df_list.append({
                    "interface": f"{chain1}-{chain2}",
                    "frame_i": f"frame_{i}",
                    "frame_j": f"frame_{j}",
                    "shared_fraction": frac
                })

        self.df_shared = pd.DataFrame(df_list)

        self.df_shared.to_csv(
            "shared_fraction.csv",
            index=False
        )

        print("Shared fraction saved → shared_fraction.csv")

    # -----------------------------
    # iRMSD (interface RMSD)
    # -----------------------------
    def compute_irmsd(self):
        """
        Compute interface RMSD (iRMSD) for all frame pairs.

        iRMSD is calculated on the set of backbone atoms that are within a cutoff
        distance (1.0 nm) between two chains. Both unweighted and weighted (by
        fraction of common atoms) versions are computed.

        Returns:
            None. Results are stored in `self.df_irmsd`.
        """
        print("Computing iRMSD...")

        interfaces = list(itertools.combinations(self.chains.keys(), 2))
        distance_cutoff = 1.0  # nm
        rows = []

        for c1, c2 in tqdm(interfaces, desc="iRMSD interfaces"):
            pair_name = f"{c1}-{c2}"

            idx1 = self.backbone_indices[c1]
            idx2 = self.backbone_indices[c2]

            pairs = np.array(np.meshgrid(idx1, idx2)).T.reshape(-1, 2)

            dist_array = md.compute_distances(self.traj, pairs)
            mask = dist_array <= distance_cutoff

            frames_idx, pairs_idx = np.where(mask)

            atom_i = pairs[pairs_idx, 0]
            atom_j = pairs[pairs_idx, 1]

            df_tmp = pd.DataFrame({
                "pair_chain": pair_name,
                "frame": frames_idx,
                "atom_i": atom_i,
                "atom_j": atom_j
            })

            rows.append(df_tmp)

        df = pd.concat(rows, ignore_index=True)
        rows = []

        # Aggregate backbone atoms per interface and frame
        for (interface, frame), group in df.groupby(["pair_chain", "frame"]):
            atoms = set(group["atom_i"]).union(set(group["atom_j"]))
            rows.append({
                "interface": interface,
                "frame": frame,
                "backbone_atoms": sorted(list(atoms))
            })

        self.df_interface_backbone = pd.DataFrame(rows)

        # Precompute common atoms and fraction between frames
        commun = {}
        fraction = {}

        for interface in self.df_interface_backbone['interface'].unique():
            df_interface = self.df_interface_backbone[
                self.df_interface_backbone['interface'] == interface
            ]

            frames = df_interface['frame'].values

            for i, j in itertools.product(frames, repeat=2):

                atoms_i = df_interface[df_interface['frame'] == i]['backbone_atoms'].values[0]
                atoms_j = df_interface[df_interface['frame'] == j]['backbone_atoms'].values[0]

                common = list(set(atoms_i).intersection(atoms_j))

                commun[(interface, i, j)] = common

                if len(atoms_i) + len(atoms_j) == 0:
                    fraction[(interface, i, j)] = 0
                else:
                    fraction[(interface, i, j)] = (
                        2 * len(common) / (len(atoms_i) + len(atoms_j))
                    )

        # Compute RMSD (both unweighted and weighted)
        rows = []

        for interface in self.df_interface_backbone['interface'].unique():
            df_interface = self.df_interface_backbone[
                self.df_interface_backbone['interface'] == interface
            ]

            frames = df_interface['frame'].tolist()

            for i, j in itertools.product(frames, repeat=2):

                common_atoms = commun.get((interface, i, j), [])

                if len(common_atoms) == 0:
                    continue

                traj_i = self.slice_frame[i]
                traj_j = self.slice_frame[j]

                traj_i.superpose(traj_j, atom_indices=common_atoms)

                rmsd_unweighted = md.rmsd(
                    traj_i, traj_j, atom_indices=common_atoms
                )[0] * 10  # convert to Angstroms

                loss = -np.log2(fraction[(interface, i, j)]**2 + 1e-15) + 1

                rmsd_weighted = rmsd_unweighted * loss

                rows.append({
                    "interface": interface,
                    "frame_i": i,
                    "frame_j": j,
                    "irmsd_unweighted": rmsd_unweighted,
                    "irmsd_weighted": rmsd_weighted
                })

        self.df_irmsd = pd.DataFrame(rows)

    # -----------------------------
    # lRMSD (ligand RMSD)
    # -----------------------------
    def compute_lrmsd(self):
        """
        Compute ligand RMSD (lRMSD) for all frame pairs.

        The longer chain is considered the receptor, the shorter the ligand.
        lRMSD is computed after aligning the receptor backbone atoms that are in
        contact. Both unweighted and weighted versions are stored.

        Returns:
            None. Results are stored in `self.df_lrmsd`.
        """
        print("Computing lRMSD...")

        df = self.df_interface_backbone.copy()

        if df["backbone_atoms"].dtype == object:
            df["backbone_atoms"] = df["backbone_atoms"].apply(
                lambda x: ast.literal_eval(x) if isinstance(x, str) else x
            )

        chain_lengths = {
            chain: len(res_range)
            for chain, res_range in self.chains.items()
        }

        def get_receptor_ligand(interface):
            c1, c2 = interface.split("-")
            if chain_lengths[c1] > chain_lengths[c2]:
                return c1, c2
            elif chain_lengths[c2] > chain_lengths[c1]:
                return c2, c1
            else:
                return sorted([c1, c2])[1], sorted([c1, c2])[0]

        receptor_atoms_per_frame = {}

        for (interface, frame), group in df.groupby(["interface", "frame"]):
            receptor, _ = get_receptor_ligand(interface)

            receptor_atoms = set(self.backbone_indices[receptor])
            atoms_interface = set().union(*group["backbone_atoms"])

            receptor_atoms_per_frame[(interface, frame)] = list(
                atoms_interface.intersection(receptor_atoms)
            )

        rows = []

        interfaces = list(itertools.combinations(self.chains.keys(), 2))

        for interface in tqdm(interfaces, desc="lRMSD"):

            interface_str = f"{interface[0]}-{interface[1]}"
            receptor, ligand = get_receptor_ligand(interface_str)

            ligand_idx = self.backbone_indices[ligand]

            for i in range(self.traj.n_frames):

                receptor_idx = receptor_atoms_per_frame.get((interface_str, i), [])

                if len(receptor_idx) == 0:
                    continue

                traj_aligned = self.traj.superpose(
                    self.traj, frame=i, atom_indices=receptor_idx
                )

                rmsd_vals = md.rmsd(
                    traj_aligned,
                    traj_aligned,
                    frame=i,
                    atom_indices=ligand_idx
                ) * 10

                for j, rmsd_unweighted in enumerate(rmsd_vals):

                    atoms_i = receptor_atoms_per_frame.get((interface_str, i), [])
                    atoms_j = receptor_atoms_per_frame.get((interface_str, j), [])

                    if len(atoms_i) + len(atoms_j) == 0:
                        frac = 0
                    else:
                        frac = 2 * len(set(atoms_i).intersection(atoms_j)) / (
                            len(atoms_i) + len(atoms_j)
                        )

                    loss = -np.log2(frac**2 + 1e-15) + 1

                    rmsd_weighted = rmsd_unweighted * loss

                    rows.append((
                        interface_str,
                        i,
                        j,
                        rmsd_unweighted,
                        rmsd_weighted
                    ))

        self.df_lrmsd = pd.DataFrame(rows, columns=[
            "interface",
            "frame_i",
            "frame_j",
            "lrmsd_unweighted",
            "lrmsd_weighted"
        ])

    # -----------------------------
    # MDQ score
    # -----------------------------
    def compute_mdq(self):
        """
        Compute a MDQ score combining shared fraction, iRMSD, and lRMSD.

        This method merges the previously computed dataframes (shared fraction,
        iRMSD, lRMSD) and calculates a weighted and unweighted score.

        The score ranges from 0 (poor similarity) to 1 (perfect similarity).

        Returns:
            None. Results are stored in `self.df_score`.
        """
        print("Computing MDQ score...")

        self.df_shared = pd.read_csv(
            "shared_fraction.csv"
        )

        def fix_frames(df):
            for col in ["frame_i", "frame_j"]:
                if col in df.columns:
                    df[col] = (
                        df[col].astype(str)
                        .str.replace("frame_", "", regex=False)
                        .astype(int)
                    )
            return df

        self.df_shared = fix_frames(self.df_shared)
        self.df_irmsd = fix_frames(self.df_irmsd)
        self.df_lrmsd = fix_frames(self.df_lrmsd)

        df = self.df_shared.merge(
            self.df_irmsd,
            on=["interface", "frame_i", "frame_j"]
        ).merge(
            self.df_lrmsd,
            on=["interface", "frame_i", "frame_j"]
        )

        df["mdq_weighted"] = (
            df["shared_fraction"]
            + 1 / (1 + (df["irmsd_weighted"] / 1.5)**2)
            + 1 / (1 + (df["lrmsd_weighted"] / 8.5)**2)
        ) / 3

        df["mdq_unweighted"] = (
            df["shared_fraction"]
            + 1 / (1 + (df["irmsd_unweighted"] / 1.5)**2)
            + 1 / (1 + (df["lrmsd_unweighted"] / 8.5)**2)
        ) / 3

        self.df_score = df

        print("MDQ score computed")

# ---------------------------------------------------------------------------------------------------------------------------------------------

class TrajectoryClusterVisualizer:
    """Cluster and visualize trajectory frames from a distance matrix."""

    def __init__(self, distance_file=None):
        if distance_file is None:
            self.create_distance_from_mdq()
            distance_file = "distance.csv"
        self.df = pd.read_csv(distance_file)
        for col in ["frame_i", "frame_j"]:
            if self.df[col].dtype == object:
                self.df[col] = self.df[col].str.replace("frame_", "", regex=False).astype(int)
        self.dbscan_labels_dict = {}
        self.hierarchical_labels_dict = {}

    # ===============================
    # UTILS
    # ===============================
    def create_distance_from_mdq(self):
        df = pd.read_csv("mdq.csv")
        df_new = df[["interface", "frame_i", "frame_j"]].copy()
        if "mdq_weighted" in df.columns:
            df_new["distance"] = 1 - df["mdq_weighted"]
        elif "mdq_weighted" in df.columns:
            df_new["distance"] = 1 - df["mdq_weighted"]
        else:
            raise ValueError("MDQ column not found")
        df_new.to_csv("distance.csv", index=False)
        print("CSV created: distance.csv")

    def _build_distance_matrix(self, df_int):
        frames = sorted(set(df_int["frame_i"]).union(set(df_int["frame_j"])))
        n = len(frames)
        frame_to_idx = {f: i for i, f in enumerate(frames)}
        dist_matrix = np.zeros((n, n))
        for _, row in df_int.iterrows():
            i, j = frame_to_idx[row["frame_i"]], frame_to_idx[row["frame_j"]]
            dist_matrix[i, j] = dist_matrix[j, i] = row["distance"]
        np.fill_diagonal(dist_matrix, 0)
        return frames, dist_matrix

    def _get_representative(self, dist_matrix, idx):
        sub = dist_matrix[np.ix_(idx, idx)]
        mean_dist = sub.mean(axis=1)
        return idx[np.argmin(mean_dist)]

    def _filter_clusters(self, labels, min_pct, total_frames):
        labels_filtered = labels.copy()
        for c in np.unique(labels):
            if c == -1:
                continue
            idx = np.where(labels == c)[0]
            if len(idx) / total_frames < min_pct:
                labels_filtered[idx] = -1
        return labels_filtered

    # ===============================
    # PLOT PROJECTIONS
    # ===============================
    def _plot_projections(self, dist_matrix, labels, frames, title):
        mds_coords = MDS(n_components=2, dissimilarity="precomputed", random_state=42).fit_transform(dist_matrix)
        umap_coords = umap.UMAP(metric="precomputed", random_state=42).fit_transform(dist_matrix)
        tsne_coords = TSNE(metric="precomputed", init="random", random_state=42).fit_transform(dist_matrix)
        coords_dict = {"MDS": mds_coords, "UMAP": umap_coords, "t-SNE": tsne_coords}
        unique_labels = np.unique(labels)
        colors = plt.cm.tab10(np.linspace(0, 1, len(unique_labels)))
        fig, axs = plt.subplots(1, 3, figsize=(18, 5))
        for j, (method_name, coords) in enumerate(coords_dict.items()):
            ax = axs[j]
            for i, label in enumerate(unique_labels):
                mask = labels == label
                idx = np.where(mask)[0]
                if label == -1:
                    ax.scatter(coords[mask, 0], coords[mask, 1], c="lightgray", alpha=0.5, s=50)
                    continue
                rep_idx = self._get_representative(dist_matrix, idx)
                normal_idx = [k for k in idx if k != rep_idx]
                ax.scatter(coords[normal_idx, 0], coords[normal_idx, 1], color=colors[i], s=50)
                ax.scatter(coords[rep_idx, 0], coords[rep_idx, 1], color=colors[i], s=150,
                           edgecolor="black", linewidth=1.5)
                for k in idx:
                    ax.text(coords[k, 0], coords[k, 1], str(frames[k]), fontsize=7)
                ax.scatter([], [], color=colors[i], label=f"Cluster {label} (n={len(idx)})")
            ax.set_title(method_name)
            ax.set_xlabel("Dim 1")
            ax.set_ylabel("Dim 2")
        axs[0].legend(loc="center left", bbox_to_anchor=(1, 0.5))
        plt.suptitle(title)
        plt.tight_layout()
        plt.show()

    # ===============================
    # DBSCAN
    # ===============================
    def run_dbscan_analysis(self, eps_dict, min_samples_pct=0.05, min_pct_dyn=0.05, total_frames_dyn):
        for interface in self.df["interface"].unique():
            df_int = self.df[self.df["interface"] == interface]
            frames, dist_matrix = self._build_distance_matrix(df_int)
            n = len(frames)
            eps = eps_dict.get(interface, 0.2)
            min_samples = max(1, int(min_samples_pct * n))
            db = DBSCAN(eps=eps, min_samples=min_samples, metric="precomputed")
            labels = db.fit_predict(dist_matrix)
            labels = self._filter_clusters(labels, min_pct_dyn, total_frames_dyn)
            self.dbscan_labels_dict[interface] = (frames, dist_matrix, labels)
            self._plot_projections(dist_matrix, labels, frames, f"DBSCAN - {interface}")

    # ===============================
    # HIERARCHICAL
    # ===============================
    def run_hierarchical_analysis(self, method="average", threshold=0.3, min_pct=0.05, total_frames):
        for interface in self.df["interface"].unique():
            df_int = self.df[self.df["interface"] == interface]
            frames, dist_matrix = self._build_distance_matrix(df_int)
            Z = linkage(squareform(dist_matrix), method=method)
            plt.figure(figsize=(10, 5))
            dendrogram(Z, labels=frames, color_threshold=threshold)
            plt.title(f"Dendrogram - {interface}")
            plt.xlabel("Frames")
            plt.ylabel("Distance")
            plt.tight_layout()
            plt.show()
            labels = fcluster(Z, t=threshold, criterion="distance")
            labels = self._relabel_clusters_sequential(labels)
            labels = self._filter_clusters(labels, min_pct, total_frames)
            self.hierarchical_labels_dict[interface] = (frames, dist_matrix, labels)
            self._plot_projections(dist_matrix, labels, frames, f"Hierarchical - {interface}")


    # ===============================
    # HEATMAP
    # ===============================
    def plot_distance_heatmap(self, start_frame=None, end_frame=None):

        df = self.df

        # Toutes les frames disponibles
        all_frames = sorted(set(df["frame_i"]).union(set(df["frame_j"])))

        # Sélection d'une sous-partie
        if start_frame is not None:
            all_frames = [f for f in all_frames if f >= start_frame]

        if end_frame is not None:
            all_frames = [f for f in all_frames if f <= end_frame]

        n = len(all_frames)

        # Mapping frame -> index
        frame_to_idx = {f: i for i, f in enumerate(all_frames)}

        # Colormap
        cmap = mpl.cm.get_cmap("viridis").copy()
        cmap.set_bad(color="red")

        for interface in df["interface"].unique():

            # Garder seulement les lignes dans la plage choisie
            df_int = df[
                (df["interface"] == interface) &
                (df["frame_i"].isin(all_frames)) &
                (df["frame_j"].isin(all_frames))
            ]

            # Matrice
            dist_matrix = np.full((n, n), np.nan)

            for _, row in df_int.iterrows():

                i = frame_to_idx[row["frame_i"]]
                j = frame_to_idx[row["frame_j"]]

                dist_matrix[i, j] = row["distance"]
                dist_matrix[j, i] = row["distance"]

            np.fill_diagonal(dist_matrix, 0)

            # Plot
            plt.figure(figsize=(6,5))

            sns.heatmap(
                dist_matrix,
                cmap=cmap,
                vmin=0,
                vmax=1,
                square=True,
                cbar_kws={"label": "Distance"}
            )

            plt.title(f"Interface {interface} - Heatmap")
            plt.xlabel("Frames")
            plt.ylabel("Frames")

            plt.tight_layout()
            plt.show()

    # ===============================
    # HEATMAP
    # ===============================
    def plot_distance_heatmap(self):
        df = self.df
        all_frames = sorted(set(df["frame_i"]).union(set(df["frame_j"])))
        n = len(all_frames)
        frame_to_idx = {f: i for i, f in enumerate(all_frames)}
        cmap = mpl.cm.get_cmap("viridis").copy()
        cmap.set_bad(color="red")
        for interface in df["interface"].unique():
            df_int = df[df["interface"] == interface]
            dist_matrix = np.full((n, n), np.nan)
            for _, row in df_int.iterrows():
                i = frame_to_idx[row["frame_i"]]
                j = frame_to_idx[row["frame_j"]]
                dist_matrix[i, j] = dist_matrix[j, i] = row["distance"]
            np.fill_diagonal(dist_matrix, 0)
            plt.figure(figsize=(6,5))
            sns.heatmap(dist_matrix, cmap=cmap, vmin=0, vmax=1, square=True, cbar_kws={"label": "Distance"})
            plt.title(f"Interface {interface} - Heatmap")
            plt.xlabel("Frames")
            plt.ylabel("Frames")
            plt.tight_layout()
            plt.show()

    # ===============================
    # CLUSTERMAP
    # ===============================
    def plot_distance_clustermap(self):
        df = self.df
        all_frames = sorted(set(df["frame_i"]).union(set(df["frame_j"])))
        n = len(all_frames)
        frame_to_idx = {f: i for i, f in enumerate(all_frames)}
        cmap = mpl.cm.get_cmap("viridis").copy()
        cmap.set_under("red")
        for interface in df["interface"].unique():
            df_int = df[df["interface"] == interface]
            dist_matrix = np.full((n, n), np.nan)
            for _, row in df_int.iterrows():
                i = frame_to_idx[row["frame_i"]]
                j = frame_to_idx[row["frame_j"]]
                dist_matrix[i, j] = dist_matrix[j, i] = row["distance"]
            np.fill_diagonal(dist_matrix, 0)
            dist_matrix_plot = np.where(np.isnan(dist_matrix), -1, dist_matrix)
            sns.clustermap(dist_matrix_plot, cmap=cmap, vmin=0, vmax=1,
                           figsize=(6,6), xticklabels=False, yticklabels=False,
                           cbar_kws={"label": "Distance"}, cbar_pos=(0.1, 0.8, 0.03, 0.18))
            plt.title(f"Interface {interface} - Clustermap")
            plt.show()

    # ===============================
    # GET CLUSTER FRAMES TABLE
    # ===============================
    def get_cluster_frames_table(self):
        def build_table(labels_dict):
            table = {}
            for interface, labels_info in labels_dict.items():
                frames, dist_matrix, labels = labels_info
                table[interface] = {}
                for c in np.unique(labels):
                    if c == -1:
                        continue
                    idx = np.where(labels == c)[0]
                    cluster_frames = [frames[i] for i in idx]
                    rep_idx = self._get_representative(dist_matrix, idx)
                    rep_frame = frames[rep_idx]
                    table[interface][c] = {"frames": cluster_frames, "representative": rep_frame}
            return table
        return build_table(self.dbscan_labels_dict), build_table(self.hierarchical_labels_dict)

    # ===============================
    # EXPORT REPRESENTATIVE FRAMES AS PDB
    # ===============================
    def export_representative_frames_pdb(self, traj_file, top_file, output_dir, cluster_tables):
        import mdtraj as md
        import os
        warnings.filterwarnings("ignore")
        os.makedirs(output_dir, exist_ok=True)
        traj = md.load(traj_file, top=top_file)
        for method_name, table in zip(["DBSCAN", "Hierarchical"], cluster_tables):
            for interface, clusters in table.items():
                for cluster_label, info in clusters.items():
                    rep_frame = info["representative"]
                    filename = os.path.join(
                        output_dir,
                        f"{method_name}_{interface}_cluster{cluster_label}_rep.pdb"
                    )
                    traj[rep_frame].save_pdb(filename)
                    print(f"Exported {filename}")

    # ===============================
    # SAVE CLUSTER STATISTICS CSV
    # ===============================
    def save_cluster_stats_csv(self, output_prefix="cluster_stats"):
        """
        Save statistics for each cluster (DBSCAN and Hierarchical) to CSV
        with mean, std, quartiles, and 95% confidence interval.
        """
        cluster_dicts = [("DBSCAN", self.dbscan_labels_dict),
                         ("Hierarchical", self.hierarchical_labels_dict)]
        
        for method_name, labels_dict in cluster_dicts:
            all_stats = []
            for interface, labels_info in labels_dict.items():
                frames, dist_matrix, labels = labels_info
                for c in np.unique(labels):
                    if c == -1:
                        continue
                    idx = np.where(labels == c)[0]
                    sub_dist = dist_matrix[np.ix_(idx, idx)]
                    mean_vals = sub_dist.mean(axis=1)
                    n = len(mean_vals)
                    mean = mean_vals.mean()
                    std = mean_vals.std(ddof=1)
                    q1 = np.percentile(mean_vals, 25)
                    median = np.percentile(mean_vals, 50)
                    q3 = np.percentile(mean_vals, 75)
                    # 95% CI
                    sem = std / np.sqrt(n)
                    ci95 = stats.t.interval(0.95, df=n-1, loc=mean, scale=sem)
                    
                    all_stats.append({
                        "interface": interface,
                        "cluster": c,
                        "n_frames": n,
                        "mean_distance": mean,
                        "std_distance": std,
                        "Q1": q1,
                        "median": median,
                        "Q3": q3,
                        "CI95_lower": ci95[0],
                        "CI95_upper": ci95[1]
                    })
            df_stats = pd.DataFrame(all_stats)
            filename = f"{output_prefix}_{method_name}.csv"
            df_stats.to_csv(filename, index=False)
            print(f"{method_name} stats saved to {filename}")
            

    def get_cluster_representatives(self):
        """
        Return dict:
        {method: {interface: {cluster: {"frames": [...], "representative": frame}}}}
        """

        def build(labels_dict):
            out = {}

            for interface, labels_info in labels_dict.items():
                frames, dist_matrix, labels = labels_info
                out[interface] = {}

                for c in np.unique(labels):
                    if c == -1:
                        continue

                    idx = np.where(labels == c)[0]
                    cluster_frames = [frames[i] for i in idx]

                    rep_idx = self._get_representative(dist_matrix, idx)
                    rep_frame = frames[rep_idx]

                    out[interface][c] = {
                        "frames": cluster_frames,
                        "representative": rep_frame
                    }

            return out

        return {
            "DBSCAN": build(self.dbscan_labels_dict),
            "Hierarchical": build(self.hierarchical_labels_dict)
        }
        
        
    def export_cluster_representatives_csv(self, output_file="cluster_representatives.csv"):

        rows = []

        for method_name, labels_dict in [
            ("DBSCAN", self.dbscan_labels_dict),
            ("Hierarchical", self.hierarchical_labels_dict)
        ]:

            for interface, labels_info in labels_dict.items():
    
                frames, dist_matrix, labels = labels_info

                for c in np.unique(labels):
                    if c == -1:
                        continue

                    idx = np.where(labels == c)[0]
                    rep_idx = self._get_representative(dist_matrix, idx)
                    rep_frame = frames[rep_idx]

                    rows.append({
                        "method": method_name,
                        "interface": interface,
                        "cluster": c,
                        "n_frames": len(idx),
                        "representative_frame": rep_frame,
                        "frames": ",".join(map(str, [frames[i] for i in idx]))
                    })

        df = pd.DataFrame(rows)
        df.to_csv(output_file, index=False)
        print(f"Saved {output_file}")
        
    def _relabel_clusters_sequential(self, labels):
        """
        Remap hierarchical cluster labels into:
        cluster0, cluster1, cluster2... in order of appearance
        """

        mapping = {}
        new_labels = np.zeros_like(labels)

        current = 0

        for i, l in enumerate(labels):
            if l == -1:
                new_labels[i] = -1
            else:
                if l not in mapping:
                    mapping[l] = current
                    current += 1
                new_labels[i] = mapping[l]

        return new_labels
