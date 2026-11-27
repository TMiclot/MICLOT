# MDQ: Molecular Dynamics DockQ

A Python toolkit to analyse protein–protein (or protein–chain) interfaces across molecular dynamics (MD) trajectories.

MDQ uses a symmetric, DockQ-inspired metric. **No native reference structure is needed**: all frames are compared to each other, pair by pair. The package calculates contact conservation, interface RMSD (iRMSD), ligand RMSD (lRMSD) and an MDQ similarity score. It then clusters and visualises the trajectory frames.

---

## Contents
  - [1. Fraction of Shared Contact (FSC)](#1-fraction-of-shared-contact-fsc)
  - [2. Weighted interface RMSD (iRMSD)](#2-weighted-interface-rmsd-irmsd)
  - [3. Weighted ligand RMSD (lRMSD)](#3-weighted-ligand-rmsd-lrmsd)
  - [4. MDQ score](#4-mdq-score)
  - [Score interpretation](#score-interpretation)
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

- the clustering algorithm and all its parameters;
- the software versions and any change made to `mdq.py`.
