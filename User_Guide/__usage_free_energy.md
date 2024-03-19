[User Guide home](Manual.md)
# Calculate free energy

To be able to calculate the Coulomb and Lennard-Jones energies, you must download the AMBER ff14SB and CHARMM36 force fields:
- [https://github.com/openmm/openmm/blob/master/wrappers/python/openmm/app/data/amber14/protein.ff14SB.xml](https://github.com/openmm/openmm/blob/master/wrappers/python/openmm/app/data/amber14/protein.ff14SB.xml)
- [https://github.com/openmm/openmm/blob/master/wrappers/python/openmm/app/data/charmm36.xml](https://github.com/openmm/openmm/blob/master/wrappers/python/openmm/app/data/charmm36.xml)


> [!IMPORTANT]  
> Here we use the files provided by [OpenMM](https://openmm.org/), because the data are presented in the same format for each force field.
> This XML format is easier to read and use thant the original files.
> This also makes it easier to compare the energies calculated by this software. For example if you use the same method as wrote in the [cookbook](https://openmm.github.io/openmm-cookbook/latest/notebooks/cookbook/Computing%20Interaction%20Energies.html).

> [!TIP]
> Because of their different file sizes, using AMBER ff14SB is faster than CHARMM36.
> Calculation with ff14SB run for 0.05 s while this time is 0.4 s.
> These are approximate values which may vary depending on your computer configuration.

## 1. Coulomb energy
### 
The Coulomb interaction between two atoms is given by the equation:

$$
\begin{equation}
E_{ij} = \frac{ 1 }{ 4 \pi \varepsilon_0 } \times \frac{ q_i q_j }{ \varepsilon_r d_{ij} }
\end{equation}
$$

| Term | Signification | Unit |
| ---- | ------------- | ---- |
| $\frac{ 1 }{ 4 \pi \varepsilon_0 }$ | Electric conversion factor. | $kJ.mol^{-1}.nm.e^{-2}$ |
| $q_i$ and $q_j$                     | Charges of atom *i* and *j*. They come from the chosen force field in the *<Atom charge="xxxx"* lines, where *xxxx* is the charge value. | |
| $d_{ij}$                            | Distance between atom *i* and *j*. | nm |

> [!NOTE]  
> The electric conversion factor is set as 138.935458 $kJ.mol^{-1}.nm.e^{-2}$. For more information, please refere to [GROMACS](https://www.gromacs.org/) documentation on [molecular quantities](https://manual.gromacs.org/current/reference-manual/definitions.html#md-units).

The Coulomb energy for a given amino acids pair is calculated by summing all $E_{ij}$, for all *i* atom in residue 1 and all *j* atom in residue 2.

$$
\begin{equation}
E_{Coulomb} = \sum_{i=1,j=1}^{N} E_{ij}
\end{equation}
$$


### 1.1. Coulomb with "hard" cutoff


### 1.2. Coulomb with reaction field


## 2. Lennard-Jones energy


## 3. *Hunter*
> "You should enjoy the little detours to the fullest. Because that's where you'll find things more important than what you want."
> Ging Freecs, *Hunter X Hunter*
