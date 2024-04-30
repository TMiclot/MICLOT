[User Guide home](Manual.md)
# Free energy


## 1. Between amino acids forming a pair - Molecular mecanics

This method use AMBERSB14 ([Maier *et al.*, 2015](https://pubs.acs.org/doi/10.1021/acs.jctc.5b00255)) and CHARMM36 ([Huang *et al.*, 2013](https://doi.org/10.1021/acs.jctc.5b00255)) force fields. To be able to calculate the Coulomb and Lennard-Jones energies, you must download the force fields files here:
- [https://github.com/openmm/openmm/blob/master/wrappers/python/openmm/app/data/amber14/protein.ff14SB.xml](https://github.com/openmm/openmm/blob/master/wrappers/python/openmm/app/data/amber14/protein.ff14SB.xml)
- [https://github.com/openmm/openmm/blob/master/wrappers/python/openmm/app/data/charmm36.xml](https://github.com/openmm/openmm/blob/master/wrappers/python/openmm/app/data/charmm36.xml)


> [!IMPORTANT]  
> Here we use the files provided by [OpenMM](https://openmm.org/), because the data are presented in the same format for each force field.
> This XML format is easier to read and use thant the original files.
> This also makes it easier to compare the energies calculated by this software. For example if you use the same method as wrote in the [cookbook](https://openmm.github.io/openmm-cookbook/latest/notebooks/cookbook/Computing%20Interaction%20Energies.html).


### 1.1. Coulomb energy

The Coulomb interaction between two atoms is given by the equation:

$$
\begin{equation}
E_{ij} = \frac{ 1 }{ 4 \pi \varepsilon_0 } \times \frac{ q_i q_j }{ \varepsilon_r r_{ij} }
\end{equation}
$$

| Term | Signification | Unit |
| ---- | ------------- | ---- |
| $\frac{ 1 }{ 4 \pi \varepsilon_0 }$ | Electric conversion factor. | $kJ.mol^{-1}.nm.e^{-2}$ |
| $q_i$ and $q_j$                     | Charges of atom *i* and *j*. They come from the chosen force field in the *<Atom charge="xxxx"* lines, where *xxxx* is the charge value. | |
| $r_{ij}$                            | Distance between atom *i* and *j*. | nm |
| $\varepsilon_r$                     | Solute dielectric constant. | |
| $E_{ij}$                            | Coulomb energy.| kJ/mol |

> [!NOTE]  
> - The electric conversion factor is set as 138.935458 $kJ.mol^{-1}.nm.e^{-2}$. For more information, please refere to [GROMACS](https://www.gromacs.org/) documentation on [molecular quantities](https://manual.gromacs.org/current/reference-manual/definitions.html#md-units).
> - By default, the solute dielectric constant is generaly set to 1.

The total Coulomb energy for a given amino acids pair is calculated by summing all $E_{ij}$, for all *i* atom in residue 1 and all *j* atom in residue 2.

$$
\begin{equation}
E_{Coulomb} = \sum_{i=1,j=1}^{N} E_{ij}
\end{equation}
$$


#### 1.1.1. Coulomb with "hard" cutoff

When we're interested in short-term interactions, it's better not to consider all the interactions between atoms.
In this case, it is necessary to apply a cutoff. This means that all pairs of atoms *ij* separated by a distance greater than the cutoff are not taken into account in the calculation. Finally, only pairs of atoms separated by a distance less than or equal to the cutoff will be used to calculate the Coulomb energy.

If $r_{ij} \leq$ cutoff, then $E_{ij}$ is calculated using the previous equation.
But, if $r_{ij} >$ cutoff, then $E_{ij} = 0$.


#### 1.1.2. Coulomb with cutoff, using reaction field

Another possibility for using a cutoff is to use the reaction field approximation ([Tironi *et al.*, 1995](https://doi.org/10.1063/1.469273)). 
This method considers thant everything above the cutoff is a constant dielectric environment. Coulomb's equation is rewritten as follows:

$$
\begin{equation}
E_{ij} = \frac{ 1 }{ 4 \pi \varepsilon_0 } \times \frac{q_i q_j}{\varepsilon_r} \left( \frac{1}{r_{ij}} + k_{fr} \space r_{ij}^2 - c_{rf} \right)
\end{equation}
$$

The components of the equation are calculated according to the following:

$$
\begin{align}
k_{fr} & = \frac{1}{r_{cutoff}^3} \times \frac{\varepsilon_{solvent} - \varepsilon_r}{2\varepsilon_{solvent} + \varepsilon_r} \\
c_{rf} & = \frac{1}{r_{cutoff}} \times \frac{3 \varepsilon_{solvent} }{2\varepsilon_{solvent} + \varepsilon_r}
\end{align}
$$

| Term | Signification | Unit |
| ---- | ------------- | ---- |
| $\varepsilon_{solvent}$ | Solvent dielectric constant. | |
| $r_{cutoff}$ | Cutoff distance. | nm |

> [!NOTE]  
> - By default, the solvent dielectric constant is generaly set to 78.5.
> - By default, the solute dielectric constant is generaly set to 1.



### 1.2. Lennard-Jones energy

The Lennard-Jones interaction between two atoms is given by the equation:

$$
\begin{equation}
E_{ij} = 4 \varepsilon_{ij} \left( \left( \frac{\sigma_{ij}}{r_{ij}} \right)^{12} - \left( \frac{\sigma_{ij}}{r_{ij}} \right)^6 \right)
\end{equation}
$$

The components of the equation are calculated according to Lorentz-Berthelot rules:

$$
\begin{align}
\sigma_{ij} & = \frac{\sigma_j + \sigma_j}{2} \\
\varepsilon_{ij} & = \sqrt{ \varepsilon_{i} \varepsilon_{j} }
\end{align}
$$

| Term | Signification | Unit |
| ---- | ------------- | ---- |
| $\sigma_{ij}$                           | Finite distance at which the inter-particle potential is zero, for the pair *ij*. | nm |
| $\sigma_{i}$ and $\sigma_{j}$           | Finite distance at which the inter-particle potential, for *i* or *j*, is zero. Values come from the chosen force field in the *<NonbondedForce* section. | nm |
| $\varepsilon_{ij}$                      | Depth of the potential wall for the pair *ij*. |  |
| $\varepsilon_{i}$ and $\varepsilon_{j}$ | Depth of the potential wall for *i* or *j*. Values come from the chosen force field in the *<NonbondedForce* section. |  |
| $r_{ij}$                                | Distance between atom *i* and *j*. | nm |
| $E_{ij}$                                | Lennard-Jones energy.| kJ/mol |

The total Lennard-Jones energy for a given amino acids pair is calculated by summing all $E_{ij}$, for all *i* atom in residue 1 and all *j* atom in residue 2.

$$
\begin{equation}
E_{Lennard-Jones} = \sum_{i=1,j=1}^{N} E_{ij}
\end{equation}
$$


#### 1.2.1. Lennard-Jones with "hard" cutoff

When we're interested in short-term interactions, it's better not to consider all the interactions between atoms.
In this case, it is necessary to apply a cutoff. This means that all pairs of atoms *ij* separated by a distance greater than the cutoff are not taken into account in the calculation. Finally, only pairs of atoms separated by a distance less than or equal to the cutoff will be used to calculate the Lennard-Jones energy.

If $r_{ij} \leq$ cutoff, then $E_{ij}$ is calculated using the previous equation.
But, if $r_{ij} >$ cutoff, then $E_{ij} = 0$.


#### 1.2.2. Lennard-Jones with cutoff, using switching function

Another possibility for using a cutoff is to use a [switching function](https://manual.gromacs.org/current/reference-manual/functions/nonbonded-interactions.html#modified-non-bonded-interactions). 
The switching function involves two distances $r_{switch}$ and $r_{cutoff}$, where $r_{switch} < r_{cutoff}$. Between these two distances, the energy is modified by a factor S that causes the energy to tend towards 0 when the integration distance reaches the cutoff. Because the 

If $r_{ij} < r_{switch}$, the energie is not modified. If $r_{ij} \geq r_{cutoff}$, the returned energy is 0.
But if $r_{switch} \leq r_{ij} < r_{cutoff}$ the energy is calculated as follow:

$$
\begin{equation}
E_{switch} = E_{ij} * S
\end{equation}
$$

Where the equation of *S* is:

$$
\begin{align}
S & = 1 - 6x^5 + 15x^4 -10x^3 \\
x & = \frac{r_{ij} - r_{switch}}{r_{cutoff} - r_{switch}}
\end{align}
$$

| Term | Signification | Unit |
| ---- | ------------- | ---- |
| $r_{switch}$ | Switchind distance. It muste be lower than the cutoff distance. | nm |
| $r_{cutoff}$ | Cutoff distance. | nm |
| $r_{ij}$  | Distance between atom *i* and *j*. | nm |
| $E_{ij}$ | Lennard-Jones energy.| kJ/mol |

> [!WARNING]
> The use of a switching potential increase the values in the switching region: from $r_{switch}$ to $r_{cutoff}$.
> As a result, the energy value is underestimated, but should remain acceptable.



### 1.3. Total pair energy

The total energy of a pair is given by the sum of the Coulomb energy and the Lennard-Jones energy.

$$
\begin{equation}
E_{Total} = E_{Lennard-Jones} + E_{Coulomb}
\end{equation}
$$


### 1.4. References

- Maier, J. A. et al. ff14sb: improving the accuracy of protein side chain and backbone parameters from ff99sb. *J. Chem. Theory Comput.* 11, 3696–3713 (2015). [https://pubs.acs.org/doi/10.1021/acs.jctc.5b00255](https://pubs.acs.org/doi/10.1021/acs.jctc.5b00255)
- Huang, J. & MacKerell, A. D. CHARMM36 all-atom additive protein force field: Validation based on comparison to NMR data. *J. Comput. Chem.* 34, 2135–2145 (2013). [https://doi.org/10.1021/acs.jctc.5b00255](https://doi.org/10.1021/acs.jctc.5b00255)
- Tironi, I. G., Sperb, R., Smith, P. E. & Van Gunsteren, W. F. A generalized reaction field method for molecular dynamics simulations. *The Journal of Chemical Physics* 102, 5451–5459 (1995). [https://doi.org/10.1063/1.469273](https://doi.org/10.1063/1.469273)





## 2. Protein binding: *4-distance description* method

*Hunter* method from [Potapov *et al.* (2010)](https://doi.org/10.1186/1471-2105-11-374) and [Cohen *et al.* (2009)](https://doi.org/10.1371/journal.pcbi.1000470) ****

> The favourable energies were accumulated in the $E_{lja}$ term and the repulsive ones in the $E_{ljr}$ term. To avoid excessive repulsion due to close placement of atoms during side chain optimization, the repulsive term $E_{ljr}$ was linearized at a cutoff distance $d_{ij}$ < 0.89 [46].

$$
\begin{equation}
	E =  w_{ScSc}E_{ScSc} + w_{ScMc}E_{ScMc} + w_{rot}E_{rot} + w_{lj}E_{lj}
\end{equation}
$$

With $w_{ScSc}$, $w_{ScMc}$, $w_{rot}$, $w_{lj}$ are the weights for each term.
> The final set of weights was 0.13, 0.13, 0.33, 0.41 for $w_{ScSc}$ , $w_{ScMc}$ , $w_{rot}$, $w_{lj}$ , respectively.

The components of the equation are calculated according to the following:

$$
\begin{align}
	E_{ScSc} & =  \sum_{i}^{M} \sum_{j,j>i}^{M} - \ln \frac{P_{real} (\{dist\}\mid AA) \times P_{real}(AA) }{P_{rand} (\{dist\}\mid AA) \times P_{rand}(AA) }  \\
	E_{ScMc} & =  \sum_{i}^{M} \sum_{j,j\neq i,j\neq i+1}^{M} - \ln \frac{P_{real} (\{dist\}\mid AX) \times P_{real}(AX) }{P_{rand} (\{dist\}\mid AX) \times P_{rand}(AX) } \\
	E_{rot} & = \sum_{i=1}^{M} - \ln (P_i^{rot} N_i^{rot}) \\
	E_{lj} & = \sum_{i=1}^{N} \sum_{j,j>1}^{N} \varepsilon_{ij} \left\lfloor \left( \frac{r_{ij}^{\mathrm{min}}}{d_{ij}} \right)^{12} -2 \left( \frac{r_{ij}^{\mathrm{min}} }{ d_{ij}} \right)^{6} \right\rfloor
\end{align}
$$



And the terms of the equations are defined as:
| Term                     | Calculated as                              | Definition |
| ------------------------ | ------------------------------------------ | ---------- |
| *M*                      |                                            | Number of residues. |
| $P(\{dist\}\mid AA)$     |                                            | Probability of observing the *ScSc* 4-distance combination (*{dist}*) for a given residue pair in real structures ($P_{real}$) or in random structures ($P_{rand}$). |
| $P(AA)$                  |                                            | Probability to observe a *ScSc* contact for a given residue pair in protein structures, in real structures ($P_{real}$) or in random structures ($P_{rand}$). |
| $P(\{dist\}\mid AX)$     |                                            | Probability of observing the *ScMc* 4-distance combination (*{dist}*) for a given residue pair in real structures ($P_{real}$) or in random structures ($P_{rand}$). |
| $P(AX)$                  |                                            | Probability to observe a *ScMc* contact for residue A in real structures ($P_{real}$) or in random structures ($P_{rand}$). |
| $P_i^{rot}$              |                                            | Rotamer probability as taken from the backbone-dependent rotamer library from [Dunbrack and Cohen (1997)](https://doi.org/10.1002/pro.5560060807). |
| $N_i^{rot}$              |                                            | Number of rotamers for a modelled residue. |
| $r_{ij}^{\mathrm{min}}$  | $r_i + r_j$                                | Distance between two atoms at the minimum of the potential. |
| $d_{ij}$                 |                                            | Measured distance between two atoms. |
| $\varepsilon_{ij}$       | $\sqrt{\varepsilon_{i} + \varepsilon_{j}}$ | Individual $\varepsilon_{i}$ values for each atom class were defined as in the CHARMM19 parameter set from [Neria *et al.* (1996)](https://doi.org/10.1063/1.472061). |


### References

- Potapov, V., Cohen, M., Inbar, Y. & Schreiber, G. Protein structure modelling and evaluation based on a 4-distance description of side-chain interactions. *BMC Bioinformatics* 11, 374 (2010). [https://doi.org/10.1186/1471-2105-11-374](https://doi.org/10.1186/1471-2105-11-374)
- Cohen, M., Potapov, V. & Schreiber, G. Four Distances between Pairs of Amino Acids Provide a Precise Description of their Interaction. *PLoS Comput Biol* 5, e1000470 (2009). [https://doi.org/10.1371/journal.pcbi.1000470](https://doi.org/10.1371/journal.pcbi.1000470)
- Dunbrack, R. L. & Cohen, F. E. Bayesian statistical analysis of protein side‐chain rotamer preferences. *Protein Science* 6, 1661–1681 (1997). [https://doi.org/10.1002/pro.5560060807](https://doi.org/10.1002/pro.5560060807)
- Neria, E., Fischer, S. & Karplus, M. Simulation of activation free energies in molecular systems. *The Journal of Chemical Physics* 105, 1902–1921 (1996). [https://doi.org/10.1063/1.472061](https://doi.org/10.1063/1.472061)





## 3. Protein binding: *Contacts-based* method

> [!IMPORTANT]  
> This is an implementation of the contacts-based method by [Vangone *et al.* (2015)](https://doi.org/10.7554/eLife.07454). This code is different and is not related to the one from [PRODIGY on GitHub](https://github.com/haddocking/prodigy/).


### 3.1. IC-NIS model

The IC-NIS model is based on the following equation:

$$
\begin{align}
	\Delta G & = (- 0.09459 \times IC_{charged/charged}) - (0.10007 \times IC_{charged/apolar}) + (0.19577 \times IC_{polar/polar}) - (0.22671 \times IC_{polar/apolar}) + (0.18681 \times NIS_{apolar}^{\\%}) + (0.13810 \times NIS_{charged}^{\\%}) - 15.9433

	Kd & = e^{\frac{\Delta G}{RT}}
\end{align}
$$


And the terms of the equations are defined as:
| Term                   | Definition | Unit |
| ---------------------- | ---------- | ---- |
| $\Delta G$             | Binding affinity. | kcal/mol |
| $IC_{charged/charged}$ | Number of contacts between two charged residues. |  |
| $IC_{charged/apolar}$  | Number of contacts between charged and apolar residues. |  |
| $IC_{polar/polar}$     | Number of contacts between two polar residues. |  |
| $NIS_{apolar}^{\\%}$    | % of apolar residue in the NIS (see 3.3.) | % |
| $NIS_{charged}^{\\%}$  | % of charged residue in the NIS (see 3.3.) | % |
| Kd                     | Dissociation constant. | M |
| R                      | Ideal gas constant. <br/> Value: 0.0019858775 | kcal/mol |
| T                      | Temperature | Kelvin |


### 3.2. Compute contacts

The contacts are computed using the [mdtraj.compute_contacts](https://mdtraj.org/1.9.4/api/generated/mdtraj.compute_contacts.html) command. It return the closest distance between any two heavy atoms in the residues pair.


### 3.3. NIS-Interface identification

To be able to reproduce the same values as returned by prodigy this code use [FreeSASA](https://freesasa.github.io/).
Because their is diffrences with [MDTraj](https://www.mdtraj.org/1.9.7/index.html) (used to calculate SASA in other functions).

1. Freesasa use Lee-Richards algorithm, but MDTraj use Shrake-Rupley algorithm.
2. Freesasa use the atom radii from [NACCESS](http://www.bioinf.manchester.ac.uk/naccess/) based on their atom name. On the contrary MDTRaj use radii based on the element.

This two diffrences lead to minor changes in the ASA of the residues. But this changes drastically modify the final results.

> [!WARNING]
> Based on the code of the function [analyse_nis](https://github.com/haddocking/prodigy/blob/main/src/prodigy_prot/predict_IC.py), filtering is done only on the relative SASA (also named *rASA*) of the complex (bound form) and take value equal or greater than 5%. This filtering din't seem to correspond to the NIS definition by [Kastritis *et al.*](https://doi.org/10.1016/j.jmb.2014.04.017), where the NIS is defined when the difference between the relative SASA of the bonded and unbonded forms is lower than 5% and, greater than 5% for the interface: <br/> 
> $\Delta rASA = rASA_{in \space monomer} - rASA_{in \space complex}$
> 
> Here, the relative SASA is the ASA of a residue divided by the maximum ASA of the residue as refered in NACCESS: <br/>
> $rASA = \frac{ASA_{residue \space in \space complex}}{MaxASA_{NACCESS}}$


### 3.4. Amino acid properties use to identify contact type and NIS-interface properties

Amino acid properties are use to identify contact type and NIS-interface properties. Please note that their is diffrences for: Cys, His, Trp and Tyr.

| Amino acid | Contact | Interface-NIS | Difference |
| ---------- | ------- | ------------- | --------- |
| ALA        | apolar  | apolar        |   |
| CYS        | apolar  | polar         | x |
| GLU        | charged | charged       |   |
| ASP        | charged | charged       |   |
| GLY        | apolar  | apolar        |   |
| PHE        | apolar  | apolar        |   |
| ILE        | apolar  | apolar        |   |
| HIS        | charged | polar         | x |
| LYS        | charged | charged       |   |
| MET        | apolar  | apolar        |   |
| LEU        | apolar  | apolar        |   |
| ASN        | polar   | polar         |   |
| GLN        | polar   | polar         |   |
| PRO        | apolar  | apolar        |   |
| SER        | polar   | polar         |   |
| ARG        | charged | charged       |   |
| THR        | polar   | polar         |   |
| TRP        | apolar  | polar         | x |
| VAL        | apolar  | apolar        |   |
| TYR        | apolar  | polar         | x |


### 3.5. References

- Vangone, A. & Bonvin, A. M. Contacts-based prediction of binding affinity in protein–protein complexes. *eLife* 4, e07454 (2015). [https://doi.org/10.7554/eLife.07454](https://doi.org/10.7554/eLife.07454)
- Kastritis, P. L., Rodrigues, J. P. G. L. M., Folkers, G. E., Boelens, R. & Bonvin, A. M. J. J. Proteins Feel More Than They See: Fine-Tuning of Binding Affinity by Properties of the Non-Interacting Surface. *Journal of Molecular Biology* 426, 2632–2652 (2014).
 [https://doi.org/10.1016/j.jmb.2014.04.017](https://doi.org/10.1016/j.jmb.2014.04.017)
