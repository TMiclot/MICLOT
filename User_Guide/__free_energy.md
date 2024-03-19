[User Guide home](Manual.md)
# Free energy

## 1. Between amino acids forming apair - Molecular mecanics

This method use AMBERSB14 ([Maier *et al.*, 2015](https://pubs.acs.org/doi/10.1021/acs.jctc.5b00255)) and CHARMM36 ([Huang *et al.*, 2013](https://doi.org/10.1021/acs.jctc.5b00255)) force fields. To be able to calculate the Coulomb and Lennard-Jones energies, you must download the force fields files here:
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

### 1.1. Coulomb energy
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
| $\varepsilon_r$                     | Solute dielectric constant. | |

> [!NOTE]  
> - The electric conversion factor is set as 138.935458 $kJ.mol^{-1}.nm.e^{-2}$. For more information, please refere to [GROMACS](https://www.gromacs.org/) documentation on [molecular quantities](https://manual.gromacs.org/current/reference-manual/definitions.html#md-units).
> - The solute dielectric constant is generaly set to 1.

The Coulomb energy for a given amino acids pair is calculated by summing all $E_{ij}$, for all *i* atom in residue 1 and all *j* atom in residue 2.

$$
\begin{equation}
E_{Coulomb} = \sum_{i=1,j=1}^{N} E_{ij}
\end{equation}
$$


#### 1.1.1. Coulomb with "hard" cutoff
When we're interested in short-term interactions, it's better not to consider all the interactions between atoms.
In this case, it is necessary to apply a cutoff. This means that all pairs of atoms *ij* separated by a distance greater than the cutoff are not taken into account in the calculation. Finally, only pairs of atoms separated by a distance less than or equal to the cutoff will be used to calculate the Coulomb energy.

If $d_{ij} \leq$ cutoff, then $E_{ij} = \frac{ 1 }{ 4 \pi \varepsilon_0 } \times \frac{ q_i q_j }{ \varepsilon_r d_{ij} }$.
But, if $d_{ij} >$ cutoff, then $E_{ij} = 0$

#### 1.1.2. Coulomb with reaction field


### 1.2. Lennard-Jones energy


#### 1.3. References
- Maier, J. A. et al. ff14sb: improving the accuracy of protein side chain and backbone parameters from ff99sb. *J. Chem. Theory Comput.* 11, 3696–3713 (2015). [https://pubs.acs.org/doi/10.1021/acs.jctc.5b00255](https://pubs.acs.org/doi/10.1021/acs.jctc.5b00255)
- Huang, J. & MacKerell, A. D. CHARMM36 all-atom additive protein force field: Validation based on comparison to NMR data. *J. Comput. Chem.* 34, 2135–2145 (2013). [https://doi.org/10.1021/acs.jctc.5b00255](https://doi.org/10.1021/acs.jctc.5b00255)
- Tironi, I. G., Sperb, R., Smith, P. E. & Van Gunsteren, W. F. A generalized reaction field method for molecular dynamics simulations. *The Journal of Chemical Physics* 102, 5451–5459 (1995). [https://doi.org/10.1063/1.469273](https://doi.org/10.1063/1.469273)




## 2. For the binding interface - *Hunter* method
> "You should enjoy the little detours to the fullest. Because that's where you'll find things more important than what you want."
> Ging Freecs, *Hunter X Hunter*

*Hunter* method from [Potapov *et al.* (2010)](https://doi.org/10.1186/1471-2105-11-374) and [Cohen *et al.* (2009)](https://doi.org/10.1371/journal.pcbi.1000470) ****

>The favourable energies were accumulated in the $E_{lja}$ term and the repulsive ones in the $E_{ljr}$ term. To avoid excessive repulsion due to close placement of atoms during side chain optimization, the repulsive term $E_{ljr}$ was linearized at a cutoff distance $d_{ij}$ < 0.89 [46].

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
| $P(\{dist\}\mid AA)$         |                                            | Probability of observing the *ScSc* 4-distance combination (*{dist}*) for a given residue pair in real structures ($P_{real}$) or in random structures ($P_{rand}$). |
| $P(AA)$                  |                                            | Probability to observe a *ScSc* contact for a given residue pair in protein structures, in real structures ($P_{real}$) or in random structures ($P_{rand}$). |
| $P(\{dist\}\mid AX)$         |                                            | Probability of observing the *ScMc* 4-distance combination (*{dist}*) for a given residue pair in real structures ($P_{real}$) or in random structures ($P_{rand}$). |
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
