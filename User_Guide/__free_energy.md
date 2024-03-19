[User Guide home](Manual.md)
# Free energy

## 1. Between amino acids forming apair - Molecular mecanics

This method use AMBERSB14 ([Maier *et al.*, 2015](https://pubs.acs.org/doi/10.1021/acs.jctc.5b00255)) and CHARMM36 ([Huang *et al.*, 2013](https://doi.org/10.1021/acs.jctc.5b00255)) force fields.



### 1.x. References
- Maier, J. A. et al. ff14sb: improving the accuracy of protein side chain and backbone parameters from ff99sb. *J. Chem. Theory Comput.* 11, 3696–3713 (2015). [https://pubs.acs.org/doi/10.1021/acs.jctc.5b00255](https://pubs.acs.org/doi/10.1021/acs.jctc.5b00255)
- Huang, J. & MacKerell, A. D. CHARMM36 all-atom additive protein force field: Validation based on comparison to NMR data. *J. Comput. Chem.* 34, 2135–2145 (2013). [https://doi.org/10.1021/acs.jctc.5b00255](https://doi.org/10.1021/acs.jctc.5b00255)




## 2. For the binding interface - *Hunter* method
*Hunter* method from [Potapov *et al.* (2010)](https://doi.org/10.1186/1471-2105-11-374) and [Cohen *et al.* (2009)](https://doi.org/10.1371/journal.pcbi.1000470).

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
