# To do list

## Interaction types

- [ ] Write documentation on structural parameters

- [x] C-bonds 
- [x] C5 H-bond 
- [x] hydrophobic / hydrophilic-hydrophobic repulsion 
- [x] charge clash/repulsion 
- [x] salt bridge 
- [x] H-bond
- [x] van der waals
- [x] Pi-Amino
- [x] Charge-aromatic (pi/quadrupole/intermediate-charge) 
- [x] aromatic-aromatic (Pi-Pi, Pi-quadrupole, T-stacking) 
- [x] ARG-ARG stacking
- [x] Pi H-bond (if not involve in amino-PI or charge-aromatic !)
- [x] N--> $\pi$*
- [ ] S/Se H-bond & S/Se chalogen bond


## Cys-Cys

- [x] disulfice & energy
- [x] diselenide
- [x] selenosulfide


## Energy

- [x] Coulomb & LJ
- [ ] "PRODIGY"
- [ ] "HUNTER"


## Structure properties

- [x] Sequence & Secondary structure
- [x] Protein region
- [ ] physico chemical properties (pI, aromaticity, etc) [https://biopython.org/docs/1.75/api/Bio.SeqUtils.ProtParam.html](https://biopython.org/docs/1.75/api/Bio.SeqUtils.ProtParam.html)


## Tools
- [x] MDTraj chain ID <--> PDB chain name
- [x] Minimization (openMM)
- [ ] preparation (PDB2PQR)


## Neighbors

- [ ] AA
- [ ] AA pairs


## Surface / interface shape

Steps:

- [ ] protein to mesh
- [ ] mesh to matrix
- [ ] matrix analysis


***

# outputs for clened data (NOT brut data!)

```
====================
OUTPUT COMPLEX
    pandas table exported to CSV
    column names are below
    
    voir https://biopython.org/docs/1.75/api/Bio.SeqUtils.ProtParam.html (ne semble nécessite internet pour les calcusl)
====================
number_AA_complex
number_AA_receptor
number_AA_ligand

SASA_complex
SASA_receptor
SASA_ligand
interface

pI_complex
pI_receptor
pI_ligand

molecular_weight_complex
molecular_weight_receptor
molecular_weight_ligand

charge_complex
charge_receptor
charge_ligand

instability_index_complex
instability_index_receptor
instability_index_ligand

GRAVY_complex
GRAVY_receptor
GRAVY_ligand

extinction_coefficients_complex
extinction_coefficients_receptor
extinction_coefficients_ligand

...


====================
OUTPUT RES-RES NEIGHBOR
====================

====================
OUTPUT PAIR-PAIR NEIGHBOR
====================

====================
OUTPUT RES-RES INTERACTION
    pandas table exported to CSV
    column names are below
====================
residue_A_id
residue_A_resSeq
residue_A_name
residue_A_chain_id

... idem for residue B

force_vdw_AMBER
force_lj_AMBER
force_vdw_CHARMM
force_lj_CHARMM

interaction_   ... 0/1 (0:False, 1:True)
interaction_
interaction_
interaction_
interaction_
interaction_
interaction_
interaction_
interaction_
interaction_
interaction_
interaction_
interaction_
```